#!/usr/bin/env perl

# This script takes a cDNA cluster produced by get_pangenes.pl and retrieves
# the collinearity data evidence supporting it, inferred from whole genome
# alignments, contained in mergedpairs.tsv.gz (presorted with -k1,1 -k4,4nr).
# It can print a representative cluster sequence with longest mode length (-s).
#
# Optionally it can also liftover/suggest fixes to the gene models based on the 
# pangene consensus (-f), this requires gmap (make install_pangenes)

# Copyright [2022]
# EMBL-European Bioinformatics Institute & Estacion Experimental Aula Dei-CSIC

$|=1;

use strict;
use warnings;
use Getopt::Std;
use DB_File;
use Compress::Zlib qw(compress uncompress);
use FindBin '$Bin';
use lib "$Bin/lib";
use pangeneTools qw( check_installed_features 
                     parse_sequence_FASTA_file extract_isoforms_FASTA
                     calc_median calc_mode get_outlier_cutoffs );

my $MINPAIRPECNONOUTLIERS = 0.25;
my $MINLIFTIDENTITY = 95.0;
my $GMAPARAMS = '-t 1 -2 -z sense_force -n 1 -F ';

my @FEATURES2CHECK = (
  'EXE_BEDTOOLS', 'EXE_GMAP', 'EXE_GZIP'
);

$ENV{'LC_ALL'}   = 'POSIX';

my $GZIPBIN = $ENV{'EXE_GZIP'} || 'gzip';
my $BEDTOOLSBIN = $ENV{'EXE_BEDTOOLS'} || 'bedtools';
my $GMAPBIN = $ENV{'EXE_GMAP'} || 'gmap';
$GMAPBIN .= " $GMAPARAMS ";


my ($INP_dir,$INP_clusterfile,$INP_noraw,$INP_fix) = ( '', '', 0 , 0 );
my ($INP_verbose,$INP_appendGFF,$INP_modeseq,$INP_outdir) = (0,0, '', '');
my ($cluster_list_file,$cluster_folder,$gdna_clusterfile, $genome_file);
my ($gene_id, $hom_gene_id, $homology_type, $species, $hom_species);
my ($isof_id, $overlap, $coords, $hom_coords, $full_id, $hom_full_id);
my ($line, $segment, $hom_segment, $dummy, $TSVdata, $cmd, $cDNA);
my (%opts,%TSVdb, @sorted_ids, @pairs, @segments, @ref_names);
my (%seen, %overlap, %cluster_gene_id, %fullid2id, %gene_length, %genome_coords);

getopts('hvacfnr:s:o:d:i:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-c print credits and checks installation\n";
  print "-d directory produced by get_pangenes.pl        (example: -d /path/data_pangenes/..._algMmap_,\n";
  print "                                                 genomic sequences usually one folder up)\n";
  print "-i cdna/cds .fna file as in .cluster_list file  (example: -i gene:ONIVA01G52180.cdna.fna)\n";
  print "-s append mode isoform sequence to file         (optional, example: -s isoforms.fna)\n";
  print "-r comma-sep string with reference taxa         (optional, requires -s)\n";
  print "-n do not print raw evidence                    (optional)\n";
  print "-f fix gene models and produce GFF              (optional, GFF printed to stdout by default)\n";
  print "-o folder to write GFF output                   (optional, requires -f, 1 file/species)\n";
  print "-a append GFF output                            (optional, requires -f -o)\n";  
  print "Note: reads the compressed merged TSV file in -d\n"; 
  exit(0);
}

if(defined($opts{'c'})) {

  print "\nPrimary citation:\n https://github.com/Ensembl/plant-scripts/pangenes\n";
  print "\nThis software uses external algorithms, please cite them accordingly:\n";
  print " gmap https://doi.org/10.1093/bioinformatics/bti310\n";

  # check all binaries needed by this program and print diagnostic info
  print check_installed_features(@FEATURES2CHECK);
  exit(0);
}

if(defined($opts{'d'})) { 
  $INP_dir = $opts{'d'} 
}
else{ die "# EXIT : need a -d directory\n" }

if(defined($opts{'i'})){  
  $INP_clusterfile = $opts{'i'};
  if($INP_clusterfile !~ /\.cdna\.fna$/ && $INP_clusterfile !~ /\.cds\.fna$/) {
    die "# EXIT : need a .fna cluster filename with parameter -i\n"
  } else {
    $gdna_clusterfile = $INP_clusterfile;
    $gdna_clusterfile =~ s/\.cdna\.fna/.gdna.fna/;
    $gdna_clusterfile =~ s/\.cds\.fna/.gdna.fna/;
  }
}
else{ die "# EXIT : need parameter -i\n" }

if(defined($opts{'s'})){
  $INP_modeseq = $opts{'s'};
 
  if(defined($opts{'r'})) {
    @ref_names = split(/,/,$opts{'r'}); 
    if(!@ref_names) {
      die "# EXIT: cannot parse reference names (-r), " .
        "make sure they have no blanks and are comma-separated\n";
    }
  }
}

if(defined($opts{'n'})){ 
  $INP_noraw = 1 
}

if(defined($opts{'f'})){
  $INP_fix = 1;

  if(defined($opts{'o'})){
    $INP_outdir = $opts{'o'};
    if(!-e $INP_outdir) {
      mkdir($INP_outdir);
    }

    if(defined($opts{'a'})){
      $INP_appendGFF = 1;

      opendir(OUTDIR,$INP_outdir) || 
        die "# EXIT: cannot list $INP_outdir\n";
      my @files = grep{/gff/} readdir(OUTDIR);
      closedir(OUTDIR);

      if(@files) {
        warn "# WARN: appending to GFF files in folder $INP_outdir/\n";
      }
    }
  }
}

if(defined($opts{'v'})){
  $INP_verbose = 1
}



# 1) locate .cluster_list file to check clusterfile is there
opendir(INPDIR,$INP_dir) || 
  die "# ERROR: cannot list $INP_dir , please check -d argument is a valid folder\n";
my @files = grep {/\.cluster_list/} readdir(INPDIR);
closedir(INPDIR);

if(@files) {
  $cluster_list_file = $files[0];
  $cluster_folder = (split(/\.cluster_list/,$cluster_list_file))[0]
} else {
  die "# ERROR: cannot find .cluster_list file in $INP_dir\n";
}

my $clusternameOK = 0;
open(LIST,"<","$INP_dir/$cluster_list_file") ||
  die "# ERROR: cannot read $INP_dir/$cluster_list_file, ".
    "please check -d argument is a valid folder\n";

while(<LIST>) {
  if(/$INP_clusterfile/) {
    $clusternameOK = 1;
  }
}
close(LIST);

if($clusternameOK == 0) {
  die "# ERROR: cannot find $INP_clusterfile in $INP_dir/$cluster_list_file, please correct\n";
}

# 2) parse FASTA headers of input cluster to extract gene names and check sequence lengths
#    (note that chr coords are parsed only for 1st isoform)
my ( $seq, %isof_len, %isof_seq, %isof_header, @len);

my ( $ref_geneid, $ref_fasta, $ref_isof_coords, $ref_taxon ) = 
  parse_sequence_FASTA_file( "$INP_dir/$cluster_folder/$INP_clusterfile" , 1);

foreach $gene_id (sort @$ref_geneid) {

  # sorted gene ids
  $cluster_gene_id{$gene_id} = 1;
  push(@sorted_ids, $gene_id);

  # length stats
  foreach $seq (split(/\n/,$ref_fasta->{$gene_id})) {
    if($seq =~ /^>(\S+)/) {
      $isof_id = $1;
      $isof_header{$gene_id}{$isof_id} = $seq;
      next;
    }
    $isof_len{$gene_id}{$isof_id} += length($seq);
    $isof_seq{$gene_id}{$isof_id} .= $seq;
  }
  push(@len, $isof_len{$gene_id}{$isof_id});
}

my ($median_length, $cutoff_low_length, $cutoff_high_length) =
  get_outlier_cutoffs( \@len , $INP_verbose );

# note: might be "mode2 mode1" for bimodal data
my @modes_length = calc_mode( \@len );

printf("# isoform length in cluster: median=%1.0f mode(s): %s\n\n",
  $median_length, join(',',@modes_length));

if($INP_modeseq) {

  my ($mode_gene_id, $mode_isof_id);

  foreach $gene_id (sort @$ref_geneid) {
    foreach $isof_id (keys(%{$isof_len{$gene_id}})) {
 
      next if($isof_len{$gene_id}{$isof_id} != $modes_length[0]);
 
      if(!$mode_isof_id || grep(/$ref_taxon->{$gene_id}/,@ref_names)) {
        ($mode_gene_id, $mode_isof_id) = ($gene_id, $isof_id);
      }    
    }    
  }

  open(ISOSEQ,">>",$INP_modeseq) ||
    die "# EXIT: cannot write to $INP_modeseq\n";
  print ISOSEQ "$isof_header{$mode_gene_id}{$mode_isof_id}\n$isof_seq{$mode_gene_id}{$mode_isof_id}\n";
  close(ISOSEQ);

  print "# mode isoform: $mode_gene_id $mode_isof_id [$ref_taxon->{$mode_gene_id}]".
    " (appended to $INP_modeseq)\n";
}

foreach $gene_id (sort @$ref_geneid) {
  foreach $isof_id (keys(%{$isof_len{$gene_id}})) {
    if($isof_len{$gene_id}{$isof_id} < $cutoff_low_length) {

      print "# short isoform: $isof_id $gene_id [$ref_taxon->{$gene_id}] ".
        "length=$isof_len{$gene_id}{$isof_id}\n";

    } if($isof_len{$gene_id}{$isof_id} > $cutoff_high_length) {

      print "# long isoform: $isof_id $gene_id [$ref_taxon->{$gene_id}] ".
       "length=$isof_len{$gene_id}{$isof_id}\n";
    }
  }
} 


# 2.1) parse segment gDNA file if required (1-based coordinates)
my ($ref_geneid_seg,$ref_fasta_seg,$ref_coords_seg,$ref_taxon_seg);

if(-e "$INP_dir/$cluster_folder/$gdna_clusterfile") {

  print "# parsing twin .gdna.fna cluster file $gdna_clusterfile\n\n";

  ( $ref_geneid_seg, $ref_fasta_seg, $ref_coords_seg, $ref_taxon_seg ) =
    parse_sequence_FASTA_file( "$INP_dir/$cluster_folder/$gdna_clusterfile" , 1);

  printf("\n# cluster %s genes = %d genomic segments = %d (%d taxa)\n",
    $INP_clusterfile,
    scalar(keys(%cluster_gene_id)),
    scalar(@$ref_geneid_seg),
    scalar(keys(%$ref_taxon)));

} else {

  printf("\n# cluster %s genes = %d (%d taxa)\n",
    $INP_clusterfile,
    scalar(keys(%cluster_gene_id)),
    scalar(keys(%$ref_taxon)));
}


# 3) parse compressed merged TSV file and feed Berkeley DB
my $mergedTSV = "$INP_dir/mergedpairs.tsv.gz";
my $TSVdb_file = "$INP_dir/mergedpairs.tsv.bdb";

# only first time
if(!-s $TSVdb_file) {

  print "\n# creating database (only first time)\n";

  tie(%TSVdb, 'DB_File', $TSVdb_file, 
    O_RDWR|O_CREAT, 0666, $DB_BTREE) ||
    die "# ERROR: cannot create file $TSVdb_file: $!\n";

  open(TSV, "$GZIPBIN -dc $mergedTSV |") ||
    die "# ERROR: cannot find file $mergedTSV, please check -d argument\n";

  my ($block, $prev_gene_id) = ('', '');
  while(<TSV>) {
    if(/^([^\t]+)/) {
      $gene_id = $1;

      if($gene_id ne $prev_gene_id) {

        if($block) {
          $TSVdb{$prev_gene_id} = compress($block);
        }
        
        $block = $_;
        $prev_gene_id = $gene_id;

      } else {
        $block .= $_;
      }
    }
  }   

  $TSVdb{$gene_id} = compress($block);

  close(TSV);

  print "# done\n";  

} else {

  print "# re-using database\n";

  tie(%TSVdb, 'DB_File', $TSVdb_file,
    O_RDWR, 0666, $DB_BTREE) ||
    die "# ERROR: cannot read file $TSVdb_file: $!\n";
}

# 4) parse TSV blocks for selected genes  
foreach $gene_id (@sorted_ids) {

  next if(!$TSVdb{$gene_id});

  $TSVdata = uncompress($TSVdb{$gene_id});
 
  foreach $line (split(/\n/,$TSVdata)) { 

    #LOC_Os01g01019 LOC_Os01g01019 oryza_sativa_MSU 1065 ortholog_collinear
    #gene:Osir64_01g0000020 gene:Osir64_01g0000020 oryza_sativa_ir64 1065	
    #NULL NULL NULL 100.00 1 Chr1:11217-12435(+);1:19441-20506(+)

    ( $gene_id, $segment, $species, $overlap, $homology_type, 
      $hom_gene_id, $hom_segment, $hom_species, $dummy, 
      $dummy, $dummy, $dummy, $dummy, $dummy, $coords   
    ) = split( /\t/, $line );
 
    if($homology_type eq 'ortholog_collinear') {

      # skip sequences from other clusters
      next if(!$cluster_gene_id{$gene_id} || 
        !$cluster_gene_id{$hom_gene_id});
      next if($ref_taxon->{$gene_id} ne $species || 
        $ref_taxon->{$hom_gene_id} ne $hom_species);

      # concat species & gene_id in case there are repeated gene ids
      $full_id = $species.$gene_id;
      $hom_full_id = $hom_species.$hom_gene_id;
      $fullid2id{$full_id} = $gene_id;
      $fullid2id{$hom_full_id} = $hom_gene_id;

      $cluster_gene_id{$gene_id}++;
      $seen{$full_id}++;
      $seen{$hom_full_id}++;

      # add overlap of gene pair (segments not considered)
      $overlap{$full_id} += $overlap;
      $overlap{$hom_full_id} += $overlap; 

      # compute gene length
      ($coords, $hom_coords) = split(/;/,$coords); 

      if(!$gene_length{$full_id}) {
        $genome_coords{$full_id} = $coords;
        if($coords =~ m/^\S+?:(\d+)-(\d+)\([+-]\)/) {
          $gene_length{$full_id} = 1+$2-$1;
        } else {
          die "# ERROR: cannot parse $coords $_\n";
        }
      }

      if(!$gene_length{$hom_full_id}) {
        $genome_coords{$hom_full_id} = $hom_coords;
        if($hom_coords =~ m/^\S+?:(\d+)-(\d+)\([+-]\)/) {
          $gene_length{$hom_full_id} = 1+$2-$1;
        } else { 
          die "# ERROR: cannot parse $hom_coords\n";
        }
      }

      push(@pairs, "$line\n");

    } elsif($homology_type eq 'segment_collinear') {

      next if(!$cluster_gene_id{$gene_id} && !$cluster_gene_id{$hom_gene_id});
      next if($cluster_gene_id{$gene_id} && $ref_taxon->{$gene_id} ne $species);
      next if($cluster_gene_id{$hom_gene_id} && $ref_taxon->{$hom_gene_id} ne $hom_species);

      #$cluster_gene_id{$gene_id}++; #segment-gene pairs not considered for stats

      if($segment eq 'segment') {
        push(@segments,$gene_id)
      } elsif($hom_segment eq 'segment') {
        push(@segments,$hom_gene_id)
      } 

      push(@pairs, "$line\n");
    }
  }
}

untie(%TSVdb);

# 5) print raw collinear TSV evidence
if(!%seen) { #scalar(keys(%seen)) != scalar(keys(%cluster_gene_id))) {
  die "# ERROR: cannot find collinear evidence for cluster $INP_clusterfile\n";

} elsif(!$INP_noraw) { 
  print "\n#gene_stable_id\tprotein_stable_id\tspecies\toverlap\thomology_type\t" .
    "homology_gene_stable_id\thomology_protein_stable_id\thomology_species\t".
    "overlap\tdn\tds\tgoc_score\twga_coverage\tis_high_confidence\tcoordinates\n";
  print @pairs;
}

# 6) print gene-level summary stats
my (%scores, %taxon_genes);

print "\n# gene-level stats\n";
print "#len\tpairs\toverlap\tgene_name\tspecies\n";
foreach $full_id (sort {$seen{$b} <=> $seen{$a}} (keys(%seen))){

  $gene_id = $fullid2id{$full_id};

  printf("%d\t%d\t%1.0f\t%s\t%s\n",
    $gene_length{$full_id},
    $seen{$full_id},
    $overlap{$full_id},
    $gene_id,
    $ref_taxon->{$gene_id}
  );

  push(@{ $scores{'pairs'} }, $seen{$full_id});
  push(@{ $scores{'overlap'} }, $overlap{$full_id});
  push(@{ $scores{'length'} }, $gene_length{$full_id});

  # group genes from same taxon/species
  push(@{ $taxon_genes{ $ref_taxon->{$gene_id} }}, $full_id);
}

printf("%d\t%d\t%1.0f\tmedian\tvalues\n",
  calc_median($scores{'length'}),
  calc_median($scores{'pairs'}),
  calc_median($scores{'overlap'}),
);

if(!$INP_fix) {
  exit(0);
} else {
  print "\n";
}

# 7) suggest fixes for poor gene models based on pan-gene consensus
#    (based on chr coords of 1st isoform of each gene)
my $non_outlier_pairs = 0;
my ($chr,$start,$end,$strand,$isof,$ref_lifted_model,$GFF);
my (@long_models, @split_models, @non_outliers, %seen_outlier_taxon);

# 7.1) get outlier cutfoff values 
my ($median_pairs, $cutoff_low_pairs, $cutoff_high_pairs) = 
  get_outlier_cutoffs( $scores{'pairs'} , $INP_verbose );
my ($median_len, $cutoff_low_len, $cutoff_high_len) = 
  get_outlier_cutoffs( $scores{'length'} , $INP_verbose );

# 7.2) identify outlier gene models
my %split_seen;
foreach $full_id (sort {$seen{$b} <=> $seen{$a}} (keys(%seen))){

  $gene_id = $fullid2id{$full_id};

  $genome_file = "$INP_dir/../_$ref_taxon->{$gene_id}.fna";
  if(!-e $genome_file) {
    die "# ERROR: cannot find genome file $genome_file\n";
  }
 
  # long models with too many collinear pairs
  if($seen{$full_id} > $cutoff_high_pairs &&
    ($gene_length{$full_id} > $cutoff_high_len ||
     $gene_length{$full_id} > $median_len)) {

      push(@long_models, $full_id);
      print "# long $gene_id\n" if($INP_verbose); 
  
  } elsif( # short models from same species with few collinear pairs
    $seen{$full_id} < $cutoff_low_pairs &&
    ($gene_length{$full_id} < $cutoff_low_len ||
       $gene_length{$full_id} < $median_len) &&
      scalar(@{ $taxon_genes{ $ref_taxon->{$gene_id} }}) > 1 &&
      !defined($split_seen{$ref_taxon->{$gene_id}})) {

        $split_seen{$ref_taxon->{$gene_id}} = 1;
        foreach my $id (@{ $taxon_genes{ $ref_taxon->{$gene_id} }}) {
          if($gene_length{$id} < $median_len) {
            push(@split_models, $id);
            print "# split $id\n" if($INP_verbose);
          }
        }

  } else { # non-outlier/consensus models

    # record taxa with 2+ non-outlier genes, 
    # these might be used to fix long genes (if there are enough) 
    if(!$seen_outlier_taxon{ $ref_taxon->{$gene_id} }) {
      $seen_outlier_taxon{ $ref_taxon->{$gene_id} } = 1;
      if( scalar(@{$taxon_genes{ $ref_taxon->{$gene_id} }}) > 1) {
        $non_outlier_pairs++;
      }
    }

    push(@non_outliers, $full_id);
  }
}


# 7.3) suggest model fixes in GFF format, order of priority: long > split > missing

if(!@non_outliers) {
  die "# ERROR: need non-outliers/consensus gene models to fix cluster, exit\n";
}

# open output GFF files if requested
my ($gfffh, %outfhandles);
if($INP_outdir) {
  foreach $species (keys(%taxon_genes)) {
    if($INP_appendGFF) {
      open(my $fh,'>>',"$INP_outdir/$species.patch.gff");
      $outfhandles{$species} = $fh;
    } else {
      open(my $fh,'>',"$INP_outdir/$species.patch.gff");
      $outfhandles{$species} = $fh;
    }
  }
}

if($INP_verbose) {
  printf("# long model candidates %d\n", scalar(@long_models));
  printf("# split model candidates %d\n", scalar(@split_models));
  printf("# non-outliers %d pairs %d\n",
    scalar(@non_outliers),$non_outlier_pairs);  
}

if(@long_models &&
  $non_outlier_pairs/scalar(@non_outliers) >= $MINPAIRPECNONOUTLIERS) {

  # hypothesis: a long model actually merges two single genes by mistake
  # proposed fix: liftover individual consensus models against genomic segment containing long gene, 
  # expect 2+ hits from same species on same strand
  foreach $full_id (@long_models) {

    my %lifted;

    $gene_id = $fullid2id{$full_id};

    # cut genomic segment harboring this (long) gene
    $genome_file = "$INP_dir/../_$ref_taxon->{$gene_id}.fna";
    $segment = cut_genomic_segment_bedtools(
      $genome_coords{$full_id},$genome_file,$BEDTOOLSBIN ); 

    # lift-over consensus models
    foreach $hom_full_id (@non_outliers) {
   
      $hom_gene_id = $fullid2id{$hom_full_id};

      # extract all isoforms/transcripts of this gene, 
      # and select the one with most matches in cDNA alignment
      my $best_isof_model;
      my @isofs = extract_isoforms_FASTA($ref_fasta->{$hom_gene_id}); 
      
      foreach $isof (@isofs) { 
        $cDNA = $isof; # leave gene name as only FASTA header of cDNA
        $cDNA =~ s/^>.*\n/>$hom_gene_id\n/;  

        $ref_lifted_model = liftover_gmap( "$gene_id,", 
          $segment, $cDNA, $GMAPBIN, $MINLIFTIDENTITY, $INP_verbose ); 

        if($ref_lifted_model->{'matches'} && (!$best_isof_model || 
            $ref_lifted_model->{'matches'} > $best_isof_model->{'matches'})) {
          $best_isof_model = $ref_lifted_model;
        } 
      } 

      if($best_isof_model->{'matches'}) { 

        $lifted{ $ref_taxon->{$hom_gene_id} }{ 'total' } ++;
        $lifted{ $ref_taxon->{$hom_gene_id} }{ 'matches' } += 
          $best_isof_model->{'matches'};
        $lifted{ $ref_taxon->{$hom_gene_id} }{ 'mismatches' } += 
          $best_isof_model->{'mismatches'};
        $lifted{ $ref_taxon->{$hom_gene_id} }{ 'indels' } += 
          $best_isof_model->{'indels'};

        push(@{ $lifted{ $ref_taxon->{$hom_gene_id} }{ 'GFF' } },
          $best_isof_model->{'GFF'} );
      }
    }

    # choose best (pair of short) models to replace long one
    foreach $species (sort { 
        $lifted{$b}{'total'} <=> $lifted{$a}{'total'} ||
        $lifted{$b}{'matches'} <=> $lifted{$a}{'matches'} ||
        $lifted{$a}{'mismatches'} <=> $lifted{$b}{'mismatches'} ||
        $lifted{$a}{'indels'} <=> $lifted{$b}{'indels'}
      } keys(%lifted)) {

      # at least from same species needed
      next if($lifted{ $species }{ 'total' } < 2);
 
      # sort GFF features and make sure they don't overlap
      $GFF = sort_check_overlap_genes_GFF( @{ $lifted{ $species }{'GFF'} } );

      next if(!$GFF); # mapped models actually overlap

      # actually print GFF and log
      $gfffh = $outfhandles{$ref_taxon->{$gene_id}} || *STDOUT;

      print "# long gene model: corrected $gene_id [$ref_taxon->{$gene_id}]\n";

      printf($gfffh "## replaces %s [%s] source=%s matches=%d mismatches=%d indels=%d\n",
        $gene_id,
        $ref_taxon->{$gene_id},
        $species,
        $lifted{$species}{'matches'},
        $lifted{$species}{'mismatches'},
        $lifted{$species}{'indels'});

      print $gfffh "$GFF\n";

      last; # take only best
    }
  }

} elsif(@split_models) {

  # hypothesis: the real gene is long and was split in 2+ partial genes
  # proposed fix: liftover consensus (longer) models

  foreach $species (keys(%taxon_genes)) {

    next if( scalar(@{$taxon_genes{$species}}) < 2 );

    # cut genomic segment harboring these (split) genes, 
    # make sure they're from the same chromosome/contig
    my %segment_data;
    $genome_file = "$INP_dir/../_$species.fna"; 

    foreach $full_id (@split_models) {

      $gene_id = $fullid2id{$full_id}; 

      next if($ref_taxon->{$gene_id} ne $species);

      if($genome_coords{$full_id} =~ m/^(\S+?):(\d+)-(\d+)\(([+-])\)/) {
        ($chr,$start,$end,$strand) = ($1, $2, $3, $4)
      } else {
        die "# ERROR: cannot parse genomic coords of $gene_id ($genome_coords{$gene_id})\n"
      }  

      if(!defined($segment_data{'chr'})) {
        $segment_data{'chr'}    = $chr;
        $segment_data{'start'}  = $start;
        $segment_data{'end'}    = $end;
        $segment_data{'strand'} = $strand;
        $segment_data{'models'} = 1;
        $segment_data{'genes'}  = "$gene_id";

      } else {

        next if($chr ne $segment_data{'chr'} || $strand ne $segment_data{'strand'});

        if($start < $segment_data{'start'}) {
          $segment_data{'start'} = $start;
          $segment_data{'models'}++;
          $segment_data{'genes'} .= ",$gene_id";

        } elsif($end > $segment_data{'end'}) {
          $segment_data{'end'} = $end;
          $segment_data{'models'}++;
          $segment_data{'genes'} .= ",$gene_id";
        }
      }
    }
   
    next if(!$segment_data{'models'} || $segment_data{'models'} < 2);

    $segment = cut_genomic_segment_bedtools(
      "$segment_data{'chr'}:$segment_data{'start'}-$segment_data{'end'}($segment_data{'strand'})",
      $genome_file,$BEDTOOLSBIN ); 

    # lift-over consensus models
    my %lifted;
    foreach $hom_full_id (@non_outliers) {

      $hom_gene_id = $fullid2id{$hom_full_id};

      # extract all isoforms/transcripts of this gene,
      # and select the one with most matches in cDNA alignment
      my $best_isof_model;
      my @isofs = extract_isoforms_FASTA($ref_fasta->{$hom_gene_id});

      foreach $isof (@isofs) {
        $cDNA = $isof; # leave gene name as only FASTA header of cDNA
        $cDNA =~ s/^>.*\n/>$hom_gene_id\n/;

        $ref_lifted_model = liftover_gmap( $segment_data{'genes'}, 
          $segment, $cDNA, $GMAPBIN, $MINLIFTIDENTITY, $INP_verbose );

        if($ref_lifted_model->{'matches'} && (!$best_isof_model ||
            $ref_lifted_model->{'matches'} > $best_isof_model->{'matches'})) {
          $best_isof_model = $ref_lifted_model;
        }
      }

      if($best_isof_model->{'matches'}) {
        $lifted{ $ref_taxon->{$hom_gene_id} }{ 'total' } ++;
        $lifted{ $ref_taxon->{$hom_gene_id} }{ 'matches' } +=
          $best_isof_model->{'matches'};
        $lifted{ $ref_taxon->{$hom_gene_id} }{ 'mismatches' } +=
          $best_isof_model->{'mismatches'};
        $lifted{ $ref_taxon->{$hom_gene_id} }{ 'indels' } +=
          $best_isof_model->{'indels'};

        push(@{ $lifted{ $ref_taxon->{$hom_gene_id} }{ 'GFF' } },
          $best_isof_model->{'GFF'} );
      }
    }

    # choose best (long) model to replace split ones
    foreach $hom_species (sort {
        $lifted{$b}{'matches'} <=> $lifted{$a}{'matches'} ||
        $lifted{$a}{'mismatches'} <=> $lifted{$b}{'mismatches'} ||
        $lifted{$a}{'indels'} <=> $lifted{$b}{'indels'}
      } keys(%lifted)) {

      # skip genes with multiple mappings
      next if(!$lifted{ $species }{ 'total' } ||
        $lifted{ $species }{ 'total' } > 1);

      $GFF = $lifted{ $hom_species }{'GFF'}->[0];

      # actually print GFF and log
      $gfffh = $outfhandles{$ref_taxon->{$gene_id}} || *STDOUT;

      print "# split gene model: corrected $segment_data{'genes'} [$species]\n";

      printf($gfffh "## replaces %s [%s] source=%s matches=%d mismatches=%d indels=%d\n",
        $segment_data{'genes'},
        $species,
        $hom_species,
        $lifted{$hom_species}{'matches'},
        $lifted{$hom_species}{'mismatches'},
        $lifted{$hom_species}{'indels'});

      print $gfffh "$GFF\n";

      last; # take only best
    }

    # TODO: are there intervening TEs?
  }

} elsif(@segments) {

  # hypothesis: model exists but failed to be annotated
  # proposed fix: liftover consensus models over matching genomic segment

  foreach my $segment_id (@$ref_geneid_seg){

    # skip if genes were parsed for this species/taxon
    next if($taxon_genes{ $ref_taxon_seg->{$segment_id} });

    # obtain DNA sequence of segment in plus strand
    $segment = $ref_fasta_seg->{$segment_id};
    if($ref_coords_seg->{$segment_id}[3] eq '-'){
        $segment = revcomp_fasta($segment);
    } 

    # lift-over consensus models
    my %lifted;
    foreach $hom_full_id (@non_outliers) {

      $hom_gene_id = $fullid2id{$hom_full_id};

      # extract all isoforms/transcripts of this gene,
      # and select the one with most matches in cDNA alignment
      my $best_isof_model;
      my @isofs = extract_isoforms_FASTA($ref_fasta->{$hom_gene_id});

      foreach $isof (@isofs) {
        $cDNA = $isof; # leave gene name as only FASTA header of cDNA
        $cDNA =~ s/^>.*\n/>$hom_gene_id\n/;

        $ref_lifted_model = liftover_gmap( 'missing', 
          $segment, $cDNA, $GMAPBIN, $MINLIFTIDENTITY, $INP_verbose );

        if($ref_lifted_model->{'matches'} && (!$best_isof_model ||
            $ref_lifted_model->{'matches'} > $best_isof_model->{'matches'})) {
          $best_isof_model = $ref_lifted_model;
        }
      }

      if($best_isof_model->{'matches'}) {
        $lifted{ $ref_taxon->{$hom_gene_id} }{ 'total' } ++;
        $lifted{ $ref_taxon->{$hom_gene_id} }{ 'matches' } +=
          $best_isof_model->{'matches'};
        $lifted{ $ref_taxon->{$hom_gene_id} }{ 'mismatches' } +=
          $best_isof_model->{'mismatches'};
        $lifted{ $ref_taxon->{$hom_gene_id} }{ 'indels' } +=
          $best_isof_model->{'indels'};

        push(@{ $lifted{ $ref_taxon->{$hom_gene_id} }{ 'GFF' } },
          $best_isof_model->{'GFF'} );
      }
    }

    # choose best model to replace missing one
    foreach $hom_species (sort {
        $lifted{$b}{'matches'} <=> $lifted{$a}{'matches'} ||
        $lifted{$a}{'mismatches'} <=> $lifted{$b}{'mismatches'} ||
        $lifted{$a}{'indels'} <=> $lifted{$b}{'indels'}
      } keys(%lifted)) {

      # skip genes with multiple mappings
      next if(!$lifted{ $species }{ 'total' } || 
        $lifted{ $species }{ 'total' } > 1);

      $GFF = $lifted{ $hom_species }{'GFF'}->[0];

      # actually print GFF and log
      $gfffh = $outfhandles{$ref_taxon->{$gene_id}} || *STDOUT;

      print "# missing gene model: corrected $segment_id [$ref_taxon_seg->{$segment_id}]\n";

      printf($gfffh "## adds %s [%s] source=%s matches=%d mismatches=%d indels=%d\n",
        $segment_id,
        $ref_taxon_seg->{$segment_id},
        $hom_species,
        $lifted{$hom_species}{'matches'},
        $lifted{$hom_species}{'mismatches'},
        $lifted{$hom_species}{'indels'});

      print $gfffh "$GFF\n";

      last; # take only best
    }


  }
}


################################################################################

# Cuts FASTA sequence for a passed genomic interval.
# Takes 3+1 params:
# i)   string with 1-based genomic coords (chr:start-end(strand)
# ii)  FASTA file with reference sequence
# iii) path to bedtools
# iv)  force strandness (optional boolean)
# Returns string with one sequence in FASTA format
sub cut_genomic_segment_bedtools {
  my ($genome_coords, $fastafile, $bedtools_path, $force_strand) = @_;
 
  # parse coord string 
  my ($chr, $start, $end, $strand);
  if($genome_coords =~ m/^(\S+?):(\d+)-(\d+)\(([+-])\)/) {
    ($chr,$start,$end,$strand) = ($1, $2, $3, $4);
    $start -= 1; # make it 0-based, as in BED format

  } else {
    die "# ERROR(cut_genomic_segment_bedtools): cannot parse genomic coords $genome_coords\n"
  }

  my $cmd = "echo '$chr\t$start\t$end' | ".
    "$bedtools_path getfasta -fi $fastafile -bed -";
  if($force_strand) {
    $cmd = "echo '$chr\t$start\t$end\tfrag\t0\t$strand' | ".
    "$bedtools_path getfasta -fi $fastafile -bed - -s";
  }

  my $fasta_segment = '';

  open(BEDTOOLS, "$cmd |") ||
    die " ERROR(cut_genomic_segment_bedtools): cannot run $cmd\n";
  while(<BEDTOOLS>) {
    $fasta_segment .= $_;
  }
  close(BEDTOOLS);

  if(!$fasta_segment) {
    die "# ERROR(cut_genomic_segment_bedtools): $cmd produced no sequence\n";
  }

  return $fasta_segment;
}


# Lifts over a gene model upon gmap mapping of a cDNA.
# Takes 5+1 params:
# i)   original gene_ids (comma separated strings)
# ii)  FASTA string of genomic sequence (target)
# iii) FASTA string of cDNA sequence (query)
# iv)  path to gmap
# v)   min identity percentage (float)
# vi)  verbose boolean, optional
# Returns ref to hash with model features, including 'GFF' and
# several scores such as 'matches', 'mismatches' or 'indels'
# Note: gene_id of lifted models if that of source cDNA
# Note: original gene_id annotated as old_locus_tag in GFF's attributes column
sub liftover_gmap {
  my ($old_gene_ids, $target_fna, $query_cdna, 
      $gmap_path, $min_identity, $verbose) = @_;

  my ($chr,$start,$end,$offset,$cmd);  
  my ($seqname,$aa, %aaseq);
  my ($identity, $match, $mismatch, $indel, $gffOK) = (0, 0, 0, 0, 0);
  my (%lifted_model);

  # check coordinates of target DNA in source, expects plus strand 
  # so that strand in the produced GFF can be applied to raw gDNA
  if($target_fna =~ m/^>(\S+?):(\d+)-(\d+)/) {
    ($chr,$start,$end) = ($1, $2, $3);
    $offset = $start;
  } else {
    die "# ERROR(liftover_gmap): cannot parse target coords ($target_fna)\n";
  }

  # 1st run: get alignment summary to parse scores
  # Note: discard stderr as it gets in teh way while parsing stdout
  $cmd = "echo '$target_fna$query_cdna' | $gmap_path -A"; 
  open(GMAP, "$cmd 2>/dev/null |") ||
    die "# ERROR(liftover_gmap): cannot run $cmd\n"; 
  while(<GMAP>) { 
    if(/^\s+Percent identity: (\S+) \((\d+) matches, (\d+) mismatches, (\d+) indels/){ 
      ($identity, $match, $mismatch, $indel) = ($1, $2, $3, $4); 
      if($identity < $min_identity) {
        ($match, $mismatch, $indel) = (0, 0, 0);
      }
    } elsif(/^aa\.([cg])\s+\d+\s(.*)/) {
      #aa.g        48  N  N  P  S  S  Q  I  T  Y  G  L  T  I  H  H  A  V 
      ($seqname,$aa) = ($1,$2);
      $aaseq{$seqname} .= $aa;
    }
    print if($verbose);
  }
  close(GMAP);

  # are there premature stop codons in the genomic sequence?
  my ($stoposg, $stoposc) = (-1,-1);
  if($aaseq{'g'} && $aaseq{'g'} =~ m/\*/) { $stoposg = $-[0] }
  if($aaseq{'c'} && $aaseq{'c'} =~ m/\*/) { $stoposc = $-[0] }
  if($stoposg < $stoposc) {
    $match = 0;
    if($verbose){
      print "# WARN(liftover_gmap): premature stop codon ($stoposg < $stoposc)\n";
    }
  }

  if($match == 0) { 
    return \%lifted_model
  } 

  # 2nd run: map query cDNA on genomic target, parse GFF and apply coord offset
  $cmd = "echo '$target_fna$query_cdna' | $gmap_path -f 2";
  
  open(GMAP, "$cmd 2>/dev/null |") || 
    die "# ERROR(liftover_gmap): cannot run $cmd\n";  
  while(<GMAP>) {

    last if(/^###/);

    if($gffOK) {
      my @GFF = split(/\t/);

      # add original chr name
      $GFF[0] = $chr;

      # add source
      $GFF[1] = 'gmap';

      # apply genomic offset to coords of lifted gene model
      $GFF[3] += $offset;
      $GFF[4] += $offset; 

      # add note saying which original gene is replaced by this
      if($GFF[2] eq 'gene') {
        $GFF[8] =~ s/\n/;/;
        foreach my $gene_id (split(/,/,$old_gene_ids)) {
          $GFF[8] .= "old_locus_tag=$gene_id;";
        } $GFF[8] .= "\n";
      }

      $lifted_model{'GFF'} .= join("\t",@GFF); # if($GFF[2] eq 'gene');
    }

    if(/^# Generated by GMAP/) {
      $gffOK = 1
    }
  }
  close(GMAP);

  $lifted_model{'identity'} = $identity;
  $lifted_model{'matches'} = $match;
  $lifted_model{'mismatches'} = $mismatch;
  $lifted_model{'indels'} = $indel;

  return \%lifted_model;
}

# Takes array with GFF string blocks, one per gene, sorts them,
# makes sure they don't overlap and return sorted GFF string.
# If two block overlap and empty string is returned
sub sort_check_overlap_genes_GFF {
  my @unsortedGFF = @_;

  my ($GFFblock, $start, $end, $key);
  my (%block, @sorted_blocks);
  my $GFFstring = '';

  # compute block start and end coords, assumes 1st line is gene
  foreach my $b (0 .. $#unsortedGFF) {
    $GFFblock = $unsortedGFF[$b];
    ($start,$end) = (split(/\t/,$GFFblock,6))[3,4];
    $key = "$start,$end";
    $block{$key} = [$start,$end,$b];  
  } 

  # sort blocks
  @sorted_blocks = sort {$block{$a}->[0]<=>$block{$b}->[0]} keys(%block);

  # make sure one block does not overlap the next one
  for $key (0 .. $#sorted_blocks-1) {
    if($block{$sorted_blocks[$key]}->[1] >= $block{$sorted_blocks[$key+1]}->[0]) {
      return $GFFstring
    }
  }

  # concat sorted, non-overlapping blocks 
  for $key (0 .. $#sorted_blocks) {
    $GFFstring .= $unsortedGFF[ $block{$sorted_blocks[$key]}->[2] ];
  } 

  return $GFFstring
}

# Takes a DNA sequence FASTA string and returns its reverse complemente string
sub revcomp_fasta {
  my ($rawseq) = @_;

  # work out header
  my ($header,$seq) = (split(/\n/,$rawseq,2));
  $header =~ s/\(-\)/(+)/;

  # reverse sequence 
  my $revcomp = $header;

  foreach my $frag (split(/\n/,reverse($seq))) {
    $frag =~ tr/ATGCNXatgcn/TACGNtacgnx/;
    $revcomp .= "$frag\n";
  }

  return $revcomp;
}
