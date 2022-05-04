#!/usr/bin/env perl

# This script takes a cDNA cluster produced by get_pangenes.pl and retrieves
# the collinearity data evidence supporting it, inferred from whole genome
# alignments, contained in mergedpairs.tsv.gz (sorted with -k1,1 -k4,4nr)
#
# Optionally it can also suggest fixes to the gene models based on the 
# pangene consensus (-f)

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
use pangeneTools qw( parse_sequence_FASTA_file calc_median get_outlier_cutoffs );

my $MINPAIRPECNONOUTLIERS = 0.50;

$ENV{'LC_ALL'}   = 'POSIX';

my $GZIPBIN = $ENV{'EXE_GZIP'} || 'gzip';
my $BEDTOOLSBIN = $ENV{'EXE_BEDTOOLS'} || 'bedtools';
my $GMAPBIN = $ENV{'EXE_GMAP'} || 'gmap-2021-12-17/bin/bin/gmap';

my ($INP_dir,$INP_clusterfile,$INP_noraw,$INP_fix) = ( '', '', 0 , 0 );
my ($cluster_list_file,$cluster_folder,$gdna_clusterfile);
my ($gene_id, $hom_gene_id, $homology_type, $species, $hom_species);
my ($overlap, $coords, $hom_coords, $full_id, $hom_full_id);
my ($line, $segment, $hom_segment, $dummy, $TSVdata, $cmd);
my (%opts,%TSVdb, @sorted_ids, @pairs, @segments);
my (%seen, %overlap, %cluster_gene_id, %fullid2id, %gene_length);

getopts('hfnd:i:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-d output directory produced by get_pangenes.pl (example: -d /path/data_pangenes/..._algMmap_,\n";
  print "                                                 genomic sequences usually one folder up)\n";
  print "-i cdna cluster as shown in .cluster_list file  (example: -i gene:ONIVA01G52180.cdna.fna)\n";
  print "-n do not print raw evidence                    (optional)\n";
  print "-f fix gene models and produce GFF              (optional)\n\n";
  print "Note: reads the compressed merged TSV file in -d\n"; 
  exit(0);
}

if(defined($opts{'d'})) { 
  $INP_dir = $opts{'d'} 
}
else{ die "# EXIT : need a -d directory\n" }

if(defined($opts{'i'})){  
  $INP_clusterfile = $opts{'i'};
  if($INP_clusterfile !~ /\.cdna\.fna$/) {
    die "# EXIT : need a .cdna.fna cluster filename with parameter -i\n"
  } else {
    $gdna_clusterfile = $INP_clusterfile;
    $gdna_clusterfile =~ s/\.cdna\.fna/.gdna.fna/;
  }
}
else{ die "# EXIT : need parameter -i\n" }

if(defined($opts{'n'})){ 
  $INP_noraw = 1 
}

if(defined($opts{'f'})){
  $INP_fix = 1
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

# 2) parse FASTA headers of input cluster to extract gene names
my ( $ref_geneid, $ref_fasta, $ref_coords, $ref_taxon ) = 
  parse_sequence_FASTA_file( "$INP_dir/$cluster_folder/$INP_clusterfile" , 1);

foreach $gene_id (sort @$ref_geneid) {
  $cluster_gene_id{$gene_id} = 1;
  push(@sorted_ids, $gene_id);
} 

# 2.1) parse segment gDNA file if required
if(-e "$INP_dir/$cluster_folder/$gdna_clusterfile") {

  print "# parsing twin .gdna.fna cluster file $gdna_clusterfile\n\n";

  my ( $ref_geneid_seg, $ref_fasta_seg, $ref_coords_seg, $ref_taxon_seg ) =
    parse_sequence_FASTA_file( "$INP_dir/$cluster_folder/$gdna_clusterfile" , 1);

  printf("\n# cluster %s genes = %d genomic segments = %d\n",
    $INP_clusterfile,
    scalar(keys(%cluster_gene_id)),
    scalar(@$ref_geneid_seg));

} else {

  printf("\n# cluster %s genes = %d\n",
    $INP_clusterfile,
    scalar(keys(%cluster_gene_id)));
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

      next if(!$cluster_gene_id{$gene_id} || !$cluster_gene_id{$hom_gene_id});
      next if($ref_taxon->{$gene_id} ne $species || $ref_taxon->{$hom_gene_id} ne $hom_species);

      # concat species & gene_id in case there are repeated gene ids
      $full_id = $species.$gene_id;
      $hom_full_id = $hom_species.$hom_gene_id;
      $fullid2id{$full_id} = $gene_id;
      $fullid2id{$hom_full_id} = $hom_gene_id;

      $cluster_gene_id{$gene_id}++;
      $seen{$full_id}++;
      $seen{$hom_full_id}++;

      # overlap only considered for gene pairs (not segments)
      $overlap{$full_id} += $overlap;
      $overlap{$hom_full_id} += $overlap; 

      # compute gene length
      ($coords, $hom_coords) = split(/;/,$coords); 

      if(!$gene_length{$full_id}) {
        if($coords =~ m/^\S+?:(\d+)-(\d+)\([+-]\)/) {
          $gene_length{$full_id} = 1+$2-$1;
        } else {
          die "# ERROR: cannot parse $coords $_\n";
        }
      }

      if(!$gene_length{$hom_full_id}) {
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

      # segment-gene pairs not considered for stats 
      #$cluster_gene_id{$gene_id}++;

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

# 6) print summary stats
my (%scores, %taxon_obs);
print "\n#length\tpairs\tgene_overlap\tgene\tspecies\n";
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

  $taxon_obs{ $ref_taxon->{$gene_id} }++;
}

printf("%d\t%d\t%1.0f\tmedian\tvalues\n",
  calc_median($scores{'length'}),
  calc_median($scores{'pairs'}),
  calc_median($scores{'overlap'}),
);

if(!$INP_fix) {
  exit(0);
}

# 7) suggest fixes for poor gene models based on pan-gene consensus

my $non_outlier_pairs = 0;
my %seen_outlier_taxon;
my (@long_models, @split_models, @non_outliers);

# 7.1) get outlier cutfoff values
my ($cutoff_low_pairs, $cutoff_high_pairs) = get_outlier_cutoffs( $scores{'pairs'} );
my ($cutoff_low_len, $cutoff_high_len) = get_outlier_cutoffs( $scores{'length'} );

# 7.2) identify outlier gene models
foreach $full_id (sort {$seen{$b} <=> $seen{$a}} (keys(%seen))){

  $gene_id = $fullid2id{$full_id};

  # long models with too many collinear pairs
  if($seen{$full_id} > $cutoff_high_pairs && 
    $gene_length{$full_id} > $cutoff_high_len) {
  
      push(@long_models, $gene_id)
  
  } elsif( # short models from same species with few collinear pairs
    $seen{$full_id} < $cutoff_low_pairs &&
    $gene_length{$full_id} < $cutoff_low_len &&
    $taxon_obs{ $ref_taxon->{$gene_id} } > 1) {

    push(@split_models, $gene_id)

  } else {

    if(!$seen_outlier_taxon{ $ref_taxon->{$gene_id} }) {
      $seen_outlier_taxon{ $ref_taxon->{$gene_id} } = 1;
      if( $taxon_obs{ $ref_taxon->{$gene_id} } > 1) {
        $non_outlier_pairs++;
      }
    }

    push(@non_outliers, $gene_id);
  }
}


# 7.3) suggest model fixes, in order of priority: long > split > missing

if(!@non_outliers) {
  die "# ERROR: need non-outliers/consensus gene models to fix cluster, exit\n";
}

if(@long_models &&
  $non_outlier_pairs/scalar(@non_outliers) >= $MINPAIRPECNONOUTLIERS) {

  # hypothesis: a long model actually merges two single genes by mistake
  # proposed fix: liftover individual consensus models, 2+ hits on same strand
  foreach $gene_id (@long_models) {

    # check whether genome sequence file is available
    my $genome_file = "$INP_dir/../_$ref_taxon->{$gene_id}.fna";
    if(!-e $genome_file) {
      die "# ERROR: cannot find genome sequence file $genome_file\n";  
    }

    my $segment = cut_genomic_segment_bedtools( 
      $ref_coords->{$gene_id}[0],
      $ref_coords->{$gene_id}[1],
      $ref_coords->{$gene_id}[2],
      $ref_coords->{$gene_id}[3],
      $genome_file,
      $BEDTOOLSBIN 
    );

    foreach $gene_id (@non_outliers) {

      my $ref_lifted_model = liftover_gmap( 
        $segment,   
        $ref_fasta->{$gene_id},
        $GMAPBIN );

      exit;
    }

  }

} elsif(@split_models) {

  # hypothesis: the real gene is long and was split in two partial genes
  # proposed fix: i) liftover consensus (longer) models, ii) are there intervening TEs?

  # print "S $gene_id $ref_taxon->{$gene_id}\n";

} elsif(@segments) {

  # hypothesis: model exists but failed to be annotated
  # proposed fix: liftover consensus models

}

# Cuts FASTA sequence for a passed genomic interval.
# Takes 6 params:
# i)   chr name
# ii)  start coord
# iii) end coord
# iv)  strand
# v)  FASTA file with reference sequence
# vi)   path to bedtools
# Returns string with one sequence in FASTA format
sub cut_genomic_segment_bedtools {
  my ($chr, $start, $end, $strand, $fastafile, $bedtools_path) = @_;
  
  my $cmd = "echo '$chr\t$start\t$end\tfrag\t0\t$strand' | ".
    "$bedtools_path getfasta -fi $fastafile -bed - -s";

  my $fasta_segment = '';

  open(BEDTOOLS, "$cmd |") ||
    die " ERROR(cut_genomic_segment_bedtools): cannot run $cmd\n";
  while(<BEDTOOLS>) {
    $fasta_segment .= $_;
  }
  close(BEDTOOLS);

  if(!$fasta_segment) {
    die " ERROR(cut_genomic_segment_bedtools): $cmd produced no sequence\n";
  }

  return $fasta_segment;
}

# Lifts over a gene model upon gmap mapping of a cDNA.
# Takes 3 params:
# i)   FASTA string of genomic sequence (target)
# ii)  FASTA string of cDNA sequence (query)
# iii) path to gmap
# Returns ref to hash with model features, including 'GFF' and
# several scores such as 'matches', 'mismatches' or 'indels'
sub liftover_gmap {
  my ($target_fna, $query_cdna, $gmap_path) = @_;
  my ($chr,$start,$end,$strand,$offset,$cmd);  
  my ($match, $mismatch, $indel, $gffOK) = (0, 0, 0, 0);
  my (%lifted_model);
 
  # check coordinates of target DNA in source 
  if($target_fna =~ m/^>(\S+?):(\d+)-(\d+)\([+-]\)/) {
    ($chr,$start,$end,$strand) = ($1, $2, $3, $4);
    $offset = $start - 1;
  } else {
    die " ERROR(liftover_gmap): cannot parse target coords ($target_fna)\n";
  }

  # 1st run: get alignment summary to parse scores
  $cmd = "echo '$target_fna$query_cdna' | $gmap_path -t 1 -2 -S -z sense_force";
  open(GMAP, "$cmd 2>&1 |") ||
    die " ERROR(liftover_gmap): cannot run $cmd\n";
  while(<GMAP>) { 
    if(/^\s+Percent identity: \S+ \((\d+) matches, (\d+) mismatches, (\d+) indels/){ 
      ($match, $mismatch, $indel) = ($1, $2, $3); 
    }
  }
  close(GMAP);

  # 2nd run: map query cDNA on genomic target, parse GFF and apply coord offset
  $cmd = "echo '$target_fna$query_cdna' | $gmap_path -t 1 -2 -f 2 -z sense_force";
  
  open(GMAP, "$cmd 2>&1 |") || 
    die " ERROR(liftover_gmap): cannot run $cmd\n";
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

      $lifted_model{'GFF'} = join("\t",@GFF);
    }

    if(/^# Generated by GMAP/) {
      $gffOK = 1
    }
  }
  close(GMAP);

  $lifted_model{'matches'} = $match;
  $lifted_model{'mismatches'} = $mismatch;
  $lifted_model{'indels'} = $indel;

  return \%lifted_model;
}
