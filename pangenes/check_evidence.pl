#!/usr/bin/env perl

# This script takes a cDNA/CDS cluster produced by get_pangenes.pl and retrieves
# the collinearity data evidence supporting it, inferred from whole genome
# alignments, often including cases with low coverage and thus not in cluster. 
# Evidence is contained in file mergedpairs.tsv.gz (presorted with -k1,1 -k4,4nr).
#
# It can output a representative cluster sequence with mode length (-s).
# When CDS clusters are analyzed a check for internal stop codons is carried out.
#
# It can produce code to plot a genomic context sketch of a cluster (-P).
#
# Optionally it can also liftover/suggest fixes to the gene models based on the 
# pangene consensus (-f), this requires gmap (make install_pangenes)

# Copyright [2022-23]
# EMBL-European Bioinformatics Institute & Estacion Experimental Aula Dei-CSIC

$|=1;

use strict;
use warnings;
use FindBin '$Bin';
use lib "$Bin/../lib";
use lib "$Bin/lib";
use Getopt::Std;
use DB_File;
use Compress::Zlib qw(compress uncompress);
use pangeneTools qw( check_installed_features feature_is_installed 
                     parse_sequence_FASTA_file extract_isoforms_FASTA
                     calc_median calc_mode get_outlier_cutoffs );

my @standard_stop_codons = qw( TAG TAA TGA );

my $MINPAIRPECNONOUTLIERS = 0.25;
my $MINLIFTIDENTITY = 90.0;
my $MINFIXOVERLAP   = 0.75; # min overlap of mapped genes to correct long/split models
my $MAXSEGMENTSIZE  = 100_000;
my $GMAPARAMS       = '-t 1 -2 -z sense_force -n 1 -F ';

# plotting global settings
my $PLOTDPI         = 300;
my $PLOTFLANKING    = 100; 
my $PLOTGENEWIDTH   = 500; 
my $PLOTNEIGH       = 1; # number of neighbors each side for -P, 
                         # must contain reference genome/annotation

my @FEATURES2CHECK = (
  'EXE_BEDTOOLS', 'EXE_GMAP', 'EXE_GZIP'
);

$ENV{'LC_ALL'}   = 'POSIX';

my $GZIPBIN = $ENV{'EXE_GZIP'} || 'gzip';
my $BEDTOOLSBIN = $ENV{'EXE_BEDTOOLS'} || 'bedtools';
my $GMAPBIN = $ENV{'EXE_GMAP'} || 'gmap';
$GMAPBIN .= " $GMAPARAMS ";

my ($INP_dir, $INP_clusterfile, $INP_noraw, $INP_fix) = ('','',0,0);
my ($INP_verbose, $INP_appendGFF, $INP_modeseq, $INP_outdir) = (0,0,'','');
my ($INP_partial, $INP_plot_code, $INP_mode_stats) = (0, 0, 0);
my ($isCDS, $seq, $gfffh, $CDSok, $outputGFF, %badCDS, %outlier_isoform) = ( 0 );
my ($cluster_list_file,$cluster_folder,$gdna_clusterfile, $genome_file);
my ($gene_id, $hom_gene_id, $homology_type, $species, $hom_species);
my ($isof_id, $overlap, $coords, $hom_coords, $full_id, $hom_full_id);
my ($line, $segment, $hom_segment, $dummy, $TSVdata, $cmd, $cDNA);
my ($chr, $start, $end, $strand, $len);
my (%isof_len, %isof_seq, %isof_header, %taxa, %outfhandles, @len);
my (%opts,%TSVdb, @sorted_ids, @pairs, @segments, @ref_names);
my (%seen, %overlap, %cluster_gene_id, %fullid2id, %gene_length);
my (%genome_coords, %scores, %taxon_genes, %taxon_segments);
my $TAB_matrix_file = 'pangene_matrix_genes.tr.tab';
my $BED_matrix_file = 'pangene_matrix.tr.bed';

getopts('hPvpacfnmr:s:o:d:i:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-c print credits and checks installation\n";
  print "-d directory produced by get_pangenes.pl        (example: -d /path/data_pangenes/..._algMmap_,\n";
  print "                                                 genomic sequences usually one folder up)\n";
  print "-i cdna/cds .fna file as in .cluster_list file  (example: -i gene:ONIVA01G52180.cdna.fna)\n";
  print "-s append mode isoform sequence to file         (optional, example: -s isoforms.fna)\n";
  print "-r CSV string with reference taxa for mode      (optional, prefers mode from ref, requires -s,\n";
  print "                                                 example: oryza_sativa1,oryza_sativa2)\n";
  print "-m print all sequences that match mode          (optional, overriden with -s -r)\n"; 
  print "-n do not print raw evidence                    (optional)\n";
  print "-f fix gene models and produce GFF              (optional, GFF printed to stdout by default)\n";
  print "-p allow partial lifted-over CDSs               (optional, by default only multiples of 3)\n";  
  print "-o folder to write GFF output                   (optional, requires -f, 1 file/species)\n";
  print "-a append GFF output                            (optional, requires -f -o)\n"; 
  print "-P make python code to plot cluster context     (optional, requires _split_ results dir, overriden with -f)\n";
  print "-v verbose, prints intermediate results         (optional, useful to see GMAP alignments)\n";
  print "\nNote: reads the compressed merged TSV file in -d\n"; 
  exit(0);
}

if(defined($opts{'c'})) {
  print "\nPrimary citation:\n https://doi.org/10.1186/s13059-023-03071-z\n";
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

    if($INP_clusterfile =~ /\.cds\.fna/) { 
      $isCDS = 1
    }
  }
}
else{ die "# EXIT : need parameter -i\n" }

if(defined($opts{'s'})){
  $INP_modeseq = $opts{'s'};
 
  if(defined($opts{'r'})) {
    @ref_names = split(/[,;]/,$opts{'r'}); 
    if(!@ref_names) {
      die "# EXIT: cannot parse reference names (-r), " .
        "make sure they have no blanks and are comma-separated\n";
    }
  }
} elsif(defined($opts{'m'})){
  $INP_mode_stats = 1
}

if(defined($opts{'n'})){ 
  $INP_noraw = 1 
}

if(defined($opts{'f'})){
  $INP_fix = 1;

  check_installed_features('EXE_GMAP');
  if(!feature_is_installed('EXE_GMAP')) {
    die "# EXIT : cannot find gmap binary, ".
      "see dependency instructions, can be installed with make install_gmap\n";
  }

  if(defined($opts{'p'})){
    $INP_partial = 1
  }

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
} elsif(defined($opts{'P'})) {

  # check required files are in place 	
  if(!-s "$INP_dir/$BED_matrix_file") {
    die "# WARN: cannot find $INP_dir/$BED_matrix_file; ".
      "please re-run get_pangenes.pl with option -s\n";

  } if(!-s "$INP_dir/$TAB_matrix_file") {
    die "# WARN: cannot find $INP_dir/$TAB_matrix_file; ".
      "please re-run get_pangenes.pl with option -s\n";

  }else {
    $INP_plot_code = 1    
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

my ( $ref_geneid, $ref_fasta, $ref_isof_coords, $ref_taxon ) = 
  parse_sequence_FASTA_file( "$INP_dir/$cluster_folder/$INP_clusterfile" , 1);

print "\n# sequence-level stats\n";

foreach $gene_id (@$ref_geneid) {

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

  # check for internal stop codons
  $CDSok = 1; # default (cDNA)
  if($isCDS) {
    $CDSok = no_premature_stops( $isof_seq{$gene_id}{$isof_id}, 
      \@standard_stop_codons, $isCDS, $INP_verbose); 
  }

  if($CDSok == 1) {
    push(@len, $isof_len{$gene_id}{$isof_id})

  } elsif($CDSok == 3) {
    print "# WARN: $gene_id $isof_id [$ref_taxon->{$gene_id}] ".
      "CDS length not multiple of 3, skip it\n";
    $full_id = $ref_taxon->{$gene_id}.$gene_id;
    $badCDS{$full_id} = 1;

  } else {
    print "# WARN: $gene_id $isof_id [$ref_taxon->{$gene_id}] ".
      "contains internal stop codons, skip it\n";

    $full_id = $ref_taxon->{$gene_id}.$gene_id;
    $badCDS{$full_id} = 1;
  }

  # taxa stats
  $taxa{ $ref_taxon->{$gene_id} }++;
}

my ($median_length, $cutoff_low_length, $cutoff_high_length) =
  get_outlier_cutoffs( \@len , $INP_verbose );

# note: if length distribution is bimodal there will be 2 modes
my @modes_length = calc_mode( \@len );

if(@modes_length) {
  printf("\n# isoform length in cluster: median=%1.0f mode(s): %s\n\n",
    $median_length, join(',',@modes_length));
} else {
  print "\n# isoform length in cluster: median=NA mode(s): NA\n\n";
}

if($INP_modeseq) {
  my ($mode_gene_id, $mode_isof_id);

  foreach $gene_id (@$ref_geneid) {
    foreach $isof_id (keys(%{$isof_len{$gene_id}})) {
 
      next if($badCDS{ $ref_taxon->{$gene_id}.$gene_id } || 
        $isof_len{$gene_id}{$isof_id} != $modes_length[0]);
 
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
    " (append to $INP_modeseq)\n";

} elsif($INP_mode_stats) {
  foreach $gene_id (@$ref_geneid) {
    foreach $isof_id (keys(%{$isof_len{$gene_id}})) {

      next if($badCDS{ $ref_taxon->{$gene_id}.$gene_id } ||
        $isof_len{$gene_id}{$isof_id} != $modes_length[0]);

      print "# mode isoform: $gene_id $isof_id [$ref_taxon->{$gene_id}]\n";
    }
  }
  print "\n";
}

foreach $gene_id (@$ref_geneid) {
  foreach $isof_id (keys(%{$isof_len{$gene_id}})) {
    if($isof_len{$gene_id}{$isof_id} < $cutoff_low_length) {

      print "# short isoform: $isof_id $gene_id [$ref_taxon->{$gene_id}] ".
        "length=$isof_len{$gene_id}{$isof_id}\n";

      # this isoform should not be used for lifting-over
      $outlier_isoform{$ref_taxon->{$gene_id}}{$isof_id} = 1;

    } elsif($isof_len{$gene_id}{$isof_id} > $cutoff_high_length) {

      print "# long isoform: $isof_id $gene_id [$ref_taxon->{$gene_id}] ".
       "length=$isof_len{$gene_id}{$isof_id}\n";

      $outlier_isoform{$ref_taxon->{$gene_id}}{$isof_id} = 1;	   
    }
  }
} 


# 2.1) parse segment gDNA cluster file if available (1-based coordinates),
#      might contain sequences with low coverage in input TSV
my ($ref_geneid_seg,$ref_fasta_seg,$ref_coords_seg,$ref_taxon_seg);

if(-e "$INP_dir/$cluster_folder/$gdna_clusterfile") {

  print "# parsing twin .gdna.fna cluster file $gdna_clusterfile\n\n";

  ( $ref_geneid_seg, $ref_fasta_seg, $ref_coords_seg, $ref_taxon_seg ) =
    parse_sequence_FASTA_file( "$INP_dir/$cluster_folder/$gdna_clusterfile" , 1);

  foreach $gene_id (@$ref_geneid_seg) {
    $taxon_segments{ $ref_taxon_seg->{$gene_id} } = 1;
  }

  printf("\n# cluster %s genes = %d genomic segments = %d (%d taxa)\n",
    $INP_clusterfile,
    scalar(keys(%cluster_gene_id)),
    scalar(@$ref_geneid_seg),
    scalar(keys(%taxa)));

} else {

  printf("\n# cluster %s genes = %d (%d taxa)\n",
    $INP_clusterfile,
    scalar(keys(%cluster_gene_id)),
    scalar(keys(%taxa)));
}


# 3) parse compressed merged TSV file and feed Berkeley DB, 
#    only 1st time this script is called
my $mergedTSVgz  = "$INP_dir/mergedpairs.tsv.gz";
my $TSVdb_file   = "$INP_dir/mergedpairs.tsv.bdb";

# only first time
if(!-s $TSVdb_file) {

  print "\n# creating database (might take long, first time)\n";

  if(!-s $mergedTSVgz) {
    die "# ERROR: cannot find $mergedTSVgz, please check -d argument\n";
  } 

  tie(%TSVdb, 'DB_File', $TSVdb_file, 
    O_RDWR|O_CREAT, 0666, $DB_BTREE) ||
    die "# ERROR: cannot create file $TSVdb_file: $!\n";

  open(TSV, "$GZIPBIN -dc $mergedTSVgz |") ||
    die "# ERROR: cannot uncompress $mergedTSVgz\n";

  my ($block, $prev_gene_id) = ('', '');
  while(<TSV>) {
    if(/^([^\t]+)/) {
      $gene_id = $1;

      if($gene_id ne $prev_gene_id) {

        if($block) {
          $TSVdb{$prev_gene_id} = compress($block);
        }
       
        # start new block 
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

  $TSVdata = uncompress($TSVdb{$gene_id}); #print ">> $TSVdata\n\n";
 
  foreach $line (split(/\n/,$TSVdata)) { 

    #LOC_Os01g01019 LOC_Os01g01019 oryza_sativa_MSU 1065 ortholog_collinear
    #gene:Osir64_01g0000020 gene:Osir64_01g0000020 oryza_sativa_ir64 1065	
    #NULL NULL NULL 100.00 1 Chr1:11217-12435(+);1:19441-20506(+)

    ( $gene_id, $segment, $species, $overlap, $homology_type, 
      $hom_gene_id, $hom_segment, $hom_species, $dummy, 
      $dummy, $dummy, $dummy, $dummy, $dummy, $coords   
    ) = split( /\t/, $line );
 
    if($homology_type eq 'ortholog_collinear') {

      # skip sequences not included in cluster,
      # note these might be collinear with not enough overlap
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

      # gather segment evidence as well (to be printed)
      push(@pairs, "$line\n");
    }
  }
}

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


# 7) produce python code to plot cluster of interest 
#    and up to $PLOTNEIGH neighbors
if($INP_plot_code) {

  print "\n# write code for plotting cluster genomic context\n";

  my ($l, $up, $dw, $cl, $col, $totalup, $totaldw);
  my ($region_size, $gene_size, $color, $taxon);
  my ($total_slots, $slot, $total_genes, $clusterOK);
  my ($clname, $shape, $clslot, $gene_in_slot);
  my (@BED, @plot_clusters_BED, @tabtaxa);
  my (%plot_coords, %plot_blocks, %plot_tracks);
  my $plot_scriptfile = $INP_clusterfile .'.plot.py';
  my $plot_file       = $INP_clusterfile .'.plot.png';
  my $plot_logfile    = $INP_clusterfile .'.plot.log.tsv';
 
  # 7.1) get taxa order
  open(TABMAT,"<","$INP_dir/$TAB_matrix_file") ||
    die "# ERROR: cannot read $INP_dir/$TAB_matrix_file\n";
  while(<TABMAT>) {
    # source:...algMmap_/MorexV3 MorexV3 Morex HOR10350 .. BarkeBaRT2v18
    if(/^source:/) {
      my @data = split(/\t/,$_);
      shift(@data);
      pop(@data);
      push(@tabtaxa, @data);
      last;
    }
  }
  close(TABMAT);

  printf("# taxon order (%d): %s\n\n",
    scalar(@tabtaxa),
    join(',',@tabtaxa)) if($INP_verbose);

  # 7.2) find neighbor clusters
  open(BEDMAT,"<","$INP_dir/$BED_matrix_file") ||
    die "# ERROR: cannot read $INP_dir/$BED_matrix_file\n";
  @BED = <BEDMAT>;
  close(BEDMAT); 

  $clusterOK = 0;
  foreach $l (0 .. $#BED) {
    my @data = split(/\t/,$BED[$l]);

    # find cluster of interest, cluster name matches 1st gene name in $ref_geneid,
    # this might happen twice if cluster split in reference annotation/genome
    if($data[3] eq $ref_geneid->[0]){

      $clusterOK++; # > 1 for split ref genes, these appease in consecutive rows  

      # recall line of input cluster
      if($clusterOK == 1) {
        $cl = $l
      }	

      # get upstream clusters
      if($clusterOK == 1) {
        $up = $l;
        $totalup = 0;
        while($up >= 0 && $totalup < $PLOTNEIGH){ 
          $up--;
          if($BED[$up] =~ m/^#/){ # non-reference cluster
            next;
	  } else {
            $totalup++;
	  }	
        }
      }

      # get downstream clusters
      $dw = $l;
      $totaldw = 0;
      while($dw <= $#BED && $totaldw < $PLOTNEIGH){
        $dw++;
        if($BED[$dw] =~ m/^#/){ # non-reference cluster
          next;
        } else {
          $totaldw++;
        }
      }
    }
  }

  print "# neighbor/slot indexes: $up < $cl < $dw\n" if($INP_verbose);
  foreach $l ($up .. $dw) {
    push(@plot_clusters_BED,$BED[$l]);
    if($l == $cl) {
      $clslot = scalar(@plot_clusters_BED)
    }
    print "# $BED[$l]\n" if($INP_verbose);
  }
  @BED=();

  # 7.3) parse BED line, FASTA clusters and extract coords
  $total_slots = 0; 
  foreach $l (@plot_clusters_BED) {

    #chr2H NA NA Horvu_MOREX_2H01G436800 1 0 NA Horvu_MOREX_2H01G436800 NA ...
    #chr2H  4 12 gene:HORVU.MOREX.r3.2HG0166240 22 + gene:HORVU.MOREX.r3.2HG0166240 ...

    if($clusterOK > 1 && $total_slots == $clslot) {
      $clusterOK-- # merge slots only as many times
    } else {	    
      $total_slots++;  
    } 

    my @data = split(/\s/,$l); 

    # queue genes from each species to track list
    foreach $col (6 .. $#data) {

      # track1/taxon1: slot1, slot2, .. total_slots
      # track2/taxon2: slot1, slot2, .. total_slots
      # Note: there might 1+ genes on same slot
	    
      my @slot_genes;	    
      $taxon = $tabtaxa[$col-6];

      # actually add genes to this slot  
      foreach $gene_id (split(/,/,$data[$col])) { 	      
        if(!grep(/^$gene_id$/,@{ $plot_tracks{ $taxon }{ $total_slots } })) {
          push(@{ $plot_tracks{ $taxon }{ $total_slots } }, $gene_id);
          push(@slot_genes, $gene_id);  
        }
      }
   
      # get coordinates of these genes from TSV (some genes lack WGA evidence though) 
      #foreach $gene_id (@slot_genes) {
      #  next if(!$TSVdb{$gene_id});
      #  $TSVdata = uncompress($TSVdb{$gene_id}); 
      #  foreach $line (split(/\n/,$TSVdata)) {
      #    ( $gene_id, $segment, $species, $overlap, $homology_type,
      #      $hom_gene_id, $hom_segment, $hom_species, $dummy,
      #      $dummy, $dummy, $dummy, $dummy, $dummy, $coords
      #    ) = split( /\t/, $line );
      #    if($segment ne 'segment' && $coords =~ m/^([^:]+):(\d+)-(\d+)\(([+-])\)/) {
      #     $plot_coords{$species}{$gene_id} = [$1, $2, $3, $4];
      #    }
      #    if($hom_segment ne 'segment' && $coords =~ m/;([^:]+):(\d+)-(\d+)\(([+-])\)/) {
      #      $plot_coords{$hom_species}{$hom_gene_id} = [$1, $2, $3, $4];
      #    } } }
    }

    # work out cluster name and parse genomic coords & strand of genes 
    $clname = $data[3] . '.cdna.fna'; # might not exist, particularly -t all was used
    if(-s "$INP_dir/$cluster_folder/$clname") {
      my ( $ref_geneid, $ref_fasta, $ref_isof_coords, $ref_taxon ) =
        parse_sequence_FASTA_file( "$INP_dir/$cluster_folder/$clname" , 1);
   
      foreach $gene_id (@$ref_geneid) {
          $taxon = $ref_taxon->{$gene_id};
          $plot_coords{$taxon}{$gene_id} = $ref_isof_coords->{$gene_id};
      }
    }
  } 

  # 7.4) write plotting script, track by track, and log
  open(PLOTSCRIPT,">",$plot_scriptfile) ||
    die "# ERROR: cannot create $plot_scriptfile\n";

  open(PLOTLOG,">",$plot_logfile) ||
    die "# ERROR: cannot create $plot_logfile\n";
    
  # add script headers
  print PLOTSCRIPT "from pygenomeviz import GenomeViz\n\n";
  print PLOTSCRIPT "# Python code automatically generated by script check_evidence.pl, \n";
  print PLOTSCRIPT "# see docs and examples at: https://github.com/Ensembl/plant-scripts/tree/master/pangenes\n";
  print PLOTSCRIPT "# It should be run as: python3 $plot_scriptfile\n";
  print PLOTSCRIPT "# If you use the resulting image please cite https://pypi.org/project/pygenomeviz\n\n";  
  print PLOTSCRIPT "genome_list = (\n";

  # compute plot dimensions
  $region_size = (2 * $PLOTFLANKING) + ($PLOTGENEWIDTH * $total_slots);

  foreach $taxon (@tabtaxa) {

    my $total_genes_track = 0;
    my $gene_coords = '';

    print PLOTSCRIPT '{"name": "'. $taxon .'", "size":' . $region_size . ', "gene_list": (';
    
    print PLOTLOG "$taxon";

    # add genes for this taxon, slot by slot
    foreach $slot (1 .. $total_slots) {
	   
       # will be used to shrink gene arrows in the same slot	    
       $total_genes = scalar(@{ $plot_tracks{ $taxon }{ $slot } }); 

       $gene_in_slot = 0;
       my @gene_coords_log;      
       foreach $gene_id (@{ $plot_tracks{ $taxon }{ $slot } }) {
         
         if($gene_id eq 'NA') {
           push(@gene_coords_log,'NA');
           next; 
         }
	
         $shape = 'bigarrow';

	 # compute length of genes in this slot
         $len = int($PLOTGENEWIDTH / $total_genes); 	 

	 # compute gene start
	 $start = 1 + $PLOTFLANKING + (($slot-1) * $PLOTGENEWIDTH) + ($gene_in_slot * $len); 

	 # compute gene end
	 $end = $start + ($len - 1);

	 # get gene strands
         if(defined($plot_coords{$taxon}{$gene_id})) { 
           $strand = 1;
           if($plot_coords{$taxon}{$gene_id}[3] eq '-') {
             $strand = -1;
	   }

           $coords = sprintf( "%s:%d-%d(%s)",
             $plot_coords{$taxon}{$gene_id}[0],
             $plot_coords{$taxon}{$gene_id}[1],
             $plot_coords{$taxon}{$gene_id}[2],
             $plot_coords{$taxon}{$gene_id}[3]);

         } else {
           $shape = 'box';		 
           $coords = 'unk:0-0(?)';
           $strand = 1;

           print "# WARN: Cannot get strand of $gene_id , will plot as box, " .
             "re-run get_pangenes.pl with options -t 0 -s\n";
         }		 

         # choose color 	 
         $color = 'white';
	 if($slot == $clslot) {
           $color = 'tab:green'
         }

         # assign label
         $l = '';
         if($slot == $clslot) {
           $l = $gene_id;
         }	 

         # store and log data for this gene	 
	 push(@gene_coords_log, "$gene_id:$coords");
	 $gene_coords .= "( $start, $end, $strand, '$l', '$shape', '$color' ), ";

         $gene_in_slot++; 	   
       }

       print PLOTLOG "\t".join(',',@gene_coords_log);
    } 

    print PLOTLOG "\n";

    print PLOTSCRIPT $gene_coords;

    print PLOTSCRIPT " )},\n";
  }
  print PLOTSCRIPT ")\n";

  print PLOTSCRIPT <<"ENDOFCODE";
gv = GenomeViz(fig_track_height=0.3, feature_track_ratio=1.0)
ngenomes=0
for genome in genome_list:
  name, size, gene_list = genome["name"], genome["size"], genome["gene_list"]
  if(ngenomes == 0):
    track = gv.add_feature_track(name, size)
  else:
    track = gv.add_feature_track(name, size, linewidth=0)
  ngenomes = ngenomes + 1
  for idx, gene in enumerate(gene_list, 1):
    start, end, strand, glabel, gstyle, color = gene 
    track.add_feature(start, end, strand, label=glabel, plotstyle=gstyle, facecolor=color, linewidth=1, labelrotation=0, labelsize=10)
ENDOFCODE

  print PLOTSCRIPT 'gv.savefig(savefile="'.$plot_file.'",dpi='.$PLOTDPI. ')'."\n"; 

  close(PLOTSCRIPT);
  close(PLOTLOG);

  print "# log file: $plot_logfile\n\n";
  print "# plotting script file: $plot_scriptfile\n\n";

  print "# install if required: pip install pygenomeviz\n"; 
  print "# see other installation options at https://pypi.org/project/pygenomeviz\n\n";
  print "# run it as: python3 $plot_scriptfile\n\n"; 
  print "# will produce: $plot_file\n";
} ## done plotting

untie(%TSVdb);

if(!$INP_fix) {
  exit(0);

} else {
  print "\n";
}




print "# FIX PARAMETERS:\n# -p $INP_partial " .
  "\$MINPAIRPECNONOUTLIERS=$MINPAIRPECNONOUTLIERS \$MINLIFTIDENTITY=$MINLIFTIDENTITY " .
  "\$MINFIXOVERLAP=$MINFIXOVERLAP \$MAXSEGMENTSIZE=$MAXSEGMENTSIZE\n\n";

# 8) suggest fixes for poor gene models based on pan-gene consensus
#    (based on chr coords of 1st isoform of each gene)
my $non_outlier_pairs = 0;
my ($isof,$outisof,$ref_lifted_model);
my ($GFF, $GFFstart, $GFFend);
my (@long_models, @split_models, @non_outliers);
my (@candidate_nonoutliers, %split_seen, %seen_nonoutlier_taxon);

# 8.1) get outlier cutfoff values 
my ($median_pairs, $cutoff_low_pairs, $cutoff_high_pairs) = 
  get_outlier_cutoffs( $scores{'pairs'} , $INP_verbose );
my ($median_len, $cutoff_low_len, $cutoff_high_len) = 
  get_outlier_cutoffs( $scores{'length'} , $INP_verbose );

# 8.2) identify outlier (long/split) models and non-outlier/consensus ones
foreach $full_id (sort {$seen{$b} <=> $seen{$a}} (keys(%seen))){

  # skip genes with internal stop codons in CDS sequences
  next if($badCDS{$full_id});

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

      # Note: ids might include other genes from same taxon 
      # that might have ~consensus length
      foreach my $id (@{ $taxon_genes{ $ref_taxon->{$gene_id} }}) {
        if($gene_length{$id} < $median_len) {

          # only taxa with 2+ split models actually considered, see below
          push(@split_models, $id);
          print "# split $id\n" if($INP_verbose);
        }
      }

  } else { 
    push(@candidate_nonoutliers, $full_id);
  }
}

# taxa with 2+ non-outlier genes might be used to fix long genes,
# others can still be used to fix split genes
foreach $full_id (@candidate_nonoutliers) {

  # this solves the problem of genes parsed in order 
  # in the previous loop, a gene in a pair of split models
  # might initially be a non-outlier candidate
  next if(grep(/^$full_id$/, @split_models));

  $seen_nonoutlier_taxon{ $ref_taxon->{ $fullid2id{$full_id} } }++;

  push(@non_outliers, $full_id);
} 

foreach $species (keys(%seen_nonoutlier_taxon)) {
  if($seen_nonoutlier_taxon{ $species } > 1) {
    $non_outlier_pairs++;
  }
}

# 8.3) suggest model fixes in GFF format, order of priority: long > split > missing

if(!@non_outliers) {
  die "# ERROR: need non-outliers/consensus gene models to fix cluster, exit\n";
}

# open output GFF files if requested
if($INP_outdir) {
  foreach $species (keys(%taxon_genes), keys(%taxon_segments)) {

    next if($outfhandles{$species}); #print ">> $species\n";

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
  # proposed fix: liftover individual consensus models on genomic segment of long gene, 
  # expect 2+ hits from same species on same strand
  foreach $full_id (@long_models) {

    my %lifted;

    $gene_id = $fullid2id{$full_id};

    # check segment size
    if($genome_coords{$full_id} =~ m/^(\S+?):(\d+)-(\d+)\(([+-])\)/) {
      ($chr,$start,$end,$strand) = ($1, $2, $3, $4);
      if($end - $start > $MAXSEGMENTSIZE ) {
        print "# skip segment $chr:$start-$end($strand) [$species] (long gene, too long)\n";
        next;
      }  
    }

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
      
      ISOF: foreach $isof (@isofs) { 

        # skip long/short isoforms, outliers are computed from gene intervals
        foreach $outisof (keys(%{ $outlier_isoform{$ref_taxon->{$hom_gene_id}} })) {
          next ISOF if($isof =~ /$outisof/);
        }

        $cDNA = $isof; # gene name only in FASTA header of cDNA
        $cDNA =~ s/^>.*\n/>$hom_gene_id\n/;  

        $ref_lifted_model = liftover_gmap( "$gene_id,", 
          $segment, $cDNA, $GMAPBIN, $MINLIFTIDENTITY, 
          \@standard_stop_codons, $INP_partial, $INP_verbose ); 

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

      # at least 2 from same species needed
      next if($lifted{ $species }{ 'total' } < 2);

      # sort GFF features and make sure they don't overlap
      $GFF = sort_check_overlap_genes_GFF( @{ $lifted{ $species }{'GFF'} } );
      
      # check overlap of lifted gene pair
      ($GFFstart, $GFFend) = get_gene_coords_GFF( $GFF );
      if($GFFend-$GFFstart < ($end-$start) * $MINFIXOVERLAP) {
        printf("# gene pair overlap with long not enough: %s,%s (%d) %s,%s (%d)\n",
          $start,$end,$end-$start,
          $GFFstart,$GFFend,$GFFend-$GFFstart) if($INP_verbose);
        next;
      }

      next if(!$GFF); # mapped models actually overlap

      # actually print GFF and log
      $gfffh = $outfhandles{$ref_taxon->{$gene_id}} || *STDOUT;

      print "# long gene model: corrected $gene_id [$ref_taxon->{$gene_id}]\n";

      $outputGFF = 
        sprintf("## replaces %s [%s] source=%s matches=%d mismatches=%d indels=%d\n",
          $gene_id,
          $ref_taxon->{$gene_id},
          $species,
          $lifted{$species}{'matches'},
          $lifted{$species}{'mismatches'},
          $lifted{$species}{'indels'});

      $outputGFF .= "$GFF\n";

      # print to GFF in one operation
      print $gfffh $outputGFF;

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

        next if($chr ne $segment_data{'chr'} || 
          $strand ne $segment_data{'strand'});

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
   
    next if( !$segment_data{'models'} || 
      $segment_data{'models'} < 2);

    if($segment_data{'end'}-$segment_data{'start'} > $MAXSEGMENTSIZE ) {

      print "# skip segment " . 
        "$segment_data{'chr'}:$segment_data{'start'}-$segment_data{'end'}($segment_data{'strand'}) " .
        "[$species] (split gene, too long)\n";
      next;
    }

    $segment = cut_genomic_segment_bedtools(
      "$segment_data{'chr'}:$segment_data{'start'}-$segment_data{'end'}($segment_data{'strand'})",
      $genome_file, $BEDTOOLSBIN ); 

    # lift-over consensus models
    my %lifted;
    foreach $hom_full_id (@non_outliers) {

      $hom_gene_id = $fullid2id{$hom_full_id};

      # extract all isoforms/transcripts of this gene,
      # and select the one with most matches in cDNA alignment
      my $best_isof_model;
      my @isofs = extract_isoforms_FASTA($ref_fasta->{$hom_gene_id});

      ISOF: foreach $isof (@isofs) {

        # skip long/short isoforms, outliers are computed from gene intervals
        foreach $outisof (keys(%{ $outlier_isoform{$ref_taxon->{$hom_gene_id}} })) {
          next ISOF if($isof =~ /$outisof/);
        }

        $cDNA = $isof; # gene name only in FASTA header of cDNA
        $cDNA =~ s/^>.*\n/>$hom_gene_id\n/;

        $ref_lifted_model = liftover_gmap( $segment_data{'genes'}, 
          $segment, $cDNA, $GMAPBIN, $MINLIFTIDENTITY, 
          \@standard_stop_codons, $INP_partial, $INP_verbose );

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

    # choose best (consensus) model to replace split ones
    foreach $hom_species (sort {
        $lifted{$b}{'matches'} <=> $lifted{$a}{'matches'} ||
        $lifted{$a}{'mismatches'} <=> $lifted{$b}{'mismatches'} ||
        $lifted{$a}{'indels'} <=> $lifted{$b}{'indels'}
      } keys(%lifted)) {

      # skip genes with multiple mappings
      next if(!$lifted{ $hom_species }{ 'total' } ||
        $lifted{ $hom_species }{ 'total' } > 1);

      $GFF = $lifted{ $hom_species }{'GFF'}->[0];

      # check overlap of lifted gene is enough
      ($GFFstart, $GFFend) = get_gene_coords_GFF( $GFF );
      if($GFFend-$GFFstart < ($segment_data{'end'}-$segment_data{'start'}) * $MINFIXOVERLAP) {
        printf("# consensus gene overlap with split segment not enough: %s,%s (%d) %s,%s (%d)\n",
          $segment_data{'start'},$segment_data{'end'},$segment_data{'end'}-$segment_data{'start'},
          $GFFstart,$GFFend,$GFFend-$GFFstart) if($INP_verbose);
        next;
      }

      # actually print GFF and log
      $gfffh = $outfhandles{ $species } || *STDOUT;

      print "# split gene model: corrected $segment_data{'genes'} [$species]\n";

      $outputGFF =
        sprintf("## replaces %s [%s] source=%s matches=%d mismatches=%d indels=%d\n",
          $segment_data{'genes'},
          $species,
          $hom_species,
          $lifted{$hom_species}{'matches'},
          $lifted{$hom_species}{'mismatches'},
          $lifted{$hom_species}{'indels'});

      $outputGFF .= "$GFF\n";

      print $gfffh $outputGFF;

      last; # take only best
    }

    # TODO: are there intervening TEs?
  }

} elsif($ref_geneid_seg && scalar(@$ref_geneid_seg)) {

  # hypothesis: model potentially exists but failed to be annotated
  # proposed fix: liftover consensus models (isoforms) over matching genomic segment

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

      ISOF: foreach $isof (@isofs) {

        # skip long/short isoforms, outliers are computed from gene intervals
        foreach $outisof (keys(%{ $outlier_isoform{$ref_taxon->{$hom_gene_id}} })) {
          next ISOF if($isof =~ /$outisof/);
        }

        $cDNA = $isof; # gene name only in FASTA header of cDNA
        $cDNA =~ s/^>.*\n/>$hom_gene_id\n/;

        $ref_lifted_model = liftover_gmap( 'missing', 
          $segment, $cDNA, $GMAPBIN, $MINLIFTIDENTITY, 
          \@standard_stop_codons, $INP_partial, $INP_verbose ); 

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
          $best_isof_model->{'GFF'} ); #print "$best_isof_model->{'GFF'}\n";
      }
    }

    # choose best model to replace missing one
    foreach $hom_species (sort {
        $lifted{$b}{'matches'} <=> $lifted{$a}{'matches'} ||
        $lifted{$a}{'mismatches'} <=> $lifted{$b}{'mismatches'} ||
        $lifted{$a}{'indels'} <=> $lifted{$b}{'indels'}
      } keys(%lifted)) {

      # skip genes with multiple mappings
      next if(!$lifted{ $hom_species }{ 'total' } || 
        $lifted{ $hom_species }{ 'total' } > 1);

      $GFF = $lifted{ $hom_species }{'GFF'}->[0];

      # actually print GFF and log
      $gfffh = $outfhandles{ $ref_taxon_seg->{$segment_id} } || *STDOUT;

      print "# missing gene model: corrected $segment_id [$ref_taxon_seg->{$segment_id}]\n";

      $outputGFF = 
        sprintf("## adds %s [%s] source=%s matches=%d mismatches=%d indels=%d\n",
          $segment_id,
          $ref_taxon_seg->{$segment_id},
          $hom_species,
          $lifted{$hom_species}{'matches'},
          $lifted{$hom_species}{'mismatches'},
          $lifted{$hom_species}{'indels'});

      $outputGFF .= "$GFF\n";
   
      print $gfffh $outputGFF;

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


# Lifts over a gene model upon gmap mapping of a query (cDNA) sequence.
# Note: requires a full (Met2stop) CDS as long as that in the cDNA
# Takes 7+1 params:
# i)    original gene_ids (comma separated strings)
# ii)   FASTA string of genomic sequence (target)
# iii)  FASTA string of cDNA sequence (query)
# iv)   path to gmap
# v)    min identity percentage (float)
# vi)   ref to list with stop codon strings
# vii)  allow partial CDS boolean
# viii) verbose boolean, optional
# Returns ref to hash with model features, including 'GFF' and
# several scores such as 'matches', 'mismatches' or 'indels'
# Note: gene_id of lifted models if that of source cDNA
# Note: original gene_id annotated as old_locus_tag in GFF's attributes column
sub liftover_gmap {
  my ($old_gene_ids, $target_fna, $query_cdna, 
      $gmap_path, $min_identity, $ref_stop_codons, $partialCDSOK, $verbose) = @_;

  my ($chr,$start,$end,$offset,$cmd,$CDSseq,$CDScheck);  
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
  
  # 1st run: get alignment summary to parse scores, prints pretty alignment if verbose
  # Note: discard stderr as it gets in the way while parsing stdout
  $cmd = "echo '$target_fna$query_cdna' | $gmap_path -A"; 
  open(GMAP, "$cmd 2>/dev/null |") ||
    die "# ERROR(liftover_gmap): cannot run $cmd\n"; 
  while(<GMAP>) { 
    if(/^\s+Percent identity: (\S+) \((\d+) matches, (\d+) mismatches, (\d+) indels/){ 
      ($identity, $match, $mismatch, $indel) = ($1, $2, $3, $4); 
      if($identity < $min_identity) {
        ($match, $mismatch, $indel) = (0, 0, 0);
      }
    }

    # fails is there are large indels
    #elsif(/^aa\.([cg])\s+\d+\s(.*)/) {
    #  #aa.g        48  N  N  P  S  S  Q  I  T  Y  G  L  T  I  H  H  A  V 
    #  ($seqname,$aa) = ($1,$2);
    #  $aaseq{$seqname} .= $aa }

    print if($verbose);
  }
  close(GMAP); 

  if($match == 0) {
    return \%lifted_model
  }
 
  # 2nd run: are there premature stop codons in the genomic sequence?
  # Note: this causes problems as CDS coords in GFF go until the last aligned base
  $cmd = "echo '$target_fna$query_cdna' | $gmap_path -3";
  open(GMAP, "$cmd 2>/dev/null |") ||
    die "# ERROR(liftover_gmap): cannot run $cmd\n";
  while(<GMAP>) {
    if($. == 2) {
      chomp;
      if(/^(.+)/) {
        $CDSseq = uc($1); # might contain ... and spaces
        $CDSseq =~ s/[^a-zA-Z]//g;
      
        $CDScheck = no_premature_stops( $CDSseq, $ref_stop_codons, $isCDS, $verbose);
        if($partialCDSOK == 0 && $CDScheck != 1 ||
          $partialCDSOK == 1 && $CDScheck < 0){
          return \%lifted_model;
        }
      }
    }
  }
  close(GMAP); 

  # 3nd run: map query cDNA on genomic target, parse GFF and apply coord offset
  $cmd = "echo '$target_fna$query_cdna' | $gmap_path -f 2";
  
  open(GMAP, "$cmd 2>/dev/null |") || 
    die "# ERROR(liftover_gmap): cannot run $cmd\n";  
  while(<GMAP>) {

    last if(/^###/);

    if($gffOK) {

      #print if($verbose);

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
# If two blocks overlap an empty string is returned
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

# Takes sorted GFF string and parses gene features.
# If there are several genes takes start from 1st and end from last.
# Returns:
# i)  start coord (integer)
# ii) end coord (integer)
sub get_gene_coords_GFF {
  my ($GFFstring) = @_;

  my ($gffstart, $gffend, $start, $end) = (-1, -1);
  
  foreach my $line (split(/\n/, $GFFstring)) {
    if($line =~ m/^[^\t]+\t[^\t]+\tgene\t(\d+)\t(\d+)\t/) { 
      ($start, $end) = ($1, $2);  
      if($gffstart == -1){ $gffstart = $start }
      if($end > $gffend) {
        $gffend = $end
      }
    }
  }
  
  return ($gffstart, $gffend);
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

# Takes:
# i)   string with CDS nucleotide sequence string
# ii)  ref to list with stop codons
# iii) optional, boolean to check sequence length is multiple of 3
# iv)  optional, boolean to request verbose output
# Returns:
# 1: if no internal stop codons found and CDS length is multiple of 3,
# 3: if CDS length is NOT multiple of 3 (only if param iii == 1)
# 0: rest of cases
sub no_premature_stops {
  my ($seq, $ref_stop_codons, $check_length, $verbose) = @_;
  my ($codon, $stop);

  if(defined($check_length) && $check_length == 1 && length($seq) % 3) {
    printf("# WARN(no_premature_stops): CDS length (%d) not multiple of 3\n",
      length($seq)) if($verbose);
    return 3
  }
  
  while($seq =~ /(\w{3})/g) {
    $codon = $1;
    foreach $stop (@$ref_stop_codons) {
      if($stop eq $codon && $+[0] < length($seq)) {
        print "# WARN(no_premature_stops): premature stop codon $+[0] < ".
          length($seq)."\n" if($verbose);
        return 0; 
      }
    }
  }

  return 1
}


