#!/usr/bin/env perl

# This script takes a cDNA/CDS cluster produced by get_pangenes.pl and 
# produces a quality control report 

# Copyright [2023]
# EMBL-European Bioinformatics Institute & Estacion Experimental Aula Dei-CSIC

$|=1;

use strict;
use warnings;
use Getopt::Std;
use File::Temp qw/ tempfile /;
use FindBin '$Bin';
use lib "$Bin/lib";
use pangeneTools qw( check_installed_features feature_is_installed 
                     parse_sequence_FASTA_file calc_stdev calc_mode );

my @FEATURES2CHECK = (
  'EXE_CLUSTALO', 'EXE_ALISTAT', 'EXE_GREP'
);

my ($INP_dir, $INP_clusterfile, $INP_first_isof, $INP_noheader, $INP_outdir) = ('','',0,0,'');
my ($isCDS, $ispep, $seq, $n_isof, $occup, $SE_len, $SE_exons) = ( 0, 0 );
my ($updir, $n_exons, $gff_file, $SE_dist, $max_dist, $c);
my ($cluster_list_file,$cluster_folder, $gene_id, $isof_id);
my ($msa_filename, $dist_filename, $fhmsa, $fhdist, $cmd);
my ($sites, $Ca, $Cr_max, $Cr_min, $Cc_max, $Cc_min, $Cij_max, $Cij_min);
my (%opts, %isof_len, %isof_seq, %isof_header, %isof_order);
my (%taxa, @len, @exons, @dist, @mode_len, @mode_exons, @mode_dist);

getopts('hnIco:d:i:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-c print credits and checks installation\n";
  print "-d directory produced by get_pangenes.pl        (example: -d /path/data_pangenes/..._algMmap_,\n";
  print "                                                 genomic sequences usually one folder up)\n";
  print "-i cdna/cds .fna/.faa file as in .cluster_list  (example: -i gene:ONIVA01G52180.cdna.fna)\n";
  print "-I take 1st isoform only                        (optional, by default takes all)\n";
  print "-o folder to write output files                 (optional, MSA files removed by default)\n";
  print "-n do not print header in text report           (optional)\n";
  exit(0);
}

if(defined($opts{'c'})) {
  print "\nPrimary citation:\n https://github.com/Ensembl/plant-scripts/pangenes\n";
  print "\nThis software uses external algorithms, please cite them accordingly:\n";
  print " clustal-omega https://doi.org/10.1002%2Fpro.3290\n";
  print " AliStat https://doi.org/10.1093/nargab/lqaa024\n";

  # check all binaries needed by this program and print diagnostic info
  print check_installed_features(@FEATURES2CHECK);
  exit(0);
}

if(defined($opts{'d'})) { 
  $INP_dir = $opts{'d'};
  $updir = $INP_dir . '/..'; 
}
else{ die "# EXIT : need a -d directory\n" }

if(defined($opts{'i'})){  
  $INP_clusterfile = $opts{'i'};
  if($INP_clusterfile !~ /\.cdna\.fna$/ && 
    $INP_clusterfile !~ /\.cds\.fna$/ &&
    $INP_clusterfile !~ /\.cds\.faa$/) {
    die "# EXIT : need a .fna/.faa cluster filename with parameter -i\n"

  } else {
    if($INP_clusterfile =~ /\.cds\.f/) {
      $isCDS = 1
    }

    if($INP_clusterfile =~ /\.cds\.faa/) { 
      $ispep = 1
    }
  }
}
else{ die "# EXIT : need parameter -i\n" }

if(defined($opts{'I'})){ 
  $INP_first_isof = 1 
}

if(defined($opts{'o'})){
  $INP_outdir = $opts{'o'};
  if(!-e $INP_outdir) {
    mkdir($INP_outdir);
  }
}

if(defined($opts{'n'})){
  $INP_noheader = 1
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

# 2) parse FASTA file, extract gene names and sequence lengths
my ( $ref_geneid, $ref_fasta, $ref_isof_coords, $ref_taxon ) = 
  parse_sequence_FASTA_file( "$INP_dir/$cluster_folder/$INP_clusterfile" , 1);

foreach $gene_id (@$ref_geneid) {

  $n_isof = 0;
  foreach $seq (split(/\n/,$ref_fasta->{$gene_id})) {

    if($seq =~ /^>(\S+)/) {
      $n_isof++;
      $isof_id = $1;
      $isof_header{$gene_id}{$isof_id} = $seq;
      $isof_order{$gene_id}{$isof_id} = $n_isof;
      next;
    }
    $isof_len{$gene_id}{$isof_id} += length($seq);
    $isof_seq{$gene_id}{$isof_id} .= $seq;
  }
}

# 3) print selected isoform sequence(s) to temp file and work out basic stats 
my ($fh, $filename) = tempfile( 'tempfasXXXXX', UNLINK => 1);

foreach $gene_id (@$ref_geneid) {
  foreach $isof_id (keys(%{$isof_len{$gene_id}})) {

    next if($INP_first_isof == 1 && $isof_order{$gene_id}{$isof_id} != 1);
    
    $taxa{ $ref_taxon->{$gene_id} }++;
    push(@len, $isof_len{$gene_id}{$isof_id});

    # find GFF file & get number of exons
    $n_exons = 0;
    $gff_file = $updir . "/_$ref_taxon->{$gene_id}.gff";
    open(GREP, "$ENV{'EXE_GREP'} '$isof_id;' $gff_file |");
    while(<GREP>) {
      #1	NAM	exon	2575663	2575953 ...
      my @data = split(/\t/,$_);

      if($isCDS == 0 && $data[2] eq 'exon') {
        $n_exons++

      } elsif($isCDS == 1 && $data[2] eq 'CDS') {
        $n_exons++
      }
    }
    close(GREP);  
    push(@exons, $n_exons);

    # actually print to temp file
    print $fh "$isof_header{$gene_id}{$isof_id}\n$isof_seq{$gene_id}{$isof_id}\n";
  }
}

$occup = scalar(keys(%taxa));
$n_isof = scalar(@len); # recompute in case inly 1st isoform taken

$SE_len = sprintf("%1.1f", calc_stdev( \@len ) / sqrt($n_isof));
@mode_len = calc_mode( \@len ); 
$SE_exons = sprintf("%1.1f", calc_stdev( \@exons ) / sqrt($n_isof));
@mode_exons = calc_mode( \@exons );


# 4) compute multiple sequence alignment (MSA), distance matrix & MSA report

if($INP_outdir ne '') {
  $msa_filename = "$INP_outdir/$INP_clusterfile";
  $msa_filename =~ s/\.(f[na]a)$/.aln.$1/;
  $dist_filename = "$INP_outdir/$INP_clusterfile";
  $dist_filename =~ s/\.(f[na]a)$/.dist.$1/;
} else {
  ($fhmsa, $msa_filename) = tempfile( 'tempmsaXXXXX', UNLINK => 1);
  ($fhdist, $dist_filename) = tempfile( 'tempdistXXXXX', UNLINK => 1);
}

$cmd = "$ENV{'EXE_CLUSTALO'} --force --full -i $filename -o $msa_filename --distmat-out=$dist_filename 2>&1";
system($cmd); 
if ( $? != 0 ) {
  die "# ERROR: failed running clustal-omega ($cmd)\n";
} elsif ( !-s $msa_filename ) {
  die "# ERROR: failed generating $msa_filename file ($cmd)\n";
}

# parse MSA distances
$max_dist = -1;
open(DIST,"<",$dist_filename) ||
  die "# ERROR: cannot read $dist_filename\n";
while(<DIST>) {
  chomp;
  my @data = split(/\s+/,$_);
  next if($#data < 1);
  foreach $c (1 .. $#data) { 
    push(@dist, $data[$c]);
    if($data[$c] > $max_dist){ $max_dist = $data[$c] } 
  }  
}
close(DIST);

$SE_dist = sprintf("%1.1f", calc_stdev( \@dist ) / $n_isof);
@mode_dist = calc_mode( \@dist ); 

# MSA report
$cmd = "$ENV{'EXE_ALISTAT'} $msa_filename 1 -b";
if($ispep) {
  $cmd = "$ENV{'EXE_ALISTAT'} $msa_filename 6 -b";
}

open(ALISTAT,"$cmd |") ||
  die "# ERROR: cannot run $cmd\n";
while(<ALISTAT>) {
  #sequences, #sites, Ca, Cr_max, Cr_min, Cc_max, Cc_min, Cij_max, Cij_min
  if(/^$msa_filename/) {
    chomp;
    my @data = split(/,\s+/,$_);
    ($sites, $Ca, $Cr_max, $Cr_min, $Cc_max, $Cc_min, $Cij_max, $Cij_min) = @data[2 .. $#data];
  }
}
close(ALISTAT);


# 5) finally print summary in one line
if($INP_noheader == 0) {
  print "file\t1stisof\toccup\tseqs\tmode_len\tSE_len\tmode_exons\tSE_exons\t" .
    "mode_dist\tmax_dist\tSE_dist\tsites\tCa\tCr_max\tCr_min\tCc_max\tCc_min\tCij_max\tCij_min\n";
}

printf(
  "%s\t%d\t%d\t%d\t%d\t%1.1f\t" .
    "%d\t%1.1f\t%1.6f\t%1.6f\t%1.6f\t%d\t" .
    "%1.6f\t%1.6f\t%1.6f\t%1.6f\t%1.6f\t%1.6f\t%1.6f\n",

  $INP_clusterfile,
  $INP_first_isof,  
  $occup,
  $n_isof,
  $mode_len[0],
  $SE_len,

  $mode_exons[0],
  $SE_exons,
  $mode_dist[0],
  $max_dist,
  $SE_dist,
  $sites, 

  $Ca, 
  $Cr_max, 
  $Cr_min, 
  $Cc_max, 
  $Cc_min, 
  $Cij_max, 
  $Cij_min);
