#!/usr/bin/env perl

# This script takes a cDNA/CDS cluster produced by get_pangenes.pl and 
# produces a quality control report 

# Copyright [2023]
# EMBL-European Bioinformatics Institute & Estacion Experimental Aula Dei-CSIC

$|=1;

use strict;
use warnings;
use Getopt::Std;
use FindBin '$Bin';
use lib "$Bin/lib";
use pangeneTools qw( check_installed_features feature_is_installed 
                     parse_sequence_FASTA_file calc_stdev );

my @FEATURES2CHECK = (
  'EXE_CLUSTALO', 'EXE_ALISTAT'
);

my ($INP_dir, $INP_clusterfile, $INP_first_isof, $INP_outdir) = ('','',0,'');
my ($isCDS, $ispep, $seq) = ( 0, 0 );
my ($cluster_list_file,$cluster_folder, $gene_id, $isof_id);
my (%opts, %isof_len, %isof_seq, %isof_header, %isof_order, @len);

getopts('hIco:d:i:', \%opts);

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
  $INP_dir = $opts{'d'} 
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

  my $n_isof = 0;
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


    last if($INP_first_isof == 1);
  }
}

# 3) print selected isoform sequence(s) to temp file and work out basic stats 
#  open(ISOSEQ,">>",$INP_modeseq) ||
#    die "# EXIT: cannot write to $INP_modeseq\n";

foreach $gene_id (@$ref_geneid) {
  foreach $isof_id (keys(%{$isof_len{$gene_id}})) {

    next if($INP_first_isof == 1 && $isof_order{$gene_id}{$isof_id} != 1);
    
    print "$isof_header{$gene_id}{$isof_id}\n"; #$isof_seq{$gene_id}{$isof_id}\n";
  }
}


#  close(ISOSEQ);





# read all or 1st isoform

# basic stats: occup, seqs, length, exons

# MSA & distance


# MSA report
