#!/usr/bin/env perl

# This script matches input sequences to clusters produced by get_pangenes.pl
# It creates a sequence index with nucleotide sequences from clusters and 
# uses GMAP tp scan them. Can use cdna [default ]or cds sequences.

# Copyright [2023]
# EMBL-European Bioinformatics Institute & Estacion Experimental Aula Dei-CSIC

$|=1;

use strict;
use warnings;
use Getopt::Std;
use File::Temp qw/ tempfile /;
use File::Basename;
use FindBin '$Bin';
use lib "$Bin/lib";
use pangeneTools qw( check_installed_features feature_is_installed 
                     parse_sequence_FASTA_file );

my @FEATURES2CHECK = (
  'EXE_GMAP'
);

my $GMAPBIN = $ENV{'EXE_GMAP'};
my $GMAPBUILDBIN = $ENV{'EXE_GMAP_BUILD'};

my ($INP_dir, $INP_seqfile, $INP_isCDS, $INP_idfile) = ('','',0,'');
my ($cluster_regex, $isCDS, $gmapdb_path, $gmapdb, $seq) = ( '.cdna.fna', 0 );
my ($cluster_file, $gene_id, $isof_id, $cmd, $index_file);
my (%opts);

getopts('hcCI:d:s:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-c print credits and checks installation\n";
  print "-d directory with pangene clusters.pl             (required, should contain *.cdna.fna or *.cds.fna files)\n";
  print "-s nucleotide sequence file in FASTA format       (required, example: -s transcripts.fna)\n";
  print "-C use CDS sequences                              (optional, cDNA sequences are scanned by default)\n";
  print "-I TSV file matching cluster names to pangene ids (optional, example: -I cluster2id.tsv)\n\n";
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
  $INP_dir = $opts{'d'};
  $gmapdb = basename($INP_dir) . '.gmap';
  $gmapdb_path = dirname($INP_dir); 
}
else{ die "# EXIT : need a -d directory\n" }

if(defined($opts{'s'})){  
  $INP_seqfile = $opts{'s'};
} 
else{ die "# EXIT : need parameter -s\n" }

if(defined($opts{'C'})){
  $INP_isCDS = 1;
  $cluster_regex = '.cds.fna';
}

if(defined($opts{'I'})){ 
  $INP_idfile = $opts{'I'};  
}

# 1) parse cluster files and write sequences to temp file
opendir(INPDIR,$INP_dir) || 
  die "# ERROR: cannot list $INP_dir , please check -d argument is a valid folder\n";
my @files = grep {/$cluster_regex$/} readdir(INPDIR);
closedir(INPDIR);

if(scalar(@files) == 0) {
  die "# ERROR: cannot find any valid clusters in $INP_dir\n";
}

my ($fh, $filename) = tempfile( 'tempgmapfnaXXXXX', UNLINK => 1);

foreach $cluster_file (@files) {

  my ( $ref_geneid, $ref_fasta, $ref_isof_coords, $ref_taxon ) = 
    parse_sequence_FASTA_file( "$INP_dir/$cluster_file" , 1);

  # prepend cluster name and print this sequence 
  foreach $gene_id (@$ref_geneid) {
    $seq = $ref_fasta->{$gene_id};
    $seq =~ s/^>/>$cluster_file-/;   
    print $fh $seq;
  }
}

## 2) index temp sequences for GMAP, skip if done previously
$cmd = "$GMAPBUILDBIN -d $gmapdb -D $gmapdb_path $filename 2>&1";
$index_file = "$gmapdb.ref153positions";

if(!-d "$gmapdb_path/$gmapdb" || 
  !-e "$gmapdb_path/$gmapdb/$index_file") {

  print "# building GMAP index $gmapdb_path/$gmapdb\n\n";
  system($cmd);
  if ( $? != 0 ) {
    die "# ERROR: failed running $GMAPBUILDBIN ($cmd)\n";
  } else {
    print "# created GMAP index $gmapdb_path/$gmapdb\n\n";
  }
} else {
  print "# re-using GMAP index $gmapdb_path/$gmapdb\n\n";
}

## 3) map input sequences to cluster sequences
$cmd = "$GMAPBIN -S -d $gmapdb -D $gmapdb_path $INP_seqfile 2>&1";
open(GMAP, "$cmd |") ||
  die "# ERROR: cannot run $cmd\n"; 
while(<GMAP>) { 
	#if(/^\s+Percent identity: (\S+) \((\d+) matches, (\d+) mismatches, (\d+) indels/){ 
	#($identity, $match, $mismatch, $indel) = ($1, $2, $3, $4); 
	#if($identity < $min_identity) {
	#($match, $mismatch, $indel) = (0, 0, 0);
	#}
	#}

    # fails is there are large indels
    #elsif(/^aa\.([cg])\s+\d+\s(.*)/) {
    #  #aa.g        48  N  N  P  S  S  Q  I  T  Y  G  L  T  I  H  H  A  V 
    #  ($seqname,$aa) = ($1,$2);
    #  $aaseq{$seqname} .= $aa }
    print;
}
close(GMAP);





#foreach $gene_id (@$ref_geneid) {
#  foreach $isof_id (keys(%{$isof_len{$gene_id}})) {

#    $taxa{ $ref_taxon->{$gene_id} }++;

    # find GFF file & get number of exons
    #$n_exons = 0;
    #$gff_file = $updir . "/_$ref_taxon->{$gene_id}.gff";
    #open(GREP, "$ENV{'EXE_GREP'} '$isof_id;' $gff_file |");
    #while(<GREP>) {
    #  #1	NAM	exon	2575663	2575953 ...
    #  my @data = split(/\t/,$_);

    #  if($isCDS == 0 && $data[2] eq 'exon') {
    #    $n_exons++

    #  } elsif($isCDS == 1 && $data[2] eq 'CDS') {
    #    $n_exons++
    #  }
    #}
    #close(GREP);  
    #push(@exons, $n_exons);

    # actually print to temp file
    #print $fh "$isof_header{$gene_id}{$isof_id}\n$isof_seq{$gene_id}{$isof_id}\n";
    #  }
    #}

#$occup = scalar(keys(%taxa));
#$n_isof = scalar(@len); # recompute in case inly 1st isoform taken


# 4) compute multiple sequence alignment (MSA), distance matrix & MSA report

#if($INP_outdir ne '') {
#  $msa_filename = "$INP_outdir/$INP_clusterfile";
#  $msa_filename =~ s/\.(f[na]a)$/.aln.$1/;
#  $dist_filename = "$INP_outdir/$INP_clusterfile";
#  $dist_filename =~ s/\.(f[na]a)$/.dist.$1/;
#} else {
#  ($fhmsa, $msa_filename) = tempfile( 'tempmsaXXXXX', UNLINK => 1);
#  ($fhdist, $dist_filename) = tempfile( 'tempdistXXXXX', UNLINK => 1);
#}

#$cmd = "$ENV{'EXE_CLUSTALO'} --force --full -i $filename -o $msa_filename --distmat-out=$dist_filename 2>&1";
#system($cmd); 
#if ( $? != 0 ) {
#  die "# ERROR: failed running clustal-omega ($cmd)\n";
#} elsif ( !-s $msa_filename ) {
#  die "# ERROR: failed generating $msa_filename file ($cmd)\n";
#}

# parse MSA distances
#$max_dist = -1;
#open(DIST,"<",$dist_filename) ||
#  die "# ERROR: cannot read $dist_filename\n";
#while(<DIST>) {
#  chomp;
#  my @data = split(/\s+/,$_);
#  next if($#data < 1);
#  foreach $c (1 .. $#data) { 
#    push(@dist, $data[$c]);
#    if($data[$c] > $max_dist){ $max_dist = $data[$c] } 
#  }  
#}
#close(DIST);


# MSA report
#$cmd = "$ENV{'EXE_ALISTAT'} $msa_filename 1 -b";
#if($ispep) {
#  $cmd = "$ENV{'EXE_ALISTAT'} $msa_filename 6 -b";
#}

