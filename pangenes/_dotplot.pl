#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw(dirname);
use FindBin '$Bin';
use lib "$Bin/lib";
use pangeneTools qw( read_FAI_regex2hash );

$|=1;

# Parses pairwise collinear TSV files made with _collinear_genes.pl to 
# produce PAF files that can be used to produce a dotplot of matched gene
# models with R package pafr [https://cran.r-project.org/package=pafr]

# Copyright [2022] 
# EMBL-European Bioinformatics Institute & Estacion Experimental de Aula Dei-CSIC

# perl _dotplot.pl _Oryza_nivara_v1.Oryza_sativa.IRGSP-1.0.algMmap.overlap0.5.tsv 

my $TSVfile = $ARGV[0] || die "# usage: $0 <TSVfile>\n";

my $resultsDIR = dirname($TSVfile);

my ($faifile, $sp1, $sp2, $species) = ('','','');

# locate FASTA index files (.fai), usually created in _cut_sequences.pl
opendir(DIR,$resultsDIR) || die "# ERROR: cannot list $resultsDIR\n";
my @faifiles = grep {/\.fai$/} readdir(DIR);
closedir(DIR);

foreach $faifile (@faifiles) {

  $species = $faifile; 
  $species =~ s/^_//;
  $species =~ s/\.fna.fai$//;

  if($TSVfile =~ m/_$species\./) {
    $sp1 = $species;
  } elsif($TSVfile =~ m/\.$species\./) {
    $sp2 = $species;
  }
} 

if($sp1 eq '' || $sp2 eq '') {
  die "# ERROR: cannot find FASTA indexes (.fai) for $TSVfile ($sp1, $sp2)\n";
}


#$index_fasta1 = dirname($chrfasta1)."/".basename($chrfasta1).".fai";

#    die "# ERROR: failed mapping $sp2 genes in WGA alignment";
