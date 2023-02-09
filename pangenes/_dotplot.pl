#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw( dirname );
use FindBin '$Bin';
use lib "$Bin/lib";
use pangeneTools qw( read_FAI_regex2hash );

$|=1;

# Parses pairwise collinear TSV files made with _collinear_genes.pl to 
# produce PAF files that can be used to produce a dotplot of matched gene
# models with R package pafr [https://cran.r-project.org/package=pafr]
# Note: contigs < $MINCONTIGSIZE are ignored for clarity

# Copyright [2022-23] 
# EMBL-European Bioinformatics Institute & Estacion Experimental de Aula Dei-CSIC

# perl _dotplot.pl _Oryza_nivara_v1.Oryza_sativa.IRGSP-1.0.algMmap.overlap0.5.tsv 

my $MINCONTIGSIZE = 100_000;
my $DUMMYQUAL = 60;
my $VERBOSE = 0;

my $TSVfile = $ARGV[0] || die "# usage: $0 <TSVfile>\n";

my $resultsDIR = dirname($TSVfile);

my ($sp1filename, $sp2filename) = ('', '');
my ($faifile, $sp1, $sp2, $species, $chr, $len) = ('','','');
my (%file,%size);

## locate FASTA index files (.fai), usually created in _cut_sequences.pl
opendir(DIR,$resultsDIR) || die "# ERROR: cannot list $resultsDIR\n";
my @faifiles = grep {/\.fai$/} readdir(DIR);
closedir(DIR);

foreach $faifile (@faifiles) {

  $species = $faifile; 
  $species =~ s/^_//;
  $species =~ s/\.fna.fai$//;

  if($TSVfile =~ m/_$species\./) {
    $sp1 = $species;
    $sp1filename = $sp1;
    $file{ $sp1 } = "$resultsDIR/$faifile";
  } elsif($TSVfile =~ m/\.$species\./) {
    $sp2 = $species;
    $sp2filename = $sp2;
    $file{ $sp2 } = "$resultsDIR/$faifile";
  }
} 

## parse contig sizes from FASTA index files
if($sp1 eq '' || $sp2 eq '') {
  die "# ERROR: cannot find FASTA indexes (.fai) for $TSVfile ($sp1, $sp2)\n";
} else {
  for $species ($sp1, $sp2) {
    my $ref_bed = read_FAI_regex2hash( $file{$species}  );
    foreach $chr (keys(%$ref_bed)) {
      $len = (split(/\t/,$ref_bed->{$chr}))[2];
      $size{$species}{$chr} = $len;
	  print "# $species $chr $len\n" if($VERBOSE);
    }
  }
}

## parse TSV file and convert it to PAF format so that pafr can take it,
## see https://dwinter.github.io/pafr/articles/Introduction_to_pafr.html

my $outPAFfile = $TSVfile;
$outPAFfile =~ s/\.tsv$/.genes.paf/;

open(PAF,">",$outPAFfile) || die "# ERROR: cannot create $outPAFfile\n";

my (@data,$chr1,$start1,$end1,$chr2,$start2,$end2);
open(TSV,"<",$TSVfile) || die "# ERROR: cannot read $TSVfile\n";
while(<TSV>) {
 
  @data = split(/\t/,$_);

  # OsMH63_01G000010 OsMH63_01G000010 oryza_sativa_mh63 .. ortholog_collinear .. oryza_sativa .. 1:7538-15379(+);1:2902-10817(+)
  # Note that after
  # https://github.com/Ensembl/plant-scripts/blob/f9c9e4e71fbb8e46f84b0609a6dfc1dd5930bacf/pangenes/_collinear_genes.pl#L1774
  # in segments species order might change
  # Os01g0100466 Os01g0100466 oryza_sativa .. segment_collinear .. oryza_sativa_mh63:1:17370-18541(-) .. 1:12807-13978(-);1:17370-18541(-)

  if($data[14] =~ m/(\S+)?:(\d+)-(\d+)\([+-]\);(\S+)?:(\d+)-(\d+)\([+-]\)/) {
   
    $sp1 = $data[2];
    $sp2 = $data[7];
    ($chr1,$start1,$end1,$chr2,$start2,$end2) = ($1,$2,$3,$4,$5,$6);

    # skip unknown contigs ie unplaced_chrUn, might occur with -s
    next if(!$size{$sp1}{$chr1} || !$size{$sp2}{$chr2});

    next if($size{$sp1}{$chr1} < $MINCONTIGSIZE || 
      $size{$sp2}{$chr2} < $MINCONTIGSIZE);

    if($sp1 eq $sp1filename) {
      printf( PAF "%s\t%d\t%d\t%d\t+\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
        $chr1,$size{$sp1}{$chr1},$start1-1,$end1,
        $chr2,$size{$sp2}{$chr2},$start2-1,$end2,
        $data[3], # overlap instead of matching bases in the mapping
        $data[3], # overlap instead of bases, including gaps, in the mapping
        $DUMMYQUAL);
	} else {
      printf( PAF "%s\t%d\t%d\t%d\t+\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
	    $chr2,$size{$sp2}{$chr2},$start2-1,$end2,
        $chr1,$size{$sp1}{$chr1},$start1-1,$end1,
        $data[3], # overlap instead of matching bases in the mapping
        $data[3], # overlap instead of bases, including gaps, in the mapping
        $DUMMYQUAL);
	}
  }
}
close(TSV);

close(PAF);

print "\n# \$MINCONTIGSIZE = $MINCONTIGSIZE\n";
print "\n# PAF file: $outPAFfile\n\n";

print "# Make a dotplot of aligned models coords with the following R script:\n";

print<<EOF;

  #https://dwinter.github.io/pafr/articles/Introduction_to_pafr.html
  #install.packages(devtools)
  #devtools::install_github("dwinter/pafr")

  library(pafr, quietly=TRUE)

  pafile = "$outPAFfile"
  ali <- read_paf(pafile)

  dotplot(ali, label_seqs = TRUE, xlab='$sp1', ylab='$sp2')

  #if chr/contig coverage wanted instead
  #plot_coverage(ali) + scale_fill_brewer()

EOF

