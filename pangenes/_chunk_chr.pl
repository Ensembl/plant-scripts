#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

# Takes a GFF & FASTA pair of files and produces a new pair of files with 
# the original chromosomes/contigs split in chunks of contiguous genes.
# A new chunk is created when the next gene on the current chr is further 
# than MAXGENEDIST bp. A chunk is a genomic block containing at least one gene,
# usually also ending with a gene.
#
# Not used anymore, legacy only.

# Copyright [2022-24] 
# EMBL-European Bioinformatics Institute & Estacion Experimental de Aula Dei-CSIC

# perl _chunk_chr.pl -sp oryza_sativa -fa Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
#   -gf Oryza_sativa.IRGSP-1.0.51.gff3 

my $BEDTOOLSEXE = 'bedtools';

my $MAXGENEDIST = 500_000; 

my %main_gff_feats = (
  'gene' => 1,
  'ncRNA_gene' => 1
);

my %skip_gff_feats = (
  'chromosome' => 1,
  'scaffold' => 1
);

my ( $help, $sp1, $fasta1, $bedtools_path, $cmd, $bed) = (0, 0);
my ( $maxdist, $gff1, $outpath ) = ($MAXGENEDIST, '', '');

GetOptions(
  "help|?"       => \$help,
  "sp|species=s" => \$sp1,
  "fa|fasta=s"   => \$fasta1,
  "gf|gff=s"     => \$gff1,
  "d|maxdist=i"  => \$maxdist,
  "o|outpath=s"  => \$outpath,
  "B|bedtools=s"  => \$bedtools_path
) || help_message();

sub help_message {
  print "\nusage: $0 [options]\n\n"
    . "-sp binomial/trinomial species name (required, example: -sp oryza_sativa, used to name outfiles)\n"
    . "-fa genome FASTA filename           (required, example: -fa oryza_sativa.fna)\n"
    . "-gf GFF filename                    (required, example: -gf oryza_sativa.RAPDB.gff)\n"
    . "-d  max distance (bp) tp next gene  (optional, example: -d $MAXGENEDIST)\n"
    . "-o  path to output folder           (optional, default current folder)\n"
    . "-B  path to bedtools binary          (optional, default: -B bedtools)\n\n"
}

if($help || (!$sp1 || !$fasta1 || !$gff1)){ 
  help_message();
  exit(0);
}  

if(!-s $fasta1 || !-s $gff1){
  print "# ERROR: please make sure all input files exist and have content\n";
  exit(-1);
} 

if($maxdist < 1){ 
  print "# ERROR: distance must be positive\n";
  exit(-1);
}

# check binaries
if(!$bedtools_path) {
  $bedtools_path = $BEDTOOLSEXE
}
if(`$bedtools_path` !~ 'sage') {
  print "# ERROR: cannot find binary file $bedtools_path , exit\n";
  exit(-1)
}

print "\n# $0 -sp $sp1 -fa $fasta1 -gf $gff1 -d $maxdist -o $outpath -B $bedtools_path\n\n";



# set output filenames
my $chunkfnafile = "$sp1.chunk$maxdist.fna";
my $chunkgfffile = "$sp1.chunk$maxdist.gff";
my $chunkbedfile = "$sp1.chunk$maxdist.bed";
if(-e $outpath) {
  $chunkfnafile = "$outpath/$chunkfnafile";
  $chunkgfffile = "$outpath/$chunkgfffile";
  $chunkbedfile = "$outpath/$chunkbedfile";
}

my ($ref_chrs, $ref_chunk_genes, $ref_chunks) = 
  chunk_GFF($gff1, $maxdist, \%main_gff_feats, \%skip_gff_feats);

if(scalar(keys(%$ref_chunks)) == 0) {
  die "# ERROR: cannot chunk GFF file ($gff1)\n";
}

open(GFFCHUNK,">",$chunkgfffile) ||
  die "# ERROR: cannot open chunk GFF file ($chunkgfffile)\n";

open(FNACHUNK,">",$chunkfnafile) ||
  die "# ERROR: cannot open chunk FASTA file ($chunkfnafile)\n";

open(BEDCHUNK,">",$chunkbedfile) ||
  die "# ERROR: cannot open chunk BED file ($chunkbedfile)\n";

my $total_chunks = 0;
foreach my $chr (@$ref_chrs) {
  foreach my $chunk (sort {$a<=>$b} keys(%{$ref_chunks->{$chr}})) {

    # print transformed gene models to chunked GFF file
    print GFFCHUNK $ref_chunk_genes->{$chr}{$chunk};

    # print sequence to chunked FASTA file
    $bed = sprintf("%s\t%d\t%d", 
      $chr,
      $ref_chunks->{$chr}{$chunk}{'start'}-1, #0-based
      $ref_chunks->{$chr}{$chunk}{'end'});

    $cmd = "echo '$bed' | $bedtools_path getfasta -fi $fasta1 -bed stdin";
    open(BEDTOOLS,"$cmd |") ||
      die "# ERROR: cannot run bedtools ($cmd)\n";
    while(<BEDTOOLS>) {
      if(/^>/) { 
        print FNACHUNK ">$chr\.chunk$chunk\n";
      } else {
        print FNACHUNK;
      }
    }
    close(BEDTOOLS);

    # log
    print BEDCHUNK "$bed\t$chr\.chunk$chunk\n";

    $total_chunks++ 
  }
}

close(BEDCHUNK);
close(FNACHUNK);
close(GFFCHUNK);

print "# chunked GFF file: $chunkgfffile\n";
print "# chunked FASTA file: $chunkfnafile\n";
print "# chunked BED file: $chunkbedfile\n";


printf("\n# total chr/contigs=%d total chunks=%d\n",
  scalar(@$ref_chrs), 
  $total_chunks);


###############################

# Parses GFF file and finds chunks of contiguous genes.
# Returns: 
# i)   ref to hash mapping chunk ID to translated gene models in chunks
# ii)  ref to hash mapping chunk ID to chunk 1-based coordinates in original FASTA
# iii) ref to list with chr names in same order as input
sub chunk_GFF {

  my ($gff_file, $maxdist, $ref_main_gff, $ref_skip_gff) = @_;

  my ($chr, $start, $end, $chunk_start, $chunk_end, $gff_line);
  my ($num_chunk, $dist, $offset, $prev_end) = (1, 0, 0, 0);
  my (%chunk_genes, %chunk, @chrs);

  open(GFF,"<",$gff_file) || 
    die "# ERROR(chunk_GFF): cannot read $gff_file\n";
  while(<GFF>){

    next if(/^#/ || /^$/);

    my @gff = split(/\t/,$_);
    ($chr, $start, $end) = @gff[0,3,4];

    next if($ref_skip_gff->{ $gff[2] });

    if($ref_main_gff->{ $gff[2] }) { 

      # new chunk with new chr
      if($num_chunk > 1 && !grep(/^$chr$/,@chrs)) {
        push(@chrs, $chr);
        $prev_end = 0;
        $num_chunk++;

      } elsif($prev_end > 0) { # new chunk if previous gene too far

        $dist = $start-$prev_end;
        if($dist > $maxdist) {
          $num_chunk++;
        }
      }

      # set chunk start and offset, 
      # used to cut sequence & to transform gene coords
      if(!defined($chunk{$chr}{$num_chunk}{'start'})) { 
        $chunk{$chr}{$num_chunk}{'start'} = $start;
        $chunk{$chr}{$num_chunk}{'offset'} = $start-1;
      }

      # update chunk last coord with every new gene
      $chunk{$chr}{$num_chunk}{'end'} = $end;

      # save end coord for next iteration
      $prev_end = $end; 
    } 

    # transform coords relative to current chunk
    $gff[0] = "$chr.chunk$num_chunk"; 
    $gff[3] -= $chunk{$chr}{$num_chunk}{'offset'};
    $gff[4] -= $chunk{$chr}{$num_chunk}{'offset'};
    $gff_line = join("\t",@gff); 

    $chunk_genes{$chr}{$num_chunk} .= $gff_line;
  }
  close(GFF);

  return (\@chrs, \%chunk_genes, \%chunk);
}

