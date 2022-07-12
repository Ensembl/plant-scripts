#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use File::Copy qw(cp);

# Takes a GFF & FASTA pair of files and produces a pair of files with 
# the original chromosomes/contigs split in chunks of contiguous genes.
# A new chunk is created when the next gene on the current chr is further 
# than MAXGENEDIST

# Copyright [2022] 
# EMBL-European Bioinformatics Institute & Estacion Experimental de Aula Dei-CSIC

# perl _cut_sequences.pl -sp oryza_sativa -fa Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
#   -gf Oryza_sativa.IRGSP-1.0.51.gff3 

my $MAXGENEDIST = 500_000; 

my ( $help, $sp1, $fasta1) = (0, 0);
my ( $maxdist, $gff1, $outpath ) = ($MAXGENEDIST, '', '');

GetOptions(
  "help|?"       => \$help,
  "sp|species=s" => \$sp1,
  "fa|fasta=s"   => \$fasta1,
  "gf|gff=s"     => \$gff1,
  "d|maxdist=i"  => \$maxdist,
  "o|outpath=s"  => \$outpath
) || help_message();

sub help_message {
  print "\nusage: $0 [options]\n\n"
    . "-sp binomial/trinomial species name (required, example: -sp oryza_sativa, used to name outfiles)\n"
    . "-fa genome FASTA filename           (required, example: -fa oryza_sativa.fna)\n"
    . "-gf GFF filename                    (required, example: -gf oryza_sativa.RAPDB.gff)\n"
    . "-d  max distance (bp) tp next gene  (optional, example: -d $MAXGENEDIST)\n"
    . "-o  path to output folder           (optional, default current folder)\n\n"
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

print "\n# $0 -sp $sp1 -fa $fasta1 -gf $gff1 -d $maxdist -o $outpath\n\n";

# set output filenames
my $chunkfnafile = "$sp1.chunk$maxdist.fna";
my $chunkgfffile = "$sp1.chunk$maxdist.gff";
if($outpath) {
  $chunkfnafile = "$outpath/$chunkfnafile";
  $chunkgfffile = "$outpath/$chunkgfffile";
}

my ($ref_chrs, $ref_chunk_genes, $ref_chunks) = chunk_GFF($gff1, $maxdist);

if(scalar(keys(%$ref_chunks)) == 0) {
  die "# ERROR: cannot chunk GFF file ($gff1)\n";
}

open(GFFCHUNK,">",$chunkgfffile) ||
  die "# ERROR: cannot open chunk GFF file ($chunkgfffile)\n";

open(FNACHUNK,">",$chunkfnafile) ||
  die "# ERROR: cannot open chunk FASTA file ($chunkfnafile)\n";

my $total_chunks = 0;
foreach my $chr (@$ref_chrs) {
  foreach my $chunk (sort {$a<=>$b} keys(%{$ref_chunks->{$chr}})) {


    # print transformed gene models to chunked GFF file
    print GFFCHUNK $ref_chunk_genes->{$chr}{$chunk};

    # print sequence to chunked FASTA file

    $total_chunks++ 
  }
}

close(FNACHUNK);
close(GFFCHUNK);

print "# chunked GFF file: $chunkgfffile\n";
print "# chunked FASTA file: $chunkfnafile\n";

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

  my ($gff_file, $maxdist) = @_;

  my ($chr, $start, $end, $chunk_start, $chunk_end, $gff_line);
  my ($num_chunk, $dist, $offset, $prev_end) = (1, 0, 0, 0);
  my (%chunk_genes, %chunk, @chrs);

  open(GFF,"<",$gff_file) || 
    die "# ERROR(chunk_GFF): cannot read $gff_file\n";
  while(<GFF>){
    my @gff = split(/\t/,$_);
    ($chr, $start, $end) = @gff[0,3,4];

    next if($gff[2] eq 'chromosome');

    push(@chrs, $chr) if not grep(/^$chr$/,@chrs);;

    if($gff[2] eq 'gene') {

      # start new chunk if previous gene too far
      if($prev_end > 0) {
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

      #$chr .= ".chunk$num_chunk";
      #print "$chr|$chunk{$chr}{$num_chunk}{'start'} $_\n";

      # save end coord for next iteration
      $prev_end = $end; 
    } 

    # transform coords relative to current chunk
    #print "# $chr $num_chunk $chunk{$chr}{$num_chunk}{'offset'}\n";  
      
    $gff[0] = "$chr.chunk$num_chunk";
    $gff[3] -= $chunk{$chr}{$num_chunk}{'offset'};
    $gff[4] -= $chunk{$chr}{$num_chunk}{'offset'};
    $gff_line = join("\t",@gff); 
    
    #print "$num_chunk $gff_line" if($gff[2] eq 'gene');

    $chunk_genes{$chr}{$num_chunk} .= $gff_line;
  }
  close(GFF);

  return (\@chrs, \%chunk_genes, \%chunk);
}

