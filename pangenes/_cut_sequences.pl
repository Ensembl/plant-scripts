#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

# Takes a GFF & FASTA files and produces FASTA files with 
# CDS nucl & pep sequences of the 1st transcript found
# Note: also creates a FASTA index file (.fai)
#
# Uses external software: gffread [https://f1000research.com/articles/9-304/v2]

# Copyright [2021-22] 
# EMBL-European Bioinformatics Institute & Estacion Experimental de Aula Dei-CSIC

# perl _cut_sequences.pl -sp oryza_sativa -fa Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
#   -gf Oryza_sativa.IRGSP-1.0.51.gff3

my $GFFREADEXE = 'gffread'; # v0.12.7

my ( $help, $nored, $gffreadpath, $sp1, $fasta1) = (0, 0);
my ( $minlen, $gff1, $patchgff1, $tname, $outpath ) = (0, '', '');

GetOptions(
  "help|?"       => \$help,
  "sp|species=s" => \$sp1,
  "fa|fasta=s"   => \$fasta1,
  "gf|gff=s"     => \$gff1,
  "pt|patch=s"   => \$patchgff1,
  "l|minlen=i"   => \$minlen,
  "nr|n"         => \$nored,
  "p|path=s"     => \$gffreadpath,
  "o|outpath=s"  => \$outpath
) || help_message();

sub help_message {
  print "\nusage: $0 [options]\n\n"
    . "-sp binomial/trinomial species name (required, example: -sp oryza_sativa, used to name outfiles)\n"
    . "-fa genome FASTA filename           (required, example: -fa oryza_sativa.fna)\n"
    . "-gf GFF filename                    (required, example: -gf oryza_sativa.RAPDB.gff)\n"
    . "-pt GFF filename with patched       (optional, example: -pt oryza_sativa.RAPDB.patch.gff)\n"
    . "-l  min length (bp) of features     (optional, example: -l 100)\n"
    . "-nr remove redundancy in seq names  (optional, ie 'gene:ONIVA01G00100')\n"
    . "-p  path to gffread binary          (optional, default: $GFFREADEXE)\n"
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

if($patchgff1 && !-e $patchgff1){
  print "# ERROR: please make sure patch GFF file exists\n";
  exit(-2);
}

if(!$gffreadpath){
  $gffreadpath = $GFFREADEXE;
}

if($minlen < 1){ 
  $minlen = 0 
}

print "\n# $0 -sp $sp1 -fa $fasta1 -gf $gff1 -pt $patchgff1 " .
  "-l $minlen -nr $nored -path $gffreadpath\n\n";

# set output filenames
my $cdnafile = "$sp1.cdna.fna";
my $cdsfile  = "$sp1.cds.fna";
my $pepfile  = "$sp1.cds.faa";
if($patchgff1){
  $cdnafile = "$sp1.patch.cdna.fna";
  $cdsfile  = "$sp1.patch.cds.fna";
  $pepfile  = "$sp1.patch.cds.faa";
}
if($outpath) {
  $cdnafile = "$outpath/$cdnafile";
  $cdsfile  = "$outpath/$cdsfile";
  $pepfile  = "$outpath/$pepfile";
}

# only bother if not empty
if(-s $patchgff1) {
  my $patched_gff_filename = patch_gff($gff1, $patchgff1);
  $gff1 = $patched_gff_filename;
}

my ($ref_names, $ref_coords) = parse_genes($gff1);

my $num_cdna = parse_gffread($gffreadpath,$fasta1,$gff1,$cdnafile,
                 'cdna',$minlen,$nored,$ref_names,$ref_coords);
my $num_cds  = parse_gffread($gffreadpath,$fasta1,$gff1,$cdsfile,
                 'cds',$minlen,$nored,$ref_names,$ref_coords);
my $num_pep  = parse_gffread($gffreadpath,$fasta1,$gff1,$pepfile,
                 'pep',$minlen,$nored,$ref_names,$ref_coords);

if(scalar(keys(%$ref_names)) == 0) {
  die "# ERROR: cannot parse Parent IDs of mRNA/transcripts, please check GFF format ($gff1)\n";	
}

if($num_cdna) {
  print "# $cdnafile n=$num_cdna\n";
} else {
  die "# ERROR: cannot extract cDNA sequences, please check GFF format and/or chr names ($gff1)\n";
}

if($num_cds) {
  print "# $cdsfile n=$num_cds\n"
} else {
  die "# ERROR: cannot extract CDS sequences, please check GFF format and/or chr names ($gff1)\n";
}

print "# $pepfile n=$num_pep\n";

unlink($patched_gff_filename);

###############################

# Runs gffread, parses its stdout and saves output in FASTA file.
# Returns number of sequences printed out.
sub parse_gffread {

  my ($gffreadexe,$fasta_file,$gff_file,$outfile,
    $seqtype,$minlen,$remove_red,$ref_tr2gene,$ref_tr2coords) = @_;

  my ($params,$mrnaid,$geneid,$coords);

  if($seqtype eq 'cds'){
    $params = '-x - ';
  } elsif($seqtype eq 'pep'){
    $params = '-y - ';
  } else {
    $params = '-w - '; # cDNA, default
  }

  if($minlen > 0) {
    $params .= " -l $minlen ";
  }

  my $num_seqs = 0;
  open(OUT,">",$outfile) || 
    die "# ERROR(parse_gffread): cannot create $outfile\n";

  open(GFFREAD,"$gffreadexe $params -g $fasta_file $gff_file |") ||
    die "# ERROR(parse_gffread): cannot run $gffreadexe\n"; 

  while(<GFFREAD>){
    if(/^>(\S+)/){
      $mrnaid = $1;
      $geneid = $ref_tr2gene->{$mrnaid} || '';
      $coords = $ref_tr2coords->{$mrnaid} || '';

      # remove redundant bits
      if($remove_red){
        $mrnaid =~ s/transcript://;
        $geneid =~ s/gene://;
      }

      print OUT ">$mrnaid $geneid $coords [$sp1]\n";
      $num_seqs++;
    } else {
      print OUT;
    }
  }
  close(GFFREAD);

  close(OUT);

  # do not remove index (used later)
  # unlink($fasta_file.'.fai'); 

  return $num_seqs
}

# Reads in GFF file to parse gene names as parent IDs of transcripts.
# Returns: 
# i) ref to hash mapping transcript ID -> gene ID
# ii) ref to hash mapping transcript ID -> gene coords
sub parse_genes {

  my ($gff_file) = @_;

  my ($mrnaid,$geneid,$coord,%names,%coords);

  open(GFF,"<",$gff_file) || die "# ERROR(parse_genenames): cannot read $gff_file\n";
  while(<GFF>){
    my @F = split(/\t/,$_);

    next if(scalar(@F)<9 || ($F[2] ne "mRNA" && $F[2] ne "transcript"));

    # take only genes where ID can be parsed
    if($F[8] =~ /ID=([^;]+).*?Parent=([^;]+)/){ 
 
      $mrnaid = $1;
      $geneid = $2;
      chomp $geneid;
 
      $coord = "$F[0]:$F[3]-$F[4]($F[6])";
      $names{$mrnaid} = $geneid;
      $coords{$mrnaid} = $coord;
    }
  }
  close(GFF);

  return (\%names,\%coords);
}

# Takes two GFF files (original and patch) and produces a new GFF file
# that includes patched gene models. Patched models match the original
# ones by means of 'old_locus_tag' tags in gene features. Two params:
# i) GFF filename
# ii) GFF filename with selected gene model patches
# Returns:
# i) output GFF filename 
# ii) integer with number of patched gene models
sub patch_gff {

  my ($gff_file, $patchfile) = @_;

  my $patched_gff_filename = $gff_file . '.patched';
  my ($total_patched,$patch, $gene_id) = (0);
  my (%depr_gene_id);

  # read in GFF patches
  open(PATCH,'<',$patchfile) ||
    die "# ERROR(patch_gff): cannot read $patchfile\n";
  while(<PATCH>) {
    $patch .= $_;

    my @gffdata = split(/\t/,$_);
    if($gffdata[2] && $gffdata[2] eq 'gene') {
      $total_patched++;
      while($gffdata[8] =~ m/old_locus_tag=([^;]+)/g) {
        $gene_id = $1;
        $depr_gene_id{ $gene_id } = 1;
      }
    }
  }
  close(PATCH);

  printf("# total patched: %d deprecated: %d\n\n",
    $total_patched,scalar(keys(%depr_gene_id)));

  # read in original GFF, apply patches and save output
  open(PATCHED,">",$patched_gff_filename) ||
    die "# ERROR(patch_gff): cannot create $patched_gff_filename\n";

  my $geneOK = 1;
  open(GFF,'<',$gff_file) ||
    die "# ERROR(patch_gff): cannot read $gff_file\n";
  while(<GFF>) {

    my @gffdata = split(/\t/,$_);
    if($gffdata[2] && $gffdata[2] eq 'gene') {
      if($gffdata[8] =~ m/ID=([^;]+)/) {
        $gene_id = $1; 
        if($depr_gene_id{ $gene_id }){
          $geneOK = 0; 
        } else { $geneOK = 1 }
      } else { $geneOK = 1 }

      print PATCHED if($geneOK);

    } else {
      print PATCHED if($geneOK)
    }
  }
  close(GFF);

  # apply patches 
  print PATCHED $patch;
  close(PATCHED);

  return $patched_gff_filename;
}
