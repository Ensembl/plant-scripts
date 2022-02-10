#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename qw(basename dirname);

$|=1;

# Takes two FASTA files with genome sequences and 2 matching GFF files with 
# annotated gene models.
# Produces a TSV file with pairs of collinear genes in a format compatible 
# with Ensembl Compara

# Copyright [2021-22] 
# EMBL-European Bioinformatics Institute & Estacion Experimental de Aula Dei-CSIC

# Uses external software:
# minimap2 [https://academic.oup.com/bioinformatics/article/34/18/3094/4994778]
# bedtools [https://academic.oup.com/bioinformatics/article/26/6/841/244688]
# samtools [https://academic.oup.com/bioinformatics/article/25/16/2078/204688]
# wfmash   [https://github.com/ekg/wfmash]

# perl get_collinear_genes.pl -sp1 oryza_sativa \
#  -fa1 Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
#  -gf1 Oryza_sativa.IRGSP-1.0.51.gff3 -sp2 oryza_nivara -fa2 \
#  Oryza_nivara.Oryza_nivara_v1.0.dna.toplevel.fa \
#  -gf2 Oryza_nivara.Oryza_nivara_v1.0.51.gff3 -r

# collinear | Osativa vs Onivara | Athaliana vs Ahalleri
# 2.17      |     24502          |     10637
# 2.22      |     18770          |      8231

my $MINIMAP2EXE = 'minimap2'; # 2.17-r941 more senstive across species than v2.22
my $MINIMAPTYPE = '-x asm20'; # https://github.com/lh3/minimap2/issues/225
my $MINIMAPPARS = "--secondary=no --cs $MINIMAPTYPE ".
                  "-r1k,5k"; # https://github.com/lh3/minimap2/issues/813, 2949 -> 2956
my $WFMASHEXE   = 'wfmash';                  # v0.7.0
my $WFMASHPARS  = '-p 95 -s 3000';
my $BEDTOOLSEXE = 'bedtools';                # v2.30.0
my $BEDINTSCPAR = '-wo -f XXX -F XXX -e';    # XXX to be replaced with [0-1]
my $SAMTOOLSEXE = 'samtools';

# might be useful to compute sequence identity for ANI
#my $BLASTNEXE   = 'blastn';   # 2.2.30+

my $THREADS      = 4;
my $SORTLIMITRAM = "500M";                               # buffer size
my $SORTBIN      = "sort --buffer-size=$SORTLIMITRAM";
$ENV{'LC_ALL'} = 'POSIX';
my $GZIPBIN      = "gzip";

my $MINMASKLEN   = 100_000;   # mask longer (intergenic, repetitive) fragments with -H
my $GENEMARGIN   = 5000;      # do not mask gene margins 
my $DUMMYSCORE   = 9999;

# while parsing PAF
my $MINQUAL    = 50;          # works well with both mapping algorithms
my $MINALNLEN  = 100;
my $SAMESTRAND = 1;
my $MINOVERLAP = 0.50;
my $VERBOSE    = 0;           # values > 1

# Work out gene names from transcripts':
# 1) remove suffix after . or -
# 2) t's in transcript name to be replaced with g' in gene names
# Example: 10t0100300.1 -> Os10g0100300
my $TRANSCRIPT2GENE = 0;

my ( $help, $do_sequence_check, $reuse, $noheader, $repetitive) = (0, 0, 0, 0, 0);
my ($dowfmash, $split_chr_regex, $tmpdir ) = ( 0, '', '' );
my ( $sp1, $fasta1, $gff1, $sp2, $fasta2, $gff2, $index_fasta1 ) = 
  ( '', '', '', '', '', '', '');
my ( $chr, $chrfasta1, $chrfasta2, $splitPAF, $ref_chr_pairs, $cmd );
my ( $indexonly, $minoverlap, $qual, $alg, $outfilename ) =
  ( 0, $MINOVERLAP, $MINQUAL, 'minimap2' );
my ( $minimap_path, $wfmash_path, $bedtools_path, $samtools_path, $threads ) =
  ( $MINIMAP2EXE, $WFMASHEXE, $BEDTOOLSEXE, $SAMTOOLSEXE, $THREADS );

GetOptions(
    "help|?"         => \$help,
    "sp1|species1=s" => \$sp1,
    "fa1|fasta1=s"   => \$fasta1,
    "gf1|gff1=s"     => \$gff1,
    "sp2|species2=s" => \$sp2,
    "fa2|fasta2=s"   => \$fasta2,
    "gf2|gff2=s"     => \$gff2,
    "out|outfile=s"  => \$outfilename,
    "ovl|overlap=f"  => \$minoverlap,
    "q|quality=i"    => \$qual,
    "s|split=s"      => \$split_chr_regex,
    "c|check"        => \$do_sequence_check,
    "r|reuse"        => \$reuse,
    "i|index"        => \$indexonly,
    "wf|wfmash"      => \$dowfmash,
    "M|minipath=s"   => \$minimap_path,
    "W|wfpath=s"     => \$wfmash_path,
    "B|btpath=s"     => \$bedtools_path,
    "S|stpath=s"     => \$samtools_path,
    "T|tmpath=s"     => \$tmpdir,
    "t|threads=i"    => \$threads,
    "H|highrep"      => \$repetitive,
    "a|add"          => \$noheader
) || help_message();

sub help_message {
    print "\nusage: $0 [options]\n\n"
      . "-sp1 binomial/trinomial species name    (required, example: -sp1 oryza_sativa)\n"
      . "-fa1 genome FASTA [.gz] filename        (required, example: -fa1 oryza_sativa.fna, use bgzip with -wf)\n"
      . "-gf1 GFF [.gz] filename                 (required, example: -gf1 oryza_sativa.RAPDB.gff)\n"
      . "-sp2 binomial/trinomial species name    (required, example: -sp2 oryza_nivara_OGE)\n"
      . "-fa2 genome FASTA [.gz] filename        (required, example: -fa2 oryza_nivara.fna)\n"
      . "-gf2 GFF [.gz] filename                 (required, example: -gf2 oryza_nivara.OGE.gff)\n"
      . "-out output filename (TSV format)       (optional, by default built from input, example: -out rice.tsv)\n"
      . "-ovl min overlap of genes               (optional, default: -ovl $MINOVERLAP)\n"
      . '-s   split genome in chrs               (optional, requires regex to match chr names ie: -s \'^\d+$\')'. "\n"
      . "-wf  use wfmash aligner                 (optional, requires samtools ; by default minimap2 is used)\n"
      . "-q   min mapping quality, minimap2 only (optional, default: -q $MINQUAL)\n"
      . "-M   path to minimap2 binary            (optional, default: -M $MINIMAP2EXE)\n"
      . "-W   path to wfmash binary              (optional, default: -W $WFMASHEXE)\n"
      . "-B   path to bedtools binary            (optional, default: -B $BEDTOOLSEXE)\n"
      . "-S   path to samtools binary            (optional, default: -S $SAMTOOLSEXE)\n"
      . "-T   path for temporary files           (optional, default current folder)\n"
      . "-t   CPU threads to use                 (optional, default: -t $THREADS)\n"
      . "-H   highly repetitive genome           (optional, masks intergenes >= $MINMASKLEN & tweaks minimap2)\n"
#. "-c   check sequences of collinear genes (optional)\n"
      . "-add concat TSV output with no header   (optional, example: -add, requires -out)\n"
      . "-r   re-use previous results & index    (optional)\n"
      . "-i   make index & tmp files, dont align (optional, requires -sp1 -fa1 -gf1)\n";
}

if ( $help 
    || (!$sp1 || !$fasta1 || !$gff1 || !$sp2 || !$fasta2 || !$gff2) ) {
    help_message();
    exit(0);
}

if ( !-s $fasta1 || !-s $gff1 || !-s $fasta2 || !-s $gff2 ) {
    print "# ERROR: please make sure all input files exist\n";
    exit(0);
}

if ( !$indexonly && $sp1 eq $sp2 ) {
    print "# ERROR: please make sure -sp1 and sp2 are different\n";
    exit(0);
}

if ( $minoverlap && ( $minoverlap < 0 || $minoverlap > 1 ) ) {
    print "# ERROR: option -ovl requires values [0,1]\n";
    exit(0);
}

if ( $qual && ( $qual < 0 || $qual > 255 ) ) {
    print "# ERROR: option -q requires values [0,255]\n";
    exit(0);
}

# add actual value to dummy values in param string
$BEDINTSCPAR =~ s/XXX/$minoverlap/g;

if ( $noheader && !$outfilename ) {
    print "# ERROR: option -add requires -out filename to concat results\n";
    exit(0);
}

# check algorithm, reset min mapping quality if wfmash
# https://github.com/ekg/wfmash/issues/96
if ($dowfmash) {
    $alg = 'wfmash';
    $qual = 1;
} elsif($repetitive) {
    #see https://github.com/lh3/minimap2/issues/813
    $MINIMAPPARS .= " -f100 ";
}

if($tmpdir ne '' && $tmpdir !~ /\/$/){
    $tmpdir .= '/'
}

# set default outfile
if ( !$outfilename ) {
    $outfilename = ucfirst($alg) . ".homologies.$sp1.$sp2.overlap$minoverlap.tsv";

    if ($split_chr_regex ne '') {
        $outfilename = ucfirst($alg) . ".homologies.$sp1.$sp2.overlap$minoverlap.split.tsv";
    }
}

print "\n# $0 -sp1 $sp1 -fa1 $fasta1 -gf1 $gff1 "
  . "-sp2 $sp2 -fa2 $fasta2 -gf2 $gff2 -out $outfilename -a $noheader "
  . "-ovl $minoverlap -q $qual -wf $dowfmash -c $do_sequence_check "
  . "-s '$split_chr_regex' -M $minimap_path -W $wfmash_path -B $bedtools_path "
  . "-T $tmpdir -t $threads -i $indexonly -r $reuse -H $repetitive\n\n";

# check binaries
if(`$bedtools_path` !~ 'sage') {
	print "# ERROR: cannot find binary file $bedtools_path , exit\n";	
	exit(-1)
} 

if ($dowfmash) {
    if(`$wfmash_path` !~ 'OPTIONS') {
        print "# ERROR: cannot find binary file $wfmash_path , exit\n";
        exit(-2)
    } elsif (`$samtools_path 2>&1` !~ 'Usage') {
        print "# ERROR: cannot find binary file $samtools_path , exit\n";
        exit(-3)
    }
} else {
    if(`$minimap_path 2>&1` !~ 'Usage') {
        print "# ERROR: cannot find binary file $minimap_path , exit\n";
        exit(-4)
    }
}

print "# mapping parameters:\n";
if ($dowfmash) {
    print "# \$WFMASHPARS: $WFMASHPARS\n\n";
}
else {
    print "# \$MINIMAPPARS: $MINIMAPPARS\n\n";
}


## 1) Parse GFFs and produce BED files with gene coords

my $geneBEDfile1 = $tmpdir . "_$sp1.gene.bed";
my $geneBEDfile2 = $tmpdir . "_$sp2.gene.bed";

if($reuse && -s $geneBEDfile1 && check_BED_format($geneBEDfile1)) {
    print "# re-using $geneBEDfile1\n";
} else {
    my ( $num_genes1, $mean_gene_len1 ) = parse_genes_GFF( $gff1, $geneBEDfile1 );
    printf( "# %d genes parsed in %s mean length=%d\n",
        $num_genes1, $gff1, $mean_gene_len1 );
}

if(!$indexonly) {
    if($reuse && -s $geneBEDfile2 && check_BED_format($geneBEDfile2)) {
        print "# re-using $geneBEDfile2\n";
    } else {
        my ( $num_genes2, $mean_gene_len2 ) = parse_genes_GFF( $gff2, $geneBEDfile2 );
        printf( "# %d genes parsed in %s mean length=%d\n",
            $num_genes2, $gff2, $mean_gene_len2 );
    }
}

# 1.1) if requested (long, repetitive genomes) mask long intergenes of $sp1
if($repetitive) {
    my $masked_fasta1 = $tmpdir . "_$sp1.mask.fna";
    my $fasta_length1 = $tmpdir . "_$sp1.tsv";
    my $maskBEDfile1  = $tmpdir . "_$sp1.mask.bed";

    if($reuse && -s $maskBEDfile1) {
        print "# re-using $maskBEDfile1\n";
    } else {
        my ($total_masked1, $median_length1) = mask_intergenic_regions(
            $fasta1,$geneBEDfile1,
            $masked_fasta1, $fasta_length1, $maskBEDfile1,
            $MINMASKLEN,$GENEMARGIN,$bedtools_path);

        printf("# %s bases masked=%d median intergene length=%d\n",
            $sp1, $total_masked1, $median_length1 );
    }
    $fasta1 = $masked_fasta1;

    my $masked_fasta2 = $tmpdir . "_$sp2.mask.fna";
    my $fasta_length2 = $tmpdir . "_$sp2.tsv";
    my $maskBEDfile2  = $tmpdir . "_$sp2.mask.bed";

    if(!$indexonly) {
        if($reuse && -s $maskBEDfile2) {
            print "# re-using $maskBEDfile2\n";
        } else {
            my ($total_masked2, $median_length2) = mask_intergenic_regions(
                $fasta2,$geneBEDfile2,
                $masked_fasta2, $fasta_length2, $maskBEDfile2,
                $MINMASKLEN,$GENEMARGIN,$bedtools_path);

            printf("# %s bases masked=%d median intergene=%d\n",
                $sp1, $total_masked2, $median_length2 );
        }
        $fasta2 = $masked_fasta2;
    }
}

## 2) align genome1 vs genome2 (WGA)
# Note: masking not recommended, see https://github.com/lh3/minimap2/issues/654

# split genome assemblies if required, 1/chr plus 'unplaced'
# Note: reduces complexity (good for large/polyploid genomes) but misses translocations
if ($split_chr_regex ne '') {
    print "\n# splitting sequences with regex\n";
    $ref_chr_pairs = 
        split_genome_sequences_per_chr($tmpdir, $fasta1, $fasta2, 
            $split_chr_regex, $indexonly, $reuse);
}
else {
    # single reference by default
    $ref_chr_pairs->{'all'} = [ $fasta1, $fasta2 ]
}

my $PAFfile = $tmpdir . "_$sp2.$sp1.$alg.paf";
if ($split_chr_regex ne '') {
    $PAFfile = $tmpdir . "_$sp2.$sp1.$alg.split.paf";
} 
if($repetitive) {
    $PAFfile =~ s/\.paf$/.highrep.paf/; 
}

if($indexonly) {
    print "# indexing genome with $alg\n\n";
} else {
    print "# computing pairwise genome alignment with $alg\n\n";
}

if ( $reuse && -s $PAFfile ) {
    print "# re-using $PAFfile\n";
}
else {

    unlink($PAFfile);

    my (@sorted_chrs,@WGAoutfiles);

    # sort chromosomes
    foreach $chr (keys(%$ref_chr_pairs)) {
        if($chr ne 'unplaced' && $chr ne 'all') {
            push(@sorted_chrs,$chr);
        }
    }

    @sorted_chrs = sort @sorted_chrs; # {$a<=>$b} not always numeric
    if($ref_chr_pairs->{'unplaced'}){ push(@sorted_chrs,'unplaced') }
    elsif($ref_chr_pairs->{'all'}){ push(@sorted_chrs,'all') }

    foreach $chr (@sorted_chrs) {

        $chrfasta1 = $ref_chr_pairs->{$chr}[0];
        $chrfasta2 = $ref_chr_pairs->{$chr}[1];
        $splitPAF =  $tmpdir . "_$sp2.$sp1.$alg.split.$chr.paf";
        if($repetitive) {
            $splitPAF =~ s/\.paf$/.highrep.paf/;
        }#print ">$chr $chrfasta1 $chrfasta2 $splitPAF\n"; 

        if ( $reuse && -s $splitPAF ) {
            print "# re-using $splitPAF\n";
            push(@WGAoutfiles, $splitPAF);
            next;
        }

        # make empty PAF file when chr files are empty
        if(!$indexonly && (!-s $chrfasta1 || !-s $chrfasta2)) {
            open(EMPTYPAF, ">", $splitPAF); 
            close(EMPTYPAF);
            push(@WGAoutfiles, $splitPAF);
            next;
        }

        if ($dowfmash) {

            # usually created in _cut_sequences.pl
            $index_fasta1 = dirname($chrfasta1)."/".basename($chrfasta1).".fai";

            if ( $reuse && -s $index_fasta1 ) {
                if ($split_chr_regex ne '') {
                    printf("# re-using $index_fasta1, make sure same regex was used\n");
                } else {
                    print "# re-using $index_fasta1\n";
                }
            } else {
                $cmd = "$samtools_path faidx $chrfasta1 -o $index_fasta1 2>&1";				
                system($cmd);
                if ( $? != 0 ) {
                    die "# ERROR: failed running samtools ($cmd)\n";
                }
                elsif ( !-s $index_fasta1 ) {
                    die "# ERROR: failed generating $index_fasta1 file ($cmd)\n";
                }
            } 

            next if($indexonly);

            $cmd = "$wfmash_path $WFMASHPARS -t $threads $chrfasta1 $chrfasta2 > $splitPAF";
            system($cmd);
            sleep(2);
            if ( $? != 0 ) {
                die "# ERROR: failed running wfmash (probably ran out of memory, $cmd)\n";
            } elsif ( !-e $splitPAF ) {
                die "# ERROR: failed generating $splitPAF file ($cmd)\n";
            } else {
                push(@WGAoutfiles, $splitPAF);
            }
        }
        else {    # default minimap2 index & alignment

            $index_fasta1 = $tmpdir . "_$sp1.$chr.mmi";
            if($repetitive) {
                $index_fasta1 =~ s/\.mmi$/.highrep.mmi/;
            }

            if ( $reuse && -s $index_fasta1 ) {
                if ($split_chr_regex ne '') {
                    printf("# re-using $index_fasta1, make sure same regex was used\n");
                } else {
                    print "# re-using $index_fasta1\n";
                }
            } else {
                $cmd = "$minimap_path $MINIMAPTYPE -t $threads -d $index_fasta1 $chrfasta1";
                system($cmd);
                if ( $? != 0 ) {
                    die "# ERROR: failed running minimap2 (probably ran out of memory, $cmd)\n";
                }
                elsif ( !-s $index_fasta1 ) {
                    die "# ERROR: failed generating $index_fasta1 file ($cmd)\n";
                }
            }

            next if($indexonly);

            $cmd = "$minimap_path $MINIMAPPARS -t $threads $index_fasta1 $chrfasta2 -o $splitPAF";
            system($cmd);
            sleep(2);
            if ( $? != 0 ) {
                die "# ERROR: failed running minimap2 (probably ran out of memory)\n";
            }
            elsif ( !-e $splitPAF ) {
                die "# ERROR: failed generating $splitPAF file ($cmd)\n";
            }
            else {
                push(@WGAoutfiles, $splitPAF);
            }
        }
    }

    if(!$indexonly) {
        # merge (split) PAF files
        $cmd = "cat "; # assumes cat is available
        foreach $splitPAF (@WGAoutfiles) {
            $cmd .= "$splitPAF ";
        }
        $cmd .= " > $PAFfile";

        system("$cmd");
        if( $? != 0 ) {
            die "# ERROR: failed merging split PAF files ($cmd)\n";
        } else {
            print("# WGA finished\n\n");
            unlink(@WGAoutfiles);
        }
    }
}

if($indexonly) {
    print "# indexes created, terminate\n";
    exit(0);
}

## 3) produce BED-like file of sp2-to-sp1 coords 10 columns
my ($cigar,@tmpBEDfiles);

my $wgaBEDfile    = $tmpdir . "_$sp2.$sp1.$alg.bed";

# used to find matching sp2 segments for unpaired sp1 genes
my $wgaBEDfilerev = $tmpdir . "_$sp1.$sp2.$alg.bed";

open( PAF, "<", $PAFfile )    || die "# ERROR: cannot read $PAFfile\n";

open( BED, ">", $wgaBEDfile ) || die "# ERROR: cannot create $wgaBEDfile\n";
open( BEDREV, ">", $wgaBEDfilerev ) || die "# ERROR: cannot create $wgaBEDfilerev\n";

while (<PAF>) {
    #Note: $F[4] in PAF conveys whether query & ref are on the same strand or not
    #Pt	1345257 118 7420 - 1 42845077 35836986 35837288	300 302	60 ... cs:Z::216*tc:69*ga:15
    my @F = split( /\t/, $_ );
    
    print BED
      "$F[0]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[7]\t$F[8]\t$F[9]\t$F[11]\t$F[$#F]";

    # reverse CIGAR after reversing the alignment
    # Note: as we only care about coords, don't do anything with substitutions
    $cigar = $F[$#F];
    $cigar =~ tr/ID\+\-/DI\-\+/;

    print BEDREV
      "$F[5]\t$F[7]\t$F[8]\t$F[4]\t$F[0]\t$F[2]\t$F[3]\t$F[9]\t$F[11]\t$cigar"
}

close(BED);
close(BEDREV);

close(PAF);

push(@tmpBEDfiles, $wgaBEDfile, $wgaBEDfilerev);

## 4) intersect gene positions with WGA, sort by gene > cDNA ovlp > genomic matches

my $sp2wgaBEDfile = $tmpdir . "_$sp2.gene.$sp1.$alg.intersect.overlap$minoverlap.bed";
my $sp2wgaBEDfile_sorted = $tmpdir . "_$sp2.gene.$sp1.$alg.intersect.overlap$minoverlap.sort.bed";

my $sp1wgaBEDfile = $tmpdir . "_$sp1.gene.$sp2.$alg.intersect.overlap$minoverlap.bed";
my $sp1wgaBEDfile_sorted = $tmpdir . "_$sp1.gene.$sp2.$alg.intersect.overlap$minoverlap.sort.bed";

$cmd = "$bedtools_path intersect -a $geneBEDfile2 -b $wgaBEDfile " .
         "$BEDINTSCPAR > $sp2wgaBEDfile";

system("$cmd");
sleep(2); #latency issues
if ( $? != 0 ) {
    die "# ERROR: failed running bedtools (WGA, $cmd)\n";
}
elsif ( !-s $sp2wgaBEDfile ) {
    die "# ERROR: failed generating $sp2wgaBEDfile file ($cmd)\n";
}

$cmd = "$SORTBIN -k4,4 -k5,5nr -k14,14nr $sp2wgaBEDfile > $sp2wgaBEDfile_sorted";
system("$cmd");
if ( $? != 0 ) {
    die "# ERROR: failed sorting (WGA, $cmd)\n";
}
elsif ( !-s $sp2wgaBEDfile_sorted ) {
    die "# ERROR: failed generating $sp2wgaBEDfile_sorted file ($cmd)\n";
}

# now with reversed WGA alignment
$cmd = "$bedtools_path intersect -a $geneBEDfile1 -b $wgaBEDfilerev " .
         "$BEDINTSCPAR > $sp1wgaBEDfile";

system("$cmd");
sleep(2); #latency issues
if ( $? != 0 ) {
    die "# ERROR: failed running bedtools (WGArev, $cmd)\n";
}
elsif ( !-s $sp1wgaBEDfile ) {
    die "# ERROR: failed generating $sp1wgaBEDfile file ($cmd)\n";
}

$cmd = "$SORTBIN -k4,4 -k5,5nr -k14,14nr $sp1wgaBEDfile > $sp1wgaBEDfile_sorted";
system("$cmd");
if ( $? != 0 ) {
    die "# ERROR: failed sorting (WGArev, $cmd)\n";
}
elsif ( !-s $sp1wgaBEDfile_sorted ) {
    die "# ERROR: failed generating $sp2wgaBEDfile_sorted file ($cmd)\n";
}

push(@tmpBEDfiles, $sp2wgaBEDfile, $sp2wgaBEDfile_sorted);
push(@tmpBEDfiles, $sp1wgaBEDfile, $sp1wgaBEDfile_sorted);

# compute coords of mapped genes 
my $geneBEDfile2mapped = $tmpdir . "_$sp2.$sp1.$alg.gene.mapped.bed";
my $geneBEDfile1mapped = $tmpdir . "_$sp1.$sp2.$alg.gene.mapped.bed";

my ( $ref_matched, $ref_unmatched ) =
  query2ref_coords( $sp2wgaBEDfile, $geneBEDfile2mapped,
    $qual, $MINALNLEN, $SAMESTRAND, $VERBOSE );

printf( "# %d genes mapped in %s (%d unmapped)\n",
    scalar(@$ref_matched), $geneBEDfile2mapped, scalar(@$ref_unmatched) );

if ( scalar(@$ref_matched) == 0 ) {
    die "# ERROR: failed mapping $sp2 genes in WGA alignment";
} 

# now with reversed WGA alignment, to find matching sp2 segments for unpaired sp1 genes

my ( $ref_matched1, $ref_unmatched1 ) =
  query2ref_coords( $sp1wgaBEDfile, $geneBEDfile1mapped,
    $qual, $MINALNLEN, $SAMESTRAND, $VERBOSE );

printf( "# %d genes mapped in %s (reverse, %d unmapped)\n",
    scalar(@$ref_matched1), $geneBEDfile1mapped, scalar(@$ref_unmatched1) );

if ( scalar(@$ref_matched1) == 0 ) {
    die "# ERROR: failed mapping $sp1 genes in WGA alignment";
}

## 5) produce list of pairs of collinear genes & genomic segments

my $gene_intersectBEDfile = 
  $tmpdir . "_$sp1.$sp2.$alg.gene.intersect.overlap$minoverlap.bed";

my $segment_intersectBEDfile =
  $tmpdir . "_$sp1.$sp2.$alg.segment.intersect.overlap$minoverlap.bed";

my $segment_intersectBEDfile1 =
  $tmpdir . "_$sp2.$sp1.$alg.segment.intersect.overlap$minoverlap.bed";

my $intersectBEDfile_sorted =
  $tmpdir . "_$sp1.$sp2.$alg.intersect.overlap$minoverlap.sort.bed";

$cmd = "$bedtools_path intersect -a $geneBEDfile1 -b $geneBEDfile2mapped " .
         "$BEDINTSCPAR -s > $gene_intersectBEDfile";

system($cmd);
sleep(2);
if ( $? != 0 ) {
    die "# ERROR: failed running bedtools (genes, $cmd)\n";
}
elsif ( !-s $gene_intersectBEDfile ) {
    die "# ERROR: failed generating $gene_intersectBEDfile file ($cmd)\n";
}

# now add collinear gene-genomic segments pairs
my $num_segments = genes_mapped2segments( $geneBEDfile2, $geneBEDfile2mapped, 
	$gene_intersectBEDfile, $segment_intersectBEDfile );

my $num_segments1 = genes_mapped2segments( $geneBEDfile1, $geneBEDfile1mapped,
        $gene_intersectBEDfile, $segment_intersectBEDfile1, 1 );

$cmd = "$SORTBIN -k1,1 -k2,2n $gene_intersectBEDfile $segment_intersectBEDfile ".
           "$segment_intersectBEDfile1 > $intersectBEDfile_sorted";
system($cmd);
if ( $? != 0 ) {
    die "# ERROR: failed sorting (genes, $cmd)\n";
}
elsif ( !-s $intersectBEDfile_sorted ) {
    die "# ERROR: failed generating $intersectBEDfile_sorted file ($cmd)\n";
}

push(@tmpBEDfiles, $segment_intersectBEDfile, $segment_intersectBEDfile1);
#push(@tmpBEDfiles, $gene_intersectBEDfile, $intersectBEDfile_sorted);

my $num_pairs = bed2compara( $intersectBEDfile_sorted, $outfilename, $sp1, $sp2,
    $noheader, $TRANSCRIPT2GENE );

printf( "# %d collinear gene/segment pairs\n", $num_pairs );

if($num_pairs > 0 && -s $outfilename) {
    print "# TSV file: $outfilename\n";

    unlink(@tmpBEDfiles);
} 



#################################

# Takes string with BED filename (produced with sub parse_genes_GFF)
# and returns 1 if format is OK
sub check_BED_format {

    my ( $bedfile ) = @_;
    my $formatOK = 1;    

    open(BED, "<", $bedfile) || die "# ERROR(check_BED_format): cannot read $bedfile\n";
    while(<BED>) {
         if($_ !~ /^[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+/) {
             $formatOK = 0;
             last;     
         }
    }
    close(BED);

    return $formatOK
}


# Takes i) input GFF filename ii) output BED filename and
# returns i) number and ii) mean length of genes parsed.
# Uses global $DUMMYSCORE as dummy scores, as opposed to cases where
# genes/transcript are actually mapped to genome
sub parse_genes_GFF {

    my ( $gff_file, $outBEDfile ) = @_;
    my $num_genes  = 0;
    my $genelength = 0;
    my ( $geneid, $magic );

    # check input file format and open it accordingly
    open( INFILE, $gff_file )
      || die "# ERROR(parse_genes_GFF): cannot read $gff_file, exit\n";
    sysread( INFILE, $magic, 2 );
    close(INFILE);

    if ( $gff_file =~ /\.gz$/ || $magic eq "\x1f\x8b" ) {  # GZIP compressed input{
        if ( !open( GFF, "$GZIPBIN -dc $gff_file |" ) ) {
            die "# ERROR(parse_genes_GFF): cannot read GZIP compressed $gff_file $!\n"
              . "# please check gzip is installed\n";
        }
    } else {
        open( GFF, "<", $gff_file )
          || die "# ERROR(parse_genes_GFF): cannot read $gff_file\n";
    }

    open( BED, ">", $outBEDfile )
      || die "# ERROR(parse_genes_GFF): cannot create $outBEDfile\n";

    while (<GFF>) {
        my @F = split( /\t/, $_ );

        next if ( scalar(@F) < 9 || $F[2] ne "gene" );

        # take only genes where ID can be parsed
        if ( $F[8] =~ /ID=([^;]+)/ ) {

            $geneid = $1;
            chomp($geneid);
            #$geneid =~ s/gene://;    # remove redundant bits

            printf( BED "%s\t%d\t%d\t%s\t%s\t%s\n",
                $F[0],
                $F[3] - 1,           # 0-based
                $F[4],
                $geneid,
                $DUMMYSCORE,
                $F[6]
            );

            $num_genes++;
            $genelength += ( $F[4] - $F[3] ) + 1;
        }
    }

    close(BED);
    close(GFF);

    return ( $num_genes, sprintf( "%1.0f", $genelength / $num_genes ) );
}

# Takes 2 strings:
# 1) name of FASTA file
# 2) regex to match chromosome names, applied to first non-blank token (\S+)
# Returns ref to hash with chr and/or 'unplaced' as keys and sequences as value
sub read_FASTA_regex2hash {
    my ($fastafile,$regex) = @_;

    my ($magic,$seqname,$seq);
    my %fasta;

    # check input file format and open it accordingly
    open(INFILE,$fastafile) || die "# ERROR(read_FASTA_regex2hash): cannot read $fastafile, exit\n";
    sysread(INFILE,$magic,2);
    close(INFILE);

    if($fastafile =~ /\.gz$/ || $magic eq "\x1f\x8b") { 
        if(!open(FASTA,"gzip -dc $fastafile |")) {
            die "# ERROR(read_FASTA_regex2hash): cannot read GZIP compressed $fastafile $!\n"
                ."# please check gzip is installed\n";
        }
    } elsif($fastafile =~ /\.bz2$/ || $magic eq "BZ") {
        if(!open(FASTA,"bzip2 -dc $fastafile |")){
            die "# ERROR(read_FASTA_regex2hash): cannot read BZIP2 compressed $fastafile $!\n"
                ."# please check bzip2 is installed\n";
        }
    } else{ open(FASTA,"<$fastafile") || 
        die "# ERROR(read_FASTA_regex2hash): cannot read $fastafile $!\n"; }

    while (<FASTA>) {
        next if(/^\s*$/ || /^#/);
        if(/^>/) { # header
            if(/^>\s*(\S+)/) { 
                $seqname = $1;
                if($seqname !~ m/^$regex$/) { $seqname = 'unplaced_'.$seqname }
            }
        } else {
            $fasta{$seqname} .= $_;
        }
    }

    close(FASTA);

    return \%fasta;
}

# Takes 6 params:
# 1) path to write files to
# 2) name of FASTA file (ref)
# 3) name of FASTA file (query)
# 4) regex to match chromosome names
# 5) indexing job (boolean)
# 6) reuse (boolean)
# Returns ref to hash with chr and/or 'unplaced' as keys and two FASTA files as value (ref, query)
# Note: 'unplaced' might hold genuine unplaced sequences but also non-shared chr names
# Note: creates FASTA files with prefix _ 
# Note: if doindex is true only creates files for $fastafile1
sub split_genome_sequences_per_chr {

    my ($path,$fastafile1,$fastafile2,$regex,$doindex,$reuse) = @_;

    my ($chr,$chrfasta1,$chrfasta2,$summaryfile);
    my (%shared_chrs,%chr_pairs);

    my $ref_fasta1 = read_FASTA_regex2hash($fastafile1,$regex);
    my $ref_fasta2 = read_FASTA_regex2hash($fastafile2,$regex);

    # check chr names found in both files
    foreach $chr (keys(%$ref_fasta1)){
        if(defined($ref_fasta2->{$chr}) && $chr !~ 'unplaced'){
            $shared_chrs{$chr} = 1;  
        } 
    } #print join ',', keys(%$ref_fasta1); 

    # write chr-specific FASTA files
    foreach $chr (keys(%shared_chrs)) { #print ">$chr\n";
        $chrfasta1 = $path . "_".basename($fastafile1).".$chr.fna";
        $chrfasta2 = $path . "_".basename($fastafile2).".$chr.fna";

        if(!$reuse || !-e $chrfasta1) {
            open( CHRFASTA1, ">$chrfasta1") ||
                die "# ERROR(split_genome_sequences_per_chr): cannot write $chrfasta1\n";
            print CHRFASTA1 ">$chr\n";
            print CHRFASTA1 $ref_fasta1->{$chr};
            close(CHRFASTA1);
        } 

        if(!$doindex && (!$reuse || !-e $chrfasta2)) {
            open( CHRFASTA2, ">$chrfasta2") ||
                die "# ERROR(split_genome_sequences_per_chr): cannot write $chrfasta2\n";
            print CHRFASTA2 ">$chr\n";
            print CHRFASTA2 $ref_fasta2->{$chr};
            close(CHRFASTA2);
        }

        $chr_pairs{$chr} = [$chrfasta1,$chrfasta2];
    }

    # write unplaced and/or not-shared chr names
    $chrfasta1 = $path . "_".basename($fastafile1).".unplaced.fna";
    $chrfasta2 = $path . "_".basename($fastafile2).".unplaced.fna";

    if(!$reuse || !-e $chrfasta1) {
        open( CHRFASTA1, ">$chrfasta1") ||
            die "# ERROR(split_genome_sequences_per_chr): cannot write $chrfasta1\n";
        foreach $chr (keys(%$ref_fasta1)){ 
            next if(defined($shared_chrs{$chr}));
            print CHRFASTA1 ">$chr\n"; 
            print CHRFASTA1 $ref_fasta1->{$chr};   
        }
        close(CHRFASTA1);
    }

    if(!$doindex && (!$reuse || !-e $chrfasta2)) {
        open( CHRFASTA2, ">$chrfasta2") ||
            die "# ERROR(split_genome_sequences_per_chr): cannot write $chrfasta2\n";
        foreach $chr (keys(%$ref_fasta2)){
            next if(defined($shared_chrs{$chr}));
            print CHRFASTA2 ">$chr\n";
            print CHRFASTA2 $ref_fasta2->{$chr};
        }
        close(CHRFASTA2);
    }

    $chr_pairs{'unplaced'} = [$chrfasta1,$chrfasta2];

    return \%chr_pairs;
}

# Takes 7 params:
# 1) name of FASTA file, must be uncompressed
# 2) name of gene BED file 
# 3) name of FASTA output file 
# 4) name of TSV output file with chr length
# 5) name of BED output file with long intergenic masked regions
# 6) min length of regions to mask (integer)
# 7) do not mask with this number of bases next to genes (integer)
# Returns number of masked bases 
sub mask_intergenic_regions {

    my ( $fasta, $geneBED, $out_fasta, $out_len, $out_bed, 
         $minlen, $margin, $bedtoolsEXE ) = @_;

    my $total_masked = 0;
    my ($chr,$seq,$cmd,$start,$end,$len);
    my @intergenes;

    # compute chr lengths and save to TSV file
    my $ref_fasta = read_FASTA_regex2hash($fasta,'\S+');
    
    open(TSV,">$out_len") || 
        die "# ERROR(mask_intergenic_regions): cannot write $out_len\n"; 

    # note sort matches sort -k1,1 with LC_ALL=POSIX
    foreach $chr (sort {$a cmp $b} keys(%$ref_fasta)) {
        $seq = $ref_fasta->{$chr};
        $seq =~ s/\n//g;
        printf(TSV "%s\t%d\n",$chr,length($seq));
    }

    close(TSV);

    # compute complement of gene space, adding $margin bases either side,
    # and write output BED file

    open(BED,">$out_bed") ||
        die "# ERROR(mask_intergenic_regions): cannot write $out_bed\n";

    $cmd = "$bedtoolsEXE complement -i $geneBED -g $out_len";
    open(BEDTOOLS,"$cmd |") || 
        die "# ERROR(mask_intergenic_regions): cannot run $cmd\n";

    while(<BEDTOOLS>) {
        my @bedata = split(/\t/);
        $start = $bedata[1]+$margin; 
        $end   = $bedata[2]-$margin; 
        $len   = $end-$start; 
        if($len > $minlen){ 
            print BED "$bedata[0]\t$start\t$end\n";
            $total_masked += $len;  
        } 
        push(@intergenes,$len) if($len > 0);
    }

    close(BEDTOOLS);
    close(BED); 
   
    # mask remaining intervals
    $cmd = "$bedtoolsEXE maskfasta -fi $fasta -bed $out_bed -fo $out_fasta";
    system($cmd);
    sleep(2);
    if ( $? != 0 ) {
        die "# ERROR(mask_intergenic_regions): failed running bedtools (intergenes, $cmd)\n";
    }
    elsif ( !-s $out_fasta ) {    
        die "# ERROR(mask_intergenic_regions): failed generating $out_fasta file ($cmd)\n";
    }

    return ($total_masked, calc_median(\@intergenes));
}

# Takes i) input BED intersect filename ii) output BED filename.
# Parses sorted BED intersect -wo output and writes to BED file
# features (cDNA/transcripts) mapped on reference genome. Note:
# features might be unsorted.
# Returns i) ref to list of matched genes and ii) ref to list of unmatched genes
# Note: able to parse cs::Z (minimap2) and cg::Z (wfmash) strings
# Note: takes first match of each cDNA only
# example input:
# 1 4848 20752 ONIVA01G00010 9999     + 1 3331 33993       + 6 26020714 26051403 29819 60 cs:Z::303*ag:30*ga... 15904
# 1 104921 116326 ONIVA01G00100 9999  + 1 103118 152580    + 1 1132 47408 45875 60 cs:Z::70*tc:...              11405
# Chr1 2903 10817 LOC_Os01g01010 9999 + Chr1 1000 10053455 + 1 1000 10053455 10052455 60 cs:Z::10052455          7914
# Chr1 2903 10817 LOC_Os01g01010 9999 + Chr1 896 705000    + 1 949 705000 704051 41 cg:Z:51=53I704000=           7914
# <--             (c)DNA/gene       --> <- (q)uery genome -> <-- (r)eference genome                         -->  ovlp
sub query2ref_coords {

    my ( $infile, $outfile, $minqual, $minalnlen, $samestrand, $verbose ) = @_;

    my ( $cchr, $cstart, $cend, $cname, $cmatch, $cstrand );
    my ( $qchr, $qstart, $qend, $cigartype, $bedline );
    my ( $WGAstrand, $rchr,      $rstart,       $rend );
    my ( $rmatch,    $rmapqual,  $SAMPAFtag,    $overlap, $done, $strand );
    my ( $SAMqcoord, $SAMrcoord, $feat,         $coordr );
    my ( $deltaq,    $deltar,    $start_deltar, $end_deltar );
    my ( %ref_coords, @matched, @unmatched, @segments );

    open( BED, "<", $infile )
      || die "# ERROR(query2ref_coords): cannot read $infile\n";
    while (<BED>) {

#cDNA format
#1 98773 99875 ONIVA01G00080.1 258 + 1 98032 101175 - 6 27346427 27348975 2375 60 cs:Z::29-ggt.. 15904
#GFF/gene format
#1  4848 20752 ONIVA01G00010 9999  + 1 3331 33993 + 6 26020714 26051403 29819 60 cs:Z::303*ag.. 15904
#1 104921 116326 ONIVA01G00100 9999  + 1 103118 152580    + 1 1132 47408 45875 60 cs:Z::70*tc:... 11405

        (
            $cchr,      $cstart, $cend,   $cname,  $cmatch,
            $cstrand,   $qchr,   $qstart, $qend,   $WGAstrand,
            $rchr,      $rstart, $rend,   $rmatch, $rmapqual,
            $SAMPAFtag, $overlap
        ) = split( /\t/, $_ );

        # skip mappings where query and ref segment on different strands if required
        next if ( $WGAstrand eq '-' && $samestrand == 1 );

        # skip poor query-to-ref WGA scores
        next if ( $rmapqual < $minqual );

        # take 1st mapping only
        next if ( defined( $ref_coords{$cname} ) );

        #next if($cname ne 'ONIVA01G00100'); # debug
        #next if($cname ne 'ONIVA01G00150');
        #next if($cname ne 'LOC_Os01g01010');
        #next if($cname ne 'LOC_Os11g34300');

        # make sure ref chr is taken
        $ref_coords{$cname}{'chr'} = $rchr;

        # correct offset for ref assembly
        # Note: this requires parsing SAM/PAF tag
        # https://github.com/lh3/minimap2#paftools
        # cs:Z::303*ag:32*ga:27+ctattcta*ag*ca:20*ag:3*ga:18*tc*ga:76-tc

        # default start coords (end in - strand),
        # in case 3'cDNA not included in WGA segment
        $SAMqcoord = $qstart;
        if ( $WGAstrand eq '+' ) {
            $SAMrcoord = $rstart;
            $ref_coords{$cname}{'start'} = $rstart;
        }
        else {
            $SAMrcoord = $rend;
            $ref_coords{$cname}{'end'} = $rend;
        }

        print "# $SAMqcoord $SAMrcoord $cname ".scalar(@matched)."\n"
            if ( $verbose > 1 );

        # check CIGAR type
        if ( $SAMPAFtag =~ m/(c\w):Z:{1,2}(\S+)/ ) {
            ( $cigartype, $SAMPAFtag ) = ( $1, $2 );
        }
        else {
            print "# ERROR(query2ref_coords): unsupported CIGAR string $SAMPAFtag\n";
            return ( \@matched, \@unmatched );
        }

        # split CIGAR string into individual feature tags
        if ( $cigartype eq 'cs' ) {
            @segments = split( /:/, $SAMPAFtag );
        }
        else {
            while ( $SAMPAFtag =~ m/(\d+[MIDNSHP=X])/g ) {
                push( @segments, $1 );
            }
        }

        # loop along features updating coords
        $done = 0;
        foreach $feat (@segments) {

            ( $deltaq, $deltar, $coordr ) =
              _parseCIGARfeature( $feat, $SAMqcoord, $SAMrcoord );

            # check if current position in alignment matches cDNA/gene coords

            # start coords (end in - strand)
            if ( $SAMqcoord < $cstart
                && ( $SAMqcoord + $deltaq ) >= $cstart ) {

                # refine delta to match exactly the start (end for strand -)
                ( $deltaq, $deltar, $coordr ) =
                  _parseCIGARfeature( $feat, $SAMqcoord, $SAMrcoord, $cstart );

                $start_deltar = -1;
                if ( $coordr > -1 ) {
                    $start_deltar = $coordr - $SAMrcoord;
                }

                print ">$deltaq $deltar $start_deltar $cstart\n"
                  if ( $verbose > 1 );

                if ( $start_deltar < 0 || $start_deltar > $deltar ) {
                    $start_deltar = $deltar;
                }

                if ( $WGAstrand eq '+' ) {
                    $ref_coords{$cname}{'start'} = $SAMrcoord + $start_deltar;
                }
                else {
                    $ref_coords{$cname}{'end'} = $SAMrcoord + $start_deltar;
                }
            }

            # end coords (start in -strand), actually copied out of loop
            if ( $SAMqcoord < $cend && ( $SAMqcoord + $deltaq ) >= $cend ) {

                # refine delta to match exactly the end (start for strand -)
                ( $deltaq, $deltar, $coordr ) =
                  _parseCIGARfeature( $feat, $SAMqcoord, $SAMrcoord, $cend );

                $end_deltar = -1;
                if ( $coordr > -1 ) {
                    $end_deltar = $coordr - $SAMrcoord;
                }

                print "<$deltaq $deltar $end_deltar $cend\n"
                  if ( $verbose > 1 );

                if ( $end_deltar < 0 || $end_deltar > $deltar ) {
                    $end_deltar = $deltar;
                }

                $done = 1;
            }

            # update coords with deltas
            $SAMqcoord += $deltaq;

            if ($done) {
                $deltar = $end_deltar;
            }

            if ( $WGAstrand eq '+' ) {
                $SAMrcoord += $deltar;
            }
            else {
                $SAMrcoord -= $deltar;
            }

            print "$SAMqcoord $SAMrcoord $feat $deltaq $deltar\n"
              if ( $verbose > 1 );

            last if ($done);
        }

        if ( $WGAstrand eq '+' ) {
            $ref_coords{$cname}{'end'} = $SAMrcoord;
        }
        else {
            $ref_coords{$cname}{'start'} = $SAMrcoord;
        }

        # work out strand of mapped cDNA/transcript/gene
        $strand = $cstrand;    # original strand
        if ( $WGAstrand eq '-' ) {
            if ( $cstrand eq '+' ) {
                $strand = '-';
            }
            else {
                $strand = '+';
            }
        }

        # print coords in ref frame only if segment is long enough
        $overlap = 1 + $ref_coords{$cname}{'end'} - $ref_coords{$cname}{'start'};
        $bedline = sprintf("%s\t%d\t%d\t%s\t%d\t%s\n",
            $ref_coords{$cname}{'chr'},
            $ref_coords{$cname}{'start'},
            $ref_coords{$cname}{'end'},
            $cname,
            $overlap,
            $strand);

        if ( $overlap >= $minalnlen ) {
            push(@matched, $bedline);
        }
        else {
            # store also short segments, useful to call unmapped regions
            push( @unmatched, $bedline );
        }
    }

    close(BED);

    # printed unsorted BED records
    open( OUTBED, ">", $outfile )
      || die "# ERROR(query2ref_coords): cannot create $outfile\n";

    foreach $feat (@matched) {
        print OUTBED $feat;
    }

    close(OUTBED);    

    return ( \@matched, \@unmatched );
}

# Takes 3 scalars:
# 1) CS/CG tag feature
# 2) start query coordinate
# 3) start reference coordinate
# 4 optional) target query coordinate
# Returns 3 scalars:
# 1) the increment (delta) in query (q) coords after adding the new CS feature
# 2) the increment in reference (r) coords
# 3) the reference coord matching to the optional target query coord, -1 by default
# Note: CS tags are CIGAR-like strings produced by minimap2 with flag --cs;
# CS tags start with preffix 'cs:Z::' and can be :-split into features.
# Note: CG tags are CIGAR-like strings produced by wfmash, with preffix 'cg::Z:'
sub _parseCIGARfeature {

    my ( $feat, $Qstart, $Rstart, $opt_query_coord ) = @_;

    my ( $deltaq, $deltar, $totsnps, $rcoord ) = ( 0, 0, 0, -1 );
    my ( $delta, $offset, $res, $qcoord, $op );

    # example CS string (coords are 0-based):
    # cs:Z::303*ag:30*ga:32*ga:27+ctattcta*ag*ca:20*ag:3*ga:18*tc:39*ga:76-tc
    # first features: 303*ag 30*ga .. 27+ctattcta .. 76-tc

    # example CG string:
    # cg:Z:25I235=1X57=1X33=1X20=1X47=
    # first features: 25I 235= 1X ..

    # CG type, see https://samtools.github.io/hts-specs/SAMv1.pdf
    if ( $feat =~ m/(\d+)([MIDNSHP=X])/ ) {

        ( $delta, $op ) = ( $1, $2 );

        if ( $op eq 'M' || $op eq '=' || $op eq 'X' ) {
            $deltar += $delta;
            $deltaq += $delta;

            # look up query coord (optional)
            if ( defined($opt_query_coord) ) {
                $offset = $opt_query_coord - $Qstart;
                if ( $offset <= $delta ) {
                    $rcoord = $Rstart + $offset;
                }
            }

        }
        elsif ( $op eq 'I' ) {

            # look up query coord (optional)
            if ( defined($opt_query_coord) && $rcoord == -1 ) {
                $offset = $opt_query_coord - $Qstart;
                if ( $offset <= $delta + $deltaq ) {
                    $rcoord = $Rstart + $deltar;
                }
            }

            $deltaq += $delta;

        }
        elsif ( $op eq 'D' ) {

            # look up query coord (optional)
            if ( defined($opt_query_coord) && $rcoord == -1 ) {
                $offset = $opt_query_coord - $Qstart;
                if ( $offset <= $deltaq ) {
                    $rcoord = $Rstart + $deltar + $offset;
                }
            }

            $deltar += $delta;

        }
        else {
            print "# ERROR(_parseCIGARfeature): unsupported CIGAR operation $feat\n";
        }

    }
    else {    # CS type

        # identical segment
        if ( $feat =~ m/^(\d+)/ ) {
            $delta = $1;
            $deltar += $delta;
            $deltaq += $delta;

            # look up query coord (optional)
            if ( defined($opt_query_coord) ) {
                $offset = $opt_query_coord - $Qstart;
                if ( $offset <= $delta ) {
                    $rcoord = $Rstart + $offset;
                }
            }
        }

        # insertion (+) / deletion (-)
        if ( $feat =~ m/([\+\-])([A-Za-z]+)/ ) {
            $delta = length($2);

            if ( $1 eq '+' ) {

                # look up query coord (optional)
                if ( defined($opt_query_coord) && $rcoord == -1 ) {
                    $offset = $opt_query_coord - $Qstart;
                    if ( $offset <= $delta + $deltaq ) {
                        $rcoord = $Rstart + $deltar;
                    }
                }

                $deltaq += $delta;

            }
            else {

                # look up query coord (optional)
                if ( defined($opt_query_coord) && $rcoord == -1 ) {
                    $offset = $opt_query_coord - $Qstart;
                    if ( $offset <= $deltaq ) {
                        $rcoord = $Rstart + $deltar + $offset;
                    }
                }

                $deltar += $delta;
            }
        }

        # SNPs
        while ( $feat =~ m/\*\w{2}/g ) {

            # look up query coord (optional)
            if ( defined($opt_query_coord) && $rcoord == -1 ) {
                $qcoord = $Qstart + $totsnps + $deltaq;
                if ( $qcoord == $opt_query_coord && $rcoord == -1 ) {
                    $rcoord = $Rstart + $totsnps + $deltar;
                }
            }

            $totsnps++;
        }
        $deltar += $totsnps;
        $deltaq += $totsnps;
    }

    return ( $deltaq, $deltar, $rcoord );
}

# Takes four file paths and optionally a boolean:
# i)   BED filename with gene model coordinates of species2 (6 cols)
# ii)  BED filename with species1 genes mapped on species2 space (6 cols)
# iii) BED filename with sp1,sp2 pairs of collinear genes (13 columns)
# iv)  BED output filename (13 columns)
# v)   boolean, columns from sp2 should be printed first
# Returns number of collinear genomic segments found
sub genes_mapped2segments {

    my ($geneBEDfile, $mappedBEDfile, $pairedBEDfile, $outBEDfile, $invert) = @_;

    my $num_segments = 0;
    my ($coords,$geneid,$len,$strand);
    my (@genes,%orig_coords,%paired);
    
    # find out which genes are paired based on WGA (sp2)
    open(PAIRBED,"<",$pairedBEDfile)
        || die "# ERROR(gene2segments): cannot read $pairedBEDfile\n";

    while(<PAIRBED>) {
        #1 30219   36442   gene:BGIOSGA002569 9999 + 1 29700 39038 gene:ONIVA01G00100  9339  +  6223
        my @data = split(/\t/,$_);
        $paired{$data[9]} = 1; # sp2
    }
    close(PAIRBED);

    # read gene model coordinates (sp2)
    open(GENEBED,"<",$geneBEDfile)
        || die "# ERROR(gene2segments): cannot read $geneBEDfile\n";

    while(<GENEBED>) {
        #1    4847    20752   gene:ONIVA01G00010      9999    +
        if(/^\S+\t\d+\t\d+\t(\S+)\t/) {
            $geneid = $1;
            next if($paired{$geneid});
            chomp;
            $orig_coords{$geneid} = $_;
        }
    }
    close(GENEBED); 

    # find sp2 genes that have mapped coords on sp1 but are not paired, print to outfile
    open(OUTBED,">",$outBEDfile) 
        || die "# ERROR(gene2segments): cannot create $outBEDfile\n";

    open(MAPBED,"<",$mappedBEDfile)
        || die "# ERROR(gene2segments): cannot read $mappedBEDfile\n";

    while(<MAPBED>) {
        #1       29700   39038   gene:ONIVA01G00100      9339    +
        if(/^(\S+\t\d+\t\d+)\t(\S+)\t(\d+)\t([+-])/) {
            ($coords, $geneid, $len, $strand) = ($1, $2, $3, $4);
            
            next if($paired{$geneid});

            # actually print to BED coordinates of genes mapped to (unannotated) genomic segments
            #1 217360 222398 segment 5039 + 1 155040 165322 gene:ONIVA01G00180 9999 + 9999

            if($invert) {
                print OUTBED "$orig_coords{$geneid}\t$coords\tsegment\t$len\t$strand\t9999\n";
            } else {
                print OUTBED "$coords\tsegment\t$len\t$strand\t$orig_coords{$geneid}\t9999\n";
            }
           
            $num_segments++;
        }
    }
    close(MAPBED);

    close(OUTBED);

    return $num_segments;
}


# Takes i) input BED intersect filename ii) output TSV filename and
# returns iii) number of collinear gene pairs & iv) number of collinear segments
# Parses bedtools intersect file and produces Ensembl Compara-like TSV file.
# Example input with explicit strand:
# 1 2983 10815 Os01g0100100 9999 + 1 2890 12378 ONIVA01G00100 9489 + 7832
# Expected TSV output:
# gene_stable_id protein_stable_id species overlap homology_type homology_gene_stable_id
# homology_protein_stable_id homology_species overlap dn ds goc_score wga_coverage
# is_high_confidence coordinates
# Note: homology_type can be 'ortholog_collinear' or 'segment_collinear'
sub bed2compara {

    my ( $infile, $TSVfile, $sp1, $sp2, $noheader, $workout_gene_names ) = @_;

    my ( $gene1, $gene2, $coords1, $coords2, $coords, $homoltype );
    my $num_pairs = 0;

    # parse input and produce TSV output
    open( BEDINT, "<", $infile )
      || die "# ERROR(bed2compara): cannot read $infile\n";

    if ($noheader) {
        open( TSV, ">>", $TSVfile )
          || die "# ERROR(bed2compara): cannot re-open $TSVfile\n";
    }
    else {
        open( TSV, ">", $TSVfile )
          || die "# ERROR(bed2compara): cannot create $TSVfile\n";

        print TSV
          "gene_stable_id\tprotein_stable_id\tspecies\toverlap\thomology_type\thomology_gene_stable_id\t"
          . "homology_protein_stable_id\thomology_species\toverlap\tdn\tds\tgoc_score\twga_coverage\t"
          . "is_high_confidence\tcoordinates\n";
    }

    while (<BEDINT>) {
        my @data = split( /\t/, $_ );

        # concat genome/graph coords
        $coords1 = "$data[0]:$data[1]-$data[2]";
        $coords2 = "$data[6]:$data[7]-$data[8]";
        $coords = "$coords1;$coords2";

        # format gene names
        $gene1 = $data[3];
        $gene2 = $data[9];

        if ($workout_gene_names) {
            $gene1 =~ s/[\.-]\d+$//;
            $gene2 =~ s/[\.-]\d+$//;
            $gene1 =~ s/t/g/;
            $gene2 =~ s/t/g/;
        }

        # check homology type
        if($data[3] eq 'segment') {
            $homoltype = 'segment_collinear';
            $gene1 = $coords1;
        } elsif($data[9] eq 'segment') {
            $homoltype = 'segment_collinear';
            $gene2 = $coords2;
        } else {
            $homoltype = 'ortholog_collinear'
        }

        printf( TSV
            "%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%d\tNULL\tNULL\tNULL\t%1.2f\t%d\t%s\n",
            $gene1,
            $data[3],
            $sp1,
            $data[12],
            $homoltype,
            $gene2,
            $data[9],
            $sp2,
            $data[12],
            100,    # max WGA score
            1,      # high confidence
            $coords
        );

        $num_pairs++;
    }

    close(TSV);
    close(BEDINT);

    return $num_pairs;
}

sub calc_median {
    
    my ($dataref) = @_;

    my $mid = int(scalar(@$dataref)/2);
    my @sorted = sort {$a<=>$b} (@$dataref);

    if(scalar(@sorted) % 2) { 
        return $sorted[ $mid ] 
    }
    else { 
        return sprintf("%1.0f",($sorted[$mid-1] + $sorted[$mid])/2) 
    }
}
