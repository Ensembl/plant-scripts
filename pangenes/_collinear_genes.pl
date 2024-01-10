#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename qw(basename dirname);
use FindBin '$Bin';
use lib "$Bin/lib";
use pangeneTools qw(calc_median N50 read_FAI_regex2hash);

$|=1;

# Takes two FASTA files with genome sequences and 2 matching GFF files with 
# annotated gene models.
# Produces a TSV file with pairs of collinear genes in a format similar to
# Ensembl Compara's, genomic coordinates are 1-based

# Copyright [2021-24] 
# EMBL-European Bioinformatics Institute & Estacion Experimental de Aula Dei-CSIC

# Uses external software:
# minimap2 [https://academic.oup.com/bioinformatics/article/34/18/3094/4994778]
# bedtools [https://academic.oup.com/bioinformatics/article/26/6/841/244688]
# samtools [https://academic.oup.com/bioinformatics/article/25/16/2078/204688]
# GSAlign  [https://doi.org/10.1186/s12864-020-6569-1]
# wfmash   [https://github.com/ekg/wfmash]

# perl _collinear_genes.pl -sp1 oryza_sativa \
#  -fa1 Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
#  -gf1 Oryza_sativa.IRGSP-1.0.51.gff3 -sp2 oryza_nivara -fa2 \
#  Oryza_nivara.Oryza_nivara_v1.0.dna.toplevel.fa \
#  -gf2 Oryza_nivara.Oryza_nivara_v1.0.51.gff3 -r

# collinear | Osativa vs Onivara | Athaliana vs Ahalleri
# 2.17      |     24502          |     10637
# 2.22      |     18770          |      8231
# 2.24      |     25227          |         ?

#hardcode different minimap
#$ENV{'EXE_MINIMAP'} = '~/soft/minimap2-2.24_x64-linux/minimap2';

my $MINIMAP2EXE = 'minimap2';  
my $MINIMAPTYPE = '-x asm20'; # https://github.com/lh3/minimap2/issues/225
my $MINIMAPPARS = "--secondary=no --cs $MINIMAPTYPE " .
                  "-r1k,5k"; # https://github.com/lh3/minimap2/issues/813, 
                  # 2949 -> 2956 2.17
                  # 2956 -> 2951 , 25277 -> 25209 2.24 

my $WFMASHEXE   = $ENV{'EXE_WFMASH'} || 'wfmash';   # v0.8.1-25-g1344b9e
my $WFMASHPARS  = '-p 80 -s 1000';   # -s: median rice gene 2362, barley 1323
                                     # -p: ~ asm20
                                     # -p 90 -s 1000 -> 1652
                                     # -p 80 -s 1000 -> 1793
                                     # -p 80 -s 2000 -> 1528

my $GSALIGNPATH = './';
my $GSAINDXEXE  = 'bwt_index';
my $GSALIGNEXE  = 'GSAlign';
my $GSALIGNPARS = '-sen -no_vcf -fmt 1';

my $BEDTOOLSEXE = 'bedtools';  # v2.30.0
my $BEDINTSCPAR = '-wo -f XXX -F XXX -e';    # XXX to be replaced with [0-1]
my $SAMTOOLSEXE = 'samtools';

my $THREADS      = 4;
my $SORTBIN      = $ENV{'EXE_SORT'} || 'sort';
my $SORTPARS     = "--buffer-size=1G "; 
$ENV{'LC_ALL'}   = 'POSIX';
my $GZIPBIN      = $ENV{'EXE_GZIP'} || 'gzip';
my $BZIP2BIN     = $ENV{'EXE_BZIP2'} ||'bzip2';

my $MINMASKLEN   = 1000_000;  # mask longer (intergenic, repetitive) fragments with -H
my $GENEMARGIN   = 5000;      # do not mask gene margins 
my $DUMMYSCORE   = 9999;

# while parsing PAF
my $MINQUAL    = 50;          # works well with minimap2
my $MINALNLEN  = 100;         # min alignment length when transforming gene coords on WGA
my $MINOVERLAP = 0.50;
my $VERBOSE    = 0;           # values > 1

my ( $help, $do_sequence_check, $reuse, $noheader, $repetitive) = (0, 0, 0, 0, 0);
my ($dowfmash, $dogsalign, $patch, $split_chr_regex, $tmpdir ) = (0, 0, 0, '', '');
my ( $sp1, $fasta1, $gff1, $sp2, $fasta2, $gff2, $index_fasta1 ) = 
  ('', '', '', '', '', '', '');
my ( $fasta1orig, $fasta2orig ) = ('', ''); 
my ( $chr, $chrfasta1, $chrfasta2, $splitPAF, $ref_chr_pairs, $cmd, $gene );
my ( $indexonly, $no_inversions, $minoverlap, $qual, $alg, $outANIfile, $outfilename ) =
  ( 0, 0, $MINOVERLAP, $MINQUAL, 'minimap2', '' );
my ( $minimap_path, $wfmash_path, $gsalign_path, $bedtools_path, $samtools_path ) =
  ( $MINIMAP2EXE, $WFMASHEXE, $GSALIGNPATH, $BEDTOOLSEXE, $SAMTOOLSEXE );
my $threads = $THREADS;

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
    "p|patch"        => \$patch,
    "q|quality=i"    => \$qual,
    "s|split=s"      => \$split_chr_regex,
    "c|check"        => \$do_sequence_check,
    "r|reuse"        => \$reuse,
    "i|index"        => \$indexonly,
    "n|noinvs"       => \$no_inversions,	
    "wf|wfmash"      => \$dowfmash,
    "gs|gsalign"     => \$dogsalign,
    "M|minipath=s"   => \$minimap_path,
    "W|wfpath=s"     => \$wfmash_path,
    "G|gspath=s"     => \$gsalign_path,
    "B|btpath=s"     => \$bedtools_path,
    "S|stpath=s"     => \$samtools_path,
    "T|tmpath=s"     => \$tmpdir,
    "t|threads=i"    => \$threads,
    "H|highrep"      => \$repetitive,
    "A|ANI=s"        => \$outANIfile,
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
      . "-p   use patched gene models            (optional, forces recalculation of gene overlaps)\n"
      . "-ovl min overlap of genes               (optional, default: -ovl $MINOVERLAP)\n"
      . '-s   split genome in chrs               (optional, requires regex to match chr names ie: -s \'^\d+$\')'. "\n"
      . "-n   dont map genes in inversions       (optional, by default all genes are mapped on WGAs on both strands\n"	  
      . "-wf  use wfmash aligner                 (optional, requires samtools ; by default minimap2 is used)\n"
      . "-gs  use GSAlign aligner                (optional, by default minimap2 is used)\n"
      . "-q   min mapping quality, minimap2 only (optional, default: -q $MINQUAL)\n"
      . "-M   path to minimap2 binary            (optional, default: -M $MINIMAP2EXE)\n"
      . "-W   path to wfmash binary              (optional, default: -W $WFMASHEXE)\n"
      . "-G   path to GSAlign bin/               (optional, default: -G $GSALIGNPATH)\n"
      . "-B   path to bedtools binary            (optional, default: -B $BEDTOOLSEXE)\n"
      . "-S   path to samtools binary            (optional, default: -S $SAMTOOLSEXE)\n"
      . "-T   path for temporary files           (optional, default current folder)\n"
      . "-t   CPU threads to use                 (optional, default: -t $THREADS)\n"
      . "-H   highly repetitive genome           (optional, masks intergenes >= $MINMASKLEN & tweaks minimap2)\n"
      . "-A   output ANI filename                (optional, requires -gs, example: out.ani)\n"
#. "-c   check sequences of collinear genes (optional)\n"
      . "-add concat TSV output with no header   (optional, example: -add, requires -out)\n"
      . "-r   re-use previous results & index    (optional, partially overriden by -p)\n"
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
} elsif($dogsalign) {
    $alg = 'gsalign';
} elsif($repetitive) {
    #see https://github.com/lh3/minimap2/issues/813
    $MINIMAPPARS .= " -f100 ";

    # see also ideas in
    # https://github.com/lh3/minimap2/issues/354
}

if(!$dogsalign){ 
    $outANIfile = ''
}

if($tmpdir ne '' && $tmpdir !~ /\/$/){
    $tmpdir .= '/'
}

# set default outfile
if ( !$outfilename ) {
    $outfilename = ucfirst($alg) . ".homologies.$sp1.$sp2.overlap$minoverlap.tsv";
    if($patch) {
        $outfilename = ucfirst($alg) . ".homologies.$sp1.$sp2.overlap$minoverlap.patch.tsv";
    }

    if ($split_chr_regex ne '') {
        $outfilename = ucfirst($alg) . ".homologies.$sp1.$sp2.overlap$minoverlap.split.tsv";
        if($patch) {
            $outfilename = ucfirst($alg) . ".homologies.$sp1.$sp2.overlap$minoverlap.split.patch.tsv";
        }
    }
}

print "\n# $0 -sp1 $sp1 -fa1 $fasta1 -gf1 $gff1 "
  . "-sp2 $sp2 -fa2 $fasta2 -gf2 $gff2 -out $outfilename -p $patch -a $noheader "
  . "-ovl $minoverlap -q $qual -wf $dowfmash -gs $dogsalign -A $outANIfile -c $do_sequence_check "
  . "-s '$split_chr_regex' -M $minimap_path -W $wfmash_path -G $gsalign_path -B $bedtools_path "
  . "-T $tmpdir -t $threads -i $indexonly -r $reuse -H $repetitive -n $no_inversions\n\n";

# save names of original input FASTA files as $fasta1/$fasta2 might change with -H
$fasta1orig = $fasta1;
$fasta2orig = $fasta2;

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
} elsif($dogsalign) {

    $GSAINDXEXE = "$gsalign_path/$GSAINDXEXE";
    $GSALIGNEXE = "$gsalign_path/$GSALIGNEXE";  
    
    if(`$GSALIGNEXE 2>&1` !~ 'Usage') {
        print "# ERROR: cannot find binary file $GSALIGNEXE , exit\n";
        exit(-4)
    }
} else {
    if(!`$minimap_path 2>&1` || `$minimap_path 2>&1` !~ 'Usage') {
        print "# ERROR: cannot find binary file $minimap_path , exit\n";
        exit(-5)
    }
}

print "# mapping parameters:\n";
if($repetitive) {
  print "# \$MINMASKLEN: $MINMASKLEN\n";
}

if ($dowfmash) {
    print "# \$WFMASHPARS: $WFMASHPARS\n\n";
} elsif ($dogsalign) {
    print "# \$GSALIGNPARS: $GSALIGNPARS\n\n";
} else {
    print "# \$MINIMAPPARS: $MINIMAPPARS\n\n";
}


## 1) Parse GFFs and produce BED files with gene coords

my $geneBEDfile1 = $tmpdir . "_$sp1.gene.bed";
my $geneBEDfile2 = $tmpdir . "_$sp2.gene.bed";
if($patch) {
    $geneBEDfile1 = $tmpdir . "_$sp1.gene.patch.bed";
    $geneBEDfile2 = $tmpdir . "_$sp2.gene.patch.bed";
}

if($reuse && -s $geneBEDfile1 && check_BED_format($geneBEDfile1) && 
    (!$patch || !$indexonly) ) {
    print "# re-using $geneBEDfile1\n";
} else {
    my ( $num_genes1, $mean_gene_len1 ) = parse_genes_GFF( $gff1, $geneBEDfile1 );
    if(check_BED_format($geneBEDfile1)) {
        printf( "# %d genes parsed in %s mean length=%d\n",
            $num_genes1, $gff1, $mean_gene_len1 );
    } else {
        die "# ERROR: failed parsing genes from $gff1\n";    
    }
}

if(!$indexonly) {
 
    # $geneBEDfile2 will be re-used if possible even with $patch if !$indexonly,
    # as patched BED files are generated earlier with $indexonly

    if($reuse && -s $geneBEDfile2 && check_BED_format($geneBEDfile2)) {
        print "# re-using $geneBEDfile2\n";
    } else {
        my ( $num_genes2, $mean_gene_len2 ) = parse_genes_GFF( $gff2, $geneBEDfile2 );
        if(check_BED_format($geneBEDfile2)) {
            printf( "# %d genes parsed in %s mean length=%d\n",
                $num_genes2, $gff2, $mean_gene_len2 );
        } else {
            die "# ERROR: failed parsing genes from $gff2\n";
        }
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
        $fasta2orig = $fasta2;	
        $fasta2 = $masked_fasta2;
    }

} 

## 2) align genome1 vs genome2 (WGA)

# split genome assemblies if required, 1/chr plus 'unplaced'
# Note: reduces complexity (good for large/polyploid genomes) but misses translocations
if ($split_chr_regex ne '') {
    print "\n# splitting sequences with regex\n";

    # bedtools approach (7m to 4m in wheat), requires .fai index
    if(-e $fasta1.'.fai' && -e $fasta2.'.fai') {
        $ref_chr_pairs =
            split_genome_sequences_per_chr_bedtools($tmpdir, $fasta1, $fasta2,
                $split_chr_regex, $bedtools_path, $indexonly, $reuse);

    } else { 

        $ref_chr_pairs = 
            split_genome_sequences_per_chr($tmpdir, $fasta1, $fasta2, 
                $split_chr_regex, $indexonly, $reuse); 
    }
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
    unlink($outANIfile) if($dogsalign && $outANIfile);

    my (@sorted_chrs,@WGAoutfiles);

    # sort chromosomes
    foreach $chr (keys(%$ref_chr_pairs)) {
        if($chr ne 'unplaced' && $chr ne 'all') {
            push(@sorted_chrs,$chr);
        }
    }

    @sorted_chrs = sort @sorted_chrs; # {$a<=>$b} not always numeric
    if($ref_chr_pairs->{'unplaced'}) { 
        push(@sorted_chrs,'unplaced') 

    } elsif($ref_chr_pairs->{'all'}){ 
        push(@sorted_chrs,'all') 
    }

    foreach $chr (@sorted_chrs) {

        $chrfasta1 = $ref_chr_pairs->{$chr}[0];
        $chrfasta2 = $ref_chr_pairs->{$chr}[1];
        $splitPAF =  $tmpdir . "_$sp2.$sp1.$alg.split.$chr.paf";
        if($repetitive) {
            $splitPAF =~ s/\.paf$/.highrep.paf/;
        }

        if ( $reuse && -s $splitPAF ) {
            print "# re-using $splitPAF\n";
            push(@WGAoutfiles, $splitPAF);
            next;
        }

        # create empty PAF file when chr files are empty, and move to next chr
        if(!-s $chrfasta1 || !-s $chrfasta2) {
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
            print "# $cmd\n";
            system($cmd);
            sleep(2);
            if ( $? != 0 ) {
                die "# ERROR: failed running wfmash (probably ran out of memory, $cmd , $?)\n";
            } elsif ( !-e $splitPAF ) {
                die "# ERROR: failed generating $splitPAF file ($cmd)\n";
            } else {
                push(@WGAoutfiles, $splitPAF);
            }
        } elsif ($dogsalign) {

            my $preffix = $tmpdir . "_$sp1.$chr";
            my $splitMAFpreffix = $tmpdir . "_$sp2.$sp1.$alg.split.$chr";
            my $splitMAF = "$splitMAFpreffix.maf";
            $index_fasta1 = "$preffix.sa";

            if ( $reuse && -s $index_fasta1 ) {
                if ($split_chr_regex ne '') {
                    printf("# re-using $index_fasta1, make sure same regex was used\n");
                } else {
                    print "# re-using $index_fasta1\n";
                }
            } else {
                $cmd = "$GSAINDXEXE $chrfasta1 $preffix 2>&1";
                system($cmd);
                if ( $? != 0 ) {
                    die "# ERROR: failed running bwt_index ($cmd, $?)\n";
                }
                elsif ( !-s $index_fasta1 ) {
                    die "# ERROR: failed generating $index_fasta1 file ($cmd)\n";
                }
            }

            next if($indexonly);

            $cmd = "$GSALIGNEXE $GSALIGNPARS -t $threads -i $preffix -q $chrfasta2 ".
                "-o $splitMAFpreffix 2> $splitMAFpreffix.log";
            print "# $cmd\n";
            system($cmd);
            sleep(2);
            if ( $? != 0 ) {
                die "# ERROR: failed running GSAlign (probably ran out of memory, $cmd, $?)\n";
            } elsif ( !-e $splitMAF ) {
                die "# ERROR: failed generating $splitMAF file ($cmd)\n";
            } else {
 
                # print ANI estimate
                open(ANI,">>",$outANIfile) ||
                    die "# ERROR: cannot create $outANIfile\n";

                open(GSALOG,"<","$splitMAFpreffix.log") ||
                    die "# ERROR: cannot find previous log $splitMAFpreffix.log\n";
                while(<GSALOG>) {
                    if(/alignment length=(\d+)\)\s+ANI=([^\%]+)/) {
                        print ANI "$chr\t$2\t$1\n";
                        last;
                    }
                }
                close(GSALOG);

                close(ANI);

                # convert MAF to PAF alignment format
                my $num_align = simpleMAF2PAF($splitMAF,$splitPAF);
                if($num_align) {
                    print("# simpleMAF2PAF : %d alignments\n",$num_align); 
                    push(@WGAoutfiles, $splitPAF);
                    system("$GZIPBIN -f $splitMAF");
                } else {
                    die "# ERROR: failed converting $splitMAF file\n"; 
                }
            }

        } else {    # default minimap2 index & alignment

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
                print "# $cmd\n";
                system($cmd);
                if ( $? != 0 ) {
                    # see https://github.com/lh3/minimap2/issues/755
                    die "# ERROR: failed running minimap2 (probably ran out of memory, ".
                        "or chr too large, $cmd, $?)\n";
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
                die "# ERROR: failed running minimap2 (probably ran out of memory, $cmd, $?)\n";
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
my ($cigar,@tmpBEDfiles, @block_length);

my $wgaBEDfile    = $tmpdir . "_$sp2.$sp1.$alg.bed";

# make up also reverse alignment, 
# used to find matching sp2 segments for unpaired sp1 genes
my $wgaBEDfilerev = $tmpdir . "_$sp1.$sp2.$alg.rev.bed";

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
      "$F[5]\t$F[7]\t$F[8]\t$F[4]\t$F[0]\t$F[2]\t$F[3]\t$F[9]\t$F[11]\t$cigar";

    # record aligned block length
    push(@block_length, $F[10]);
}

close(BED);
close(BEDREV);

close(PAF);

push(@tmpBEDfiles, $wgaBEDfile, $wgaBEDfilerev);

printf("# WGA blocks: N50 %d median %1.0f\n",
    N50(\@block_length),
    calc_median(\@block_length));
@block_length = ();


## 4) intersect gene positions with WGA, sort by gene > cDNA ovlp > genomic matches

my $sp2wgaBEDfile = $tmpdir . "_$sp2.gene.$sp1.$alg.intersect.overlap$minoverlap.bed";
my $sp2wgaBEDfile_sorted = $tmpdir . "_$sp2.gene.$sp1.$alg.intersect.overlap$minoverlap.sort.bed";

my $sp1wgaBEDfile = $tmpdir . "_$sp1.gene.$sp2.$alg.intersect.overlap$minoverlap.rev.bed";
my $sp1wgaBEDfile_sorted = $tmpdir . "_$sp1.gene.$sp2.$alg.intersect.overlap$minoverlap.sort.rev.bed";

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

$cmd = "$SORTBIN $SORTPARS -k4,4 -k5,5nr -k14,14nr $sp2wgaBEDfile > $sp2wgaBEDfile_sorted";
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

$cmd = "$SORTBIN $SORTPARS -k4,4 -k5,5nr -k14,14nr $sp1wgaBEDfile > $sp1wgaBEDfile_sorted";
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
my $geneBEDfile1mapped = $tmpdir . "_$sp1.$sp2.$alg.gene.mapped.rev.bed";

my ( $ref_matched, $ref_unmatched, $perc_blocks_3genes ) =
  query2ref_coords( "$fasta1orig.fai", $sp2wgaBEDfile_sorted, $geneBEDfile2mapped,
    $qual, $MINALNLEN, $no_inversions, $VERBOSE ); 

printf( "# %d genes mapped (%1.1f%% in 3+blocks) in %s (%d unmapped)\n\n",
    scalar(@$ref_matched), $perc_blocks_3genes,
    $geneBEDfile2mapped, scalar(@$ref_unmatched) );

if ( scalar(@$ref_matched) == 0 ) {
    die "# ERROR: failed mapping $sp2 genes in WGA alignment";
} else {
    foreach $gene (@$ref_unmatched) {
        print "# unmapped: $gene\n"
    }
}


# now with reversed WGA alignment, to find matching sp2 segments for unpaired sp1 genes

my ( $ref_matched1, $ref_unmatched1, $perc_blocks_3genes1 ) =
  query2ref_coords( "$fasta2orig.fai", $sp1wgaBEDfile_sorted, $geneBEDfile1mapped,
    $qual, $MINALNLEN, $no_inversions, $VERBOSE );

printf( "# %d genes mapped (%1.1f%% in 3+blocks) in %s (reverse, %d unmapped)\n\n",
    scalar(@$ref_matched1), $perc_blocks_3genes1,
    $geneBEDfile1mapped, scalar(@$ref_unmatched1) );

if ( scalar(@$ref_matched1) == 0 ) {
    die "# ERROR: failed mapping $sp1 genes in WGA alignment";
} else {
    foreach $gene (@$ref_unmatched1) {
        print "# unmapped: $gene\n"
    }
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

unlink($intersectBEDfile_sorted) if(-e $intersectBEDfile_sorted);

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
# Note: genes in such pairs can actually be $hom_gene_id
my $num_segments2 = genes_mapped2segments( $geneBEDfile2, $geneBEDfile2mapped, 
	$gene_intersectBEDfile, $segment_intersectBEDfile, 1 );

my $num_segments1 = genes_mapped2segments( $geneBEDfile1, $geneBEDfile1mapped,
        $gene_intersectBEDfile, $segment_intersectBEDfile1, 1 );

$cmd = "$SORTBIN $SORTPARS -k1,1 -k2,2n $gene_intersectBEDfile $segment_intersectBEDfile ".
           "$segment_intersectBEDfile1 > $intersectBEDfile_sorted";

system($cmd);
if ( $? != 0 ) {
    die "# ERROR: failed sorting (genes, $cmd)\n";
}
elsif ( !-s $intersectBEDfile_sorted ) {
    die "# ERROR: failed generating $intersectBEDfile_sorted file ($cmd)\n";
}

push(@tmpBEDfiles, $segment_intersectBEDfile, $segment_intersectBEDfile1);
push(@tmpBEDfiles, $gene_intersectBEDfile, $intersectBEDfile_sorted);

my ($num_pairs, $num_segments, $hits_per_gene) = 
    bed2compara( $intersectBEDfile_sorted, $geneBEDfile1, $geneBEDfile2,
        $outfilename, $sp1, $sp2, $noheader);

printf( "# %d collinear gene pairs , %d collinear segments, %1.3f hits/gene\n", 
    $num_pairs, $num_segments, $hits_per_gene );

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
        if(!open(FASTA,"$GZIPBIN -dc $fastafile |")) {
            die "# ERROR(read_FASTA_regex2hash): cannot read GZIP compressed $fastafile $!\n"
                ."# please check gzip is installed\n";
        }
    } elsif($fastafile =~ /\.bz2$/ || $magic eq "BZ") {
        if(!open(FASTA,"$BZIP2BIN -dc $fastafile |")){
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

# Takes 7 params:
# 1) path to write files to
# 2) name of FASTA file (ref)
# 3) name of FASTA file (query)
# 4) regex to match chromosome names
# 5) bedtools path
# 6) indexing job (boolean)
# 7) reuse (boolean)
# Returns ref to hash with chr and/or 'unplaced' as keys and two FASTA files as value (ref, query)
# Note: 'unplaced' might hold genuine unplaced sequences but also non-shared chr names
# Note: creates FASTA & BED files with prefix _
# Note: if doindex is true only creates files for $fastafile1
sub split_genome_sequences_per_chr_bedtools {
    my ($path,$fastafile1,$fastafile2,$regex,$bedtools_path,$doindex,$reuse) = @_;

    my ($faifile1, $faifile2) = ( $fastafile1.'.fai' , $fastafile2.'.fai' );
    my ($chr,$chrfasta1,$chrfasta2,$unplacedbed1,$unplacedbed2,$cmd);
    my (%shared_chrs,%chr_pairs);

    my $ref_bed1 = read_FAI_regex2hash($faifile1,$regex);
    my $ref_bed2 = read_FAI_regex2hash($faifile2,$regex);

    # check chr names found in both files
    foreach $chr (keys(%$ref_bed1)){
        if(defined($ref_bed2->{$chr}) && $chr !~ 'unplaced'){
            $shared_chrs{$chr} = 1;
        }
    } 

    # write chr-specific FASTA files
    foreach $chr (keys(%shared_chrs)) { 
        $chrfasta1 = $path . "_".basename($fastafile1).".$chr.fna";
        $chrfasta2 = $path . "_".basename($fastafile2).".$chr.fna";

        if(!$reuse || !-e $chrfasta1) {
            $cmd = "echo '$ref_bed1->{$chr}' | $bedtools_path getfasta -fi $fastafile1 -bed stdin | " .
                " perl -lne 's/^>(\\S+?):\\d+-\\d+/>\$1/; print' > $chrfasta1";
            system($cmd); #print($cmd);
            sleep(2);
            if ( $? != 0 ) {
                die "# ERROR(split_genome_sequences_per_chr_bedtools): failed running bedtools (chr, $cmd)\n";
            } elsif ( !-s $chrfasta1 ) {
                die "# ERROR(split_genome_sequences_per_chr_bedtools): failed generating $chrfasta1 file ($cmd)\n";
            }
        } 

        if(!$doindex && (!$reuse || !-e $chrfasta2)) {
            $cmd = "echo '$ref_bed2->{$chr}' | $bedtools_path getfasta -fi $fastafile2 -bed stdin | " .
                " perl -lne 's/^>(\\S+?):\\d+-\\d+/>\$1/; print' > $chrfasta2";
            system($cmd);
            sleep(2);
            if ( $? != 0 ) {
                die "# ERROR(split_genome_sequences_per_chr_bedtools): failed running bedtools (chr, $cmd)\n";
            } elsif ( !-s $chrfasta1 ) {
                die "# ERROR(split_genome_sequences_per_chr_bedtools): failed generating $chrfasta2 file ($cmd)\n";
            }
        }

        $chr_pairs{$chr} = [$chrfasta1,$chrfasta2];
    }

    # write unplaced and/or not-shared chr names
    $chrfasta1 = $path . "_".basename($fastafile1).".unplaced.fna";
    $chrfasta2 = $path . "_".basename($fastafile2).".unplaced.fna";

    # as there might be many unplaced contigs, these are better written to file
    $unplacedbed1 = $path . "_".basename($fastafile1).".unplaced.bed";
    $unplacedbed2 = $path . "_".basename($fastafile2).".unplaced.bed";

    if($ref_bed1->{'unplaced'} && (!$reuse || !-e $chrfasta1)) {

        open( UBED1, ">", $unplacedbed1) ||
            die "# ERROR(split_genome_sequences_per_chr_bedtools): cannot write $unplacedbed1\n";
        print UBED1 $ref_bed1->{'unplaced'};
        close(UBED1);        

        $cmd = "$bedtools_path getfasta -fi $fastafile1 -bed $unplacedbed1 | " .
            " perl -lne 's/^>(\\S+?):\\d+-\\d+/>\$1/; print' > $chrfasta1";
        system($cmd); 
        sleep(2);
        if ( $? != 0 ) {
            die "# ERROR(split_genome_sequences_per_chr_bedtools): failed running bedtools (chr, $cmd)\n";
        } elsif ( !-s $chrfasta1 ) {
            die "# ERROR(split_genome_sequences_per_chr_bedtools): failed generating $chrfasta1 file ($cmd)\n";
        }
    }

    if(!$doindex && $ref_bed2->{'unplaced'} && (!$reuse || !-e $chrfasta2)) {

        open( UBED2, ">", $unplacedbed2) ||
            die "# ERROR(split_genome_sequences_per_chr_bedtools): cannot write $unplacedbed2\n";
        print UBED2 $ref_bed2->{'unplaced'};
        close(UBED2);

        $cmd = "$bedtools_path getfasta -fi $fastafile2 -bed $unplacedbed2 | " .
            " perl -lne 's/^>(\\S+?):\\d+-\\d+/>\$1/; print' > $chrfasta2";
        system($cmd); 
        sleep(2);
        if ( $? != 0 ) {
            die "# ERROR(split_genome_sequences_per_chr_bedtools): failed running bedtools (chr, $cmd)\n";
        } elsif ( !-s $chrfasta1 ) {
            die "# ERROR(split_genome_sequences_per_chr_bedtools): failed generating $chrfasta2 file ($cmd)\n";
        }
    }

    $chr_pairs{'unplaced'} = [$chrfasta1,$chrfasta2];

    return \%chr_pairs;
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

    my ($chr,$chrfasta1,$chrfasta2);
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

# Takes 
# i) reference FAI filename (string)
# ii) input BED intersect filename (string)
# iii) output BED filename (string)
# iv) min quality score (real)
# v) min alignment length (natural)
# vi) same strand only (boolean) 
# vii) verbose, optional (boolean)
#
# Parses sorted BED intersect -wo output and writes to BED file
# features (cDNA/transcripts) mapped on reference genome. Note:
# features might be unsorted.
#
# Returns 
# i) ref to list of BED-like lines of matched genes (coordinates in reference space)
# ii) ref to list of BED-like lines of unmatched genes
# iii) % genes in WGA blocks of at least 3 genes (float)
#
# Note: able to parse cs::Z (minimap2) and cg::Z (CGAlign, wfmash) strings
# Note: takes first match of each cDNA/gene only
# Note: discards partially/poorly mapped genes 
#
# example input (infile), all coords are 0-based:
# 1 4848 20752 ONIVA01G00010 9999     + 1 3331 33993       + 6 26020714 26051403 29819 60 cs:Z::303*ag:30*ga... 15904
# 1 104921 116326 ONIVA01G00100 9999  + 1 103118 152580    + 1 1132 47408 45875 60 cs:Z::70*tc:...              11405
# Chr1 2903 10817 LOC_Os01g01010 9999 + Chr1 1000 10053455 + 1 1000 10053455 10052455 60 cs:Z::10052455          7914
# Chr1 2903 10817 LOC_Os01g01010 9999 + Chr1 896 705000    + 1 949 705000 704051 41 cg:Z:51=53I704000=           7914
# <--             (c)DNA/gene       --> <- (q)uery genome -> <-- (r)eference genome                         -->  ovlp
#
# Note: ovlp is the actual intersection overlap computed by bedtools intersect
#
sub query2ref_coords {

    my ( $refaifile, $infile, $outfile, $minqual, $minalnlen, $samestrand, $verbose ) = @_;

    my ( $cchr, $cstart, $cend, $cname, $cmatch, $cstrand );
    my ( $qchr, $qstart, $qend, $cigartype, $bedline );
    my ( $WGAstrand, $rchr,      $rstart,       $rend );
    my ( $rmatch,    $rmapqual,  $SAMPAFtag,    $overlap, $done, $strand );
    my ( $SAMqcoord, $SAMrcoord, $feat,         $coordr );
    my ( $deltaq,    $deltar,    $start_deltar, $end_deltar );
    my ( %ref_coords, %genes_per_block, %matched_gene, %unmatched, %ref_max_length);
    my ( @matched, @filt_unmatched, @segments );

    # parse reference chr lengths
    my $ref_chr_bed = read_FAI_regex2hash($refaifile);
    foreach $rchr (keys(%$ref_chr_bed)){
      chomp($ref_chr_bed->{$rchr});
      $ref_max_length{$rchr} = (split(/\t/,$ref_chr_bed->{$rchr}))[2];	
    }

    # parse input file
    open( BED, "<", $infile )
      || die "# ERROR(query2ref_coords): cannot read $infile\n";
    while (<BED>) {

        chomp;

        (
            $cchr, $cstart, $cend, $cname, $cmatch, $cstrand,   # cDNA/gene, [cstart-cend] to be found within [qstart-qend]
			$qchr, $qstart, $qend, $WGAstrand,                  # query in WGA 
            $rchr, $rstart, $rend, $rmatch,                     # reference in WGA
			$rmapqual, $SAMPAFtag,                              # WGA details
			$overlap                                            # overlap of cDNA/gene with query in WGA

        ) = split( /\t/, $_ ); 

        # take only 1st mapping passing QC
        next if ( defined($ref_coords{$cname}) );

        # skip mappings where query and ref segment on different strands if requested
        if ( $WGAstrand eq '-' && $samestrand == 1 ) {
            if(!$unmatched{$cname}) {
			    $unmatched{$cname} =  
                    "[different strand] $cchr\t$cstart\t$cend\t$cname\t$cmatch\t" .
                    "$cstrand\t$qchr\t$qstart\t$qend\t$WGAstrand\t$rchr" .
                    "$rstart\t$rend\t$rmatch\t$rmapqual\t$overlap";
            }
            next;
        }

        # skip poor query-to-ref WGA scores
        if ( $rmapqual < $minqual ) {
		    if(!$unmatched{$cname}) {
                $unmatched{$cname} =
                    "[quality $rmapqual < $minqual] $cchr\t$cstart\t$cend\t$cname\t$cmatch\t" .
                    "$cstrand\t$qchr\t$qstart\t$qend\t$WGAstrand\t$rchr" .
                    "$rstart\t$rend\t$rmatch\t$rmapqual\t$overlap";
            }
            next;
        }

        #next if($cname ne 'ONIVA01G00100'); # debug
		#next if($cname ne 'gene:OsIR64_12g0001370'); 
		#next if($cname ne 'gene:OsIR64_05g0012650'); 
        #next if($cname ne 'LOC_Os01g13840'); print "$_\n"; # + strand
	    #next if($cname ne 'LOC_Os10g29730'); print "$_\n"; # - strand 

        # record <genes per alignment block
        $genes_per_block{"$qchr:$qstart-$qend"}++;

        # make sure ref chr is taken
        $ref_coords{$cname}{'chr'} = $rchr;

        ## correct offset for ref assembly

        # Note: this requires parsing SAM/PAF tag
        # https://github.com/lh3/minimap2#paftools
        # cs:Z::303*ag:32*ga:27+ctattcta*ag*ca:20*ag:3*ga:18*tc*ga:76-tc
        # or
        # cg:Z:2310M12I805M2I225M13I895M2I3402M2I1184M1I18

        # default start coords (end in - strand),
        # in case 3'cDNA not included in WGA segment (ref_coords),
        # scan is always qstart -> qend 
        $SAMqcoord = $qstart;
        if ( $WGAstrand eq '+' ) {
            # q ------------>			
            # r ------------>
            $SAMrcoord = $rstart;
            $ref_coords{$cname}{'start'} = $rstart;
        }
        else {
            # q ------------>
            # r <------------
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
            return ( \@matched, \@filt_unmatched );
        }

        # split CIGAR string into individual feature tags
        if ( $cigartype eq 'cs' ) {

            if ( $WGAstrand eq '+' ) {
                @segments = split( /:/, $SAMPAFtag )
            } else {
                @segments = reverse split( /:/, $SAMPAFtag )
            }
        }
        else {
            while ( $SAMPAFtag =~ m/(\d+[MIDNSHP=X])/g ) {
                push( @segments, $1 );
            }

            if ( $WGAstrand ne '+' ) {
                @segments = reverse @segments;
            }
        }

        # loop along PAF alignment features updating coords, 
        # qcoord overlaps are used as exit condition
        $done = 0;
        foreach $feat (@segments) {

            ( $deltaq, $deltar, $coordr ) =
              _parseCIGARfeature( $feat, $SAMqcoord, $SAMrcoord );

            # check if current position in query alignment overlaps cDNA/gene coords

            # start coords (end in - strand)
            if ( $SAMqcoord < $cstart
                && ( $SAMqcoord + $deltaq ) >= $cstart ) {

                # refine delta to match exactly the start (end for strand -)
                ( $deltaq, $deltar, $coordr ) =
                  _parseCIGARfeature( $feat, $SAMqcoord, $SAMrcoord, $cstart ); 
                  print "ss $deltaq, $deltar, $coordr : $feat, $SAMqcoord, $SAMrcoord, $cstart\n"
                    if ( $verbose > 1 );				  

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
                    $ref_coords{$cname}{'end'} = $SAMrcoord - $start_deltar;
                }
            }

            # end coords (start in -strand), actually copied out of loop
            if ( $SAMqcoord < $cend && ( $SAMqcoord + $deltaq ) >= $cend ) {

                # refine delta to match exactly the end (start for strand -)
                ( $deltaq, $deltar, $coordr ) =
                  _parseCIGARfeature( $feat, $SAMqcoord, $SAMrcoord, $cend ); 
                  print "ee $deltaq, $deltar, $coordr : $feat, $SAMqcoord, $SAMrcoord, $cend\n"
                    if ( $verbose > 1 );

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
			
            # CGATCGATAAATAGAGTAG---GAATAGCA
            # ||||||   ||||||||||   |||| |||
            # CGATCG---AATAGAGTAGGTCGAATtGCA
            # ..6-ata..                       deltaq=9  deltar=6
            #          ....10+gtc...          deltaq=10 deltar=13
            #                       .4*at     deltaq=5  deltar=5
            #                            ..3  deltaq=3  deltar=3

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

        # put together BED line
        $bedline = sprintf("%s\t%d\t%d\t%s\t%d\t%s\n",
            $ref_coords{$cname}{'chr'},
            $ref_coords{$cname}{'start'},
            $ref_coords{$cname}{'end'},
            $cname,
            $overlap,
            $strand);

        # compute actual overlap between cDNA/CDS and ref genome
        $overlap = 1 + $ref_coords{$cname}{'end'} - $ref_coords{$cname}{'start'};

        # check match is within reference chr bounds
        if($ref_coords{$cname}{'start'} < 0 ||
            (defined($ref_max_length{$ref_coords{$cname}{'chr'}}) && 
            $ref_coords{$cname}{'end'} > $ref_max_length{$ref_coords{$cname}{'chr'}})) {
            $unmatched{$cname} = "[overlap outside chr] $bedline";

        } elsif ( $overlap >= $minalnlen ) { # if overlapping segment is long enough
            push(@matched, $bedline);
            $matched_gene{ $cname } = 1;
            print "$bedline\n" if($verbose > 1)
        }
        else {
            # skip gene models with short overlap with WGA segments
            $unmatched{$cname} = "[overlap $overlap < $minalnlen] $bedline";
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

    # compute % matched genes in WGA blocks of at least 3 genes
    my $totgenes = 0;
    foreach my $block 
        (sort {$genes_per_block{$b}<=>$genes_per_block{$a}} 
        keys(%genes_per_block)) {

        last if($genes_per_block{$block} < 3);
        $totgenes += $genes_per_block{$block};   
    } 

    # remove matched genes from raw list of unmatched
    foreach $cname (keys(%unmatched)) {
        next if( $matched_gene{ $cname } ); 
        push(@filt_unmatched, $unmatched{$cname});
    }

    return ( \@matched, \@filt_unmatched, 100*$totgenes / scalar(@matched) );
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

    # q ------------> (strand not considered)
    # r ------------>

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

# Takes 2 params:
# i)  BED filaname 
# ii) (optional) ref to hash with names (4th column) to skip
# Returns a hash ref with names (keys) pointing to full BED lines (values)
sub _BED2hash {

    my ($BEDfile, $ref_skip_names) = @_;

    my ($name,%coords);

    open(BED,"<",$BEDfile)
        || die "# ERROR(_BED2hash): cannot read $BEDfile\n";

    while(<BED>) {
        #1    4847    20752   gene:ONIVA01G00010      9999    +
        if(/^\S+\t\d+\t\d+\t(\S+)\t/) {

            $name = $1;
            next if($ref_skip_names && $ref_skip_names->{$name});

            chomp;
            $coords{$name} = $_;
        }
    }
    close(BED);

    return \%coords;
}


# Takes four file paths and optionally a boolean:
# i)   BED filename with gene model coordinates of species2 (6 cols)
# ii)  BED filename with species1 genes mapped on species2 space (6 cols)
# iii) BED filename with sp1,sp2 pairs of collinear genes (13 columns)
# iv)  BED output filename (13 columns)
# v)   (optional) boolean, columns from sp2 should be printed first
# Returns number of genomic segments in sp1 collinear to genes in sp2
sub genes_mapped2segments {

    my ($geneBEDfile, $mappedBEDfile, $pairedBEDfile, $outBEDfile, $invert) = @_;

    my $num_segments = 0;
    my ($chr,$sta,$end,$geneid,$len,$strand,$overlap);
    my (@genes,$ref_orig_coords,%paired);
    
    # find out which genes are paired based on WGA (sp2)
    open(PAIRBED,"<",$pairedBEDfile)
        || die "# ERROR(gene2segments): cannot read $pairedBEDfile\n";

    while(<PAIRBED>) {
        #1 30219   36442   gene:BGIOSGA002569 9999 + 1 29700 39038 gene:ONIVA01G00100  9339  +  6223
        my @data = split(/\t/,$_);
        $paired{$data[3]} = 1; # sp1
        $paired{$data[9]} = 1; # sp2
    }
    close(PAIRBED);

    # read gene model coordinates (sp2)
    $ref_orig_coords = _BED2hash($geneBEDfile,\%paired);

    # find sp2 genes that have mapped coords on sp1 but are not paired, print to outfile
    open(OUTBED,">",$outBEDfile) 
        || die "# ERROR(gene2segments): cannot create $outBEDfile\n";

    open(MAPBED,"<",$mappedBEDfile)
        || die "# ERROR(gene2segments): cannot read $mappedBEDfile\n";

    while(<MAPBED>) {
        #1       29700   39038   gene:ONIVA01G00100      9339    +
        if(/^(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\d+)\t([+-])/) {
            ($chr, $sta, $end, $geneid, $len, $strand) = ($1, $2, $3, $4, $5, $6);
            
            next if($paired{$geneid});

            # actually print to BED coordinates of genes mapped to (unannotated) genomic segments
            #1 217360 222398 segment 5039 + 1 155040 165322 gene:ONIVA01G00180 9999 + 9999

            $overlap = $end-$sta;

            if($invert) {
                print OUTBED "$ref_orig_coords->{$geneid}\t$chr\t$sta\t$end\tsegment\t$len\t$strand\t$overlap\n";
            } else {
                print OUTBED "$chr\t$sta\t$end\tsegment\t$len\t$strand\t$ref_orig_coords->{$geneid}\t$overlap\n";
            }
           
            $num_segments++;
        }
    }
    close(MAPBED);

    close(OUTBED);

    return $num_segments;
}


# Parses bedtools intersect file and BED files with original gene coords 
# to produce Ensembl Compara-like TSV file, with 1-based coordinates
# Takes 7 parameters:
# i)    BED intersect filename (see format below)
# ii)   BED file with sp1 gene coords
# iii)  BED file with sp2 gene coords
# iv)   TSV output filename (Compara-like format) 
# v)    string with name of species1
# vi)   string with name of species2
# vii)  boolean to remove header from TSV
# TODO) boolean to remove redundant parts from gene names (experimental, see below)
# Returns:  
# i) number of collinear gene pairs 
# ii) number of collinear segments
# iii) float with average hits per gene (on same species)
#
# Columns in TSV output format, homology_type = ['ortholog_collinear','segment_collinear']:
# gene_stable_id protein_stable_id species overlap homology_type homology_gene_stable_id
# homology_protein_stable_id homology_species overlap dn ds goc_score wga_coverage
# is_high_confidence coordinates
sub bed2compara {

    my ( $infile, $geneBEDfile1, $geneBEDfile2, $TSVfile, $sp1, $sp2, $noheader ) = @_;

    my ( $gene1, $gene2, $coords1, $coords2, $coords, $homoltype );
    my ($num_pairs, $num_segments) = (0, 0);
    my ($num_matched_genes, $num_hits_genes, %hit) = (0, 0);
    my ($ref_orig_coords1,$ref_orig_coords2, $realsp1, $realsp2);

    # read gene model coordinates 
    $ref_orig_coords1 = _BED2hash($geneBEDfile1);
    $ref_orig_coords2 = _BED2hash($geneBEDfile2);

    # parse intersection BED and produce TSV output
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
        #1 116435  120177  ONIV..  9999 +  1 116779  119742  gene:Os01g0100400       2964    + 2963
        #1 116866  118035  segment 1170 -  1 12807   13978   gene:Os01g0100466       9999    - 1169
        #1 160018  166571  ..00200 9999 -  1 191201  197773  segment                 6573    - 6572

        # old behavior:
		#
        #0 invert false
        #1  4843    11631   segment 6789    -   6   26022262    26029047    gene:Os06g0639700   9999    -   6788
		# ->
        #Oryza_nivara:1:4843-11631(-)   segment Oryza_nivara    6788    segment_collinear   gene:Os06g0639700   \
        #  gene:Os06g0639700   Oryza_sativa    6788    NULL    NULL    NULL    100.00  1   1:4843-11631(-);6:26022262-26029047(-)

        #1 invert true 
        #6  26022262    26029047    gene:Os06g0639700   9999    -   1   4843    11631   segment 6789    -   6788
        # ->
        #gene:Os06g0639700  gene:Os06g0639700   Oryza_nivara    6788    segment_collinear   Oryza_sativa:1:4843-11631(-) \
        #  segment Oryza_sativa    6788    NULL    NULL    NULL    100.00  1   NA;1:4843-11631(-)

        # behavior after setting invert true in all cases:
		#gene:Os06g0639700	gene:Os06g0639700	Oryza_sativa	6788	segment_collinear	Oryza_nivara:1:4843-11631(-) \
        #  segment	Oryza_nivara	6788	NULL	NULL	NULL	100.00	1	6:26022262-26029047(-);1:4843-11631(-)

        my @data = split( /\t/, $_ );

        if(scalar(@data) < 13) {
            print "# WARN(bed2compara): skip short line ($. $infile)\n";
            next;
        }

        # format gene names
        $gene1 = $data[3];
        $gene2 = $data[9];

        # init species names		
        $realsp1 = $sp1;
        $realsp2 = $sp2;		
  
        # workout genomic coordinates (1-based)
        if($data[3] ne 'segment') {
            if( $ref_orig_coords1->{$gene1} &&
                $ref_orig_coords1->{$gene1} =~ m/^(\S+)\t(\d+)\t(\d+)\t\S+\t\S+\t(\S+)/) {
                $coords1 = "$1:$2-$3($4)"
            } elsif( $ref_orig_coords2->{$gene1} &&
                $ref_orig_coords2->{$gene1} =~ m/^(\S+)\t(\d+)\t(\d+)\t\S+\t\S+\t(\S+)/) {
        
                # new behavior, gene1 might actually be from sp2 (segments)   				
                $coords1 = "$1:$2-$3($4)";
                $realsp1 = $sp2;
                $realsp2 = $sp1;

            } else { $coords1 = 'NA' }
        } else {
            $coords1 = "$data[0]:$data[1]-$data[2]($data[5])";
        }

        if($data[9] ne 'segment') { 
            if( $ref_orig_coords2->{$gene2} &&
                $ref_orig_coords2->{$gene2} =~ m/^(\S+)\t(\d+)\t(\d+)\t\S+\t\S+\t(\S+)/) {
                $coords2 = "$1:$2-$3($4)"
            } else { $coords2 = 'NA' }
        } else {
            $coords2 = "$data[6]:$data[7]-$data[8]($data[11])";
        }
 
        $coords = "$coords1;$coords2";

        # check homology type and work out genomic coordinates
        if($data[3] eq 'segment') {
            $homoltype = 'segment_collinear';
            $gene1 = "$realsp1:$coords1";
            $num_segments++;
        } elsif($data[9] eq 'segment') {
            $homoltype = 'segment_collinear';
            $gene2 = "$realsp2:$coords2";
            $num_segments++;
        } else {
            $homoltype = 'ortholog_collinear';
            $num_pairs++;
        }

        printf( TSV
            "%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%d\tNULL\tNULL\tNULL\t%1.2f\t%d\t%s\n",
            $gene1,
            $data[3],
            $realsp1,
            $data[12],
            $homoltype,
            $gene2,
            $data[9],
            $realsp2,
            $data[12],
            100,    # max WGA score
            1,      # high confidence
            $coords
        );

        # compile match stats
        if($homoltype ne 'segment_collinear') {
            if(!$hit{$gene1}{$sp2}) { 
                $num_matched_genes++ 
            }
            $hit{$gene1}{$sp2}=1;
            $num_hits_genes++;
        }
    }

    close(TSV);
    close(BEDINT);

    return ($num_pairs, $num_segments, $num_hits_genes/$num_matched_genes);
}

# Takes up to 3 parameters:
# i)   filename of MAF alignment
# ii)  filename of output PAF alignment
# iii) boolean to request quality control (optional)
# Converts MAF format (https://genome.ucsc.edu/FAQ/FAQformat.html#format5) 
# to simplified 3-state PAF (M,D,I, http://samtools.github.io/hts-specs/SAMv1.pdf). 
# Coordinates of strand - alignments are reversed to match the PAF style of minimap2.
# PAF output includes 12 standard columns + cg:Z: CIGAR (13th column).
# Returns number of alignments parsed
sub simpleMAF2PAF {

    my ($MAFfile, $outPAFfile, $QC) = @_;
   
    my ($isQuery,$cigar,$segment_len,$indels,$aligned);
    my ($chr, $start, $len, $strand, $chrlen, $seq);
    my ($qchr, $qstart, $qlen, $qstrand, $qchrlen, $qseq);   
    my (@gaps,%gap,$last,$pos,$nextpos,$g,$score);
    my $total_alignments = 0;

    open(MAF,"<",$MAFfile) || 
        die "# ERROR(convertMAF2PAF): cannot read $MAFfile\n";

    open(PAF,">",$outPAFfile) || 
        die "# ERROR(convertMAF2PAF): cannot create $MAFfile\n";

    while(<MAF>) {

        #a score=18842
        #s ref.1 17370551 22818 + 44309235 CACGAAGTGAGTTTTTGCCATGAACG...
        #s qry.1 27077898 22753 - 44350042 CACGAAGTGAGTTTTTGCCATGAACG...
        #
        #a score=18776
        #s ref.1 35900090 18888 + 44309235 AGGTATAAATGGA--GAGAGAGAGCA...
        #s qry.1 36094320 18894 + 44350042 AGGTATAAATGGAGAGAGAGAGAGCA...

        if(/^a score=(\S+)/) {
            $score = $1;
            ($isQuery, $aligned, $indels) = (0,0,0);
            $total_alignments++;

        } elsif(/^s \w{3}\.(\S+)\s+(\d+)\s+(\d+)\s+(\S)\s+(\d+)\s+(\S+)/) {
            
            if($isQuery) { # query, 2nd in pair
                $isQuery++;
                ($qchr, $qstart, $qlen, $qstrand, $qchrlen, $qseq) = ($1,$2,$3,$4,$5,$6);

                # reverse alignment so that output PAF format matches that of minimap2
                if($qstrand eq '-') {
                    $qstart = $qchrlen - $qstart;
                    # Note: if actual sequence wanted should be rev complemented
                }

            } else { # reference, 1st in pair, always + strand (as opposed to minimap2)
                $isQuery++;
                ($chr, $start, $len, $strand, $chrlen, $seq) = ($1,$2,$3,$4,$5,$6);
            } 

            # actually parse aligned sequences
            if($isQuery == 2) {

                # gapless
                if($qlen == $len && $qseq !~ '-' && $seq !~ '-') {
                    $cigar = "cg:Z:$qlen".'M';
                    $aligned = $qlen;   
                    $indels = 0; 
                } else { # gapped alignments

                    @gaps = ();
                    %gap = ();
                    $cigar = 'cg:Z:';

                    # find gap positions (D & I),
                    # rest are matched (M) segments
                    while($qseq =~ m/(\-+)/g) {
                        #     start    type  totalgaps 
                        $gap{$-[0]} = ['I', length($1)] # insertion to reference
                    } 

                    while($seq =~ m/(\-+)/g) {
                        $gap{$-[0]} = ['D', length($1)] # deletion in reference
                    }        

                    # sort 0-based positions and get bounds 
                    $last = length($seq) - 1;
                    @gaps = sort {$a<=>$b} keys(%gap); #print join(',',@gaps)."\n";

                    # add segment before first gap
                    if($gaps[0] > 0) {
                        $cigar .= $gaps[0] .'M';
                        $aligned += $gaps[0];   
                    }

                    # add gaps and matched segment after
                    foreach $g (0 .. $#gaps-1) {

                        # gap
                        $pos = $gaps[$g];
                        $cigar .= "$gap{$pos}->[1]$gap{$pos}->[0]";
                        $indels += $gap{$pos}->[1];
                    
                        # segment after
                        $pos += $gap{$pos}->[1]; 
                        $nextpos = $gaps[$g+1];
                        $segment_len = $nextpos - $pos;
                        $cigar .= $segment_len.'M';
                        $aligned += $segment_len;
                    }

                    # last gap
                    $pos = $gaps[$#gaps]; 
                    $cigar .= "$gap{$pos}->[1]$gap{$pos}->[0]";
                    $indels += $gap{$pos}->[1];

                    # add segment after last gap
                    $pos = $gaps[$#gaps] + $gap{$pos}->[1]; #print "$pos $last\n";
                    if($pos <= $last) {
                        $segment_len = ($last-$pos)+1;
                        $cigar .= $segment_len.'M';
                        $aligned += $segment_len;
                    }

                    # QC
                    if($QC && $aligned + $indels != $last+1){
                        print "# ERROR: alignment length does not match CIGAR segments\n$_\n";
                    }
                }

                # print PAF
                #https://lh3.github.io/minimap2/minimap2.html#10
                #1	string	Query sequence name
                #2	int	Query sequence length
                #3	int	Query start coordinate (0-based)
                #4 	int	Query end coordinate (0-based)
                #5	char	+ if query/target on the same strand; - if opposite
                #6	string	Target sequence name
                #7	int	Target sequence length
                #8	int	Target start coordinate on the original strand
                #9	int	Target end coordinate on the original strand
                #10	int	Number of matching bases in the mapping
                #11	int	Number bases, including gaps, in the mapping
                #12	int	Mapping quality (0-255 with 255 for missing)
                # +
                #13     string  CIGAR string

                printf(PAF "%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n",
                    $qchr, 
                    $qchrlen,
                    ($qstrand eq $strand ? $qstart       : $qstart-$qlen),
                    ($qstrand eq $strand ? $qstart+$qlen : $qstart+1),
                    ($qstrand eq $strand ? '+'           : '-'),
                    $chr,
                    $chrlen,
                    $start,
                    $start+$len,
                    $aligned,
                    ($qlen > $len ? $qlen : $len),
                    $score,
                    $cigar
                );
            }
        }
    }

    close(PAF);

    close(MAF);

    return $total_alignments;
}

