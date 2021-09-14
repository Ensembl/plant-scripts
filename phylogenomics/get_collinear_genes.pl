#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

# Takes two FASTA files with genome sequences and 2 matching GFF files with annotated gene models.
# Produces TSV file with pairs of collinear genes in a format compatible with Ensembl Compara.
# Uses external software:
# minimap2 [https://academic.oup.com/bioinformatics/article/34/18/3094/4994778]
# bedtools [https://academic.oup.com/bioinformatics/article/26/6/841/244688]
# wfmash   [https://github.com/ekg/wfmash]

#perl get_collinear_genes.pl -sp1 oryza_sativa -fa1 Oryza_sativa.IRGSP-1.0.dna.toplevel.fa -gf1 Oryza_sativa.IRGSP-1.0.51.gff3 -al1 IRGSP -sp2 oryza_nivara -fa2 Oryza_nivara.Oryza_nivara_v1.0.dna.toplevel.fa -gf2 Oryza_nivara.Oryza_nivara_v1.0.51.gff3 -al2 OGE -r

my $MINIMAP2EXE = './minimap2'; # 2.22
my $MINIMAPTYPE = '-x asm20';
my $MINIMAPPARS = "--secondary=no --cs --cap-kalloc=1g $MINIMAPTYPE";
my $WFMASHEXE   = './wfmash'; # v0.7.0
my $WFMASHPARS  = '-s 2000 -t 3'; # works for small genomes such as rice
my $BEDTOOLSEXE = 'bedtools'; # v2.30.0
my $BEDINTSCPAR = '-wo -f XXX -F XXX -e'; # XXX to be replaced with [0-1]

#my $BLASTNEXE   = 'blastn';   # 2.2.30+

my $SORTLIMITRAM = "500M"; # buffer size
my $SORTBIN      = "sort --buffer-size=$SORTLIMITRAM";
my $GZIPBIN      = "gzip";

my $DUMMYSCORE = 9999;

# while parsing PAF
my $MINQUAL    = 60; 
my $MINALNLEN  = 100;
my $SAMESTRAND = 1;
my $MINOVERLAP = 0.50;
my $VERBOSE    = 0; # values > 1

# Work out gene names from transcripts':
# 1) remove suffix after . or -
# 2) t's in transcript name to be replaced with g' in gene names
# Example: 10t0100300.1 -> Os10g0100300
my $TRANSCRIPT2GENE = 1;

my ( $help, $do_sequence_check, $reuse, $noheader, $dowfmash ) = (0,0,0,0,0);
my ( $sp1, $fasta1, $gff1, $sp2, $fasta2, $gff2, $label1, $label2 );
my ( $minoverlap, $qual, $alg, $outfilename ) = ($MINOVERLAP, $MINQUAL, 'minimap2');

GetOptions(
	"help|?"         => \$help,
	"sp1|species1=s" => \$sp1,
	"fa1|fasta1=s"   => \$fasta1,
	"gf1|gff1=s"     => \$gff1,
	"al1|label1=s"   => \$label1,
	"sp2|species2=s" => \$sp2,
	"fa2|fasta2=s"   => \$fasta2,
	"gf2|gff2=s"     => \$gff2,
	"al2|label2=s"   => \$label2,
	"out|outfile=s"  => \$outfilename,
	"ovl|overlap=f"  => \$minoverlap,
	"q|quality=i"    => \$qual,
	"c|check"        => \$do_sequence_check,
	"r|reuse"        => \$reuse,
	"wf|wfmash"      => \$dowfmash,
	"a|add"          => \$noheader
) || help_message();

sub help_message {
	print "\nusage: $0 [options]\n\n"
		. "-sp1 binomial/trinomial species name    (required, example: -sp1 oryza_sativa)\n"
		. "-fa1 genome FASTA [.gz] filename        (required, example: -fa1 oryza_sativa.fna)\n"
		. "-gf1 GFF [.gz] filename                 (required, example: -gf1 oryza_sativa.RAPDB.gff)\n"
		. "-al1 annotation label                   (required, example: -al1 RAPDB)\n"
		. "-sp2 binomial/trinomial species name    (required, example: -sp2 oryza_nivara)\n"
		. "-fa2 genome FASTA [.gz] filename        (required, example: -fa2 oryza_nivara.fna)\n"
		. "-gf2 GFF [.gz] filename                 (required, example: -gf2 oryza_nivara.OGE.gff)\n"
		. "-al2 annotation label                   (required, example: -al2 OGE)\n"
		. "-out output filename (TSV format)       (optional, by default built from input, example: -out rice.tsv)\n"
		. "-ovl min overlap of genes               (optional, default: -ovl $MINOVERLAP)\n" 
		. "-wf  use wfmash aligner                 (optional, by default minimap2 is used)\n"
		. "-q   min mapping quality                (optional, default: -q $MINQUAL)\n"
		#. "-c   check sequences of collinear genes (optional)\n"
		. "-add concat TSV output with no header   (optional, example: -add, requires -out)\n"
		. "-r   re-use previous minimap2 results   (optional)\n";
}

if($help || 
	(!$sp1 || !$fasta1 || !$gff1 || !$label1 ||
		!$sp2 || !$fasta2 || !$gff2 || !$label2)){ 
	help_message();
	exit(0);
}

if($minoverlap && ($minoverlap < 0 || $minoverlap > 1)) {
	print "# ERROR: option -ovl requires values [0,1]\n";
    exit(0);
} 

if($qual && ($qual < 0 || $qual > 255)) {
    print "# ERROR: option -q requires values [0,255]\n";
    exit(0);
}

# add actual value to dummy values in param string
$BEDINTSCPAR =~ s/XXX/$minoverlap/g; 

if($noheader && !$outfilename) {
	print "# ERROR: option -add requires -out filename to concat results\n";
	exit(0);
}

# check algorithm
if($dowfmash){
	$alg = 'wfmash';
}

# set default outfile
if(!$outfilename) {
	$outfilename = ucfirst($alg).".homologies.$sp1.$label1.$sp2.$label2.overlap$minoverlap.tsv";
}

print "\n# $0 -sp1 $sp1 -fa1 $fasta1 -gf1 $gff1 -al1 $label1 ".
	"-sp2 $sp2 -fa2 $fasta2 -gf2 $gff2 -al2 $label2 -out $outfilename ".
	"-ovl $minoverlap -q $qual -wf $dowfmash -c $do_sequence_check -r $reuse\n\n";

print "# mapping parameters:\n";
if($dowfmash){
	print "# \$WFMASHPARS: $WFMASHPARS\n\n"
} else {
	print "# \$MINIMAPPARS: $MINIMAPPARS\n\n";
}

## 1) align genome1 vs genome2 with minimap2 (WGA)
## Note: masking not recommended, see https://github.com/lh3/minimap2/issues/654

my $PAFfile = "_$sp2.$label2.$sp1.$label1.$alg.paf";

print "# computing pairwise genome alignment with $alg\n\n";

if($reuse && -s $PAFfile){
	print "# re-using $PAFfile\n";
} else {

	if($dowfmash) {

		system("$WFMASHEXE $WFMASHPARS $fasta1 $fasta2 > $PAFfile");
        if($? != 0){
            die "# ERROR: failed running wfmash (probably ran out of memory)\n";
        }
        elsif(!-s $PAFfile){
            die "# ERROR: failed generating $PAFfile file (wfmash)\n";
        }
        else{
            print("# wfmash finished\n\n");
        }

	} else { # default minimap2 index & alignment

		my $index_fasta1 = "_$sp1.$label1.mmi";

		if($reuse && -s $index_fasta1){
			print "# re-using $index_fasta1\n";
		} else {
			system("$MINIMAP2EXE $MINIMAPTYPE -d $index_fasta1 $fasta1 2>&1");
			if($? != 0){
				die "# ERROR: failed running minimap2 (probably ran out of memory)\n";
			}
			elsif(!-s $index_fasta1){
        		die "# ERROR: failed generating $index_fasta1 file (minimap2)\n";
    		}
		}

		system("$MINIMAP2EXE $MINIMAPPARS $index_fasta1 $fasta2 -o $PAFfile 2>&1");
		if($? != 0){
			die "# ERROR: failed running minimap2 (probably ran out of memory)\n";
		}
		elsif(!-s $PAFfile){
			die "# ERROR: failed generating $PAFfile file (minimap2)\n";
		}
		else{
			print("# minimap2 finished\n\n");
		}
	}
}

## 2) produce BED-like file of sp2-to-sp1 coords 10 columns
## Note: $F[4] in PAF conveys whether query & ref are on the same strand or not

my $wgaBEDfile = "_$sp2.$label2.$sp1.$label1.$alg.bed";

open(PAF,"<",$PAFfile) || die "# ERROR: cannot read $PAFfile\n";
open(BED,">",$wgaBEDfile) || die "# ERROR: cannot create $wgaBEDfile\n";

while(<PAF>){
	my @F = split(/\t/,$_);
	print BED "$F[0]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[7]\t$F[8]\t$F[9]\t$F[11]\t$F[$#F]";
}

close(BED);
close(PAF);

## 3) Parse GFFs and produce BED files with gene coords

my $geneBEDfile1 = "_$sp1.$label1.gene.bed";
my $geneBEDfile2 = "_$sp2.$label2.gene.bed";

my ($num_genes1, $mean_gene_len1) = parse_genes_GFF($gff1,$geneBEDfile1);
printf("# %d genes parsed in %s mean length=%d\n",
	$num_genes1,$gff1,$mean_gene_len1);

my ($num_genes2, $mean_gene_len2) = parse_genes_GFF($gff2,$geneBEDfile2);
printf("# %d genes parsed in %s mean length=%\n",
	$num_genes2,$gff2,$mean_gene_len2);

## 4) intersect gene positions with WGA, sort by gene > cDNA ovlp > genomic matches

my $sp2wgaBEDfile = "_$sp2.$label2.gene.$sp1.$alg.intersect.overlap$minoverlap.bed";

system("$BEDTOOLSEXE intersect -a $geneBEDfile2 -b $wgaBEDfile $BEDINTSCPAR | ".
	"$SORTBIN -k4,4 -k5,5nr -k14,14nr > $sp2wgaBEDfile");
if($? != 0){
	die "# ERROR: failed running bedtools (WGA)\n";
}
elsif(!-s $sp2wgaBEDfile){
	die "# ERROR: failed generating $sp2wgaBEDfile file (bedtools)\n";
}

my $geneBEDfile2mapped = "_$sp2.$label2.$alg.gene.mapped.bed";

my ($num_matched, @unmatched) = 
	query2ref_coords($sp2wgaBEDfile, $geneBEDfile2mapped, 
		$qual, $MINALNLEN, $SAMESTRAND, $VERBOSE);

printf("# %d genes mapped in %s (%d unmapped)\n",
	$num_matched,$geneBEDfile2mapped,scalar(@unmatched));

## 5) produce list of pairs of collinear genes

my $gene_intersectBEDfile = "_$sp1.$label1.$sp2.$label2.$alg.gene.intersect.overlap$minoverlap.bed";

system("$BEDTOOLSEXE intersect -a $geneBEDfile1 -b $geneBEDfile2mapped $BEDINTSCPAR -s | ".
    "$SORTBIN -k4,4 -k13,13nr > $gene_intersectBEDfile");
if($? != 0){
    die "# ERROR: failed running bedtools (genes)\n";
}
elsif(!-s $gene_intersectBEDfile){
    die "# ERROR: failed generating $gene_intersectBEDfile file (bedtools)\n";
}

my $num_pairs = bed2compara($gene_intersectBEDfile, $outfilename, $sp1, $sp2,
	$noheader, $TRANSCRIPT2GENE);

printf("# %d collinear gene pairs\n",$num_pairs); 
print "# TSV file: $outfilename\n";


## sequence check (TODO)

#bedtools getfasta -fi Oryza_sativa.IRGSP-1.0.dna.toplevel.fa -bed test4.bed > test4.ref.fna
#bedtools getfasta -fi Oryza_nivara.Oryza_nivara_v1.0.dna.toplevel.fa \
#	-bed Oryza_nivara.gene.bed -nameOnly > test4.gene.fna

#makeblastdb -dbtype nucl -in test4.ref.fna
#blastn -query test4.gene.fna -db test4.ref.fna -outfmt 6 > test4.blastn

#perl -lane 'print "$F[3]\t$F[0]:$F[1]-$F[2]"' test4.bed > test4.bed4blastn
#grep -f test4.bed4blastn test4.blastn | cut -f 1 | sort -u | wc -l
#31020

#31020/31080 = 0.998 of transferred genomic segments confirmed by BLASTN

#################################

# Takes i) input GFF filename ii) output BED filename and 
# returns i) number and ii) mean length of genes parsed.
# Uses global $DUMMYSCORE as dummy scores, as opposed to cases where
# genes/transcript are actually mapped to genome
sub parse_genes_GFF {

	my ($gff_file, $outBEDfile) = @_;
	my $num_genes = 0;
	my $genelength = 0;
	my ($geneid, $magic);
	
	# check input file format and open it accordingly
	open(INFILE,$gff_file) || die "# ERROR(parse_genes_GFF): cannot read $gff_file, exit\n";
	sysread(INFILE,$magic,2);
	close(INFILE);

	if($gff_file =~ /\.gz$/ || $magic eq "\x1f\x8b"){ # GZIP compressed input{
		if(!open(GFF, "$GZIPBIN -dc $gff_file |")){
			die "# ERROR(parse_genes_GFF): cannot read GZIP compressed $gff_file $!\n"
				."# please check gzip is installed\n";
		}	
	} else {
		open(GFF,"<",$gff_file) || die "# ERROR(parse_genes_GFF): cannot read $gff_file\n";
	}	
	
	open(BED,">",$outBEDfile) || die "# ERROR(parse_genes_GFF): cannot create $outBEDfile\n";

	while(<GFF>){
		my @F = split(/\t/,$_);

		next if(scalar(@F)<9 || $F[2] ne "gene");

		# take only genes where ID can be parsed
		if($F[8] =~ /ID=([^;]+)/){ 
			
			$geneid = $1;
			$geneid =~ s/gene://; # remove redundant bits

			print BED "$F[0]\t$F[3]\t$F[4]\t$geneid\t$DUMMYSCORE\t$F[6]\n"; 
			$num_genes++;
			$genelength += ($F[4] - $F[3]) + 1;
		}
	}

	close(BED);
	close(GFF);

	return ($num_genes, sprintf("%1.0f",$genelength/$num_genes));
}


# Takes i) input BED intersect filename ii) output BED filename.
# Parses sorted BED intersect -wo output and writes to BED file 
# features (cDNA/transcripts) mapped on reference genome.
# Returns i) number of matched genes and ii) list of unmatched genes
# Note: able to parse cs::Z (minimap2) and cg::Z (wfmash) strings
# Note: takes first match of each cDNA only
# example input:
# 1 4848 20752 ONIVA01G00010 9999     + 1 3331 33993       + 6 26020714 26051403 29819 60 cs:Z::303*ag:30*ga... 15904
# 1 104921 116326 ONIVA01G00100 9999  + 1 103118 152580    + 1 1132 47408 45875 60 cs:Z::70*tc:...              11405
# Chr1 2903 10817 LOC_Os01g01010 9999 + Chr1 1000 10053455 + 1 1000 10053455 10052455 60 cs:Z::10052455          7914
# Chr1 2903 10817 LOC_Os01g01010 9999 + Chr1 896 705000    + 1 949 705000 704051 41 cg:Z:51=53I704000=           7914
# <--             (c)DNA/gene       --> <- (q)uery genome -> <-- (r)eference genome                         -->  ovlp
sub query2ref_coords {

	my ($infile, $outfile, $minqual, $minalnlen, $samestrand, $verbose) = @_;

	my ($cchr,$cstart,$cend,$cname,$cmatch,$cstrand);
	my ($qchr,$qstart,$qend,$cigartype);
	my ($WGAstrand,$rchr,$rstart,$rend);
	my ($rmatch,$rmapqual,$SAMPAFtag,$overlap,$done,$strand);
	my ($SAMqcoord,$SAMrcoord,$feat,$coordr);
	my ($deltaq,$deltar,$start_deltar,$end_deltar);
	my $num_matched = 0;
	my (%ref_coords,@unmatched,@segments);

	open(BED,"<",$infile) || die "# ERROR(query2ref_coords): cannot read $infile\n";

	open(OUTBED,">",$outfile) || die "# ERROR(query2ref_coords): cannot create $outfile\n";

	while(<BED>){

    #cDNA format
	#1 98773 99875 ONIVA01G00080.1 258 + 1 98032 101175 - 6 27346427 27348975 2375 60 cs:Z::29-ggt.. 15904
    #GFF/gene format
	#1  4848 20752 ONIVA01G00010 9999  + 1 3331 33993 + 6 26020714 26051403 29819 60 cs:Z::303*ag.. 15904
	#1 104921 116326 ONIVA01G00100 9999  + 1 103118 152580    + 1 1132 47408 45875 60 cs:Z::70*tc:... 11405

		($cchr,$cstart,$cend,$cname,$cmatch,$cstrand,
			$qchr,$qstart,$qend,$WGAstrand,
			$rchr,$rstart,$rend,$rmatch,$rmapqual,
			$SAMPAFtag,$overlap) = split(/\t/,$_);

		# skip mappings where query and ref segment on different strands if required
		next if($WGAstrand eq '-' && $samestrand == 1);

		# skip poor query-to-ref WGA scores
		next if($rmapqual < $minqual);

		# take 1st mapping only
		next if(defined($ref_coords{$cname}));

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
		if($WGAstrand eq '+') {
			$SAMrcoord = $rstart; 
			$ref_coords{$cname}{'start'} = $rstart;
		} else {
			$SAMrcoord = $rend;
			$ref_coords{$cname}{'end'} = $rend;
		}

		print "# $SAMqcoord $SAMrcoord $cname $num_matched\n" if($verbose > 1);

		# check CIGAR type
		if($SAMPAFtag =~ m/(c\w):Z:{1,2}(\S+)/){
			($cigartype,$SAMPAFtag) = ($1,$2);
		} else {
			print "# ERROR(query2ref_coords): unsupported CIGAR string $SAMPAFtag\n";
			return ($num_matched, @unmatched);	
		}

		# split CIGAR string into individual feature tags
		if($cigartype eq 'cs'){
			@segments = split(/:/,$SAMPAFtag);
		} else {
			while($SAMPAFtag =~ m/(\d+[MIDNSHP=X])/g){
				push(@segments,$1);
			}
		}

		# loop along features updating coords
		$done = 0;
		foreach $feat (@segments){

            ($deltaq,$deltar,$coordr) = _parseCIGARfeature($feat,$SAMqcoord,$SAMrcoord);

			## check if current position in alignment matches cDNA/gene coords  

			# start coords (end in - strand)
			if($SAMqcoord < $cstart && 
				($SAMqcoord + $deltaq) >= $cstart){

				# refine delta to match exactly the start (end for strand -)
				($deltaq,$deltar,$coordr) = _parseCIGARfeature($feat,$SAMqcoord,$SAMrcoord,$cstart);
				$start_deltar = -1;
				if($coordr > -1) {
					$start_deltar = $coordr - $SAMrcoord;
				}

				print ">$deltaq $deltar $start_deltar $cstart\n" if($verbose > 1);

				if($start_deltar < 0 || $start_deltar > $deltar){ $start_deltar = $deltar }

				if($WGAstrand eq '+') {
					$ref_coords{$cname}{'start'} = $SAMrcoord + $start_deltar;
				} else {
					$ref_coords{$cname}{'end'} = $SAMrcoord + $start_deltar;
				}	
			} 

			# end coords (start in -strand), actually copied out of loop
			if($SAMqcoord < $cend &&
                ($SAMqcoord + $deltaq) >= $cend){

				# refine delta to match exactly the end (start for strand -)
				($deltaq,$deltar,$coordr) = _parseCIGARfeature($feat,$SAMqcoord,$SAMrcoord,$cend);
				$end_deltar = -1;
                if($coordr > -1) {
                    $end_deltar = $coordr - $SAMrcoord;
                }

				print "<$deltaq $deltar $end_deltar $cend\n" if($verbose > 1);

				if($end_deltar < 0 || $end_deltar > $deltar){ $end_deltar = $deltar }

				$done = 1;
			}

			# update coords with deltas
	        $SAMqcoord += $deltaq;

			if($done) {
				$deltar = $end_deltar;
			}

	        if($WGAstrand eq '+') {
	            $SAMrcoord += $deltar;
	        } else {
	            $SAMrcoord -= $deltar;
	        }

			print "$SAMqcoord $SAMrcoord $feat $deltaq $deltar\n" if($verbose > 1);

			last if($done);
		}

		if($WGAstrand eq '+') {
		    $ref_coords{$cname}{'end'} = $SAMrcoord;
		} else {
			$ref_coords{$cname}{'start'} = $SAMrcoord;
		}

		# work out strand of mapped cDNA/transcript/gene
		$strand = $cstrand; # original strand 
		if($WGAstrand eq '-') {
			if($cstrand eq '+'){ 
				$strand = '-';
			} else {
				$strand = '+';
			}
		}

		# print coords in ref frame only if segment is long enough
		if(1+$ref_coords{$cname}{'end'}-$ref_coords{$cname}{'start'} >= $minalnlen){
			printf(OUTBED "%s\t%d\t%d\t%s\t%d\t%s\n",
				$ref_coords{$cname}{'chr'},
				$ref_coords{$cname}{'start'},
				$ref_coords{$cname}{'end'},
				$cname,
				1+$ref_coords{$cname}{'end'}-$ref_coords{$cname}{'start'},			
				$strand);
			$num_matched++;
		} else {
			# shell genes cannot be matched in ref genome
			push(@unmatched,$cname);
		}
	}

	close(OUTBED);
	close(BED);

	return ($num_matched, @unmatched);
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

	my ($feat,$Qstart,$Rstart,$opt_query_coord) = @_;

	my ($deltaq,$deltar,$totsnps,$rcoord) = (0,0,0,-1);
	my ($delta,$offset,$res,$qcoord,$op);

	# example CS string (coords are 0-based):
	# cs:Z::303*ag:30*ga:32*ga:27+ctattcta*ag*ca:20*ag:3*ga:18*tc:39*ga:76-tc
	# first features: 303*ag 30*ga .. 27+ctattcta .. 76-tc
	
	# example CG string:
	# cg:Z:25I235=1X57=1X33=1X20=1X47=
	# first features: 25I 235= 1X ..

	# CG type, see https://samtools.github.io/hts-specs/SAMv1.pdf
	if($feat =~ m/(\d+)([MIDNSHP=X])/){

		($delta,$op) = ($1,$2);

		if($op eq 'M' || $op eq '=' || $op eq 'X') {
			$deltar += $delta;
            $deltaq += $delta;

            # look up query coord (optional)
            if(defined($opt_query_coord)) {
                $offset = $opt_query_coord - $Qstart;
                if($offset <= $delta) {
                    $rcoord = $Rstart + $offset;
                }
            }

		} elsif($op eq 'I') {
			
			# look up query coord (optional)
			if(defined($opt_query_coord) && $rcoord == -1) {
				$offset = $opt_query_coord - $Qstart;
				if($offset <= $delta + $deltaq) {
					$rcoord = $Rstart + $deltar;
				}
			}

			$deltaq += $delta;

		} elsif($op eq 'D') {

			# look up query coord (optional)
			if(defined($opt_query_coord) && $rcoord == -1) {		
				$offset = $opt_query_coord - $Qstart;
				if($offset <= $deltaq) {	
					$rcoord = $Rstart + $deltar + $offset;
				}
			}

			$deltar += $delta;

		} else {
			print "# ERROR(_parseCIGARfeature): unsupported CIGAR operation $feat\n";
		}

	} else { # CS type

		# identical segment
		if($feat =~ m/^(\d+)/) {
			$delta = $1;	
			$deltar += $delta;
			$deltaq += $delta;

			# look up query coord (optional)
			if(defined($opt_query_coord)) {
				$offset = $opt_query_coord - $Qstart;
				if($offset <= $delta) {	
					$rcoord = $Rstart + $offset;
				}
			}
		}

		# insertion (+) / deletion (-)
		if($feat =~ m/([\+\-])([A-Za-z]+)/) {
			$delta = length($2);
                
			if($1 eq '+') {
                    
				# look up query coord (optional)
				if(defined($opt_query_coord) && $rcoord == -1) {
					$offset = $opt_query_coord - $Qstart;
					if($offset <= $delta + $deltaq) { 
						$rcoord	= $Rstart + $deltar;
					}
				}
					
				$deltaq += $delta;

			} else {
			
				# look up query coord (optional)
				if(defined($opt_query_coord) && $rcoord == -1) {	
					$offset = $opt_query_coord - $Qstart;
					if($offset <= $deltaq) {	
						$rcoord = $Rstart + $deltar + $offset;
					}	
				}
	
				$deltar += $delta;
			}
		}

		# SNPs
		while($feat =~ m/\*\w{2}/g){

			# look up query coord (optional)
			if(defined($opt_query_coord) && $rcoord == -1) {
				$qcoord = $Qstart + $totsnps + $deltaq;
				if($qcoord == $opt_query_coord && $rcoord == -1) {
					$rcoord = $Rstart + $totsnps + $deltar;
				}
			}

			$totsnps++;
		}
		$deltar += $totsnps;
		$deltaq += $totsnps;
	}


	return ($deltaq, $deltar, $rcoord);
}

# Takes i) input BED intersect filename ii) output TSV filename and
# returns i) number of collinear gene pairs.
# Parse bedtools intersect file and produces Ensembl Compara-like TSV file.
# Example input with explicit strand:
# 1 2983 10815 Os01g0100100 9999 + 1 2890 12378 ONIVA01G00100 9489 + 7832
# Expected TSV output:
# gene_stable_id protein_stable_id species overlap homology_type homology_gene_stable_id 
# homology_protein_stable_id homology_species overlap dn ds goc_score wga_coverage 
# is_high_confidence coordinates
sub bed2compara {
	
	my ($infile, $TSVfile, $sp1, $sp2, $noheader, $workout_gene_names) = @_;
	my ($gene1,$gene2,$coords);
	my $num_pairs = 0;

	# parse input and produce TSV output
	open(BEDINT,"<",$infile) || die "# ERROR(bed2compara): cannot read $infile\n";

	if($noheader) {
		open(TSV,">>",$TSVfile) || die "# ERROR(bed2compara): cannot re-open $TSVfile\n";
	} else {
		open(TSV,">",$TSVfile) || die "# ERROR(bed2compara): cannot create $TSVfile\n";
	
        print TSV "gene_stable_id\tprotein_stable_id\tspecies\toverlap\thomology_type\thomology_gene_stable_id\t".
            "homology_protein_stable_id\thomology_species\toverlap\tdn\tds\tgoc_score\twga_coverage\t".
            "is_high_confidence\tcoordinates\n";
    }

	while(<BEDINT>){
		my @data = split(/\t/,$_);

		# concat genome/graph coords
		$coords = "$data[0]:$data[1]-$data[2]";   # 1st genome
		$coords .= ";$data[6]:$data[7]-$data[8]"; # 2nd genome

		# format gene names
		$gene1 = $data[3];
		$gene2 = $data[9];

		if($workout_gene_names) {
			$gene1 =~ s/[\.-]\d+$//; 
    		$gene2 =~ s/[\.-]\d+$//;
        	$gene1 =~ s/t/g/; 
			$gene2 =~ s/t/g/;
    	}

		printf(TSV "%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%d\tNULL\tNULL\tNULL\t%1.2f\t%d\t%s\n",
			$gene1,
			$data[3],
			$sp1,
			$data[12],
			'ortholog_collinear',
			$gene2,
        	$data[9],
			$sp2,
			$data[12],
			100, # max WGA score
			1, # high confidence
			$coords);

		$num_pairs++;
	}

	close(TSV);
	close(BEDINT);

	return $num_pairs;
}	



