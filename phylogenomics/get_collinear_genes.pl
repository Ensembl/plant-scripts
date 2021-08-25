#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

# Takes two FASTA files with genome sequences and 2 matching GFF files with annotated gene models.
# Produces TSV file with pairs of collinear genes in a format compatible with Ensembl Compara.
# Uses external software:
# minimap2 [https://academic.oup.com/bioinformatics/article/34/18/3094/4994778]
# bedtools [https://academic.oup.com/bioinformatics/article/26/6/841/244688]

#perl get_collinear_genes.pl -sp1 oryza_sativa -fa1 Oryza_sativa.IRGSP-1.0.dna.toplevel.fa -gf1 Oryza_sativa.IRGSP-1.0.51.gff3 -sp2 oryza_nivara -fa2 Oryza_nivara.Oryza_nivara_v1.0.dna.toplevel.fa -gf2 Oryza_nivara.Oryza_nivara_v1.0.51.gff3 -o kk -r


my $MINIMAP2EXE = 'minimap2'; # tested with 2.17-r941
my $MINIMAPPARS = '--cs -x asm20';
my $BEDTOOLSEXE = 'bedtools'; # tested with v2.30.0
my $BEDINTSCPAR = '-wo -f 0.5 -F 0.5 -e';
my $BLASTNEXE   = 'blastn';   # tested with 2.2.30+

our $SORTLIMITRAM = "500M"; # buffer size
our $SORTBIN      = "sort --buffer-size=$SORTLIMITRAM";

my $DUMMYSCORE = 9999;

# while parsing PAF
my $MINQUAL    = 60; 
my $MINALNLEN  = 100;
my $SAMESTRAND = 1;
my $VERBOSE    = 0;

# Work out gene names from transcripts':
# 1) remove suffix after . or -
# 2) t's in transcript name to be replaced with g' in gene names
# Example: 10t0100300.1 -> Os10g0100300
my $TRANSCRIPT2GENE = 1;

my ( $help, $do_sequence_check, $reuse, $noheader, $outfile ) = (0,0,0,0,'');
my ( $sp1, $fasta1, $gff1, $sp2, $fasta2, $gff2  );

GetOptions(
	"help|?"         => \$help,
	"sp1|species1=s" => \$sp1,
	"fa1|fasta1=s"   => \$fasta1,
	"gf1|gff1=s"     => \$gff1,
	"sp2|species2=s" => \$sp2,
	"fa2|fasta2=s"   => \$fasta2,
	"gf2|gff2=s"     => \$gff2,
	"o|outfile=s"    => \$outfile,
	"c|check"        => \$do_sequence_check,
	"r|reuse"        => \$reuse,
	"n|noheader"     => \$noheader
) || help_message();

sub help_message {
	print "\nusage: $0 [options]\n\n"
		. "-sp1 binomial/trinomial species name    (required, example: -sp1 oryza_sativa)\n"
		. "-fa1 genome FASTA filename              (required, example: -fa1 oryza_sativa.fna)\n"
		. "-gf1 GFF filename                       (required, example: -gf1 oryza_sativa.gff)\n"
		. "-sp2 binomial/trinomial species name    (required, example: -sp2 oryza_nivara)\n"
		. "-fa2 genome FASTA filename              (required, example: -fa2 oryza_nivara.fna)\n"
		. "-gf2 GFF filename                       (required, example: -gf2 oryza_nivara.gff)\n"
		. "-o output TSV filename                  (required, example: -o oryza_sativa.oryza_nivara.tsv)\n"
		. "-c compare sequences of collinear genes (optional)\n"
		. "-n no header in TSV file                (optional, useful to concatenate files)\n"
		. "-r re-use previous minimap2 results     (optional)\n";
}

if($help || 
	(!$sp1 || !$fasta1 || !$gff1 || 
		!$sp2 || !$fasta2 || !$gff2 || !$outfile)){ 
	help_message();
	exit(0);
}

print "\n# $0 -sp1 $sp1 -fa1 $fasta1 -gf1 $gff1 ".
	"-sp2 $sp2 -fa2 $fasta2 -gf2 $gff2 -o $outfile -c $do_sequence_check -r $reuse\n\n";

## 1) align genome1 vs genome2 with minimap2 (WGA)
## Note: no masking required, see https://github.com/lh3/minimap2/issues/654

my $PAFfile = "$sp2.$sp1.minimap.paf";

print "# computing pairwise genome alignment with minimap2\n\n";

if($reuse && -s $PAFfile){
	print "# re-using $PAFfile\n";
} else {
	system("$MINIMAP2EXE $MINIMAPPARS $fasta1 $fasta2 -o $PAFfile 2>&1");
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

## 2) produce BED-like file of sp2-to-sp1 coords 10 columns
## Note: $F[4] in PAF conveys whether query & ref are on the same strand or not

my $wgaBEDfile = "$sp2.$sp1.minimap.bed";

open(PAF,"<",$PAFfile) || die "# ERROR: cannot read $PAFfile\n";
open(BED,">",$wgaBEDfile) || die "# ERROR: cannot create $wgaBEDfile\n";

while(<PAF>){
	my @F = split(/\t/,$_);
	print BED "$F[0]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[7]\t$F[8]\t$F[9]\t$F[11]\t$F[$#F]";
}

close(BED);
close(PAF);

## 3) Parse GFFs and produce BED files with gene coords
## Note: add dummy match score=9999

my $geneBEDfile1 = "$sp1.gene.bed";
my $geneBEDfile2 = "$sp2.gene.bed";

my $num_genes1 = parse_genes_GFF($gff1,$geneBEDfile1);
printf("# %d genes parsed in %s\n",$num_genes1,$gff1);

my $num_genes2 = parse_genes_GFF($gff2,$geneBEDfile2);
printf("# %d genes parsed in %s\n",$num_genes2,$gff2);

## 4) intersect gene positions with WGA, sort by gene > cDNA ovlp > genomic matches

my $sp2wgaBEDfile = "$sp2.gene.$sp2.$sp1.minimap.intersect.bed";

system("$BEDTOOLSEXE intersect -a $geneBEDfile2 -b $wgaBEDfile $BEDINTSCPAR | ".
	"$SORTBIN -k4,4 -k5,5nr -k14,14nr > $sp2wgaBEDfile");
if($? != 0){
	die "# ERROR: failed running bedtools\n";
}
elsif(!-s $sp2wgaBEDfile){
	die "# ERROR: failed generating $sp2wgaBEDfile file (bedtools)\n";
}

my $geneBEDfile2mapped = "$sp2.gene.mapped.bed";

my ($num_matched, @unmatched) = 
	query2ref_coords($sp2wgaBEDfile, $geneBEDfile2mapped, 
		$MINQUAL, $MINALNLEN, $SAMESTRAND, $VERBOSE);

printf("# %d genes mapped in %s (%d unmapped)\n",
	$num_matched,$geneBEDfile2mapped,scalar(@unmatched));

## 5) produce list of pairs of collinear genes

my $gene_intersectBEDfile = "$sp1.$sp2.gene.intersect.bed";

system("$BEDTOOLSEXE intersect -a $geneBEDfile1 -b $geneBEDfile2mapped $BEDINTSCPAR -s | ".
    "$SORTBIN -k4,4 -k13,13nr > $gene_intersectBEDfile");
if($? != 0){
    die "# ERROR: failed running bedtools\n";
}
elsif(!-s $gene_intersectBEDfile){
    die "# ERROR: failed generating $gene_intersectBEDfile file (bedtools)\n";
}

my $TSVfile = "Minimap.homologies.$sp1.$sp2.tsv";
my $num_pairs = bed2compara($gene_intersectBEDfile, $TSVfile, $sp1, $sp2,
	$noheader, $TRANSCRIPT2GENE);

print "# TSV file: $TSVfile\n";


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
# returns number of genes parsed.
# Uses global $DUMMYSCORE as dummy scores, as opposed to cases where
# genes/transcript are actually mapped to genome
sub parse_genes_GFF {

	my ($gff_file, $outBEDfile) = @_;
	my $num_genes = 0;
	my $geneid;

	open(GFF,"<",$gff_file) || die "# ERROR(parse_genes_GFF): cannot read $gff_file\n";
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
		}
	}

	close(BED);
	close(GFF);

	return $num_genes;
}


# Takes i) input BED intersect filename ii) output BED filename and
# returns i) number of matched genes and ii) list of unmatched genes
# Parses sorted BED intersect -wo output and produces BED file with cDNA/transcripts mapped on ref
# Note: takes first match of each cDNA only
# example input:
# 7	17764202 17769979 ONIVA07G18210.2 1684 7 17758528 17798398 + 7 23334256 23377175 32397 60 cs:Z::23-cc:...    5777
# <--             (c)DNA/gene          --> <- (q)uery genome-> <-- (r)eference genome                        --> ovlp
sub query2ref_coords {

	my ($infile, $outfile, $minqual, $minalnlen, $samestrand, $verbose) = @_;

	my ($cchr,$cstart,$cend,$cname,$cmatch,$cstrand);
	my ($qchr,$qstart,$qend);
	my ($WGAstrand,$rchr,$rstart,$rend);
	my ($rmatch,$rmapqual,$SAMPAFtag,$overlap);
	my ($qoffset,$roffset,$length,$done,$strand);
	my ($SAMqcoord,$SAMrcoord,$feat,$deltaq,$deltar);
	my $num_matched = 0;
	my (%ref_coords,@unmatched);

	open(BED,"<",$infile) || die "# ERROR(query2ref_coords): cannot read $infile\n";

	open(OUTBED,">",$outfile) || die "# ERROR(query2ref_coords): cannot create $outfile\n";

	while(<BED>){

    #cDNA format
	#1 98773 99875 ONIVA01G00080.1 258 + 1 98032 101175 - 6 27346427 27348975 2375 60 cs:Z::29-ggt.. 15904
    #GFF/gene format
	#1  4848 20752 ONIVA01G00010 9999  + 1 3331 33993 + 6 26020714 26051403 29819 60 cs:Z::303*ag.. 15904
	
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

		#next if($cname ne 'ONIVA06G10510.1'); # debug
		#print;

		# estimate offset of aligned cDNA in genomic coords of query 
		$qoffset = $cstart - $qstart;
		if($qoffset < 0){ $qoffset = 0 } # genomic segment contains 5' half of cDNA only
		$length = $cend - $cstart;

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

		print "$SAMqcoord $SAMrcoord\n" if($verbose > 1);

		$SAMPAFtag =~ s/cs:Z:://;
		$done = 0;
		foreach $feat (split(/:/,$SAMPAFtag)){

			$deltaq = $deltar = 0;
		
			# identical segment
			if($feat =~ m/(\d+)/){ 
				$deltar += $1; 
				$deltaq += $1;
			} 
		
			# insertion (+) / deletion (-)
			if($feat =~ m/([\+\-])([A-Za-z]+)/){
				if($1 eq '+'){
					$deltaq += length($2);
				} else {
					$deltar += length($2);
				}
			}

			# SNPs
			while($feat =~ m/\*\w{2}/g){ 
				$deltar++; 
				$deltaq++;
        	}

			## check if current position in alignment matches cDNA coords  

			# start coords (end in - strand)
			if($SAMqcoord < ($qstart + $qoffset) && 
				($SAMqcoord + $deltaq) >= ($qstart + $qoffset)){

				if($WGAstrand eq '+') {
					$ref_coords{$cname}{'start'} = $SAMrcoord;
				} else {
					$ref_coords{$cname}{'end'} = $SAMrcoord;
				}	
				print ">$deltaq" if($verbose > 1);
			} 

			# end coords (start in -strand), actually copied out of loop
			if($SAMqcoord < ($qstart + $qoffset + $length) &&
	            $SAMqcoord + $deltaq >= ($qstart + $qoffset + $length)){	

				$done = 1;
				print "<$deltaq\n" if($verbose > 1);
			}

			# update coords with deltas
	        $SAMqcoord += $deltaq;
	        if($WGAstrand eq '+') {
	            $SAMrcoord += $deltar;
	        } else {
	            $SAMrcoord -= $deltar;
	        }

			print "$SAMqcoord $SAMrcoord $feat $deltaq ".
				"($qstart + $qoffset + $length) $deltar\n" if($verbose > 1);

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

	open(TSV,">",$TSVfile) || die "# ERROR(bed2compara): cannot create $TSVfile\n";

	# print header
    if(!$noheader){
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



