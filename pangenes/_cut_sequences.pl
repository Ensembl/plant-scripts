#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

# Takes a GFF & FASTA files and produces FASTA files with 
# CDS nucl & pep sequences of the 1st transcript found
#
# Uses external software: gffread [https://f1000research.com/articles/9-304/v2]

# Copyright [2021-22] 
# EMBL-European Bioinformatics Institute & Estacion Experimental de Aula Dei-CSIC

# perl _cut_sequences.pl -sp oryza_sativa -fa Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
#   -gf Oryza_sativa.IRGSP-1.0.51.gff3

my $GFFREADEXE = 'gffread'; # v0.12.7

my ( $help, $nored, $gffreadpath, $sp1, $fasta1) = (0, 0);
my ( $minlen, $gff1, $tname, $outpath ) = (0);

GetOptions(
	"help|?"       => \$help,
	"sp|species=s" => \$sp1,
	"fa|fasta=s"   => \$fasta1,
	"gf|gff=s"     => \$gff1,
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
		. "-l  min length (bp) of features     (optional, example: -l 100)\n"
		. "-nr remove redundancy in seq names  (optional, ie 'gene:ONIVA01G00100')\n"
		. "-p  path to gffread binary          (optional, default: $GFFREADEXE)\n"
		. "-o  path to output folder           (optional, default current folder)\n\n"
}

if($help || 
	(!$sp1 || !$fasta1 || !$gff1)){ 
	help_message();
	exit(0);
}  

if(!-s $fasta1 || !-s $gff1){
	print "# ERROR: please make sure all input files exist\n";
    exit(0);
} 

if(!$gffreadpath){
	$gffreadpath = $GFFREADEXE;
}

if($minlen < 1){ 
	$minlen = 0 
}

# set output filenames
my $cdnafile = "$sp1.cdna.fna";
my $cdsfile  = "$sp1.cds.fna";
my $pepfile  = "$sp1.cds.faa";
if($outpath) {
	$cdnafile = "$outpath/$sp1.cdna.fna";
	$cdsfile  = "$outpath/$sp1.cds.fna";
	$pepfile  = "$outpath/$sp1.cds.faa";
}

print "\n# $0 -sp $sp1 -fa $fasta1 -gf $gff1 -l $minlen -nr $nored -path $gffreadpath\n\n";

my ($ref_names, $ref_coords) = parse_genes($gff1);

my $num_cdna = parse_gffread($gffreadpath,$fasta1,$gff1,$cdnafile,
	'cdna',$minlen,$nored,$ref_names,$ref_coords);
my $num_cds  = parse_gffread($gffreadpath,$fasta1,$gff1,$cdsfile,
	'cds',$minlen,$nored,$ref_names,$ref_coords);
my $num_pep  = parse_gffread($gffreadpath,$fasta1,$gff1,$pepfile,
	'pep',$minlen,$nored,$ref_names,$ref_coords);

if(scalar(keys(%$ref_names))) {
	printf("# genes n=%d\n",scalar(keys(%$ref_names)));
} else {
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

# Parses gene names as parent IDs of transcripts.
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

			$coord = "$F[0]:$F[3]-$F[4]";

			$names{$mrnaid} = $geneid;
			$coords{$mrnaid} = $coord;
		}
	}
	close(GFF);

	return (\%names,\%coords);
}

