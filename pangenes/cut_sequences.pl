#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

# Takes a GGF and a FASTA file and produces FASTA files with CDS nucl & pep sequences of the 1st transcript found

# Uses external software:
# gffread [https://f1000research.com/articles/9-304/v2]

#perl cut_sequences.pl -sp oryza_sativa -fa Oryza_sativa.IRGSP-1.0.dna.toplevel.fa -gf Oryza_sativa.IRGSP-1.0.51.gff3

my $GFFREADEXE = 'gffread'; # v0.12.7

my ( $help, $gffreadpath, $sp1, $fasta1, $gff1, $tname, $nored ) = (0);

GetOptions(
	"help|?"       => \$help,
	"sp|species=s" => \$sp1,
	"fa|fasta=s"   => \$fasta1,
	"gf|gff=s"     => \$gff1,
	"nr|n"         => \$nored,
	"p|path=s"     => \$gffreadpath
) || help_message();

sub help_message {
	print "\nusage: $0 [options]\n\n"
		. "-sp binomial/trinomial species name (required, example: -sp oryza_sativa, used to name outfiles)\n"
		. "-fa genome FASTA filename           (required, example: -fa oryza_sativa.fna)\n"
		. "-gf GFF filename                    (required, example: -gf oryza_sativa.RAPDB.gff)\n"
		. "-nr remove redundancy in seq names  (optional, ie 'gene:ONIVA01G00100')\n"
		. "-p  path to gffread binary          (optional, default: $GFFREADEXE)\n\n"
}

if($help || 
	(!$sp1 || !$fasta1 || !$gff1)){ 
	help_message();
	exit(0);
}  

if(!-s $fasta1 || !-s $gff1){
	print "# ERROR: please make sure all input files exist\n";
    exit(0);
} 

if(!$gffreadpath){
	$gffreadpath = $GFFREADEXE;
}

# set output filenames
my $cdnafile = "$sp1.cdna.fna";
my $cdsfile  = "$sp1.cds.fna";
my $pepfile  = "$sp1.cds.faa";

print "\n# $0 -sp $sp1 -fa $fasta1 -gf $gff1 -nr $nored -path $gffreadpath\n\n";

my ($ref_names, $ref_coords) = parse_genes($gff1);

my $num_cdna = parse_gffread($gffreadpath,$fasta1,$gff1,$cdnafile,
	'cdna',$nored,$ref_names,$ref_coords);
my $num_cds  = parse_gffread($gffreadpath,$fasta1,$gff1,$cdsfile,
	'cds',$nored,$ref_names,$ref_coords);
my $num_pep  = parse_gffread($gffreadpath,$fasta1,$gff1,$pepfile,
	'pep',$nored,$ref_names,$ref_coords);

print "# $cdnafile n=$num_cdna\n";
print "# $cdsfile n=$num_cds\n";
print "# $pepfile n=$num_pep\n";

###############################

# Runs gffread, parses its stdout and saves output in FASTA FILE.
# Returns number of sequences printed out.
sub parse_gffread {

	my ($gffreadexe,$fasta_file,$gff_file,$outfile,
		$seqtype,$remove_red,$ref_tr2gene,$ref_tr2coords) = @_;

	my ($params,$mrnaid,$geneid,$coords);

	if($seqtype eq 'cds'){
		$params = '-x - ';
	} elsif($seqtype eq 'pep'){
        $params = '-y - ';
	} else {
		$params = '-w - '; # cDNA, default
	}

	my $num_seqs = 0;
	open(OUT,">",$outfile) || die "# ERROR(parse_gffread): cannot create $pepfile\n";
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

	return $num_seqs
}

# Parses gene names as parent IDs of transcripts.
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

