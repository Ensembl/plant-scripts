#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

# Takes a GGF and a FASTA file and produces FASTA files with CDS nucl & pep sequences of the 1st transcript found

# Uses external software:
# gffread [https://f1000research.com/articles/9-304/v2]

#perl cut_sequences.pl -sp oryza_sativa -fa Oryza_sativa.IRGSP-1.0.dna.toplevel.fa -gf Oryza_sativa.IRGSP-1.0.51.gff3

my $GFFREADEXE = 'gffread'; # v0.12.7

my ( $help, $gffreadpath, $sp1, $fasta1, $gff1, $tname ) = (0);

GetOptions(
	"help|?"       => \$help,
	"sp|species=s" => \$sp1,
	"fa|fasta=s"   => \$fasta1,
	"gf|gff=s"     => \$gff1,
	"p|path=s"     => \$gffreadpath
) || help_message();

sub help_message {
	print "\nusage: $0 [options]\n\n"
		. "-sp binomial/trinomial species name    (required, example: -sp oryza_sativa, used to name outfiles)\n"
		. "-fa genome FASTA [.gz] filename        (required, example: -fa oryza_sativa.fna)\n"
		. "-gf GFF [.gz] filename                 (required, example: -gf oryza_sativa.RAPDB.gff)\n"
		. "-p  path to gffread binary             (optional, default: $GFFREADEXE)\n\n"
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

print "\n# $0 -sp $sp1 -fa $fasta1 -gf $gff1 -path $gffreadpath\n\n";

my $ref_genes = parse_genenames($gff1);

my $num_cdna = parse_gffread($gffreadpath,$fasta1,$gff1,$cdnafile,'cdna',$ref_genes);
my $num_cds  = parse_gffread($gffreadpath,$fasta1,$gff1,$cdsfile,'cds',$ref_genes);
my $num_pep  = parse_gffread($gffreadpath,$fasta1,$gff1,$pepfile,'pep',$ref_genes);

print "# $cdnafile n=$num_cdna\n";
print "# $cdsfile n=$num_cds\n";
print "# $pepfile n=$num_pep\n";

###############################

# Runs gffread, parses its stdout and saves output in FASTA FILE.
# Returns number of sequences printed out.
sub parse_gffread {

	my ($gffreadexe,$fasta_file,$gff_file,$outfile,$seqtype,$ref_tr2gene) = @_;
	my ($params,$mrnaid,$geneid);

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

        	print OUT ">$mrnaid $geneid [$sp1]\n";
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
# Return ref to hash mapping transcript ID -> gene ID
sub parse_genenames {

	my ($gff_file) = @_;

	my ($mrnaid,$geneid,%names);

	open(GFF,"<",$gff_file) || die "# ERROR(parse_genenames): cannot read $gff_file\n";
	while(<GFF>){
		my @F = split(/\t/,$_);

		next if(scalar(@F)<9 || ($F[2] ne "mRNA" && $F[2] ne "transcript"));

		# take only genes where ID can be parsed
		if($F[8] =~ /ID=([^;]+).*?Parent=([^;]+)/){ 

			$mrnaid = $1;
			#$mrnaid =~ s/transcript://; # remove redundant bits
			$geneid = $2;
			#$geneid =~ s/gene://; # remove redundant bits

			$names{$mrnaid} = $geneid;
		}
	}
	close(GFF);

	return \%names;
}

