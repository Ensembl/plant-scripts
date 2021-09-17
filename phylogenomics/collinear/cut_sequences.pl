#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

# Takes a GGF and a FASTA file and produces FASTA files with CDS nucl & pep sequences of the 1st transcript found

# Uses external software:
# gffread [https://f1000research.com/articles/9-304/v2]

#perl cut_sequences.pl -sp1 oryza_sativa -fa1 Oryza_sativa.IRGSP-1.0.dna.toplevel.fa -gf1 Oryza_sativa.IRGSP-1.0.51.gff3

my $GFFREADEXE = 'gffread'; # v0.12.7

my ( $help, $gffreadpath, $sp1, $fasta1, $gff1 ) = (0);

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

open(CDNA,">",$cdnafile) || die "# ERROR: cannot create $cdnafile\n";
open(GFFREAD,"$gffreadpath -w - -g $fasta1 $gff1 |") ||
	die "# ERROR: cannot run $gffreadpath\n";
while(<GFFREAD>){
	if(/^>(\S+)/){
		print CDNA ">$1 [$sp1]\n"	
	} else {
		print CDNA;
	}
}
close(GFFREAD);
close(CDNA);

open(CDS,">",$cdsfile) || die "# ERROR: cannot create $cdsfile\n";
open(GFFREAD,"$gffreadpath -x - -g $fasta1 $gff1 |") ||
    die "# ERROR: cannot run $gffreadpath\n";
while(<GFFREAD>){
    if(/^>(\S+)/){
        print CDS ">$1 [$sp1]\n"
    } else {
        print CDS;
    }
}
close(GFFREAD);
close(CDS);

open(PEP,">",$pepfile) || die "# ERROR: cannot create $pepfile\n";
open(GFFREAD,"$gffreadpath -y - -g $fasta1 $gff1 |") ||
    die "# ERROR: cannot run $gffreadpath\n";
while(<GFFREAD>){
    if(/^>(\S+)/){
        print PEP ">$1 [$sp1]\n"
    } else {
        print PEP;
    }
}
close(GFFREAD);
close(PEP);

print "# created $cdnafile $cdsfile $pepfile\n\n";
