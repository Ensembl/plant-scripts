#!/usr/bin/env perl
use strict;
use warnings;

# Takes two FASTA files and checks whether they contain the same headers & sequences
# I use to compare assemblies from ENA and other sources, such as JGI, to see whether
# GFF files map to the same sequences
# Bruno Contreras Moreira EMBL-EBI 2019 

if(!$ARGV[1]){
	die "# usage: $0 <FASTA file1> <FASTA file2>\n";
}

my ($header,$seq,$header2,$seq2);
my (@order1,%FASTA1,%FASTA2);

# read and store 1st FASTA file
open(FASTA1,"<",$ARGV[0]) || die "# ERROR: cannot read $ARGV[0]\n";
while(<FASTA1>){
	if(/^>(\S+)/){ 
		$header = $1;
	}
	elsif(/^(\S+)/){
		chomp;
		$seq = $1;
		$FASTA1{$header} .= $seq;		
		push(@order1,$header);
	}
}
close(FASTA1);

# read and store 2nd FASTA file
open(FASTA2,"<",$ARGV[1]) || die "# ERROR: cannot read $ARGV[1]\n";
while(<FASTA2>){
        if(/^>(\S+)/){ 
		$header = $1; 
	}
        elsif(/^(\S+)/){
                chomp;
                $seq = $1;
                $FASTA2{$header} .= $seq;
        }
}
close(FASTA2);

# compare headers
foreach $header (@order1) {
	
	# check header is present in 2nd FASTA
	if(defined($FASTA2{$header})) {
		
		# check sequence is the same
		if($FASTA1{$header} ne $FASTA2{$header}){
			printf("# sequence $header is different in $ARGV[2] (%d,%d)\n",
				length($FASTA1{$header}),length($FASTA2{$header}));
		} else {
			print "# sequence $header is OK\n";
		}		
	} else {
		print "# cannot find sequence $header in $ARGV[1]\n"
	}
}
