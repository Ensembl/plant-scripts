#!/usr/bin/env perl

# takes a FASTA file with contigs from BIGD (non INSDC archive) and attempts to match them
# to sequences in a ENA sequence report asumming order is conserved. Uses sequence length as a control.
# Bruno Contreras Moreira 2019

use strict;
use warnings;

if(!$ARGV[1]){ 
	die "# usage: $0 <non INDSC FASTA> <ENA sequence report file>\n";
}

my ($fasta_gz,$ENA_repfile) = @ARGV; 
my ($synname,$enaname,$len,%length,%syn);


open(SEQRP,"<",$ENA_repfile) || die "# ERROR: cannot read $ENA_repfile\n";
while(<SEQRP>){
	#accession       sequence-name   sequence-length sequence-role   replicon-name   replicon-type   assembly-unit
	#VTWI01000001.1  Contig0 1029    unplaced-scaffold       NA      NA      Primary Assembly
	next if(/^accession/);
	my @data = split(/\t/,$_);
	$syn{ $data[0] } = $data[1];
	$length{ $data[1] } = $data[2];
}
close(SEQRP);

open(FASTAGZ,"zcat $fasta_gz |") || die "# ERROR: cannot read $fasta_gz\n";
while(<FASTAGZ>){
	#>GWHAAAP00000001.1      OriSeqID=Contig0        Len=1029
	#CGAATCCACTGGTTTCACCTCGAATCTTCACACTTTGTTTAATCAAAAATCACTCA...
	if(/^>(\S+)\s+OriSeqID=(\S+)\s+Len=(\d+)/){
		($synname,$enaname,$len) = ($1,$2,$3);
		if($len != $length{ $enaname }){
			warn "# ERROR: length of $1 in $fasta_gz does not match that in $ENA_repfile ($len != $length{ $enaname })\n";
			print "$1\t$2\t$len\tmismatch\n";
		}
		else{ 
			print "$1\t$2\t$len\tOK\n";
		}
	}
}
close(FASTAGZ);

