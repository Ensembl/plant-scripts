#!/bin/env perl
use strict;
use warnings;

# fix a GFF with no transcripts, only with gene and CDS features
# by Bruno Contreras Moreira EMBL-EBI 2018

if(!$ARGV[0]){ die "# usage: $0 <GFF>\n" }

my ($tparent, $tname, $CDS, $transcript);

open(GFF,"<",$ARGV[0]) || die "# ERROR: cannot read $ARGV[0]\n";
while(<GFF>){

	#X15901.1	EMBL	gene	82	1143	.	-	.	ID=gene-psbA;Name=psbA;gbkey=Gene;gene=psbA;gene_biotype=protein_coding
	#X15901.1	EMBL	CDS	82	1143	.	-	0	ID=cds-CAA34007.1;Parent=gene-psbA;Dbxref=GOA:P0C434,InterPro:IPR000484,InterPro:IPR005867,UniProtKB/Swiss-Prot:P0C434,NCBI_GP:CAA34007.1;Name=CAA34007.1;gbkey=CDS;gene=psbA;product=PSII...

	my @col = split(/\t/,$_);

	if(defined($col[2]) && $col[2] eq 'CDS'){
		
		$CDS = $_;
		$transcript = $_;

		if($col[8] =~ m/ID=([^;]+);Parent=([^;]+)/){

			# create a new transcript
			$tparent = $2;
			$tname = $tparent;
			$tname =~ s/gene/transcript/;	
			$transcript =~ s/\tCDS\t/\ttranscript\t/;	
			$transcript =~ s/ID=[^;]+/ID=$tname/;
			$transcript =~ s/Parent=[^;]+;.*/Parent=$tparent/;
			print $transcript;
	
			# now print the updated CDS
			$CDS =~ s/$tparent/$tname/;
			print $CDS;
		}
	
	} else {
		print;
	}

}
close(GFF);
