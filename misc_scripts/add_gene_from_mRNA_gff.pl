#!/bin/env perl
use strict;
use warnings;

# add genes to a GFF file containing only mRNA features
# by Bruno Contreras Moreira EMBL-EBI 2019

if(!$ARGV[0]){ die "# usage: $0 <GFF>\n" }

my ($tparent, $tname, $gene, $transcript);

open(GFF,"<",$ARGV[0]) || die "# ERROR: cannot read $ARGV[0]\n";
while(<GFF>){

	##gff-version 3
	#Chr01	gt4sp_v3	mRNA	13661	18457	.	+	.	ID=itb01g00040.t1;Name=itb01g00040.t1;Parent=itb01g00040
	#Chr01	gt4sp_v3	exon	13661	13835	.	+	.	Parent=itb01g00040.t1
	#Chr01	gt4sp_v3	exon	15840	16032	.	+	.	Parent=itb01g00040.t1
	#Chr01	gt4sp_v3	exon	17868	17904	.	+	.	Parent=itb01g00040.t1
	#Chr01	gt4sp_v3	exon	18177	18259	.	+	.	Parent=itb01g00040.t1
	#Chr01	gt4sp_v3	exon	18391	18457	.	+	.	Parent=itb01g00040.t1
	#Chr01	gt4sp_v3	CDS	13661	13835	.	+	0	Parent=itb01g00040.t1
	#Chr01	gt4sp_v3	CDS	15840	16032	.	+	2	Parent=itb01g00040.t1
	#Chr01	gt4sp_v3	CDS	17868	17904	.	+	1	Parent=itb01g00040.t1
	#Chr01	gt4sp_v3	CDS	18177	18259	.	+	0	Parent=itb01g00040.t1
	#Chr01	gt4sp_v3	CDS	18391	18457	.	+	1	Parent=itb01g00040.t1

	my @col = split(/\t/,$_);

	if(defined($col[2]) && $col[2] eq 'mRNA'){
		
		$gene = $_;
		$transcript = $_;

		if($col[8] =~ m/ID=([^;]+);.*?Parent=([^;\n]+)/){

			# create a new gene
			$tname = $1;
			$tparent = $2;
			$gene =~ s/\tmRNA\t/\tgene\t/;	
			$gene =~ s/ID=[^;]+;.*/ID=$tparent/;
	
			# first print the new gene feature
			print $gene;

			# now print transcript
			print $transcript;
		}
	
	} else {
		print;
	}

}
close(GFF);
