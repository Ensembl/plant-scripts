#!/bin/env perl
use strict;
use warnings;

# script to fix GFF files like the one below where a transcript strand does not match that of its gene
# by Bruno Contreras Moreira EMBL-EBI 2020

#BLBR01000003.1  ips     gene    16774145        16776582        .       +       .       ID=DRNTG_06943
# ...
#BLBR01000003.1  ips     mRNA    16774261        16775052        .       -       .       ID=DRNTG_06943.5;Parent=DRNTG_06943;
#BLBR01000003.1  ips     exon    16774261        16775052        .       -       .       ID=DRNTG_06943.5.exon1;Parent=DRNTG_06943.5
#BLBR01000003.1  ips     CDS     16774361        16774552        .       -       0       ID=DRNTG_06943.5.cds1;Parent=DRNTG_06943.5


die "# usage: $0 <GFF file>\n" if(!$ARGV[0]);

my $feat_type;
my $gene_strand = '';
open(GFF,'<',$ARGV[0]) || die "# ERROR: cannot read $ARGV[0]\n";
while(<GFF>){
	next if(/^#/);
	my @gffdata = split(/\t/,$_);

	$feat_type = $gffdata[2];	

	if($feat_type eq 'gene'){
		$gene_strand = $gffdata[6];
		print;
    }
	else { # not gene

		if($gffdata[6] ne $gene_strand) { next }
		else{ print }
	}
}
close(GFF);

