#!/bin/env perl
use strict;
use warnings;

# add transcripts to a GFF file containing gene, exon and CDS features
# by Bruno Contreras Moreira EMBL-EBI 2019

if(!$ARGV[0]){ die "# usage: $0 <GFF>\n" }

my ($old_parent, $tparent, $tname, $feature, $transcript);

open(GFF,"<",$ARGV[0]) || die "# ERROR: cannot read $ARGV[0]\n";
while(<GFF>){

	# Pav_co3988703.1.	augustus	gene	1	195	0.7	-	.	ID=Pav_co3988703.1_g010.1.br;Note=PREDICTED: disease ...
	# Pav_co3988703.1.	augustus	exon	16	195	.	-	.	ID=Pav_co3988703.1_g010.1.br:exon:1;Parent=Pav_co3988703.1_g010.1.br
	# Pav_co3988703.1.	augustus	CDS	16	195	0.7	-	0	ID=Pav_co3988703.1_g010.1.br:CDS:1;Parent=Pav_co3988703.1

	my @col = split(/\t/,$_);
	
	if(defined($col[2]) && $col[2] eq 'gene'){

		print;

		# copy gene as a seed for a new transcript
		$transcript = $_;
		
		if($col[8] =~ m/ID=([^;\n]+)/){

			# create a new transcript from copied gene
			$tparent = $1;
			$tname = $tparent . ':mrna';
			$transcript =~ s/\tgene\t/\ttranscript\t/;
			$transcript =~ s/ID=[^;\n]+/ID=$tname;Parent=$tparent/;
			print $transcript;
		}
	} elsif(defined($col[2]) && ($col[2] eq 'exon' || $col[2] eq 'CDS')){
		
		$feature = $_;

		if($col[8] =~ m/ID=[^;]+;Parent=([^;\n]+)/){

			# update parent
			$old_parent = $1;
			$feature =~ s/Parent=$old_parent/Parent=$tname/;
			print $feature;
		}
	} else {
		print;
	}

}
close(GFF);
