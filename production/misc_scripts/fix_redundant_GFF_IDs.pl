#!/bin/env perl
use strict;
use warnings;

# script to fix GFF files like the one below so that load GFF3 works
# by Bruno Contreras Moreira EMBL-EBI 2018

die "# usage: $0 <GFF file>\n" if(!$ARGV[0]);

my ($feat_type);

# check ID=names in GFF3 file to warn about gene:, mRNA:,... tags, which otherwise are added as stable_ids to db
#SL3.0ch00	maker_ITAG	gene	16480	17940	.	+	.	ID=gene:Solyc00g005000.3...
#SL3.0ch00	maker_ITAG	mRNA	16480	17940	.	+	.	ID=mRNA:Solyc00g005000.3.1...
#SC0002872	maker_ITAG	mRNA	162	2000	.	+	.	ID=mRNA:Solyc00g294237.1.1;Parent=gene:Solyc00g294237.1;
#SL3.0ch00	maker_ITAG	exon	16480	16794	.	+	.	ID=exon:Solyc00g005000.3.1.1...
#SC0003073	maker_ITAG	exon	854	1150	.	+	.	ID=exon:Solyc00g318835.1.1.1;Parent=mRNA:Solyc00g318835.1.1;
##SL3.0ch00	maker_ITAG	CDS	16480	16794	.	+	0	ID=CDS:Solyc00g005000.3.1.1...
open(GFF,'<',$ARGV[0]) || die "# ERROR: cannot read $ARGV[0]\n";
while(<GFF>){
	next if(/^#/);
	my @gffdata = split(/\t/,$_);

	if($gffdata[8] && $gffdata[8] =~ m/ID=(\w+):/){ 
		
		$feat_type = $1;	

		if($feat_type eq 'gene'){
			$gffdata[8] =~ s/gene://;					
		}
		elsif($feat_type eq 'mRNA'){
                        $gffdata[8] =~ s/ID=mRNA:/ID=/;                              
			$gffdata[8] =~ s/Parent=gene:/Parent=/;
                }
		elsif($feat_type eq 'exon'){
                        $gffdata[8] =~ s/ID=exon:(\S+?);/ID=$1.exon;/;
			$gffdata[8] =~ s/Parent=mRNA:/Parent=/;
                }
		elsif($feat_type eq 'CDS'){
                        $gffdata[8] =~ s/ID=CDS:(\S+?);/ID=$1.CDS;/;
                        $gffdata[8] =~ s/Parent=mRNA:/Parent=/;
                }
	}

	print join("\t",@gffdata);
}
close(GFF);

