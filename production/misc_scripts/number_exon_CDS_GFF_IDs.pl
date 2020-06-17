#!/bin/env perl
use strict;
use warnings;

# script to fix GFF files like the one below so that load GFF3 works
# by Bruno Contreras Moreira EMBL-EBI 2020

die "# usage: $0 <GFF file>\n" if(!$ARGV[0]);

my ($feat_type);

#BLBR01000018.1 taco_ips        gene    22416637        22421843        .       -       .       ID=DRNTG_23552;
#BLBR01000018.1 taco_ips        mRNA    22416637        22421843        .       -       .       ID=DRNTG_23552.1;Parent=DRNTG_23552;
#BR01000018.1 taco_ips        exon    22416637        22417837        .       -       .       ID=DRNTG_23552.1.exon;Parent=DRNTG_23552.1
#BLBR01000018.1 taco_ips        CDS     22417203        22417837        .       -       2       ID=DRNTG_23552.1.cds;Parent=DRNTG_23552.1
#BLBR01000018.1 taco_ips        exon    22417937        22418084        .       -       .       ID=DRNTG_23552.1.exon;Parent=DRNTG_23552.1
#BLBR01000018.1 taco_ips        CDS     22417937        22418084        .       -       0       ID=DRNTG_23552.1.cds;Parent=DRNTG_23552.1

my $count_exon = 0;
my $count_cds = 0;
open(GFF,'<',$ARGV[0]) || die "# ERROR: cannot read $ARGV[0]\n";
while(<GFF>){
	next if(/^#/);
	my @gffdata = split(/\t/,$_);

	if($gffdata[8] && $gffdata[8] =~ m/ID=[^;]+/){ 
			
		$feat_type = $gffdata[2];	

		if($feat_type eq 'mRNA' || $feat_type eq 'transcript'){
			$count_exon = 0;
			$count_cds = 0;
        }
		elsif($feat_type eq 'exon'){
			$count_exon++;
			$gffdata[8] =~ s/ID=(\S+?);/ID=$1$count_exon;/;
        }
		elsif($feat_type eq 'CDS'){
			$count_cds++;
			$gffdata[8] =~ s/ID=(\S+?);/ID=$1$count_cds;/;
		}
	}

	print join("\t",@gffdata);
}
close(GFF);

