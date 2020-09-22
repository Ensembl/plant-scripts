#!/bin/env perl
use strict;
use warnings;

# script to fix GFF files like the one below so that load GFF3 works
# by Bruno Contreras Moreira EMBL-EBI 2020

die "# usage: $0 <GFF file>\n" if(!$ARGV[0]);

my ($feat_type);

#Chr1    maker   gene    14121985        14124796        .       -       .       ID=NC1G0017580;Name=NC1G0017580
#Chr1    maker   mRNA    14121985        14124796        .       -       .       ID=NC1G0017580;Parent=NC1G0017580;...
#Chr1    maker   exon    14124473        14124796        .       -       .       ID=NC1G0017580:1;Parent=NC1G0017580
#Chr1    maker   exon    14122986        14123030        .       -       .       ID=NC1G0017580:2;Parent=NC1G0017580
#Chr1    maker   exon    14122430        14122500        .       -       .       ID=NC1G0017580:3;Parent=NC1G0017580
#Chr1    maker   exon    14121985        14122336        .       -       .       ID=NC1G0017580:4;Parent=NC1G0017580
#Chr1    maker   five_prime_UTR  14124557        14124796        .       -       .       ID=NC1G0017580:five_prime_utr;Parent=NC1G0017580
#Chr1    maker   CDS     14124473        14124556        .       -       0       ID=NC1G0017580:cds;Parent=NC1G0017580
#Chr1    maker   CDS     14122986        14123030        .       -       0       ID=NC1G0017580:cds;Parent=NC1G0017580
#Chr1    maker   CDS     14122430        14122500        .       -       0       ID=NC1G0017580:cds;Parent=NC1G0017580
#Chr1    maker   CDS     14122192        14122336        .       -       1       ID=NC1G0017580:cds;Parent=NC1G0017580
#Chr1    maker   three_prime_UTR 14121985        14122191        .       -       .       ID=NC1G0017580:three_prime_utr;Parent=NC1G0017580

my ($ID,$mrnaID,%isoform) = ('', '');

open(GFF,'<',$ARGV[0]) || die "# ERROR: cannot read $ARGV[0]\n";
while(<GFF>){

	if(/^#/){
		print;
		next;
	}

	my @gffdata = split(/\t/,$_);

	if($gffdata[8]){ 
		
		$feat_type = $gffdata[2];	

		if($feat_type eq 'gene'){
			$ID = ''; # init
			if($gffdata[8] =~ m/ID=([^;]+)/){ $ID = $1 }					
		}
		elsif($feat_type eq 'mRNA'){
			$isoform{ $ID }++;
            $mrnaID = $ID .'.'. $isoform{ $ID };
			$gffdata[8] =~ s/ID=[^;]+/ID=$mrnaID/;                              
		}
		elsif($feat_type eq 'exon'){
			$gffdata[8] =~ s/ID=[^:]+/ID=$mrnaID/;
			$gffdata[8] =~ s/Parent=[^;\n]+/Parent=$mrnaID/;
    	}
		elsif($feat_type eq 'CDS'){
			$gffdata[8] =~ s/ID=[^:]+/ID=$mrnaID/;
			$gffdata[8] =~ s/Parent=[^;\n]+/Parent=$mrnaID/;
		}
	}

	print join("\t",@gffdata);
}
close(GFF);


