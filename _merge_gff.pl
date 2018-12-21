#!/bin/env perl
use strict;
use warnings;

# This script takes two GFF3 files, one with genes,
#chr01	irgsp1_locus	gene	2983	10815	.	+	.	ID=Os01g0100100;
#
# and another with transcripts,
#chr01  irgsp1_rep      mRNA    2983    10815   .       +       .       ID=Os01t0100100-01;..Locus_id=Os01g0100100;...
#chr01  irgsp1_rep      five_prime_UTR  2983    3268    .       +       .       Parent=Os01t0100100-01
#chr01  irgsp1_rep      five_prime_UTR  3354    3448    .       +       .       Parent=Os01t0100100-01
#chr01  irgsp1_rep      CDS     3449    3616    .       +       0       Parent=Os01t0100100-01
#
# and merges them

my $PARENTtranscript = 'Locus_id';
my $FEATtranscript = 'mRNA';

my ($gff_gene_file,$gff_transcript_file);

if(!$ARGV[1]){ die "# usage: $0 <genes_GFF> <transcript_GFF>\n\n" }
else{ 
	($gff_gene_file,$gff_transcript_file) = @ARGV;
}

my (%transcript,$feature,$ID);

## read transcript GFF
$ID = 'NA';
open(GFF,'<',$gff_transcript_file) || die "# ERROR: cannot read $gff_transcript_file\n";
while(<GFF>){
	next if(/^#/); 
	my @gffcol = split(/\t/,$_);
	
	if($gffcol[2] eq $FEATtranscript){

		if($gffcol[8] =~ m/$PARENTtranscript=(\S+?);/){ $ID = $1; } 
		else{
			print "# ERROR: cannot parse $PARENTtranscript, please edit \$PARENTtranscript or ".
				"check GFF format ($gff_transcript_file)\n\n";
			exit(0);
		}
	}

	$transcript{$ID} .= $_;
}
close(GFF);

## now read gene GFF and interleave lines from transcript GFF
$ID = 'NA';
open(GFF,'<',$gff_gene_file) || die "# ERROR: cannot read $gff_gene_file\n";
while(<GFF>){
        next if(/^#/);
        my @gffcol = split(/\t/,$_);
        
        if($gffcol[2] eq 'gene'){

                if($gffcol[8] =~ m/ID=(\S+?);/){ $ID = $1; } #print $ID; }
                else{
                        print "# ERROR: cannot parse ID, please check GFF format ($gff_gene_file)\n\n";
                        exit(0);
                }

		if($transcript{$ID}){
			print $_; # print current gene
			print $transcript{$ID};
		}
		else{
			print "# ERROR: cannot find transcripts for gene ID: $ID\n\n";
                        exit(0);
		}
        }
}
close(GFF);

