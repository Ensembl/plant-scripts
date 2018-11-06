#!/usr/bin/env perl
#
# takes ENA accession and visits NCBI FTP to retrieve contig/chr synonyms from full assembly report
# 2018 Bruno Contreras Moreira EMBL-EBI
#
# Example calls: 
# $ ./get_chr_synonyms.pl GCA_000188115.3
# $ ./get_chr_synonyms.pl Viridiplantae.2018-10-09.report.tsv

use strict;
use warnings;
use Net::FTP;

# hard-coded paths, based on
# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/188/115/GCA_000188115.3_SL3.0/GCA_000188115.3_SL3.0_assembly_report.txt
my $NCBIFTPURL = 'ftp.ncbi.nlm.nih.gov';
my $GGAGENPATH = '/genomes/all/GCA/'; # absolute

my ($acc2path,$fullpath,@ENA_accessions);

## 0) check input params
if(defined($ARGV[0])){ 

	if(-s $ARGV[0]){ # input file instead of accession code

		open(INFILE,'<',$ARGV[0]) || 
			die "# ERROR: cannot read $ARGV[0]\n";
		while(<INFILE>){
			next if(/^#/); # skip commented lines
			push(@ENA_accessions,(split)[0]);
		}
		close(INFILE);
		
	} else { # single accession
		push(@ENA_accessions,$ARGV[0]);
	}
}
else{ die "# usage: $0 <valid ENA accession>\n\n£ example: $0 GCA_000188115.3\n"; }

## 1) connect to FTP site 
my $ftp = Net::FTP->new($NCBIFTPURL, Debug => 0) ||
	die "# ERROR: cannot connect to $NCBIFTPURL: $@";

$ftp->login("anonymous",'-anonymous@') || 
	die "# ERROR: cannot login ", $ftp->message();

## 2) check all ENA accessions
foreach my $ENA_accession (@ENA_accessions){

	# 2.0) check accession format
	# https://stackoverflow.com/questions/25523899/in-perl-best-way-to-insert-a-char-every-n-chars
	if($ENA_accession =~ /GCA_(\d+)\.\d+$/){
		$acc2path = $1;
                $acc2path =~ s/...\K(?=.)/\//sg;
        }
        else{ die "# EXIT : bad ENA accession: $ENA_accession lacks .version\n"; }

	# 2.1) compose folder path based on ENA accession
	$fullpath = $GGAGENPATH . $acc2path;

	# 2.2) get to appropriate remote folder
	$ftp->cwd($fullpath) || 
		die "# ERROR: cannot change working directory ", $ftp->message();

	# 2.3) check right assembly subfolder and try to copy assembly report file
	foreach my $subfolder ( $ftp->ls() ){
		if($subfolder =~ $ENA_accession) {
		
			$ftp->cwd($subfolder) || 
				die "# ERROR: cannot change working directory ", $ftp->message();

			foreach my $file ( $ftp->ls() ){

				if($file =~ $ENA_accession && $file =~ /assembly_report.txt/){

					$ftp->get($file) || 
						die "# ERROR: cannot download $file ", $ftp->message;	

					if(-s $file){
						print "# assembly report: $file\n";
					}
					last;
				}
			}
			last;
		}
	} 
}

# 3) disconnect and exit
$ftp->quit();
