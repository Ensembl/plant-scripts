#!/usr/bin/env perl

# Example of Perl client to browse RNA-seq CRAM files from FTP server

# Sequence reads from RNA-seq studies at the European Nucleotide Archive
# are regularly mapped to genome assemblies in Ensembl Plants. For each study 
# CRAM files are created with the https://www.ebi.ac.uk/fg/rnaseq/api pipeline 
# and published on FTP site ftp://ftp.ensemblgenomes.org/pub/misc_data/Track_Hubs
#
# Each study contains a separate folder for each assembly used for mapping. 
# For instance, study SRP133995 was mapped to tomato assembly SL3.0 and the 
# tracksDb.txt file therein indicates the full path to the relevant CRAM file 
# next to its metadata. 

# Copyright [2020] EMBL-European Bioinformatics Institute

use strict;
use warnings;
use Net::FTP;

## C1) Find RNA-seq CRAM files for a genome assembly

# Note: assembly name is 'assembly_default' in recipe R2
# Note: can take a few minutes

my $FTPURL = 'ftp.ensemblgenomes.org';
my $FTPDIR = '/pub/misc_data/Track_Hubs';

my $assembly_name = '';
if($ARGV[0]){ $assembly_name = $ARGV[0] }
else{ die "# usage: $0 <assembly name, ie SL3.0>\n" }

my ($study,$file,$descr,$cramfile,$subgroup);

if( my $ftp = Net::FTP->new( $FTPURL, Passive=>1, Debug=>0, Timeout=>60) ){

	$ftp->login( "anonymous", '-anonymous@' )|| 
		die "# ERROR: cannot login " . $ftp->message();
	$ftp->cwd($FTPDIR) || 
		die "# ERROR: cannot change working directory to $FTPDIR " .
		$ftp->message();

	# print header
	print "cramfile\tstudy\tassembly\tsubgroup\tdescription\n";

	# stop if test only
	if($assembly_name eq 'test'){
		exit(0);
	}

	# list all ENA studies
	foreach $study ( $ftp->ls() ) {
		
		$ftp->cwd($study);

		my @contents = $ftp->ls();

		# skip other assemblies
		if(!grep(/^$assembly_name$/, @contents)){
			$ftp->cdup();
			next;
		}

		# get description from hub.txt
		$descr = 'NA';
		$ftp->get("hub.txt");
		if(open(HUB,"<","hub.txt")){
			while(<HUB>){
				if(/^longLabel ([^;]*)/){
					$descr = $1;
				}
			}
			close(HUB);
			unlink('hub.txt');
		}
		else { warn "# WARN: cannot get $study/hub.txt\n" }

		# look for CRAM files
		foreach $file (@contents) {
			if($file eq $assembly_name) {

				$ftp->cwd($file);

				# get and parse trackDb.txt
				$ftp->get("trackDb.txt");
				if(open(TRACKDB,"<","trackDb.txt")){
				
					while(<TRACKDB>){
						if(/track/){ $subgroup = 'NA' }
						elsif(/subGroups (.*)$/){ $subgroup = $1 }
						elsif(/bigDataUrl (\S+\.cram)/){
							$cramfile = $1;

							# print this CRAM file
							print "$cramfile\t$study\t$file\t$subgroup\t$descr\n";
						}
					}
					close(HUB);
					unlink('trackDb.txt');
				}
				else { warn "# WARN: cannot get $study/$file/trackDb.txt\n" }

				$ftp->cdup();
			}
		}

		# up to study level
		$ftp->cdup();
	}

	$ftp->close()
}
