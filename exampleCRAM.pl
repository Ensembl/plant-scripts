#!/usr/bin/env perl

# Example of Perl client to browse RNA-seq CRAM files from FTP server

# Copyright [2020] EMBL-European Bioinformatics Institute

use strict;
use warnings;
use Net::FTP;

## C1) Find RNA-seq CRAM files for a genome assembly

# Note: assembly name is 'assembly_default' in recipe R2

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
	print "study\tassembly\tsubgroup\tcramfile\tdescription\n";

	# stop if test only
	if($assembly_name eq 'test'){
		exit(0);
	}

	# list all ENA studies
	foreach $study ( $ftp->ls() ) {
		
		$ftp->cwd($study);

		my @contents = $ftp->ls();

		# skip other assemblies
		next if(!grep(/^$assembly_name$/, @contents));

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
							print "$study\t$file\t$subgroup\t$cramfile\t$descr\n";
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
