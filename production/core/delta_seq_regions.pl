#!/bin/env perl
use strict;
use warnings;
use Getopt::Std;

# Checks the evolution of the number of seq_regions of a species core db
# by querying a production server and the mirror
#
# Example: perl delta_seq_regions.pl -s glycine_max -p eg-p2
#
# Bruno Contreras Moreira EMBL-EBI 2019

my (%opts,$species,$prod_server);
my $mirror_server = 'mysql-eg-mirror';

my (@old_cores,@new_cores);

getopts('hs:p:m:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)){
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-s species_name                                (required, example: -s arabidopsis_thaliana)\n";
  print "-p production db server, where data is loaded  (required, example: -p eg-p2)\n";
  print "-m mirror db server, with historical data      (optional, default: -m mysql-eg-mirror)\n\n";
  exit(0);
}

if($opts{'s'}){ $species = $opts{'s'} } 
else{ die "# EXIT : need a valid -s species_name, such as -s arabidopsis_thaliana\n" }

if($opts{'p'}){ 
	$prod_server = $opts{'p'};
}
else{ die "# EXIT : need a valid -p production db server, such as -p eg-p3-w\n" }

if($opts{'m'}){
        $mirror_server = $opts{'m'};
}

## first check mirror core dbs

# find out previous cores
open(MIRROR,"$mirror_server -e \"SHOW DATABASES;\" | grep $species\_core |");
while(<MIRROR>){
	push(@old_cores, (split)[0]);
} 
close(MIRROR);

# count seq_regions in those
foreach my $coredb (@old_cores){
	open(MIRRORCORE,"$mirror_server $coredb -e \"SELECT COUNT(*) from seq_region\" |");
	while(<MIRRORCORE>){
		next if(/^COUNT/);
        	print "$mirror_server\t$coredb\t$_";
	}
	close(MIRRORCORE);
}

## now check current production server

open(PROD,"$prod_server -e \"SHOW DATABASES;\" | grep $species\_core |");
while(<PROD>){
        push(@new_cores, (split)[0]);
}
close(PROD);	

foreach my $coredb (@new_cores){
        open(PRODCORE,"$prod_server $coredb -e \"SELECT COUNT(*) from seq_region\" |");
        while(<PRODCORE>){
                next if(/^COUNT/);
                print "$prod_server\t$coredb\t$_";
        }
        close(PRODCORE);
}

