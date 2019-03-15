#!/bin/env perl
use strict;
use warnings;
use Getopt::Std;

# This script adds sample features (location, gene, sample) to (table meta of) a core database
# 
# This might be useful after running load_GFF_hive.pl 
#
# by Bruno Contreras Moreira EMBL-EBI 2019

my ($db_name,$prod_server,$prod_server_details);
my (%opts, @dbs);

getopts('hc:s:p:d:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)){
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-d species database name or file with 1name/line (required, example: -d solanum_lycopersicum_core_42_95_3)\n";
  print "-p production db server                          (required, example: -p mysql-ens-plants-prod-1-ensrw)\n\n";
  exit(0);
}

if($opts{'p'}){ 
	$prod_server = $opts{'p'}; 
	chomp( $prod_server_details = `$prod_server details url` );
}
else{ die "# EXIT : need a valid -p production db server, such as -p mysql-ens-plants-prod-1-ensrw\n" }

if($opts{'d'}){ $db_name = $opts{'d'} }
else{ die "# EXIT : need a valid -d species db name, such as -d solanum_lycopersicum_core_42_95_3\n" }

if(-s $db_name){ # file with several dbs

	open(LIST,'<',$db_name) || die "# ERROR: cannot read $db_name\n";
	while(<LIST>){
		my $this_db_name = (split)[0];
		push(@dbs, $this_db_name);
	}
	close(LIST);
} else{ # single db_name

	push(@dbs, $db_name);
}

## loop through dbs and add sample features
foreach $db_name (@dbs){

	print "## $db_name\n";

}
