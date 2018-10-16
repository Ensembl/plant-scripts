#!/bin/env perl
use strict;
use warnings;
use Getopt::Std;

# call to health check python script from API, uses HV env variables
#
# by Bruno Contreras Moreira EMBL-EBI 2018

my ($db_name,$prod_server,$stage_server,%opts);
my ($HCSERVER,$HCSTAGING);

getopts('hs:p:d:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)){
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-d species database name                       (required, example: -d solanum_lycopersicum_core_42_95_3)\n";
  print "-p production db server                        (required, example: -p eg-p3-w)\n";
  print "-s staging db server                           (required, example: -s eg-s2)\n\n";
  exit(0);
}

if($opts{'d'}){ $db_name = $opts{'d'} }
else{ die "# EXIT : need a valid -d species db name, such as -f solanum_lycopersicum_core_42_95_3\n" }

if($opts{'p'}){ 
	$prod_server = $opts{'p'}; 
	chomp( $HCSERVER = `$prod_server details url` );
}
else{ die "# EXIT : need a valid -p production db server, such as -p eg-p3-w\n" }

if($opts{'s'}){ 
	$stage_server = $opts{'s'}; 
	chomp( $HCSTAGING = `$stage_server details url` );
}
else{ die "# EXIT : need a valid -s staging db server, such as -s eg-s2\n" }

my $hccmd="python $ENV{'ENSAPIPATH'}/ensembl-prodinf-core/ensembl_prodinf/hc_client.py ".
    "--uri $ENV{'HCENDPOINT'} ".
    "--db_uri \"$HCSERVER$db_name\" ".
    "--production_uri \"$ENV{'HCPRODUCTION'}ensembl_production\" " .
    "--staging_uri $HCSTAGING ".
    "--live_uri $ENV{'HCLIVE'} ".
    "--compara_uri \"$ENV{'HCCOMPARA_MASTER'}ensembl_compara_master\" ". 
    "--hc_groups $ENV{'HCGROUP'} ".
    "--data_files_path $ENV{'HCDATA_FILE_PATH'} ".
    "--tag $ENV{'HCTAG'} ".
    "--action submit ";

print $hccmd;
system("$hccmd");

