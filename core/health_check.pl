#!/bin/env perl
use strict;
use warnings;
use Getopt::Std;

# call to health check python script from API, uses HV env variables
# by Bruno Contreras Moreira EMBL-EBI 2018

my ($db_name,$prod_server,$stage_server,$compara_server,$hccmd,%opts);
my ($HCSERVER,$HCSTAGING,$HCCOMPARA);

my $HCTEMPLATE = "python $ENV{'ENSAPIPATH'}/ensembl-prodinf-core/ensembl_prodinf/hc_client.py ".
	"--uri $ENV{'HCENDPOINT'} ".
	"--production_uri \"$ENV{'HCPRODUCTION'}ensembl_production\" " .
	"--live_uri $ENV{'HCLIVE'} ".
	"--hc_groups $ENV{'HCGROUP'} ".
	"--data_files_path $ENV{'HCDATA_FILE_PATH'} ".
	"--action submit ";

getopts('hc:s:p:d:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)){
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-d species database name or file with 1name/line (required, example: -d solanum_lycopersicum_core_42_95_3)\n";
  print "-p production db server                          (required, example: -p eg-p3-w)\n";
  print "-c compara db server                             (required, example: -c mysql-ens-compara-prod-5)\n";
  print "-s staging db server                             (required, example: -s eg-s2)\n\n";
  exit(0);
}

if($opts{'d'}){ $db_name = $opts{'d'} }
else{ die "# EXIT : need a valid -d species db name, such as -d solanum_lycopersicum_core_42_95_3\n" }

if($opts{'p'}){ 
	$prod_server = $opts{'p'}; 
	chomp( $HCSERVER = `$prod_server details url` );
}
else{ die "# EXIT : need a valid -p production db server, such as -p eg-p3-w\n" }

if($opts{'c'}){
        $compara_server = $opts{'c'};
        chomp( $HCCOMPARA = `$compara_server details url` );
        $HCTEMPLATE .= "--compara_uri \"$HCCOMPARA";
	$HCTEMPLATE .= "ensembl_compara_master_plants\" ";
}
else{ die "# EXIT : need a valid -c compara db server, such as -c mysql-ens-compara-prod-5\n" }

if($opts{'s'}){ 
	$stage_server = $opts{'s'}; 
	chomp( $HCSTAGING = `$stage_server details url` );
	$HCTEMPLATE .= "--staging_uri $HCSTAGING ";
}
else{ die "# EXIT : need a valid -s staging db server, such as -s eg-s2\n" }

if(-s $db_name){ # file with several dbs to HC

	open(LIST,'<',$db_name) || die "# ERROR: cannot read $db_name\n";
	while(<LIST>){
		my $this_db_name = (split)[0];
	
		$hccmd = "$HCTEMPLATE --db_uri \"$HCSERVER$this_db_name\" --tag $db_name";

	        print $hccmd;
        	system("$hccmd");
	}
	close(LIST);
}
else{ # single db_name

	$hccmd = "$HCTEMPLATE --db_uri \"$HCSERVER$db_name\" --tag $db_name";

	print $hccmd;
	system("$hccmd");
}
