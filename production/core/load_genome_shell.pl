#!/bin/env perl
use strict;
use warnings;
use Getopt::Std;
use Sys::Hostname;

# This script takes an ENA assembly's GCA accession and attempts to 
# load it into a production mysql db. Creates a local folder, clones 
# E! genome loader there and also sym-links the config file there
# 
# NOTE: This script should be run in an interactive farm shell
# NOTE: ENA seems to consistently name chromosome names as integers from 1 to N
# NOTE: pseudo-chromosomes 0 are not supported and instead represented as a set of "unplaced-scaffolds"
#
# Adapted from DBolser's run_the_genome_loader_pipeline_4.sh
# by Bruno Contreras Moreira EMBL-EBI 2018-9

## See: https://github.com/Ensembl/ensembl-genomeloader
# https://github.com/Ensembl/ensembl-genomeloader/blob/master/CONFIG.md
# Dan's copy was ~/EG_Places/Devel/lib/ensembl-genomeloader/enagenome_config.xml
# https://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/Oracle+Instances

my (%opts,$GCA_accession,$species,$division,$ensembl_version,$eg_version,$nogenes);
my ($prod_server,$egserver,$config_file,$argsline,$GCA_core_acc,$GCA_version);
my ($prod_db_args,$taxonomy_db_args,$analysis_db_args,$db_name,$interpro,$xrefs,$jsonpath);
my ($species_display_name,$species_production_name,$cmd) = ('','');

getopts('hXINJ:D:P:s:d:G:v:p:c:e:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)){
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-s species_name                                (required, example: -s arabidopsis_thaliana)\n";
  print "-d Ensembl division                            (required, example: -d EnsemblPlants)\n";
  print "-v next Ensembl version                        (required, example: -v 95)\n";
  print "-G GCA accession                               (required, example: -G GCA_000188115.3)\n";
  print "-p production db server, where data is loaded  (required, example: -p eg-p3-w)\n";
  print "-e server with ensembl_production db           (required, example: -e \$egprodserver)\n";
  print "-c full-path to +/- genes ENA config file      (required, example: -c \$enaconfigng)\n";
  print "-P species production_name                     (optional, example: -P marchantia_polymorpha)\n";
  print "-D species display name                        (optional, example: -D 'Marchantia polymorpha')\n";
  print "-J save JSON dump                              (optional, example: -J ./)\n";
  print "-N don't load genes or any other features      (optional, overrides -X, -I)\n";
  print "-I load InterPro mappings                      (optional)\n";
  print "-X load cross-references other than InterPro   (optional)\n\n";
  exit(0);
}

if($opts{'s'}){ $species = $opts{'s'} } 
else{ die "# EXIT : need a valid -s species_name, such as -s arabidopsis_thaliana\n" }

if($opts{'d'}){ $division = $opts{'d'} }
else{ die "# EXIT : need a valid -d division, such as -d EnsemblPlants\n" }

if($opts{'v'}){ 
	$ensembl_version = $opts{'v'}; 
	$eg_version = $ensembl_version-53;
}
else{ die "# EXIT : need a valid -v version, such as -v 95\n" }

if($opts{'G'}){ 
	$GCA_accession = $opts{'G'};
	if($GCA_accession =~ m/(\S+?)\.(\d+)$/){ 
		($GCA_core_acc,$GCA_version) = ($1,$2) 
	}
	else{ die "# EXIT : need a valid -G GCA accession, such as -G GCA_000188115.3\n" }

	# compose db name
	$db_name = $species.'_core_'.$eg_version.'_'.$ensembl_version.'_'.$GCA_version;
}
else{ die "# EXIT : need a valid -G GCA accession, such as -G GCA_000188115.3\n" }

if($opts{'p'}){ 
	$prod_server = $opts{'p'};
	chomp( $prod_db_args = `$prod_server -details script` );
}
else{ die "# EXIT : need a valid -p production db server, such as -p eg-p3-w\n" }

if($opts{'e'}){ 
	$egserver = $opts{'e'};
	chomp( $taxonomy_db_args = `$egserver --details script_tax_` );
	chomp( $analysis_db_args = `$egserver --details script_prod_` );
}
else{ die "# EXIT : need a valid -e EG production server, such as -e \$egprodserver\n" }

if($opts{'c'}){ $config_file = $opts{'c'} }
else{ die "# EXIT : need a valid -c file, such as -c \$enaconfigng\n" }

if($opts{'P'}){
	$species_production_name = $opts{'P'};
}

if($opts{'D'}){
   $species_display_name = $opts{'D'};
}

if($opts{'N'}){ $nogenes = 1 }
else{ 
	$nogenes = 0;

	if($opts{'I'}){ $interpro = 1 }
	else{ $interpro = 0 }

	if($opts{'X'}){ $xrefs = 1 }
	else{ $xrefs = 0 }
}

if($opts{'J'}){ $jsonpath = $opts{'J'} }
else{ $jsonpath = '' }

$argsline = sprintf("%s -s %s -d %s -v %s -G %s -p %s -e %s -c %s -D '%s' -P %s -I %d -X %d -N %d -J %s",
	$0, $species, $division, $ensembl_version, $GCA_accession, 
	$prod_server, $egserver, $config_file, 
	$species_display_name,$species_production_name, $interpro, $xrefs, $nogenes, 
	$jsonpath );

print "# $argsline\n\n";

# check this is not a login farm node
my $thishost = hostname();
if($thishost =~ m/ebi-login/ || $thishost =~ m/ebi-cli/){
	print "# ERROR: you should run this job in a farm interactive shell, not a login node\n\n";
	exit(0);
}

# clone and build ensembl-genomeloader in folder called $GCA_accession
# NOTE: this avoids conflicts with other instances of GL 
# NOTE: you'll need to add your public key to the list if github SSH keys
if(!-d $GCA_accession){
	system("git clone git\@github.com:Ensembl/ensembl-genomeloader.git $GCA_accession");
	#testing fork
	#system("git clone git\@github.com:brunocontrerasmoreira/ensembl-genomeloader.git $GCA_accession");
}
else {
    chdir($GCA_accession);
    system("git status");
    system("git branch -v");
    chdir('..');
}

chdir($GCA_accession.'/genome_materializer');
$ENV{'GRADLE_OPTS'} = "-Djava.io.tmpdir=/scratch"; 
# /scratch typically larger than /tmp, to avoid "java.io.IOException: No space left on device"
system("./gradlew fatJar");
chdir('..');

# link config file
#symlink($config_file,'enagenome_config.xml');

# actually call genome loader script

print "# loading $GCA_accession ($species) on database $db_name at $prod_server\n\n";

$cmd = "perl -I ./modules ./scripts/load_genome.pl ".
	"-a $GCA_core_acc ".
	"--division $division ".
   "--config_file $config_file ".
	"$prod_db_args --dbname $db_name ".
	"$taxonomy_db_args --tax_dbname ncbi_taxonomy ".
	"$analysis_db_args --prod_dbname ensembl_production "; #_$ensembl_version ";

if($species_display_name){
	$cmd .= " --display_name '$species_display_name' ";
}

if($species_production_name){
	$cmd .= " --production_name '$species_production_name' ";
}

if($jsonpath){
	$cmd .= " --save_json $jsonpath";
}

if($nogenes){ 
	$cmd .= " --nogenes ";
} else {
	if($interpro){
		$cmd .= " --load_interpro ";
	}
	if($xrefs){
	   $cmd .= " --load_xrefs ";
	}
}	



print "$cmd\n\n";

system("time $cmd");
if($? != 0){
    die "# ERROR: failed running $cmd\n";
}
