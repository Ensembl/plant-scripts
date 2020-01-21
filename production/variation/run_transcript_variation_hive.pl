#!/bin/env perl
use strict;
use warnings;
use Getopt::Std;
use File::Basename;
use Bio::EnsEMBL::Registry;

# This script annotattes the consequences of variants overlapping transcripts
#
# NOTE: hive pipelines must be run in eb-cli nodes
#
# Uses env $ENSAPIPATH to locate ensembl-variation & ensembl-hive API
# and env $USER to create hive job names 
#
# Adapted from Dan Bolser's run_the_transcript_variation_pipeline.sh
# by B Contreras Moreira 2019-20
#
# http://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/The+transcript+variation+pipeline
#
## check user arguments ######################################################
##############################################################################

my (%opts,$species_db,$species,$ensembl_version,$prod_db_cmd);
my ($pipeline_dir,$reg_file,$hive_args,$hive_db,$hive_url,$argsline);
my ($rerun,$overwrite,$hive_path,$var_path) = (0,0);
my $hive_db_cmd = 'mysql-eg-hive-ensrw';
my $ensemblpath = $ENV{'ENSAPIPATH'};

my $max_distance = 500; 
# as in Francesc Montardit 2018 MSc thesis on brachy, prunus and A.thaliana
# The default 'max distance to transcript' seems to be about 5k
# Within this distance the variation will be called 'updstream or downstream' 
# Dan Bolser had 200, as a 'default' run for rice showed that a max_distance of
# 200 gives approximately 3 times as many up/down TVs as genic TVs

getopts('huwrp:m:s:v:R:H:P:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)){
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-s species_variation_db                        (required, example: -s oryza_sativa_variation_43_96_7)\n";
  print "-v next Ensembl version                        (required, example: -v 96)\n";
  print "-p production server command                   (required, should match -R, example: -p mysql-eg-prod-2-ensrw\n)"; 
  print "-R registry file, can be env variable          (required, example: -R \$p2panreg)\n";
  print "-P folder to put pipeline files, can be env    (required, example: -P \$tvartmp)\n";
  print "-H hive server command                         (optional, default: $hive_db_cmd)\n";
  print "-m max distance                                (optional, default: $max_distance)\n";
  print "-w over-write db (hive_force_init)             (optional, useful when a previous run failed)\n";                             
  print "-r re-run jump to beekeper.pl                  (optional, default: run init script from scratch)\n\n";
  exit(0);
}

if($opts{'v'}){
        $ensembl_version = $opts{'v'};

        # check ensembl-hive API exists for this version
        $hive_path = "$ensemblpath/ensembl-$ensembl_version/ensembl-hive/";
        if(!-d $hive_path) {
                die "# EXIT : cannot find ensembl-$ensembl_version/ensembl-hive,".
			"\n# make sure \$ENSAPIPATH is set\n";
        } 
	# now check ensembl-variation
	$var_path = "$ensemblpath/ensembl-$ensembl_version/ensembl-variation/";
        if(!-d $var_path) {
                die "# EXIT : cannot find ensembl-$ensembl_version/ensembl-variation,".
                        "\n# make sure \$ENSAPIPATH is set\n";
        }
}
else{ die "# EXIT : need a valid -v version, such as -v 95\n" }

if($opts{'s'}){
        $species_db = $opts{'s'};
	$species = (split(/_variation/,$species_db))[0]; 
        $hive_db = $ENV{'USER'}."_variation_consequence_$species";
}
else{ die "# EXIT : need a valid -s species_variation_db, such as -s oryza_sativa_variation_43_96_7\n" }

if($opts{'p'}){
        $prod_db_cmd = $opts{'p'};
}
else{ die "# EXIT : need a valid -p production server cmd\n" }

if($opts{'P'} && -d $opts{'P'}){ 
	$pipeline_dir = "$opts{'P'}/$species\_$ensembl_version";
}
else{ die "# EXIT : need a valid -P folder to put pipeline files, such as -P \$tvartmp\n" }

if($opts{'R'} && -e $opts{'R'}){ $reg_file = $opts{'R'} }
else{ die "# EXIT : need a valid -R file, such as -R \$p2panreg\n" }

if($opts{'H'}){ $hive_db_cmd = $opts{'H'} }
chomp( $hive_args = `$hive_db_cmd details script` );
if($hive_args =~ m/--host (\S+) --port (\d+) --user (\S+) --pass (\S+)/){
	$hive_args = "-hive_db_host=$1 -hive_db_port=$2 -hive_db_user=$3 -hive_db_password $4"; 	
}

chomp( $hive_url  = `$hive_db_cmd --details url` );
$hive_url .= $hive_db;

if($opts{'m'}){ $max_distance = $opts{'m'} }

if($opts{'r'}){ $rerun = 1 }

if($opts{'w'}){ $overwrite = 1 }

$argsline = sprintf("%s -s %s -v %s -p %s -R %s -H %s -P %s -m %d -w %d -r %d ",
  $0, $species_db, $ensembl_version, $prod_db_cmd, $reg_file, $hive_db_cmd, 
	$pipeline_dir, $max_distance, $overwrite, $rerun );

print "# $argsline\n\n";

## Add variation attributes to variation.attrib_type
system("$prod_db_cmd $species_db < $var_path/sql/attrib_entries.sql");

## Run init script and produce a hive_db with all tasks to be carried out
#########################################################################

my $initcmd = "init_pipeline.pl Bio::EnsEMBL::Variation::Pipeline::VariationConsequence_conf ".
    	"$hive_args ".
	"-hive_root_dir $hive_path ".
    	"-reg_file $reg_file ".
	"-species $species ".
    	"-pipeline_dir $pipeline_dir ".
	"-hive_force_init $overwrite -hive_debug_init 1 ".
	"-hive_default_max_retry_count 1 ".
	"-disambiguate_single_nucleotide_alleles 1 ".
	# these were used apparently by Dan Bolser
	"-default_lsf_options '-q production-rh74 -M  2000 -R \"rusage[mem= 2000]\"' ".
	"-highmem_lsf_options '-q production-rh74 -M 15000 -R \"rusage[mem=15000]\"' ".
	"-urgent_lsf_options  '-q production-rh74 -M  2000 -R \"rusage[mem= 2000]\"' ".
	"-long_lsf_options    '-q production-rh74 -M  2000 -R \"rusage[mem= 2000]\"' ";

print "# $initcmd\n\n";

if($rerun == 0){

	open(INITRUN,"$initcmd |") || die "# ERROR: cannot run $initcmd\n";
	while(<INITRUN>){
		print;
	}
	close(INITRUN);
}

## Send jobs to hive 
######################################################################### 

print "# hive job URL: $hive_url";

system("beekeeper.pl -url '$hive_url;reconnect_when_lost=1' -sync");
system("runWorker.pl -url '$hive_url;reconnect_when_lost=1' -reg_conf $reg_file");
system("beekeeper.pl -url '$hive_url;reconnect_when_lost=1' -reg_conf $reg_file -loop");

print "# hive job URL: $hive_url\n\n";
