#!/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# This script submits a core stats job to hive
#
# It uses env $USER to create hive job names and assumes Ensembl-version API
# is loaded in @INC / $PERL5LIB
#
# Adapted from Dan Bolser's run_the_core_stats_pipeline.sh
# by B Contreras Moreira 2018-9
#
# http://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/Core+Statistics+Pipeline
#
## check user arguments ######################################################
##############################################################################

my $hive_db_cmd = 'mysql-ens-hive-prod-2-ensrw';
my $overwrite = 0;
my ($help,$reg_file,@species,$species_cmd,$ensembl_version);
my ($hive_args,$hive_url,$hive_db);      

GetOptions(	
		"help|?" => \$help,
		"overwrite|w" => \$overwrite, 
		"version|v=s" => \$ensembl_version,
                "species|s=s" => \@species,
		"hivecmd|H=s" => \$hive_db_cmd,    
		"regfile|R=s" => \$reg_file,
) || help_message(); 

if($help){ help_message() }

sub help_message {
	print "\nusage: $0 [options]\n\n".
	"-s species_name(s)                     (required, example: -s arabidopsis_thaliana -s brachypodium_distachyon)\n".
	"-v next Ensembl version                (required, example: -v 95)\n".
	"-R registry file, can be env variable  (required, example: -R \$p2panreg)\n".
	"-H hive database command               (optional, default: $hive_db_cmd)\n".
	"-w over-write db (hive_force_init)     (optional, useful when a previous run failed)\n";
	exit(0);
}

if($ensembl_version){
	# check Ensembl API is in env
	if(!grep(/ensembl-$ensembl_version\/ensembl-hive\/modules/,@INC)){
                die "# EXIT : cannot find ensembl-$ensembl_version/ensembl-hive/modules in \$PERL5LIB / \@INC\n"
        }
}
else{ die "# EXIT : need a valid -v version, such as -v 95\n" } 

if(@species){
	foreach my $sp (@species){ 
		$species_cmd .= "--species $sp "; 
	}
} 
else{ die "# EXIT : need a valid -s species, such as -s arabidopsis_thaliana -s brachypodium_distachyon\n" }

if(!$reg_file || !-e $reg_file){ die "# EXIT : need a valid -R file, such as -R \$p2panreg\n" }

chomp( $hive_args = `$hive_db_cmd details script` );
$hive_db = $ENV{'USER'}."_core_statistics_$ensembl_version";
chomp( $hive_url  = `$hive_db_cmd --details url` );
$hive_url .= $hive_db;

## Run init script and produce a hive_db with all tasks to be carried out
#########################################################################

my $initcmd = "init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::CoreStatistics_conf ".
    	"$hive_args ".
    	"--registry $reg_file ".
    	"$species_cmd ".
	"--hive_force_init $overwrite";

print "# $initcmd\n\n";

open(INITRUN,"$initcmd |") || die "# ERROR: cannot run $initcmd\n";
while(<INITRUN>){
	print;
}
close(INITRUN);


## Send jobs to hive 
######################################################################### 

print "# hive job URL: $hive_url";

system("beekeeper.pl -url '$hive_url;reconnect_when_lost=1' -sync");
system("runWorker.pl -url '$hive_url;reconnect_when_lost=1' -reg_conf $reg_file");
system("beekeeper.pl -url '$hive_url;reconnect_when_lost=1' -reg_conf $reg_file -loop");

print "# hive job URL: $hive_url\n\n";

