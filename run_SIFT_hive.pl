#!/bin/env perl
use strict;
use warnings;
use Getopt::Std;

# This script runs SIFT for all single amino acid substitutions in protein
# sequence of a selected species.
# By looking at patterns of position-specific amino acid substitutions
# observed in BLASTP hits, SIFT predicts whether variations are likely to 
# be deletereous (rare) or not (frequent)
#
# NOTE: FULL Mode populates these tables: meta, protein_function_predictions, 
# protein_function_predictions_attrib, translation_md5,translation_mapping
#
# NOTE: hive pipelines must be run in eb-cli nodes
#
# Uses env $ENSAPIPATH to locate ensembl-hive API
# and env $USER to create hive job names 
#
# B Contreras Moreira EMBL-EBI 2019
#
# https://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/Running+the+Sift+pipeline
#
## check user arguments ######################################################
##############################################################################

my (%opts,$species,$protein_fasta_file,$ensembl_version);
my ($pipeline_dir,$reg_file,$hive_args,$hive_db,$hive_url,$argsline);
my ($rerun,$update,$overwrite,$hive_path) = (0,0,0);
my $hive_db_cmd = 'mysql-eg-hive-ensrw';
my $ensemblpath = $ENV{'ENSAPIPATH'};
my $blastbin = $ENV{'blastbin'};
my $blastdb = $ENV{'uniref90'};

getopts('huwrb:p:s:v:R:H:P:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)){
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-s species_name                                (required, example: -s arabidopsis_thaliana)\n";
  print "-v next Ensembl version                        (required, example: -v 95)\n";
  print "-R registry file, can be env variable          (required, example: -R \$p2panreg)\n";
  print "-P folder to put pipeline files, can be env    (required, example: -P \$siftmp)\n";
  print "-H hive database command                       (optional, default: $hive_db_cmd)\n";
  print "-p path to BLAST+ binaries                     (optional, default: $blastbin)\n";
  print "-b path to BLASTDB fasta file                  (optional, default: $blastdb)\n";
  print "-u update SIFT predictions [UPDATE]            (optional, default: reset and predict from scratch [FULL])\n";
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
}
else{ die "# EXIT : need a valid -v version, such as -v 95\n" }

if($opts{'s'}){
        $species = $opts{'s'};
        $hive_db = $ENV{'USER'}."_protein_function_$species";
}
else{ die "# EXIT : need a valid -s species_name, such as -s arabidopsis_thaliana\n" }

if($opts{'P'} && -d $opts{'P'}){ 
	$pipeline_dir = "$opts{'P'}/$species\_$ensembl_version";
	$protein_fasta_file = "$pipeline_dir/$species.faa"; # output by the pipeline
}
else{ die "# EXIT : need a valid -P folder to put pipeline files, such as -P \$siftmp\n" }

if($opts{'R'} && -e $opts{'R'}){ $reg_file = $opts{'R'} }
else{ die "# EXIT : need a valid -R file, such as -R \$p2panreg\n" }

if($opts{'p'}){ 
	$blastbin = $opts{'p'}; 
	if(!-e "$blastbin/blastp"){
		die "# EXIT : need a valid -p BLAST+ path\n"
	}
}

if($opts{'b'}){ 
        $blastdb = $opts{'b'}; 
        if(!-e "$blastdb.00.phr"){
                die "# EXIT : please run\n# $blastbin/makeblastdb -dbtype prot -in $blastdb\n# and re-rerun\n"
        }
}

if($opts{'H'}){ $hive_db_cmd = $opts{'H'} }
chomp( $hive_args = `$hive_db_cmd details script` );
if($hive_args =~ m/--host (\S+) --port (\d+) --user (\S+) --pass (\S+)/){
	$hive_args = "-pipeline_db -host=$1 --pipeline_db -port=$2 --pipeline_db -user=$3 -password $4"; 	
}

chomp( $hive_url  = `$hive_db_cmd --details url` );
$hive_url .= $hive_db;

if($opts{'r'}){ $rerun = 1 }

if($opts{'w'}){ $overwrite = 1 }

if($opts{'u'}){ $update = 1 }


$argsline = sprintf("%s -s %s -f %s -v %s -R %s -H %s -P %s -p %s -b %s -u %d -w %d -r %d",
  $0, $species, $protein_fasta_file,  
  $ensembl_version, $reg_file, $hive_db_cmd, $pipeline_dir, 
  $blastbin,$blastdb,
  $update, $overwrite, $rerun );

print "# $argsline\n\n";

## Check BLASTDB 


## Run init script and produce a hive_db with all tasks to be carried out
#########################################################################

my $initcmd = "init_pipeline.pl Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::ProteinFunction_conf ".
    	"$hive_args ".
	"-hive_root_dir $hive_path ".
    	"-ensembl_registry $reg_file -species $species ".
    	"-pipeline_dir $pipeline_dir -sift_dir $pipeline_dir -sift_working $pipeline_dir ".
	"-ncbi_dir $blastbin -blastdb $blastdb ".
	"--hive_force_init $overwrite ";

if(!$update){ $initcmd .= "-sift_run_type 1 " }

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
