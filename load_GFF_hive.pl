#!/bin/env perl
use strict;
use warnings;
use Getopt::Std;
use Cwd;

# This script takes a GFF3 & a peptide FASTA file and attempts to load the 
# features on top of a previously loaded ENA genome assembly
#
# This should be run after HOWTO_genome_loader.sh
# It submits related tasks to hive 
#
# It uses env $USER to create hive job names
#
# Adapted from Dan Bolser's run_the_gff_loader2.sh by B Contreras Moreira

## https://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/Load+GFF3+Pipeline

# set lib paths from env
my $ENSAPIPATH = $ENV{'ENSAPIPATH'}; 

## check user arguments ######################################################
##############################################################################

my (%opts,$species,$protein_fasta_file,$gff3_file,$gene_source,$ensembl_version);
my ($pipeline_dir,$reg_file,$hive_args,$hive_db,$hive_url);
my $hive_db_cmd = 'mysql-eg-hive-ensrw';
my $rerun = 0;

getopts('hrs:f:g:S:v:R:H:P:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)){
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-s species_name                                (required, example: -s arabidopsis_thaliana)\n";
  print "-f protein FASTA file                          (required, example: -f atha.pep.fasta)\n";
  print "-v next Ensembl version                        (required, example: -v 95)\n";
  print "-g GFF3 file                                   (required, example: -g atha.gff)\n";
  print "-R registry file, can be env variable          (required, example: -R \$p2panreg)\n";
  print "-P folder to put pipeline files, can be env    (required, example: -P \$gfftmp)\n";
  print "-S source of gene annotation, one word         (optional, default: taken from 3rd column of GFF3 file)\n";
  print "-H hive database command                       (optional, default: $hive_db_cmd)\n";
  print "-r re-run jump to beekeper.pl                  (optional, default: run init script from scratch)\n\n";
  exit(0);
}

if($opts{'s'}){ 
	$species = $opts{'s'}; 
	$hive_db = "$ENV{'USER'}_load_gff3_$species"; 
} 
else{ die "# EXIT : need a valid -s species_name, such as -s arabidopsis_thaliana\n" }

if($opts{'f'} && -e $opts{'f'}){ $protein_fasta_file = $opts{'f'} }
else{ die "# EXIT : need a valid -f file, such as -f atha.pep.fasta\n" }

if($opts{'g'} && -e $opts{'g'}){ $gff3_file = $opts{'g'} }
else{ die "# EXIT : need a valid -g file, such as -g atha.gff\n" }

if($opts{'S'}){ $gene_source = $opts{'S'} }
else{
	$gene_source = `grep -v "^#" $gff3_file | grep gene | cut -f 2 | head -1`;
	if(!$gene_source){ die "# EXIT : cannot parse annotation source from $gff3_file\n" }
}

if($opts{'v'}){
	$ensembl_version = $opts{'v'};	

	# check Ensembl API is in place, make sure it is up-to-date, and source it
	my $hive_api_dir = "$ENSAPIPATH/ensembl-$ensembl_version/ensembl-hive";
	if(-d $hive_api_dir){
		my $cwd = getcwd();
		chdir($hive_api_dir);
		system("git checkout origin/version/2.4");
		chdir($ENSAPIPATH);
		system("source ensembl-$ensembl_version/setup.sh");
		chdir($cwd);
	}
	else { die "# EXIT : cannot find $hive_api_dir, is \$ENSAPIPATH set?\n" }
}
else{ die "# EXIT : need a valid -v version, such as -v 95\n" }

if($opts{'R'} && -e $opts{'R'}){ $reg_file = $opts{'R'} }
else{ die "# EXIT : need a valid -R file, such as -R \$p2panreg\n" }

if($opts{'H'}){ $hive_db_cmd = $opts{'H'} }
$hive_args = `$hive_db_cmd details script`;
$hive_url  = `$hive_db_cmd --details url`.$hive_db;

if($opts{'P'} && -d $opts{'P'}){ $pipeline_dir = "$opts{'P'}/$species" }
else{ die "# EXIT : need a valid -P folder to put pipeline files, such as -P \$gfftmp\n" }

if($opts{'r'}){ $rerun = 1 }


## Run init script and produce a hive_db with all tasks to be carried out
#########################################################################

my $initcmd="init_pipeline.pl \
  	Bio::EnsEMBL::EGPipeline::PipeConfig::LoadGFF3_conf \
    	$hive_args \
    	--registry $reg_file \
    	--pipeline_dir $pipeline_dir \
    	--species $species \
    	--gff3_file $gff3_file \
    	--protein_fasta_file $protein_fasta_file \
    	--gene_source '$gene_source' \
	--hive_force_init 1"; # previous db instance will be overwritten

if($rerun == 0){

	open(INITRUN,"$initcmd |") || die "# ERROR: cannot run $initcmd\n";
	while(<INITRUN>){
		print;
	}
	close(INITRUN);
}

## Send jobs to hive 
######################################################################### 

system("beekeeper.pl -url $hive_url;reconnect_when_lost=1 -sync");
system("runWorker.pl -url $hive_url;reconnect_when_lost=1 -reg_conf $reg_file");
system("beekeeper.pl -url $hive_url;reconnect_when_lost=1 -reg_conf $reg_file -loop");

print "# hive job URL: $hive_url";

