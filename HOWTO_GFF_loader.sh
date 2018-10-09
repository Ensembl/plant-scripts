#!/bin/bash

# This protocol was written as follow up of HOWTO_genome_loader.sh
# Most of it was adapted from Dan Bolser's run_the_gff_loader2.sh
# Unlike the former, this pipeline uses hive
# by Bruno Contreras Moreira EMBL-EBI 2018

## https://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/Load+GFF3+Pipeline

######## manually edit as needed #############################################

# set input GFF and protein FASTA files
species=solanum_lycopersicum
dir=/nfs/production/panda/ensemblgenomes/data/Plants/Tomato/ITAG3.0_release/
protein_fasta_file=$dir/ITAG3.0_proteins.fasta
gff3_file=$dir/ITAG3.0_gene_models.gff
gene_source=$(grep -v "^#" $gff3_file | grep gene | cut -f 2 | head -1)
#gene_source=Gnomon
#gene_source=maker

# set pipeline and api paths, registry file
pipeline_dir=/hps/cstor01/nobackup/crop_genomics/Production_Pipelines
pipeline_dir=$pipeline_dir/${USER}/load_gff3/$species
ensemblapipath=/nfs/production/panda/ensemblgenomes/apis
mydevelpath=/nfs/panda/ensemblgenomes/development/$USER
registry=$mydevelpath/Registries/p2pan.reg

# manually set current Ensembl version, division and hive server 
# check https://www.ensembl.org for current version
web_ensembl_version=94
next_ensembl_version=$(echo $web_ensembl_version+1 | bc)

# hive db server and web frontend
hive_server=mysql-eg-hive-ensrw
hive_url=http://guihive.ebi.ac.uk:8080/

##############################################################################
##############################################################################

## check main arguments
if [ ! -d "$dir" ] ; then
	echo "# ERROR: cannot find $dir"
	exit 1
elif [ ! -e "$gff3_file" ] ; then
	echo "# ERROR: cannot find $gff3_file"
	exit 1
elif [ ! -s "$protein_fasta_file" ] ; then
        echo "# ERROR: cannot find $protein_fasta_file"
	exit 1
else
	echo "# input: $gff3_file $protein_fasta_file"
fi

## checkout current Ensembl API and eg-pipelines 
cd $mydevelpath
cd ensembl-$next_ensembl_version/ensembl-hive
git checkout origin/version/2.4
cd $mydevelpath
source ensembl-$next_ensembl_version/setup.sh

# get to-be-deprecated EG version and check version for the last time
ensembl_version=$(perl -MBio::EnsEMBL::ApiVersion -e "print software_version")
eg_version=$(echo $ensembl_version-53 | bc)

if [ "$ensembl_version" != "$next_ensembl_version" ]; then
        echo "ERROR: make sure web_ensembl_version is correct ($web_ensembl_version)"
        exit 1
else
        echo "About to load $species GFF with version $ensembl_version ($eg_version)"
fi

# now update eg-pipelines and add modules to $PERL5LIB
if [ ! -d eg-pipelines ]; then
    git clone git@github.com:EnsemblGenomes/eg-pipelines.git 
else
    cd eg-pipelines
    git pull
    cd ..
fi

PERL5LIB=$PERL5LIB:$mydevelpath/eg-pipelines/modules
export PERL5LIB=$(brew --prefix bioperl-169)/libexec:$PERL5LIB # adds Bio/DB/SeqFeature/Store.pm

## Run init script and produce a hive_db with all tasks to be carried out
cmd="init_pipeline.pl \
  	Bio::EnsEMBL::EGPipeline::PipeConfig::LoadGFF3_conf \
    	$($hive_server details script) \
    	--registry $registry \
    	--pipeline_dir $pipeline_dir \
    	--species $species \
    	--gff3_file $gff3_file \
    	--protein_fasta_file $protein_fasta_file \
    	--gene_source "$gene_source" \
	--hive_force_init 1" # make sure a previous db instance is overwritten

echo $cmd
time $cmd

# this might change with future pipelines updates
hive_db=${USER}_${pipeline}_${species}

url=$($hive_server --details url)$hive_db

echo "# hive job URL: $url"

url="${url};reconnect_when_lost=1"

beekeeper.pl -url ${url} -sync
runWorker.pl -url ${url} -reg_conf ${registry}
beekeeper.pl -url ${url} -reg_conf ${registry} -loop

echo "# Monitor this job at $hive_url"

