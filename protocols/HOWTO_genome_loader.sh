#!/bin/bash

# Protocol to load an ENA genome assembly into a production mysql db
# NOTE: ENA seems to consistently name chromosome names as integers from 1 to N
# NOTE: pseudo-chromosomes 0 are not supported and instead represented as a set of 
# sequences with "unplaced-scaffold" annotation
#
# Adapted from DBolser's run_the_genome_loader_pipeline_4.sh
# by Bruno Contreras Moreira EMBL-EBI 2018

## See: https://github.com/Ensembl/ensembl-genomeloader

# open interactive shell (requires ssh -X)
bshell25

######## manually edit as needed #############################################

# edit lc species name and full.\d GCA accession from ENA 
export species=solanum_lycopersicum
export gca=GCA_000188115.3
export genes_are_annotated=false # set to false if your assembly carries no gene annotation 

## manually set current Ensembl version, division, paths
# check https://www.ensembl.org for current version
# https://www.ebi.ac.uk/seqdb/confluence/pages/viewpage.action?spaceKey=EnsGen&title=MySQL+commands

web_ensembl_version=94
next_ensembl_version=$(echo $web_ensembl_version+1 | bc)
division=EnsemblPlants

ensemblapipath=/nfs/production/panda/ensemblgenomes/apis
mydevelpath=/nfs/panda/ensemblgenomes/development/$USER 

## manually set stage and production servers

stage_db=eg-s2  
prod_db=eg-p3-w # -w rewrite, ideally not the same copied over from stage

## manually set Health Check parameters

ENDPOINT=http://eg-prod-01.ebi.ac.uk:7000/hc/

SERVER=$(         $prod_db  details url)
PRODUCTION=$(     eg-pan    details url)
STAGING=$(        $stage_db  details url) # master_schema_$ensembl_version is looked up here
LIVE=$(           eg-sql    details url)
COMPARA_MASTER=$( eg-pan    details url)

GROUP=EGCoreHandover

DATA_FILE_PATH=/nfs/panda/ensembl/production/ensemblftp/data_files/

TAG=my_gl_hc_run

##############################################################################
##############################################################################

## checkout current Ensembl API 
# See: https://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/Creating+a+work+directory
# NOTE: it will complain if it was cloned earliear, you can safely skip errors such as
# "fatal: destination path 'ensembl' already exists and is not an empty directory."
# "Could not check out ensembl (release/95)"
cd $mydevelpath
$ensemblapipath/eg-utils/bin/checkout_ensembl.sh ensembl-$next_ensembl_version release/$next_ensembl_version
source ensembl-$next_ensembl_version/setup.sh

## get to-be-deprecated EG version and check version for the last time
ensembl_version=$(perl -MBio::EnsEMBL::ApiVersion -e "print software_version")
eg_version=$(echo $ensembl_version-53 | bc)

if [ "$ensembl_version" != "$next_ensembl_version" ]; then
	echo "ERROR: make sure web_ensembl_version is correct ($web_ensembl_version)"
	exit 1	
else
	echo "About to load $species $gca with version $ensembl_version ($eg_version)"
fi

## create folder for input assembly, clone and build Ensembl genome loader there
# NOTE: this avoids conflicts when different instances of GL are running at the same time
# NOTE2: I had to add my public farm node key to my github SSH keys
if [ ! -d "$gca" ]; then
    git clone git@github.com:Ensembl/ensembl-genomeloader.git "$gca"
else
    cd "$gca"
    git status
    git branch -v
    cd ..
fi

cd "$gca"/genome_materializer
time ./gradlew fatJar
cd ..

## Get Genome loader appropriate configuration
# See: https://github.com/Ensembl/ensembl-genomeloader/blob/master/CONFIG.md
# NOTE: Dan's copy was ~/EG_Places/Devel/lib/ensembl-genomeloader/enagenome_config.xml
# NOTE: might be useful comparing to https://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/Oracle+Instances
if [ "$genes_are_annotated" = true ] ; then
    cp ../enagenome_config.xml ./enagenome_config.xml
else
    cp ../enagenome_config.nogenes.xml ./enagenome_config.xml
fi


## set database name for this assembly and choose db server
db_suffix=${gca#*.}
db_name=${species}_core_${eg_version}_${ensembl_version}_${db_suffix}

echo "About to load $species $gca to database $db_name at server $prod_db"

cmd="perl \
    -I ./modules \
    ./scripts/load_genome.pl \
    \
    -a ${gca%.*}     \
    --division $division \
    \
    $($prod_db --details script) \
    --dbname $db_name
    \
    $(eg-pan     --details script_tax_) \
    --tax_dbname ncbi_taxonomy \
    \
    $(eg-pan     --details script_prod_) \
    --prod_dbname ensembl_production"

echo $cmd
time $cmd

# printout diagnostic info
if [ $? -eq 0 ]; then
    echo OK 
    hostname
    echo $species $gca
else
    echo FAIL
    hostname
    echo $species $gca
    exit
fi


### Now lets Health Check (HC) #############################################
############################################################################
# NOTE: docs at https://github.com/Ensembl/ensembl-prodinf-core/blob/master/docs/bulk_hc_submission.rst

if [ ! -d "$mydevelpath/ensembl-prodinf-core" ]; then
    git clone git@github.com:Ensembl/ensembl-prodinf-core.git "$mydevelpath/ensembl-prodinf-core"
else
    cd "$mydevelpath/ensembl-prodinf-core"
    git pull
    cd ..
fi

hccmd="python \
    $mydevelpath/ensembl-prodinf-core/ensembl_prodinf/hc_client.py \
    --uri $ENDPOINT \
    --db_uri "${SERVER}${db_name}" \
    --production_uri "${PRODUCTION}ensembl_production" \
    --staging_uri "${stage_db}\
    --live_uri $LIVE \
    --compara_uri "${COMPARA_MASTER}ensembl_compara_master" \
    --hc_groups $GROUP \
    --data_files_path $DATA_FILE_PATH \
    --tag $TAG  \
    --action submit"

echo $hccmd
$hccmd

echo "# NOTE: interactive submission of jobs might fail, you might have to re-send"
