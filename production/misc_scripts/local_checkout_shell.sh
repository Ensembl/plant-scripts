#!/bin/bash

# This script does a local checkout of Ensembl APIs using env variables

# Should be run in interactive shell such as
# bsub -q production-rh7 -M4000 -R "select[mem>4000] rusage[mem=4000]" -Is bash 

# by Bruno Contreras Moreira EMBL-EBI 2018-20

eg_api_path=$EGAPIPATH
ensembl_api_path=$ENSAPIPATH
ensembl_version=$1

if [ $# -eq 0 ]
  then
    echo "# ERROR: need ensembl_version"
    exit 1
elif [ -z "$EGAPIPATH" ]; then
	echo "# ERROR: make sure \$EGAPIPATH is set"
	exit 1
elif [ -z "$ENSAPIPATH" ]; then
	echo "# ERROR: make sure \$ENSAPIPATH is set"
	exit 1
fi


# checkout current Ensembl API 
# See: https://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/Creating+a+work+directory
# NOTE: it will complain if it was cloned earliear, you can safely skip errors such as
# "fatal: destination path 'ensembl' already exists and is not an empty directory."
# "Could not check out ensembl (release/95)"
cd $ensembl_api_path
#$eg_api_path/eg-utils/bin/checkout_ensembl.sh ensembl-$ensembl_version release/$ensembl_version
#checkout master branch by default, so that you can run it before next release branch is frozen
$eg_api_path/eg-utils/bin/checkout_or_update_ensembl.sh ensembl-$ensembl_version
cd ensembl-$ensembl_version
git ensembl --clone regulation --branch master
git ensembl --pull regulation --branch master
cd ..
source ensembl-$ensembl_version/setup.sh

# checkout hive  
cd ensembl-$ensembl_version/ensembl-hive
git checkout version/2.4
cd ../..

# update Health Check code
# NOTE: docs at https://github.com/Ensembl/ensembl-prodinf-core/blob/master/docs/bulk_hc_submission.rst
cd ensembl-prodinf-core
git pull
cd ..

# datacheck for core_stats
cd ensembl-datacheck
git pull
cd ..

# for ATAC assembly mappings, commented in case local edited copy is not synced
#cd eg-assemblyconverter
#git pull
#cd ..

# check ensembl version matches API
api_ensembl_version=$(perl -MBio::EnsEMBL::ApiVersion -e "print software_version")
if [ "$ensembl_version" != "$api_ensembl_version" ]; then
        echo "ERROR: make sure ensembl_version is correct ($ensembl_version)"
        exit 1
fi

# now update eg-pipelines and add modules to $PERL5LIB
cd eg-pipelines
git pull
cd ..

