#!/bin/bash

# Copyright [2020] EMBL-European Bioinformatics Institute

# documentation about Ensembl VEP can be found at 
# http://www.ensembl.org/info/docs/tools/vep/index.html

# set Ensembl Plants release number, check it
# at the bottom of http://plants.ensembl.org
EGRELEASE=47

RELEASE=$(( $EGRELEASE + 53)); # do not change

echo "EGRELEASE=${EGRELEASE}"
echo

## V1) Download, install and update VEP

# Fresh install
#git clone https://github.com/Ensembl/ensembl-vep.git
#cd ensembl-vep
#perl INSTALL.pl

# To update from a previous version:
#cd ensembl-vep
#git pull
#git checkout release/$RELEASE
#perl INSTALL.pl


