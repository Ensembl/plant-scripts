#!/bin/bash

# Copyright [2020] EMBL-European Bioinformatics Institute

# documentation about Ensembl VEP can be found at 
# http://www.ensembl.org/info/docs/tools/vep/index.html

# set Ensembl Plants release number, check it
# at the bottom of http://plants.ensembl.org
# EG stands for Ensembl Genomes
EGRELEASE=47

# work out Ensembl release, do not change
RELEASE=$(( EGRELEASE + 53));

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

## V2) Unpack downloaded cache file & check SIFT support 

# Note: cache downloaded in recipe F8
# Note: look for "sift  b"

SPECIES=arabidopsis_thaliana
VEPCACHE="${SPECIES}*.tar.gz*"
tar xfz $VEPCACHE
grep sift "${SPECIES}/${EGRELEASE}_*/info.txt"
echo "${SPECIES}/${EGRELEASE}_*/info.txt"

## V3) Predict effect of variants 

# See more options and examples at 
# http://www.ensembl.org/info/docs/tools/vep/script/vep_options.html
# http://www.ensembl.org/info/docs/tools/vep/script/vep_example.html

VCFILE=ensembl-vep/examples/arabidopsis_thaliana.TAIR10.vcf

VEPOPTIONS=(
	--genomes              # Ensembl Genomes, for Plants
	--species $SPECIES 
	--cache                # use local cache file, opposed to --database
	--dir_cache ./         # location of unpacked cache $SPECIES folder
	--check_existing       # co-located known variants
	--distance 5000        # max dist between variant and transcript
	--biotype              # show biotype of neighbor transcript
)

ensembl-vep/vep "${VEPOPTIONS[@]}"

## V4) Predict effect of variants for species not in Ensembl

# GFF file must be sorted and indexed with BGZIP and TABIX, see 
# http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#gff

FASTAGZFILE=      # GZIP-compressed file of genome FASTA file
GFFILE=           # gene models matching sequences in FASTAGZFILE
GZGFFILE=$GFFILE.sorted.gz

grep -v "#" $GFFILE | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > $GZGFFILE
tabix -p gff $GZGFFILE

ensembl-vep/vep -i $VCFILE -gff $GZGFFILE -fasta $FASTAGZFILE
