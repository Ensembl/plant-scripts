#!/bin/bash

# Copyright [2020-21] EMBL-European Bioinformatics Institute

# documentation about Ensembl VEP can be found at 
# http://www.ensembl.org/info/docs/tools/vep/index.html

# set Ensembl Plants release number, check it
# at the bottom of http://plants.ensembl.org
# EG stands for Ensembl Genomes

# example call: VEPATH=/path/to/ensembl-vep EGRELEASE=50 ./exampleVEP.sh

if [ -z "${EGRELEASE}" ]; then
	EGRELEASE=49
fi

# edit if needed to point to ensembl-vep
if [ -z "${VEPATH}" ]; then
	VEPATH=./
fi

# check VEP
if [ ! -f "${VEPATH}/ensembl-vep/vep" ]; then
	echo "# ERROR: Cannot find ${VEPATH}/ensembl-vep/vep not found, please set VEPATH accordingly"
    exit 1
fi


# work out Ensembl release, do not change
RELEASE=$((EGRELEASE + 53));

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

if [ ! -f ${VEPCACHE} ]; then
	echo "# ERROR: Cache file ${VEPCACHE} not found, get it with recipe F8"
	exit 1
else
	tar xfz $VEPCACHE
	pattern="${SPECIES}/${EGRELEASE}_*/info.txt"
	files=( $pattern )
	INFOFILE="${files[0]}" 
	if [ -f "${INFOFILE}" ]; then
		grep sift "${INFOFILE}"
		echo "${INFOFILE}"
	else
		echo "# ERROR: Cannot find file ${INFOFILE}, please correct/set variable EGRELEASE"
		exit 1
	fi
fi

## V3) Predict effect of variants 

# See more options and examples at 
# http://www.ensembl.org/info/docs/tools/vep/script/vep_options.html
# http://www.ensembl.org/info/docs/tools/vep/script/vep_example.html 

VCFILE="${VEPATH}/ensembl-vep/examples/arabidopsis_thaliana.TAIR10.vcf"
OUTFILE='arabidopsis_thaliana.vep.output'

VEPOPTIONS=(
	--genomes              # Ensembl Genomes, for Plants
	--species $SPECIES 
	--cache                # use local cache file, opposed to --database
	--dir_cache ./         # location of unpacked cache $SPECIES folder
	--cache_version $EGRELEASE
	--check_existing       # co-located known variants
	--distance 5000        # max dist between variant and transcript
	--biotype              # show biotype of neighbor transcript
    --sift b               # note not all species have SIFT precomputed
	--input_file $VCFILE
	--output_file $OUTFILE
)

${VEPATH}/ensembl-vep/vep "${VEPOPTIONS[@]}"

## V4) Predict effect of variants for species not in Ensembl

# GFF file must be sorted and indexed with BGZIP and TABIX, see 
# http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#gff

FASTAGZFILE=      # GZIP-compressed file of genome FASTA file
GFFILE=           # gene models matching sequences in FASTAGZFILE
GZGFFILE=$GFFILE.sorted.gz

if [[ -f $GFFILE && -f $FASTAGZFILE ]]; then
	# sort and index
	grep -v "#" $GFFILE | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > $GZGFFILE
	tabix -p gff $GZGFFILE
	
	# actually call vep
	${VEPATH}/ensembl-vep/vep -i $VCFILE -gff $GZGFFILE -fasta $FASTAGZFILE
fi
