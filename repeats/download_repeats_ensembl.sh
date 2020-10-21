#!/usr/bin/env bash

# Script to download the annotated repeated elements of a selected
# Ensembl Genomes species

# Copyright [2020] EMBL-European Bioinformatics Institute

# documentation about Ensembl schemas can be found at 
# http://www.ensembl.org/info/docs/api/index.html

if [[ $# -eq 0 ]] ; then
	echo "# usage: $0 production_name, example: arabidopsis_thaliana"
	exit 0
else
	SPECIES=$1
fi

# set server details
SERVER=mysql-eg-publicsql.ebi.ac.uk
USER=anonymous
PORT=4157

# get Ensembl Plants current release number from FTP server
# Note: wget is used, this can be modified to use alternatives ie curl
FTPSERVER="ftp://ftp.ensemblgenomes.org/pub"
DIV=plants
SUMFILE="${FTPSERVER}/${DIV}/current/summary.txt"
RELEASE=`wget --quiet -O - $SUMFILE | \
	perl -lne 'if(/Release (\d+) of Ensembl/){ print $1 }'`

# work out Ensembl Genomes release
EGRELEASE=$(( RELEASE - 53));

SPECIESCORE=$(mysql --host $SERVER --user $USER --port $PORT \
	-e "show databases" | grep "${SPECIES}_core_${EGRELEASE}_${RELEASE}")

MINLEN=90

# note these might be reduntant
#1       3       106     trf
#1       4       91      trf
#1       755     761     dust
#1       1064    1098    trf
#1       1066    1097    trf

mysql --host $SERVER --user $USER --port $PORT $SPECIESCORE -Nb -e "SELECT sr.name,r.seq_region_start,r.seq_region_end,rc.repeat_class FROM repeat_feature r JOIN seq_region sr JOIN repeat_consensus rc WHERE r.seq_region_id=sr.seq_region_id AND r.repeat_consensus_id=rc.repeat_consensus_id AND (r.seq_region_end-r.seq_region_start+1) > $MINLEN  ORDER BY sr.name,r.seq_region_start" > ${SPECIES}.repeats.bed

mysql --host $SERVER --user $USER --port $PORT $SPECIESCORE -Nb -e "SELECT sr.name,g.seq_region_start,g.seq_region_end,g.stable_id FROM gene g JOIN seq_region sr WHERE g.seq_region_id=sr.seq_region_id  ORDER BY sr.name,g.seq_region_start" > ${SPECIES}.genes.bed

~/soft/bedtools2/bin/bedtools subtract -a ${SPECIES}.repeats.bed -b ${SPECIES}.genes.bed > ${SPECIES}.repeats.curated.bed

sort -k1,1 -k2,2n ${SPECIES}.repeats.curated.bed > ${SPECIES}.repeats.sorted.bed

~/soft/bedtools2/bin/bedtools merge -c 4 -o collapse -i ${SPECIES}.repeats.sorted.bed > ${SPECIES}.repeats.merged.bed


