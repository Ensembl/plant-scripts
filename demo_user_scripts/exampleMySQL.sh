#!/bin/bash

# Your usage of the data returned by the REST service is
# subject to same conditions as laid out on the Ensembl website.
#
# Copyright [2020] EMBL-European Bioinformatics Institute

# documentation about Ensembl schemas can be found at 
# http://www.ensembl.org/info/docs/api/index.html

# set server details
SERVER=mysql-eg-publicsql.ebi.ac.uk
USER=anonymous
PORT=4157

# set Ensembl Plants release number, check it
# at the bottom of http://plants.ensembl.org
EGRELEASE=47

RELEASE=$(( $EGRELEASE + 53)); # do not change

echo "EGRELEASE=${EGRELEASE}"
echo

## 1) check currently supported Ensembl Genomes (EG) core schemas,
##  note it includes non-plants as well

mysql --host $SERVER --user $USER --port $PORT \
	-e "show databases" | grep "core_${EGRELEASE}_${RELEASE}"

## 2) count protein-coding genes of a particular species

SPECIES=arabidopsis_thaliana
SPECIESCORE=`mysql --host $SERVER --user $USER --port $PORT \
    -e "show databases" | grep "${SPECIES}_core_${EGRELEASE}_${RELEASE}"`

mysql --host $SERVER --user $USER --port $PORT \
	$SPECIESCORE -e "SELECT COUNT(*) FROM gene WHERE biotype='protein_coding'"

## 3) get Triticum aestivum homeologous genes across A,B & D subgenomes
## Compara schema is described at 
## https://m.ensembl.org/info/docs/api/compara/compara_schema.html

# find out correct method_link_species_set.name
mysql --host $SERVER --user $USER --port $PORT \
	ensembl_compara_plants_${EGRELEASE}_$RELEASE -e \
	"SELECT name from method_link_species_set" | \
	grep "T.aes" | grep homoeologues

# remove LIMIT 10 if you want the complete set
mysql --host $SERVER --user $USER --port $PORT \
    ensembl_compara_plants_${EGRELEASE}_$RELEASE <<SQL
SELECT 
	homology_member.homology_id, cigar_line, perc_cov, perc_id, 
	perc_pos, gene_member.stable_id as genes, gene_member.genome_db_id
FROM
	homology_member 
	INNER JOIN homology USING (homology_id)
	INNER JOIN method_link_species_set USING (method_link_species_set_id)
	INNER JOIN gene_member USING (gene_member_id)
WHERE method_link_species_set.name="T.aes homoeologues" LIMIT 10
SQL



