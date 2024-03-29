#!/bin/bash

# Your usage of the data returned by the public MySQL server is
# subject to same conditions as laid out on the Ensembl website.
#
# Copyright [2020-21] EMBL-European Bioinformatics Institute

# documentation about Ensembl schemas can be found at 
# http://www.ensembl.org/info/docs/api/index.html

# set server details
SERVER=mysql-eg-publicsql.ebi.ac.uk
USER=anonymous
PORT=4157

# get Ensembl Plants current release number from FTP server
# Note: wget is used, this can be modified to use alternatives ie curl
FTPSERVER="ftp://ftp.ensemblgenomes.org/pub"
DIV=plants
SUMFILE="${FTPSERVER}/${DIV}/current/summary.txt"
RELEASE=$(wget --quiet -O - $SUMFILE | \
	perl -lne 'if(/Release (\d+) of Ensembl/){ print $1 }')

# work out Ensembl Genomes release
EGRELEASE=$((RELEASE - 53));

# alternatively set other EG release number
# EGRELEASE=

echo "EGRELEASE=${EGRELEASE}"
echo

# stop here if just a test
if [[ $# -ge 0 ]] && [[ $1 = "test" ]]; then
	mysql --host $SERVER --user $USER --port $PORT -e "show databases"
	exit 0
fi

## S1) Check currently supported Ensembl Genomes (EG)/non-vertebrates core schemas,

# Note: includes non-plants as well

mysql --host $SERVER --user $USER --port $PORT \
	-e "show databases" | grep "core_${EGRELEASE}_${RELEASE}"

# The following API script can also be used:
# https://github.com/Ensembl/ensembl-metadata/blob/master/misc_scripts/get_list_databases_for_division.pl

## S2) Count protein-coding genes of a particular species

SPECIES=arabidopsis_thaliana
SPECIESCORE=$(mysql --host $SERVER --user $USER --port $PORT \
	-e "show databases" | grep "${SPECIES}_core_${EGRELEASE}_${RELEASE}")

mysql --host $SERVER --user $USER --port $PORT \
	$SPECIESCORE -e "SELECT COUNT(*) FROM gene WHERE biotype='protein_coding'"

## S3) Get stable_ids of transcripts used in Compara analyses 

# Canonical transcripts are used in the gene tree analysis,
# which usually are the longest translations with no stop codons.
# This file can be combined to that obtained in recipe F3 to
# obtain the sequences

mysql --host $SERVER --user $USER --port $PORT \
	"ensembl_compara_plants_${EGRELEASE}_${RELEASE}" \
    -e "SELECT sm.stable_id \
		FROM seq_member sm, gene_member gm, genome_db gdb \
		WHERE sm.seq_member_id = gm.canonical_member_id \
		AND sm.genome_db_id = gdb.genome_db_id \
		AND gdb.name = '$SPECIES' \
		LIMIT 10"

## S4) Get variants significantly associated to phenotypes

# Variation schema documented at 
# http://www.ensembl.org/info/docs/api/variation/variation_schema.html

SPECIESVAR=$(mysql --host $SERVER --user $USER --port $PORT \
	-e "show databases" | grep "${SPECIES}_variation_${EGRELEASE}_${RELEASE}")

mysql --host $SERVER --user $USER --port $PORT \
    $SPECIESVAR<<SQL
	SELECT f.object_id, s.name, f.seq_region_start, 
		f.seq_region_end, p.description
    FROM phenotype p 
		JOIN phenotype_feature f ON p.phenotype_id = f.phenotype_id  
		JOIN seq_region s ON f.seq_region_id = s.name
    WHERE f.type = 'Variation' AND f.is_significant=1 LIMIT 10
SQL


## S5) Get Triticum aestivum homeologous genes across A,B & D subgenomes

# Compara schema is described at 
# https://m.ensembl.org/info/docs/api/compara/compara_schema.html

# find out the correct method_link_species_set_id in this release
MLSSID=$(mysql --host $SERVER --user $USER --port $PORT \
	ensembl_compara_plants_${EGRELEASE}_$RELEASE -Nb -e \
	"SELECT method_link_species_set_id \
	FROM method_link_species_set \
	WHERE name LIKE 'T%aes%homoeologues'")

# actually retrieve the homeologues using the MLSSID retrieved above
# remove LIMIT 10 if you want the complete set
mysql --host $SERVER --user $USER --port $PORT \
    ensembl_compara_plants_${EGRELEASE}_$RELEASE <<SQL
SELECT 
	homology_member.homology_id, cigar_line, perc_cov, perc_id, 
	perc_pos, gene_member.stable_id as genes, gene_member.genome_db_id
FROM
	homology_member 
	INNER JOIN homology USING (homology_id)
	INNER JOIN gene_member USING (gene_member_id)
WHERE method_link_species_set_id = ${MLSSID}
LIMIT 10
SQL

## S6) Count the number of whole-genome alignments of all genomes

# Compara schema is described at 
# https://m.ensembl.org/info/docs/api/compara/compara_schema.html

mysql --host $SERVER --user $USER --port $PORT \
    ensembl_compara_plants_${EGRELEASE}_$RELEASE <<SQL
SELECT
    genome_db.name,
    SUM(type = "LASTZ_NET") AS n_lastz,
    SUM(type = "SYNTENY") AS n_syntenies,
    SUM(type IN ("EPO", "EPO_LOW_COVERAGE")) AS n_multiple
FROM
    genome_db
    JOIN species_set USING (genome_db_id)
    JOIN method_link_species_set USING (species_set_id)
    JOIN method_link USING (method_link_id)
WHERE
    genome_component IS NULL
    AND genome_db.name != "ancestral_sequences"
GROUP BY genome_db_id
ORDER BY genome_db.name
SQL

## S7) Extract all the mutations and consequences for a selected wheat line

# The variation name contains the name of the line where the muation is present. 
# Variation schema documented at 
# http://www.ensembl.org/info/docs/api/variation/variation_schema.html

SPECIES=triticum_aestivum
LINE="Cadenza1441"
SPECIESVAR=$(mysql --host $SERVER --user $USER --port $PORT \
	-e "show databases" | grep "${SPECIES}_variation_${EGRELEASE}_${RELEASE}")
mysql --host $SERVER --user $USER --port $PORT \
    $SPECIESVAR <<SQL
SELECT 
    seq_region.name as CHR, 
    variation_feature.seq_region_start as POS, 
    variation_feature.allele_string,  
    variation.name as ID, 
    transcript_variation.feature_stable_id, 
    transcript_variation.consequence_types, 
    transcript_variation.sift_prediction, 
    transcript_variation.sift_score,
    variation_set.name as variation_set
FROM variation 
JOIN variation_set_variation
    ON variation.variation_id =variation_set_variation.variation_id
JOIN variation_set 
    ON variation_set.variation_set_id = variation_set_variation.variation_set_id 
JOIN variation_feature 
    ON variation_feature.variation_id = variation.variation_id
JOIN transcript_variation 
    ON transcript_variation.variation_feature_id = variation_feature.variation_feature_id
JOIN seq_region
    ON variation_feature.seq_region_id = seq_region.seq_region_id
WHERE 
    variation.name LIKE "${LINE}%" 
SQL

## S8) Get FASTA of repeated sequences from selected species

# This recipe first interrogates the MySQL server (BED file) and 
# then the FTP server (genome FASTA file) and finally produces a 
# FASTA file with repeated sequences. Uses MINLEN to skip short repeats
# and skips repeats starting on some coordinate.
# Note: requires wget, sort, gzip and bedtools

MINLEN=90
SPECIES=oryza_sativa
SPECIESCORE=$(mysql --host $SERVER --user $USER --port $PORT \
	-e "show databases" | grep "${SPECIES}_core_${EGRELEASE}_${RELEASE}") 
BEDFILE=${SPECIES}.repeats.bed
FASTANAME="*${SPECIES^}*.dna.toplevel.fa.gz"
FASTA=${SPECIES}.toplevel.fasta
REPFASTA=${SPECIES}.repeats.fasta

mysql --host $SERVER --user $USER --port $PORT $SPECIESCORE -Nb -e \
	"SELECT sr.name,(r.seq_region_start-1),r.seq_region_end,rc.repeat_class \
	FROM repeat_feature r JOIN seq_region sr JOIN repeat_consensus rc \
	WHERE r.seq_region_id=sr.seq_region_id \
	AND r.repeat_consensus_id=rc.repeat_consensus_id \
	AND (r.seq_region_end-r.seq_region_start+1) > ${MINLEN}" | \
	sort -u -k1,1 -k2,2n > $BEDFILE

# similar to recipe F4
URL="${FTPSERVER}/${DIV}/current/fasta/${SPECIES}/dna/${FASTANAME}"
wget -c $URL -O- | gunzip > $FASTA
bedtools getfasta -name -fi $FASTA -bed $BEDFILE > $REPFASTA

## S9) Get GFF of repeated sequences from selected species

# This recipe first interrogates the MySQL server and produces a GFF file 
# with repeat sequences. Uses MINLEN to skip short repeats.
# Note: requires wget, sort and perl

mysql --host $SERVER --user $USER --port $PORT $SPECIESCORE -Nb -e \
	"SELECT sr.name,rc.repeat_class,'Repeat',r.seq_region_start, \
	r.seq_region_end,r.score,r.seq_region_strand,0,rc.repeat_name \
	FROM repeat_feature r JOIN seq_region sr JOIN repeat_consensus rc \
	WHERE r.seq_region_id=sr.seq_region_id \
	AND r.repeat_consensus_id=rc.repeat_consensus_id \
	AND (r.seq_region_end-r.seq_region_start+1) > 90" | sort -k1,1 -k4,4n | \
	perl -lane 'if($F[5] eq "NULL"){ $F[5]="."}; if($F[6]==1){ $F[6]="+" } else {$F[6]="-" }; $F[8]="Note=\"$F[8]\";"; print join("\t",@F)'

:
