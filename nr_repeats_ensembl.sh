#!/usr/bin/env bash

# Script to download non-redundant annotated repeated elements 
# of a selected species in Ensembl Plants

# Required binaries: wget, sort, perl, mysql, bedtools, cd-hit-est

# Copyright [2020] EMBL-European Bioinformatics Institute

# documentation about Ensembl schemas can be found at 
# http://www.ensembl.org/info/docs/api/index.html

if [[ $# -eq 0 ]] ; then
	echo "# example usage: $0 arabidopsis_thaliana [MB RAM]"
	exit 0
else
	SPECIES=$1
fi

# PARAMS
MINLEN=90
MAXDEGENPERC=10
MAXIDFRAC=0.95
DEBUG=0

# SERVER DETAILS
FTPSERVER="ftp://ftp.ensemblgenomes.org/pub"
DIV=plants
SERVER=mysql-eg-publicsql.ebi.ac.uk
USER=anonymous
PORT=4157

## 1) get Ensembl Plants current release number from FTP server
# Note: wget is used, this can be modified to use alternatives ie curl
SUMFILE="${FTPSERVER}/${DIV}/current/summary.txt"
RELEASE=`wget --quiet -O - $SUMFILE | \
	perl -lne 'if(/Release (\d+) of Ensembl/){ print $1 }'`

# work out Ensembl Genomes release
EGRELEASE=$(( RELEASE - 53));

## 2) select core db matching selected species
SPECIESCORE=$(mysql --host $SERVER --user $USER --port $PORT \
	-e "show databases" | grep "${SPECIES}_core_${EGRELEASE}_${RELEASE}")

if [ -z "$SPECIESCORE" ]; then
	echo "# ERROR: cannot find species $SPECIES"
	exit 1
else
	echo "# Ensembl core db: $SPECIESCORE";
fi

## 3) retrieve 1-based coords of repeats

# note these might be redundant/overlapping
#1       3       106     trf
#1       4       91      trf

mysql --host $SERVER --user $USER --port $PORT $SPECIESCORE -Nb -e \
	"SELECT sr.name,r.seq_region_start,r.seq_region_end,rc.repeat_class \
	FROM repeat_feature r JOIN seq_region sr JOIN repeat_consensus rc \
	WHERE r.seq_region_id=sr.seq_region_id \
	AND r.repeat_consensus_id=rc.repeat_consensus_id \
	AND (r.seq_region_end-r.seq_region_start+1) > $MINLEN" | \
	sort -u -k1,1 -k2,2n > _${SPECIES}.repeats1.bed

## 4) retrieve 1-based coords of genes
mysql --host $SERVER --user $USER --port $PORT $SPECIESCORE -Nb -e \
	"SELECT sr.name,g.seq_region_start,g.seq_region_end,g.stable_id \
	FROM gene g JOIN seq_region sr \
	WHERE g.seq_region_id=sr.seq_region_id" | \
	sort -k1,1 -k2,2n > _${SPECIES}.genes1.bed

## 5) curate repeats by substracting annotated genes and
##    convert to 0-based BED format
bedtools subtract -sorted \
	-a _${SPECIES}.repeats1.bed -b _${SPECIES}.genes1.bed | \
	perl -lane '$F[1]-=1; print join("\t",@F)' >\
	_${SPECIES}.repeats.bed

if [ ! -s  _${SPECIES}.repeats.bed ]; then
	echo "# no repeats found"
	exit 2
fi

## 6) download and uncompress genomic sequence 
FASTA="*${SPECIES^}*.dna.toplevel.fa.gz"
URL="${FTPSERVER}/${DIV}/current/fasta/${SPECIES}/dna/${FASTA}"
if [ ! -s  _${SPECIES}.toplevel.fasta ]; then
	echo "# downloading $URL"
	wget -c $URL -O- | gunzip > _${SPECIES}.toplevel.fasta
else
	echo "# re-using _${SPECIES}.toplevel.fasta"
fi

## 7) extract repeat sequences 
bedtools getfasta -name -fi _${SPECIES}.toplevel.fasta -bed _${SPECIES}.repeats.bed >\
	_${SPECIES}.repeats.fasta

## 8) eliminate degenerate (MAXDEGENPERC) repeat sequences 
cat _${SPECIES}.repeats.fasta | \
	perl -slne 'if(/^(>.*)/){$h=$1} else {$fa{$h}.=$_} END{ foreach $h (keys(%fa)){ $l=length($fa{$h}); $dg=($fa{$h}=~tr/Nn//); print "$h\n$fa{$h}" if(100*$dg/$l<=$maxdeg) }}' \
	-- -maxdeg=$MAXDEGENPERC > _${SPECIES}.repeats.nondeg.fasta

## 9) eliminate short & redundant sequences
if [[ $# -eq 2 ]] ; then
	RAM=$2
else
	RAM=1024
fi

OUTFILE=${SPECIES}.${EGRELEASE}.repeats.nr${MAXIDFRAC}.fasta

cd-hit-est -M $RAM -c $MAXIDFRAC -l $MINLEN \
	-i _${SPECIES}.repeats.nondeg.fasta \
	-o $OUTFILE

if [ -s "$OUTFILE" ]; then
	echo "# output: $OUTFILE"
else
	echo "# job failed, temp files conserved"
	exit 3
fi

## 10) clean temp files
if [ -z "$DEBUG" ] || [ "$DEBUG" -eq "0" ]; then
	echo
	echo "# removing temp files"; 
	rm _${SPECIES}.*.bed _${SPECIES}.*.fasta _${SPECIES}.*.fai ${SPECIES}.*.clstr
fi

exit 0
