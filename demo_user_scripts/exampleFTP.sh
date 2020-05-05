#!/bin/bash

# Your usage of the data returned by the REST service is
# subject to same conditions as laid out on the Ensembl website.
#
# Copyright [2020] EMBL-European Bioinformatics Institute

# set FTP server & division
SERVER=ftp://ftp.ensemblgenomes.org/pub
DIV=plants

# set Ensembl Plants release number, check it 
# at the bottom of http://plants.ensembl.org
EGRELEASE=47      

RELEASE=$(( $EGRELEASE + 53)); # do not change

# set example species
SPECIES=Brachypodium_distachyon

echo "EGRELEASE=${EGRELEASE}"
echo

# 1) download peptide sequences in FASTA format

FASTAPEP=${SPECIES}*pep.all.fa.gz
URL=${SERVER}/release-${EGRELEASE}/${DIV}/fasta/${SPECIES,,}/pep/${FASTAPEP}
echo "# downloading $URL"
wget -c $URL 

# 2) download CDS nucleotide sequences in FASTA format

FASTACDS=${SPECIES}*cds.all.fa.gz
URL=${SERVER}/release-${EGRELEASE}/${DIV}/fasta/${SPECIES,,}/cds/${FASTACDS}
echo "# downloading $URL"
wget -c $URL 

# 3) download transcripts (cDNA) in FASTA format

FASTACDNA=${SPECIES}*cdna.all.fa.gz
URL=${SERVER}/release-${EGRELEASE}/${DIV}/fasta/${SPECIES,,}/cdna/${FASTACDNA}
echo "# downloading $URL"
wget -c $URL 

# 4) download soft-masked genomic sequences

FASTASM=${SPECIES}*.dna_sm.toplevel.fa.gz
URL=${SERVER}/release-${EGRELEASE}/${DIV}/fasta/${SPECIES,,}/dna/${FASTASM}
echo "# downloading $URL"
wget -c $URL

# 5) download all homologies in a single TSV file, several GBs

TSVFILE=Compara.${RELEASE}.protein_default.homologies.tsv.gz
URL=${SERVER}/${DIV}/release-${EGRELEASE}/tsv/ensembl-compara/homologies/${TSVFILE}
echo "# downloading $URL"
wget -c $URL 

