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

OPTARG=$1

# set example species
SPECIES=Brachypodium_distachyon
SPECIESBIOMART=bdistachyon_eg_genes

echo "EGRELEASE=${EGRELEASE} OPTARG=$OPTARG"
echo

## F1) download peptide sequences in FASTA format

FASTAPEP=${SPECIES}*pep.all.fa.gz
URL=${SERVER}/release-${EGRELEASE}/${DIV}/fasta/${SPECIES,,}/pep/${FASTAPEP}
echo "# downloading $URL"
wget $OPTARG -c $URL 

## F2) download CDS nucleotide sequences in FASTA format

FASTACDS=${SPECIES}*cds.all.fa.gz
URL=${SERVER}/release-${EGRELEASE}/${DIV}/fasta/${SPECIES,,}/cds/${FASTACDS}
echo "# downloading $URL"
wget $OPTARG -c $URL 

## F3) download transcripts (cDNA) in FASTA format

FASTACDNA=${SPECIES}*cdna.all.fa.gz
URL=${SERVER}/release-${EGRELEASE}/${DIV}/fasta/${SPECIES,,}/cdna/${FASTACDNA}
echo "# downloading $URL"
wget $OPTARG -c $URL 

## F4) download soft-masked genomic sequences

FASTASM=${SPECIES}*.dna_sm.toplevel.fa.gz
URL=${SERVER}/release-${EGRELEASE}/${DIV}/fasta/${SPECIES,,}/dna/${FASTASM}
echo "# downloading $URL"
wget $OPTARG -c $URL


## F5) Upstream/downstream sequences

# You can obtained these in http://plants.ensembl.org/biomart/martview
# by constructing queries which result in URLs such as this for 
# Coffea canephora:

#http://plants.ensembl.org/biomart/martview?VIRTUALSCHEMANAME=plants_mart&ATTRIBUTES=ccanephora_eg_gene.default.sequences.ensembl_transcript_id|ccanephora_eg_gene.default.sequences.5utr|ccanephora_eg_gene.default.sequences.upstream_flank."500"&FILTERS=&VISIBLEPANEL=resultspanel


## F6) get mappings to UniProt proteins

UNIPTSV=${SPECIES}*.uniprot.tsv.gz
URL=${SERVER}/${DIV}/release-${EGRELEASE}/tsv/${SPECIES,,}//$UNIPTSV
echo "# downloading $URL"
wget $OPTARG -c $URL

## F7) get indexed, bgzipped VCF file with variants mapped

# Note: this contains all variants known to Ensembl Plants,
# no individual genotypes are conserved

VCF=${SPECIES,,}.vcf.gz*
URL=${SERVER}/${DIV}/release-${EGRELEASE}/variation/vcf/${SPECIES,,}/${VCF}
echo "# downloading $URL"
wget $OPTARG -c $URL

## F8) get precomputed VEP cache files

VEPCACHE=${SPECIES,,}*.tar.gz*
URL=${SERVER}/${DIV}/release-${EGRELEASE}/variation/vep/${VEPCACHE}
echo "# downloading $URL"
wget $OPTARG -c $URL

# Note: you can get indexed cached files instead from 
# URL=${SERVER}/${DIV}/release-${EGRELEASE}/variation/indexed_vep_cache/${VEPCACHE}

## F9) download all homologies in a single TSV file, several GBs

TSVFILE=Compara.${RELEASE}.protein_default.homologies.tsv.gz
URL=${SERVER}/${DIV}/release-${EGRELEASE}/tsv/ensembl-compara/homologies/${TSVFILE}
echo "# downloading $URL"
wget $OPTARG -c $URL

# Note: you can extract homologies from this file by parsing it
# in the command line. Example:
# zcat $TSVFILE | grep triticum_aestivum | grep oryza_sativa | grep ortholog 


## F10) download UniProt report of Ensembl Plants, 
# summarized how many protein sequences from each species
# have been annotated in SwissProt & TrEMBL

UNIPFILE=uniprot_report_EnsemblPlants.txt
URL=${SERVER}/${DIV}/release-${EGRELEASE}/$UNIPFILE
echo "# downloading $URL"
wget $OPTARG -c $URL

## F11) retrieve list of new species in current release

NEWLIST=new_genomes.txt
URL=${SERVER}/${DIV}/release-${EGRELEASE}/$NEWLIST
echo "# downloading $URL"
wget $OPTARG -c $URL

## F12) get current plant species tree (cladogram)
TREEFILE=/plants_protein-trees_default.nh
URL=${SERVER}/${DIV}/release-${EGRELEASE}/compara/species_trees/$TREEFILE
echo "# downloading $URL"
wget $OPTARG -c $URL

