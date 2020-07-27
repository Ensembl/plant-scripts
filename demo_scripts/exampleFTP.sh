#!/bin/bash

# Copyright [2020] EMBL-European Bioinformatics Institute

# The recipes below use wget to download files, 
# you could use curl instead

# set servers & division
SERVER=ftp://ftp.ensemblgenomes.org/pub
DIV=plants
BIOMARTSERVICE=http://plants.ensembl.org/biomart/martservice

# get Ensembl Plants current release number
SUMFILE=${SERVER}/${DIV}/current/summary.txt
RELEASE=`wget --quiet -O - $SUMFILE | \
	perl -lne 'if(/Release (\d+) of Ensembl/){ print $1 }'`

EGRELEASE=$(( $RELEASE - 53));    

# alternatively set other EG release number
# EGRELEASE=

OPTARG=$1

echo "EGRELEASE=${EGRELEASE} OPTARG=${OPTARG}"
echo

# set example species
SPECIES=Brachypodium_distachyon

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

# Note: this is actually a BioMart query.
# You can construct your queries at 
# http://plants.ensembl.org/biomart/martview
# and then export them as XML

BIOMARTQUERY=$(cat <<XMLQUERY
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "plants_mart" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
	<Dataset name = "bdistachyon_eg_gene" interface = "default" >
		<Filter name = "upstream_flank" value = "500"/>
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "5utr" />
	</Dataset>
</Query>
XMLQUERY
)
FASTAUP=${SPECIES}.upstream_flank500.fa
URL=${BIOMARTSERVICE}?query=$BIOMARTQUERY
echo "# downloading $URL"
wget $OPTARG -c "$URL" -O $FASTAUP

## F6) get mappings to UniProt proteins

UNIPTSV=${SPECIES}*.uniprot.tsv.gz
URL=${SERVER}/${DIV}/release-${EGRELEASE}/tsv/${SPECIES,,}//$UNIPTSV
echo "# downloading $URL"
wget $OPTARG -c $URL

## F7) get indexed, bgzipped VCF file with variants mapped

# Note: this contains all variants known to Ensembl Plants,
# individual genotypes are not necessarily conserved

VCF=${SPECIES,,}.vcf.gz*
URL=${SERVER}/${DIV}/release-${EGRELEASE}/variation/vcf/${SPECIES,,}/${VCF}
echo "# downloading $URL"
wget $OPTARG -c $URL

# wheat is an exception, as you can tell from the VCF file which EMS lines
# share a certain mutation, as in this excerpt:
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
#1A      238016  Cadenza0202.chr1A.238016        G       A       .       .       EMS-induced mutation;TSA=SNV
#1A      238016  Cadenza0230.chr1A.238016        G       A       .       .       EMS-induced mutation;TSA=SNV
#1A      238016  Cadenza1874.chr1A.238016        G       A       .       .       EMS-induced mutation;TSA=SNV
#1A      406098  Cadenza0148.chr1A.406098        T       C       .       .       EMS-induced mutation;TSA=SNV
#1A      406098  Cadenza0877.chr1A.406098        T       C       .       .       EMS-induced mutation;TSA=SNV
#1A      406098  Cadenza1340.chr1A.406098        T       C       .       .       EMS-induced mutation;TSA=SNV

## F8) get precomputed VEP cache files

SPECIES=arabidopsis_thaliana
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

# Alternatively a smaller OrthoXML-format file can be obtained
# OXMLFILE=Compara.${RELEASE}.protein_default.allhomologies.orthoxml.xml.gz
# URL=${SERVER}/${DIV}/release-${EGRELEASE}/xml/ensembl-compara/homologies/${OXMLFILE}

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

