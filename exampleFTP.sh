#!/bin/bash

# Copyright [2020] EMBL-European Bioinformatics Institute

# The recipes below use wget to download files 
# Please change these env variables to use other tools ie curl
EXE="wget"
ARGSDEF=" -c " 
ARGSTOFILE=" -O "
ARGSTDOUT=" --quiet $ARGSTOFILE - "

# set servers & division
SERVER="ftp://ftp.ensemblgenomes.org/pub"
DIV=plants
BIOMARTSERVICE="http://plants.ensembl.org/biomart/martservice"

# get Ensembl Plants current release number
SUMFILE="${SERVER}/${DIV}/current/summary.txt"
RELEASE=$($EXE $ARGSTDOUT $SUMFILE | \
	perl -lne 'if(/Release (\d+) of Ensembl/){ print $1 }')

# work out Ensembl Genomes release
EGRELEASE=$(( RELEASE - 53));

# alternatively set a different Ensembl Genomes (EG) release
# EGRELEASE=

# optional arguments, if any
OPTARG=$1

echo "EGRELEASE=${EGRELEASE} OPTARG=${OPTARG}"
echo

# set example species
SPECIES=Brachypodium_distachyon

## F1) Download peptide sequences in FASTA format

FASTAPEP="${SPECIES}*pep.all.fa.gz"
URL="${SERVER}/release-${EGRELEASE}/${DIV}/fasta/${SPECIES,,}/pep/${FASTAPEP}"
echo "# downloading $URL"
$EXE $OPTARG $ARGSDEF $URL

## F2) Download CDS nucleotide sequences in FASTA format

FASTACDS="${SPECIES}*cds.all.fa.gz"
URL="${SERVER}/release-${EGRELEASE}/${DIV}/fasta/${SPECIES,,}/cds/${FASTACDS}"
echo "# downloading $URL"
$EXE $OPTARG $ARGSDEF $URL 

## F3) Download transcripts (cDNA) in FASTA format

FASTACDNA="${SPECIES}*cdna.all.fa.gz"
URL="${SERVER}/release-${EGRELEASE}/${DIV}/fasta/${SPECIES,,}/cdna/${FASTACDNA}"
echo "# downloading $URL"
$EXE $OPTARG $ARGSDEF $URL

## F4) Download soft-masked genomic sequences

FASTASM="${SPECIES}*.dna_sm.toplevel.fa.gz"
URL="${SERVER}/release-${EGRELEASE}/${DIV}/fasta/${SPECIES,,}/dna/${FASTASM}"
echo "# downloading $URL"
$EXE $OPTARG $ARGSDEF $URL

## F5) Upstream/downstream sequences

# Note: this is actually a precompiled BioMart query.
# You can construct your queries at http://plants.ensembl.org/biomart/martview
# and export them as XML

MARTSPECIES=bdistachyon_eg_gene
BIOMARTQUERY=$(cat <<-XMLQUERY
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "plants_mart" formatter = "FASTA" header = "0" uniqueRows = "0" count = "0" datasetConfigVersion = "0.6" >
	<Dataset name = "$MARTSPECIES" interface = "default" >
		<Filter name = "chromosome_name" value = "5"/>        
		<Filter name = "upstream_flank" value = "100"/>
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "5utr" />
	</Dataset>
</Query>
XMLQUERY
)

FASTAUP="${SPECIES}.upstream_flank100.chr5.fa"
URL="${BIOMARTSERVICE}?query=$BIOMARTQUERY"
echo "# downloading $FASTAUP"
if [[ $OPTARG == "--spider" ]]; then
	echo "# skip this recipe in test"
	echo
else
	$EXE $OPTARG $ARGSDEF "$URL" $ARGSTOFILE $FASTAUP
fi

## F6) Get mappings to UniProt proteins

UNIPTSV="${SPECIES}*.uniprot.tsv.gz"
URL="${SERVER}/${DIV}/release-${EGRELEASE}/tsv/${SPECIES,,}/$UNIPTSV"
echo "# downloading $URL"
$EXE $OPTARG $ARGSDEF $URL

## F7) Get indexed, bgzipped VCF file with variants mapped

# Note: this file contains all variants known to Ensembl Plants,
# individual genotypes are not necessarily conserved

VCF="${SPECIES,,}.vcf.gz*"
URL="${SERVER}/${DIV}/release-${EGRELEASE}/variation/vcf/${SPECIES,,}/${VCF}"
echo "# downloading $URL"
$EXE $OPTARG $ARGSDEF $URL

# wheat is an exception, as you can tell from the VCF file which EMS lines
# share a certain mutation, as in this excerpt:
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
#1A      238016  Cadenza0202.chr1A.238016        G       A       .       .       EMS-induced mutation;TSA=SNV
#1A      238016  Cadenza0230.chr1A.238016        G       A       .       .       EMS-induced mutation;TSA=SNV
#1A      238016  Cadenza1874.chr1A.238016        G       A       .       .       EMS-induced mutation;TSA=SNV
#1A      406098  Cadenza0148.chr1A.406098        T       C       .       .       EMS-induced mutation;TSA=SNV
#1A      406098  Cadenza0877.chr1A.406098        T       C       .       .       EMS-induced mutation;TSA=SNV
#1A      406098  Cadenza1340.chr1A.406098        T       C       .       .       EMS-induced mutation;TSA=SNV

## F8) Get precomputed VEP cache files

SPECIES=arabidopsis_thaliana
VEPCACHE="${SPECIES,,}*.tar.gz*"
URL="${SERVER}/${DIV}/release-${EGRELEASE}/variation/vep/${VEPCACHE}"
echo "# downloading $URL"
$EXE $OPTARG $ARGSDEF $URL

# Note: you can get indexed cached files instead from 
# URL=${SERVER}/${DIV}/release-${EGRELEASE}/variation/indexed_vep_cache/${VEPCACHE}

## F9) Download all homologies in a single TSV file, several GBs

TSVFILE="Compara.${RELEASE}.protein_default.homologies.tsv.gz"
URL="${SERVER}/${DIV}/release-${EGRELEASE}/tsv/ensembl-compara/homologies/${TSVFILE}"
echo "# downloading $URL"
$EXE $OPTARG $ARGSDEF $URL

# Note: you can extract homologies from this file by parsing it
# in the command line. Example:
# zcat $TSVFILE | grep triticum_aestivum | grep oryza_sativa | grep ortholog 

# Note: homologies of each species can be retrieved from a more specific file
# SPECIES=Triticum_aestivum
#URL="${SERVER}/${DIV}/release-${EGRELEASE}/tsv/ensembl-compara/homologies/${SPECIES,,}${TSVFILE}"
#wget -c "$URL"
#zcat "$TSVFILE" | grep oryza_sativa | grep ortholog

# Note: Alternatively a smaller file in OrthoXML format can be obtained
# OXMLFILE="Compara.${RELEASE}.protein_default.allhomologies.orthoxml.xml.gz"
# URL="${SERVER}/${DIV}/release-${EGRELEASE}/xml/ensembl-compara/homologies/${OXMLFILE}"

## F10) download UniProt report of Ensembl Plants, 
# summarized how many protein sequences from each species
# have been annotated in SwissProt & TrEMBL

UNIPFILE="uniprot_report_EnsemblPlants.txt"
URL="${SERVER}/${DIV}/release-${EGRELEASE}/$UNIPFILE"
echo "# downloading $URL"
$EXE $OPTARG $ARGSDEF $URL

## F11) Retrieve list of new species in current release

NEWLIST="new_genomes.txt"
URL="${SERVER}/${DIV}/release-${EGRELEASE}/$NEWLIST"
echo "# downloading $URL"
$EXE $OPTARG $ARGSDEF $URL

## F12) Get current plant species tree (cladogram)
TREEFILE="plants_species-tree_Ensembl.nh"
URL="${SERVER}/${DIV}/release-${EGRELEASE}/compara/species_trees/$TREEFILE"
echo "# downloading $URL"
$EXE $OPTARG $ARGSDEF $URL

