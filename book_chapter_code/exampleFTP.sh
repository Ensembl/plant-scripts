
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
wget -c ${SERVER}/release-${EGRELEASE}/${DIV}/fasta/${SPECIES,,}/pep/${FASTAPEP} .

# 2) download transcripts in FASTA format



# 3) download raw and repeat-masked genomic sequence

# x) download all homologies in a single TSV file, several GBs

TSVFILE=Compara.${RELEASE}.protein_default.homologies.tsv.gz
wget -c ${SERVER}/${DIV}/release-${EGRELEASE}/tsv/ensembl-compara/homologies/${TSVFILE} .

