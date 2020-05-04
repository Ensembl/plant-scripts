
# set FTP server
SERVER=ftp://ftp.ensemblgenomes.org/pub/plants

# set Ensembl and Ensembl Genomes release numbers,
# check them out at the bottom of the Web sites
RELEASE=100       # http://www.ensembl.org
EGRELEASE=47      # http://plants.ensembl.org

# Ensembl Plants FTP server can be found at
# ftp://ftp.ensemblgenomes.org/pub/plants/current
URL=${SERVER}/current

# you could also choose a specific version
#URL=${SERVER}/release-${EGRELEASE}

SPECIES=brachypodium_distachyon

# 1) download peptide sequences in FASTA format


# 2) download transcripts in FASTA format

# 3) download raw and repeat-masked genomic sequence

# x) download all homologies in a single TSV file
# See ${URL}/tsv/ensembl-compara/homologies/README.gene_trees.tsv_dumps.txt)
TSVFILE=Compara.${RELEASE}.protein_default.homologies.tsv.gz
wget -c ${URL}/tsv/ensembl-compara/homologies/${TSVFILE} .

