
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

# set example species
SPECIES=arabidopsis_thaliana

echo "EGRELEASE=${EGRELEASE}"
echo


# 1) check currently supported species and schemas

SPECIESCORE=`mysql --host $SERVER --user $USER --port $PORT \
	-e "show databases" | grep "${SPECIES}_core_${EGRELEASE}_${RELEASE}"`

# 2) count protein-coding genes

mysql --host $SERVER --user $USER --port $PORT \
	$SPECIESCORE -e "SELECT COUNT(*) FROM gene WHERE biotype='protein_coding'"
