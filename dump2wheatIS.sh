#!/usr/bin/env bash

# usage: dump2wheatIS.sh 42
# produces core & otherfeatures TSV files for https://urgi.versailles.inra.fr/wheatis/
#
# requires $ENSAPIPATH pointing to ensembl-metadata
#
# by Bruno Contreras Moreira EMBL-EBI 2019

if [ $# -eq 0 ]; then
	echo "# ERROR: need ensembl_genomes_version"
	exit
fi

EGRELEASE=$1

PRODSERVER=mysql-ens-meta-prod-1
MIRRSERVER=mysql-eg-mirror 
PLANTCOREDBLIST=plant_list-core-$EGRELEASE.txt
PLANTOTHERDBLIST=plant_list-otherfeatures-$EGRELEASE.txt
WHEATISCOREFILE=transplant-EBI-core-$EGRELEASE.tsv
WHEATISOTHERFILE=transplant-EBI-otherfeatures-$EGRELEASE.tsv

# make lists of plant databases in selected EGRELEASE
export PERL5LIB=$PERL5LIB:$ENSAPIPATH/ensembl-metadata/modules/

perl $ENSAPIPATH/ensembl-metadata/misc_scripts/get_list_databases_for_division.pl \
	$($PRODSERVER details script) -release $EGRELEASE -division plants \
	| grep -P '_core_'  > $PLANTCOREDBLIST

perl $ENSAPIPATH/ensembl-metadata/misc_scripts/get_list_databases_for_division.pl \
        $($PRODSERVER details script) -release $EGRELEASE -division plants \
        | grep -P '_otherfeatures_'  > $PLANTOTHERDBLIST

SQL='
  -- same as used by Dan Bolser
  SELECT
    # ID is created by solr
    "Sequence feature" AS entry_type,
    "Ensembl Plants"   AS database_name,
    stable_id          AS db_id,
    gene.version       AS db_version,
    COALESCE(xref.display_label, stable_id)
                       AS name,
    COALESCE(gene.description, "")
                       AS description,
    CONCAT("http://plants.ensembl.org/",
           (SELECT meta_value FROM meta
            WHERE meta_key = "species.url"
           ), "/Gene/Summary?g=", stable_id
          )            AS url,
    (SELECT meta_value FROM meta
     WHERE meta_key = "species.scientific_name"
    )                  AS species,
    biotype            AS feature_type,
    name               AS sequence_id,
    seq_region_start   AS start_position,
    seq_region_end     AS end_position
  FROM
    seq_region
  INNER JOIN
    gene
  USING
    (seq_region_id)
  LEFT JOIN
    xref
  ON
    display_xref_id = xref_id
'

# extract metadata to TSV files
while read -r db; do
    >&2 echo $db
    $MIRRSERVER $db -Ne "$SQL"
done < $PLANTCOREDBLIST >> $WHEATISCOREFILE

while read -r db; do
    >&2 echo $db
    $MIRRSERVER $db -Ne "$SQL"
done < $PLANTOTHERDBLIST >> $WHEATISOTHERFILE

exit
