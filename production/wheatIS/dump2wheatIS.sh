#!/usr/bin/env bash

# usage: dump2wheatIS.sh 42
# produces core & otherfeatures TSV files for https://urgi.versailles.inra.fr/wheatis/
#
# requires $ENSAPIPATH pointing to ensembl-metadata
#
# by Bruno Contreras Moreira, Raphael Flores EMBL-EBI 2019-20

if [ $# -eq 0 ]; then
	echo "# ERROR: need ensembl_genomes_version"
	exit
fi

EGRELEASE=$1

PRODSERVER=mysql-ens-meta-prod-1
MIRRSERVER=mysql-ens-sta-3 
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
  -- matches https://urgi.versailles.inra.fr/wheatis/join#tsv-tabulation-separated-values
  SELECT
    # Unicity within a dataset is handled by the `name` field
    COALESCE(xref.display_label, stable_id)
                        AS name,
    "Genome annotation" AS entryType,
    "EBI"               AS node,
    "Ensembl Plants"    AS databaseName,
    CONCAT("http://plants.ensembl.org/",
           (SELECT meta_value FROM meta
            WHERE meta_key = "species.url"
           ), "/Gene/Summary?g=", stable_id
          )             AS url,
    (SELECT meta_value FROM meta
     WHERE meta_key = "species.display_name"
    )                   AS species,
    CONCAT(COALESCE(gene.description, ""),
     " feature type = ", COALESCE(biotype,"n/a")
     )                  AS description
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

# Add missing URL param for otherfeatures databases
perl -p -i -e 's/(Gene\/Summary\?)/$1db=otherfeatures;/' $WHEATISOTHERFILE

exit
