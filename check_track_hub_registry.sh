#!/usr/local/bin/bash
# usage: run without arguments, but redirect output (is large)
# some errors will happen (grep: genomes.txt: No such file or directory)
# Original by Marc Rosello

PLANTHUBPATH="/nfs/ensemblgenomes/ftp/pub/misc_data/Track_Hubs"

FILES=${PLANTHUBPATH}/*
COUNTER=0
for f in $FILES
do
  study=$(echo $f | sed "s/^.*Hubs\///")
  COUNTER=$((COUNTER + 1))
  cd $f
  genome_count=$(grep -c "trackDb" genomes.txt)
  for trackDb in $(grep "trackDb" genomes.txt | awk '{print $2}')
  do
      #PRINT STUDY ID AND GENOME COUNT
      echo -n "$study	"
      echo -n "$genome_count	" #not needed, can count recurring study ids

      #PRINT GENOME NAME OF GENOME COUNT
      genome_name=$(echo $trackDb | sed 's/\/trackDb.txt//')
      echo -n "$genome_name	"

      #PRINT SPECIES LIST OF GENOME NAME
      spec_list=$(egrep -o "scientific_name=\"[^\"]+" $trackDb | sed 's/scientific_name="//' | sort -u)
      while read -r line;do
	  echo -n "$line,"
      done <<< "$spec_list" #for loop doesn't see each line as separate iteration
      echo -n "	"
      #COUNT TRACKS
      tracksC=$(grep -c "parent " $trackDb)
      tracks=$(grep "parent " $trackDb)
      samplesC=$(grep "parent " $trackDb | sort -u | wc -l)
      echo "$tracksC tracks over $samplesC samples"
  done

#  if [ "$COUNTER" -gt 20 ];then
#      break
#  fi

done
echo ""

exit
