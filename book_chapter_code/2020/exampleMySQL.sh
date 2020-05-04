
# documentation about Ensembl schemas can be found at 
# http://www.ensembl.org/info/docs/api/index.html

mysql --host mysql-eg-publicsql.ebi.ac.uk --user anonymous \
	--port 4157 arabidopsis_thaliana_core_47_100_11 \
	-e "select * from gene"
