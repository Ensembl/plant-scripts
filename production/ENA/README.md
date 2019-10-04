
This script checks the current assemblies at ENA matching -t taxon name and compares them to the current contents of Ensembl division (default=Plants). It creates a TSV report which looks like this:

| ENA_assembly_id | scientific_name | strain | chromosomes | GenBank | genes | RefSeq | genes |
|-----------------|-----------------|--------| ------------|---------|-------|--------|-------|
|GCA_000001735.2|Arabidopsis thaliana|NA|7|GCA_000001735.2_TAIR10.1_genomic.gff.gz|33268|GCF_000001735.4_TAIR10.1_genomic.gff.gz|33467|
|GCA_000321445.1|Oryza sativa Japonica Group|Hitomebore|12|GCA_000321445.1_Osat_hitom_01_genomic.gff.gz|0|NA|0|
