
This script checks the current assemblies at ENA matching -t taxon name and compares them to the current contents of Ensembl division (default=Plants). It creates a TSV report which looks like this:

# ENA_assembly_id	scientific_name	annotated_genes	chromosomes	chromosome_accessions
GCA_000001735.2	Arabidopsis thaliana	37964	7	CP002684.1,CP002685.1,CP002686.1,CP002687.1,CP002688.1,BK010421.1,AP000423.1
GCA_000182155.3	Oryza barthii	0	12	CM002639.1,CM002643.2,CM002644.1,CM002645.1,CM002646.2,CM002647.1,CM002648.1,CM002649.1,CM002650.
1,CM002640.2,CM002641.2,CM002642.2


