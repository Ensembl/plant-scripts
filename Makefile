test_demo:
	cd demo_scripts && perl demo_test.t

clean_demo:
	cd demo_scripts && rm -f *rachypodium* && rm -f Compara*gz 
	cd demo_scripts && rm -f new_genomes.txt && rm -f uniprot_report_EnsemblPlants.txt
	cd demo_scripts && rm -f arabidopsis_thaliana*.tar.gz
	cd demo_scripts && rm -f plants_species-tree*.nh

test_phylo:
	cd phylogenomics && perl phylo_test.t
