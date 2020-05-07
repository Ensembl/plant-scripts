test_demo:
	cd demo_user_scripts && perl demo_test.t

clean_demo:
	cd demo_user_scripts && rm -f *rachypodium*gz* && rm -f Compara*gz 
	cd demo_user_scripts && rm -f new_genomes.txt && rm -f uniprot_report_EnsemblPlants.txt

test_phylo:
	cd phylogenomics && perl phylo_test.t
