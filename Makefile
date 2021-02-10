test:
	perl demo_test.t

clean:
	rm -f *rachypodium* && rm -f Compara*gz 
	rm -f new_genomes.txt && rm -f uniprot_report_EnsemblPlants.txt
	rm -f arabidopsis_thaliana*.tar.gz
	rm -f plants_species-tree*.nh
	rm -f oryza_sativa*

installR:
	Rscript install_R_deps.R

install_repeats:
	pip install --user -r files/pythonlist
	cd files && git clone https://github.com/EnsemblGenomes/Red.git && cd Red/src_2.0 && make bin && make
	cd files && git clone https://github.com/lh3/minimap2.git && cd minimap2 && make
	cd files && wget -c https://github.com/Ensembl/plant_tools/releases/download/v0.3/nrTEplantsJune2020.fna.bz2 && bunzip2 nrTEplantsJune2020.fna.bz2

clean_repeats:
	cd files && rm -rf Red minimap2 nrTEplantsJune2020.fna*
