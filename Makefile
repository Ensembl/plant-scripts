test:
	perl demo_test.t

clean:
	rm -f *rachypodium* && rm -f Compara*gz 
	rm -f new_genomes.txt && rm -f uniprot_report_EnsemblPlants.txt
	rm -f arabidopsis_thaliana*.tar.gz
	rm -f plants_species-tree*.nh
	rm -f oryza_sativa*

install:
	Rscript install_R_deps.R
	cd repeats && git clone https://github.com/EnsemblGenomes/Red.git && d Red && make bin && make
	cd repeats && git clone https://github.com/lh3/minimap2.git && cd minimap2 && make
	cd repeats wget -c https://github.com/Ensembl/plant_tools/releases/download/Jun2020/nrTEplantsJune2020.fna.bz2 && bunzip2 /nrTEplantsJune2020.fna.bz2

