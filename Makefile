test:
	perl demo_test.t

test_travis:
	perl demo_test.t travis

clean:
	rm -f *rachypodium* && rm -f Compara*gz 
	rm -f new_genomes.txt && rm -f uniprot_report_EnsemblPlants.txt
	rm -f arabidopsis_thaliana*.tar.gz
	rm -f plants_species-tree*.nh
	rm -f oryza_sativa*

install:
	sudo apt-get install -y wget mysql-client libmysqlclient-dev bedtools

install_REST:
	cpanm --local-lib lib --installdeps --notest --cpanfile lib/cpanfileREST .
	pip3 install --user requests

install_biomart_r:
	Rscript install_R_deps.R

install_ensembl:
	cpanm --local-lib lib --installdeps --notest --cpanfile lib/cpanfileEnsembl .
	cd lib && git clone https://github.com/Ensembl/ensembl.git
	cd lib && git clone https://github.com/Ensembl/ensembl-variation.git
	cd lib && git clone https://github.com/Ensembl/ensembl-funcgen.git
	cd lib && git clone https://github.com/Ensembl/ensembl-compara.git
	cd lib && git clone https://github.com/Ensembl/ensembl-metadata.git
	cd lib && git clone -b release-1-6-924 --depth 1 https://github.com/bioperl/bioperl-live.git

install_repeats:
	pip3 install --user -r lib/requirements.txt
	cd lib && git clone https://github.com/EnsemblGenomes/Red.git && cd Red/src_2.0 && make bin && make CXX=g++-10
	cd lib && git clone https://github.com/lh3/minimap2.git && cd minimap2 && make
	cd files && wget -c https://github.com/Ensembl/plant-scripts/releases/download/v0.3/nrTEplantsJune2020.fna.bz2 && bunzip2 nrTEplantsJune2020.fna.bz2

install_redat:
	cd files && wget -c ftp://ftpmips.helmholtz-muenchen.de/plants/REdat/mipsREdat_9.3p_ALL.fasta.gz && gunzip mipsREdat_9.3p_ALL.fasta.gz

test_repeats_travis:
	cd repeats && ./Red2Ensembl.py ../files/Arabidopsis_thaliana.fna.gz test_Atha_chr4 --msk_file Atha.sm.fna

test_repeats:
	cd repeats && ./Red2Ensembl.py ../files/Arabidopsis_thaliana.fna.gz test_Atha_chr4 --msk_file Atha.sm.fna && ./AnnotRedRepeats.py ../files/nrTEplantsJune2020.fna test_Atha_chr4 --bed_file test.nrTEplants.bed

uninstall_repeats:
	cd files && rm -rf nrTEplantsJune2020.fna*
	cd lib && rm -rf Red minimap2 

clean_repeats:
	cd repeats && rm -rf test_Atha_chr4 Atha.sm.fna test.nrTEplants.bed
