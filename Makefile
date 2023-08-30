
minimap2release = 2.24
gffreadrelease  = 0.12.7
gmaprelease     = 2021-12-17
clustalorelease = 1.2.4
alistatrelease  = 1.14

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
	-sudo apt install -y wget mysql-client libmysqlclient-dev libdb-dev bedtools pip cpanminus

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

install_minimap2:
	if [ ! -d "lib/minimap2" ]; then \
		cd lib && wget https://github.com/lh3/minimap2/releases/download/v${minimap2release}/minimap2-${minimap2release}.tar.bz2 && \
			tar xfj minimap2-${minimap2release}.tar.bz2 && cd minimap2-${minimap2release} && make && cd .. && \
			rm -f minimap2-${minimap2release}.tar.bz2 && ln -fs minimap2-${minimap2release} minimap2; \
	fi

install_Red:
	cd lib && git clone https://github.com/EnsemblGenomes/Red.git && cd Red/src_2.0 && make bin && make
	#in case you need to use an alternative g++ compiler
        #cd lib && git clone https://github.com/EnsemblGenomes/Red.git && cd Red/src_2.0 && make bin && make CXX=g++-10

install_repeats: install_minimap2 install_Red
	pip3 install --user -r lib/requirements.txt
	cd files && wget -c https://github.com/Ensembl/plant-scripts/releases/download/v0.3/nrTEplantsJune2020.fna.bz2 && bunzip2 nrTEplantsJune2020.fna.bz2

install_redat:
	cd files && wget -c ftp://ftpmips.helmholtz-muenchen.de/plants/REdat/mipsREdat_9.3p_ALL.fasta.gz && gunzip mipsREdat_9.3p_ALL.fasta.gz

test_repeats_travis:
	cd repeats && ./Red2Ensembl.py ../files/Arabidopsis_thaliana.fna.gz test_Atha_chr4 --msk_file Atha.sm.fna

test_repeats:
	cd repeats && ./Red2Ensembl.py ../files/Arabidopsis_thaliana.fna.gz test_Atha_chr4 --msk_file Atha.sm.fna && \
		./AnnotRedRepeats.py ../files/nrTEplantsJune2020.fna test_Atha_chr4 --bed_file test.nrTEplants.bed

uninstall_repeats:
	cd files && rm -rf nrTEplantsJune2020.fna*
	cd lib && rm -rf Red minimap2-${minimap2release} minimap2

clean_repeats:
	cd repeats && rm -rf test_Atha_chr4 Atha.sm.fna test.nrTEplants.bed

# gmap takes several minutes to compile
install_gmap: 
	cd pangenes/bin && wget http://research-pub.gene.com/gmap/src/gmap-gsnap-${gmaprelease}.tar.gz && tar xfz gmap-gsnap-${gmaprelease}.tar.gz && \
		cd gmap-${gmaprelease} && ./configure --prefix=${PWD}/pangenes/bin/gmap-${gmaprelease}/exe && \
		make && make install && cd .. && rm -rf gmap-gsnap-${gmaprelease}.tar.gz && ln -fs gmap-${gmaprelease} gmap
	
install_gffread:
	cd pangenes/bin && wget https://github.com/gpertea/gffread/releases/download/v${gffreadrelease}/gffread-${gffreadrelease}.tar.gz && \
		tar xfz gffread-${gffreadrelease}.tar.gz && cd gffread-${gffreadrelease} && make && cd .. && \
		rm -f gffread-${gffreadrelease}.tar.gz && ln -fs gffread-${gffreadrelease} gffread

install_pangenes: install_minimap2 install_gffread install_gmap
        cpanm --sudo -v --installdeps --notest --cpanfile cpanfile .
	cd files && wget -c https://github.com/Ensembl/plant-scripts/releases/download/v0.4/test_rice.tgz && tar xfz test_rice.tgz && rm -f test_rice.tgz

# see https://github.com/ekg/wfmash for other options
install_wfmash:
	-sudo apt install cmake libjemalloc-dev zlib1g-dev libgsl-dev libhts-dev
	cd pangenes/bin && git clone https://github.com/ekg/wfmash && cd wfmash && cmake -H. -Bbuild && cmake --build build -- -j 3

install_gsalign:
	cd pangenes/bin && git clone https://github.com/hsinnan75/GSAlign.git && cd GSAlign && make

install_pangenes_quality:
	cd pangenes/bin && wget http://www.clustal.org/omega/clustalo-${clustalorelease}-Ubuntu-x86_64 && \
	chmod +x clustalo-${clustalorelease}-Ubuntu-x86_64 && \
	ln -fs clustalo-${clustalorelease}-Ubuntu-x86_64 clustalo && \
	wget https://github.com/thomaskf/AliStat/archive/refs/tags/v${alistatrelease}.tar.gz && \
	tar xfz v${alistatrelease}.tar.gz && cd AliStat-${alistatrelease} && make && cd .. && \
	rm -f v${alistatrelease}.tar.gz && ln -s AliStat-${alistatrelease} AliStat

uninstall_pangenes:
	cd pangenes/bin && rm -rf gffread-${gffreadrelease} gmap-${gmaprelease} gffread wfmash GSAlign gmap \
		clustalo-${clustalorelease}-Ubuntu-x86_64 clustalo AliStat-${alistatrelease} AliStat
	cd lib && rm -rf minimap2-${minimap2release} minimap2
	cd files && rm -rf test_rice

test_pangenes:
	cd pangenes && perl get_pangenes.pl -d ../files/test_rice && \
		perl get_pangenes.pl -d ../files/test_rice -s '^\d+$$' && \
		perl get_pangenes.pl -d ../files/test_rice -H

clean_pangenes:
	cd pangenes && rm -rf test_rice_pangenes
