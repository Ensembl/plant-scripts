
These scripts can be used to call and annotate repeated sequences in genomes using the [Repeat detector](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0654-5) (Red) and the library of repeats [nrTEplants](https://github.com/Ensembl/plant_tools/tree/master/bench/repeat_libs)

## Dependencies

Binaries: 

* A clone of Red from https://github.com/brunocontrerasmoreira/Red ; the original repo is [here](https://github.com/BioinformaticsToolsmith/Red)

* A clone of minimap2 from https://github.com/lh3/minimap2

* A copy of nrTEplants

Python3 modules:

* [sqlalchemy](https://pypi.org/project/SQLAlchemy)
* [sqlalchemy_utils](https://pypi.org/project/SQLAlchemy-Utils)
* [pymysql](https://pypi.org/project/PyMySQL)

which can be installed with: 
```
git clone https://github.com/brunocontrerasmoreira/Red.git
cd Red && make bin && make
cd ..

git clone https://github.com/lh3/minimap2.git
cd minimap2 && make
cd ..

wget -c https://github.com/Ensembl/plant_tools/releases/download/Jun2020/nrTEplantsJune2020.fna.bz2
bunzip2 /nrTEplantsJune2020.fna.bz2

python3 -m pip install sqlalchemy sqlalchemy_utils PyMySQL
```

## Examples

For large genomes such as barley or wheat you will need a large amount of RAM (~20GB) to run Red:

```
# local run, saves results in folder 'Camelina_sativa' 
./Red2Ensembl.py Camelina_sativa.Cs.dna.toplevel.fa Camelina_sativa --exe /path/to/Red --cor 4 

# local run & loading repeats in core db
./Red2Ensembl.py Brachypodium_distachyon_v3.0.dna.toplevel.fa Brachypodium_distachyon \
	--exe /path/to/Red --cor 4 --host mysql-ens-plants-prod-1 --user xyz --pw XYZ \
	--port 123 --db brachypodium_distachyon_core_49_102_4
```

The repeats called by Red can be optionally annotated by similarity to sequences in an external FASTA file, such as the library nrTEplants. The script does not load the resulting annotations in a core db just yet:
```
./AnnotRedRepeats.py nrTEplantsJune2020.fna Camelina_sativa --exe /path/to/minimap2 --cor 4
```

