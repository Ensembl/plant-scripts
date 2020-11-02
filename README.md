
These scripts can be used to call and annotate repeated sequences in genomes using the [Repeat detector](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0654-5) (Red) and the library of repeats [nrTEplants](https://github.com/Ensembl/plant_tools/tree/master/bench/repeat_libs)

## Dependencies

Binaries: 

* A clone of Red from https://github.com/brunocontrerasmoreira/Red ; the original repo is [here](https://github.com/BioinformaticsToolsmith/Red)

* A clone of minimap2 from https://github.com/lh3/minimap2

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

python3 -m pip install sqlalchemy sqlalchemy_utils PyMySQL
```

## Example

This script can be used to obtain single-copy core genes present within a clade.
Example calls include:

```
perl ens_single-copy_core_genes.pl -c Brassicaceae -f Brassicaceae
perl ens_single-copy_core_genes.pl -c Brassicaceae -f Brassicaceae -t cdna -o beta_vulgaris
perl ens_single-copy_core_genes.pl -f poaceae -c 4479 -r oryza_sativa -WGA 75
perl ens_single-copy_core_genes.pl -f all -c 33090 -m all -r physcomitrella_patens
```

