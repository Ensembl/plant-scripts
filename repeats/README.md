
These scripts can be used to i) mask and ii) annotate repeated sequences in plants genomes using the 
[Repeat detector](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0654-5) (Red) 
and the curated library of repeats 
[nrTEplants](https://github.com/Ensembl/plant_tools/releases/download/v0.3/nrTEplantsJune2020.fna.bz2),
described [here](https://github.com/Ensembl/plant_tools/tree/master/bench/repeat_libs). 

Optionally, repeats can be loaded into an Ensembl core db as an analysis with logic_name='repeatdetector'. 
Annotated repeats are imported with logic_name='repeatdetector_annotated'.

## Dependencies

The following dependencies can be installed in the parent folder with:

    make install_repeats

There are two binaries, for which the GNU C++ compiler v8 is needed (g++-8):

* A clone of Red from https://github.com/EnsemblGenomes/Red (the original repo is [here](https://github.com/BioinformaticsToolsmith/Red))
* A clone of minimap2 from https://github.com/lh3/minimap2
* A copy of the [nrTEplants library](https://github.com/Ensembl/plant_tools/releases/download/v0.3/nrTEplantsJune2020.fna.bz2)

Plus three Python3 modules:

* [sqlalchemy](https://pypi.org/project/SQLAlchemy)
* [sqlalchemy_utils](https://pypi.org/project/SQLAlchemy-Utils)
* [pymysql](https://pypi.org/project/PyMySQL)


## Examples

Note that the input FASTA file can be GZIP/BZIP2 compressed; 
for large genomes you will need 20-40 GB of RAM to run Red (see below):

```
# test run, saves results in folder 'test_Atha_chr4' 
./Red2Ensembl.py ../files/Arabidopsis_thaliana.fna.gz test_Atha_chr4 --exe /path/to/Red

# real example, with several chromosomes, 4 cores and saved masked sequences
./Red2Ensembl.py Brachypodium_distachyon_v3.0.dna.toplevel.fa Brachypodium_distachyon \
	--exe /path/to/Red --cor 4 --msk_file 

# local run & loading repeats in core Ensembl db
./Red2Ensembl.py Brachypodium_distachyon_v3.0.dna.toplevel.fa Brachypodium_distachyon \
	--exe /path/to/Red --cor 4 --host pl1 --user xyz --pw XYZ \
	--port 123 --db brachypodium_distachyon_core_49_102
```

The repeats called by Red can be optionally annotated by similarity to sequences in an external FASTA file, 
such as the library nrTEplants. The script does not load the resulting annotations in a core db just yet:
```
# produce annotation report
./AnnotRedRepeats.py nrTEplantsJune2020.fna Camelina_sativa --exe /path/to/minimap2 --cor 4

# add annotated repeats to core db
./AnnotRedRepeats.py nrTEplantsJune2020.fna Camelina_sativa --exe /path/to/minimap2 --cor 4 \
    --host pl1 --user xyz --pw XYZ \
    --port 123 --db camelina_sativa_core_50_103_1
```

Annotation reports look like this:
```
# FASTA file with repeat sequences (length>90): Camelina_sativa/annot/Red_repeats.fna

# re-using previously formatted repeat library  Camelina_sativa/annot/repeat_lib.mmi
# running minimap
/path/to/minimap2 -K100M --score-N 0  -x map-ont  -t 4 Camelina_sativa/annot/repeat_lib.mmi Camelina_sativa/annot/Red_repeats.fna
# mapped repeats:  Camelina_sativa/annot/repeat_mappings.sort.paf
# Genome length: 641355730 Repeated content: 230739085 36.0% Annotated: 42633901 6.6%

Crypton 580
DIRS    222132
DNA     456403
DNA/CACTA       312
DNA/En-Spm      729412
DNA/HAT 119922
DNA/Harbinger   56149
DNA/MuDR        598270
DNA/Pogo        243186
DNA/hAT 21260
Helitron        290418
Helitron/Helitron       1673
Helitron|TRIM   1432
LARD    239264
LINE    819243
LINE/L1 80920
LTR     4254019
LTR/Copia       5839801
LTR/Gypsy       18494851
MITE    312641
Maverick        1258
MobileElement   619
Mutator 122
Other   998
Other/Centromeric       194
Other/Simple    1503
RC/Helitron     312013
Retroelement    995
SINE    65415
SINE|LARD       8392
SINE|TRIM       282
TIR     5128793
TIR/Mutator     35807
TIR/PIF-Harbinger       379
TIR/hAT 8030
TRIM    1531956
Unassigned      75795
Unclassified    2676694
nonLTR  101
rRNA    2667
```

## Runtime and RAM requirements

These data were measured on a CentOS7.9 computer using 4 cores of a Xeon E5-2620 v4 (2.10GHz) CPU.

![](../files/runtime_ram.png)

