
These scripts can be used to:
+ i) mask repeated sequences in plant genomes
+ ii) annotate repeats with the 
[Repeat detector](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0654-5) (Red) 
and 
	* the curated library of repeats 
[nrTEplants](https://github.com/Ensembl/plant-scripts/releases/download/v0.3/nrTEplantsJune2020.fna.bz2),
described in full detail [here](https://github.com/Ensembl/plant_tools/tree/master/bench/repeat_libs)
	* annotated repeats from species in Ensembl Plants, obtained with *get_repeats_ensembl.sh*

Optionally, repeats can be loaded into an Ensembl core db as an analysis with logic_name='repeatdetector'. 
Annotated repeats are imported with logic_name='repeatdetector_annotated'.

## Dependencies

The following dependencies can be installed in the parent folder with:

    make install install_repeats

There are two binaries, for which the GNU C++ compiler v8 is needed (g++-8):

* A clone of Red from https://github.com/EnsemblGenomes/Red (the original repo is [here](https://github.com/BioinformaticsToolsmith/Red))
* A clone of minimap2 from https://github.com/lh3/minimap2

Plus:

* A copy of the [nrTEplants library](https://github.com/Ensembl/plant-scripts/releases/download/v0.3/nrTEplantsJune2020.fna.bz2)

And three Python3 modules:

* [sqlalchemy](https://pypi.org/project/SQLAlchemy)
* [sqlalchemy_utils](https://pypi.org/project/SQLAlchemy-Utils)
* [pymysql](https://pypi.org/project/PyMySQL)

Note that script [get_repeats_ensembl.sh](./get_repeats_ensembl.sh) has some more dependencies listed at its header.

## Examples

Note that the input FASTA file can be GZIP/BZIP2 compressed.
The script *Red2Ensembl.py* will attempt to estimate the GB RAM needed for the input genome.

### i) Masking

```
## test run, saves results in folder 'test_Atha_chr4' 
./Red2Ensembl.py ../files/Arabidopsis_thaliana.fna.gz test_Atha_chr4 --msk_file Atha.sm.fna 

# parsing FASTA file
# genome length = 18585056 bp
...


## real example, with several chromosomes, taking 4 CPU cores 
./Red2Ensembl.py Brachypodium_distachyon_v3.0.dna.toplevel.fa Brachypodium_distachyon --cor 4 

## local run & loading repeats in core Ensembl db
./Red2Ensembl.py Brachypodium_distachyon_v3.0.dna.toplevel.fa Brachypodium_distachyon \
	--host pl1 --user xyz --pw XYZ --port 123 --db brachypodium_distachyon_core_49_102
```

### ii) Annotating masked repeated sequences

The repeats called by Red can be optionally annotated by similarity to sequences in an external FASTA file, 
such as the library **nrTEplants**. The script does not load the resulting annotations in a core db just yet:
```
## test run, re-uses folder 'test_Atha_chr4'
./AnnotRedRepeats.py ../files/nrTEplantsJune2020.fna test_Atha_chr4 --bed_file test.bed

## consider only repeats with length >= 200 bp
./AnnotRedRepeats.py ../files/nrTEplantsJune2020.fna Brachypodium_distachyon --cor 4 \
	--minlen 200

## add annotated repeats to Ensembl core db and use a different minimap2 binary
./AnnotRedRepeats.py ../files/nrTEplantsJune2020.fna Brachypodium_distachyon --exe /path/to/minimap2 --cor 4 \
    --host pl1 --user xyz --pw XYZ \
    --port 123 --db brachypodium_distachyon_core_49_102
```

Note that any FASTA file can be used to annotate the repeats. For instance, repeats annotated
in current species in Ensembl can be retrieved and used as well:
```
./get_repeats_ensembl.sh arabidopsis_thaliana

# This will produce file: arabidopsis_thaliana.repeats.nondeg.fasta

# Note this file can be highly redundant; redundancy can be eliminated with linclust,
# see https://github.com/soedinglab/MMseqs2

./AnnotRedRepeats.py arabidopsis_thaliana.repeats.nondeg.fasta test_Atha_chr4 --bed_file test.ensembl.bed
```

## Annotation summary 

If a library such as nrTEplants or any other RepBase-formatted file is used, 
an annotation report like this is produced. These are valid examples of FASTA headers:

    >TEdenovo-B-R2315-Map11:repetDB.Mar2020#TIR @Brassica_rapa [S:]
	>AT1TE94285:TAIR10_TE#DNA/MuDR @Arabidopsis_thaliana [S:]

The repeat classification is then parsed to produce a report like this:

```
# Genome length: 18585056 Repeated content: 6837303 36.8% Annotated: 2748796 14.8%

class	bp
DIRS	1212
DNA	32110
DNA/En-Spm	50044
DNA/HAT	33911
DNA/Harbinger	13879
DNA/Mariner	3935
DNA/MuDR	283157
DNA/Pogo	21954
DNA/Tc1	3467
Helitron	20670
LARD	83725
LINE	2384
LINE/L1	9898
LINE?	1235
LTR	113511
LTR/Copia	88739
LTR/Gypsy	900679
MITE	2502
Other	42920
Other/Simple	1596
RC/Helitron	766803
RathE1_cons	1188
RathE2_cons	245
RathE3_cons	196
SINE	9192
Satellite	132
TIR	79579
TIR/Mutator	364
TRIM	53086
Unclassified	126483
```

## Runtime and RAM requirements

These data were measured on a CentOS7.9 computer using 4 cores of a Xeon E5-2620 v4 (2.10GHz) CPU.

![](../files/runtime_ram.png)

