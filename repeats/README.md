
These scripts can be used to:
+ i) mask repeated sequences in plant genomes with the 
[Repeat detector](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0654-5) (Red)
+ ii) annotate repeats with https://github.com/lh3/minimap2 and  
	* the curated library of repeats 
[nrTEplants](https://github.com/Ensembl/plant-scripts/releases/download/v0.3/nrTEplantsJune2020.fna.bz2),
described in full detail [here](https://github.com/Ensembl/plant_tools/tree/master/bench/repeat_libs)
	* annotated repeats from species in Ensembl Plants, obtained with *get_repeats_ensembl.sh*

Optionally, repeats and annotated repeats can be loaded into an [Ensembl core db](https://www.ensembl.org/info/docs/api/core/core_schema.html) 
as new analyses with default logic names are 'repeatdetector' and 'repeatdetector_annotated'. 

## Dependencies

The following dependencies can be installed in the parent folder with:

    make install install_repeats

There are two required binaries, for which version 10 of the GNU C++ compiler (actually any g++ >= 8 should work, please edit the [Makefile](../Makefile) accordingly):

* A clone of Red from https://github.com/EnsemblGenomes/Red (the original repo is [here](https://github.com/BioinformaticsToolsmith/Red))
* A clone of minimap2 from https://github.com/lh3/minimap2

Plus:

* A copy of the [nrTEplants library](https://github.com/Ensembl/plant-scripts/releases/download/v0.3/nrTEplantsJune2020.fna.bz2)

And three Python3 modules:

* [sqlalchemy](https://pypi.org/project/SQLAlchemy)
* [sqlalchemy_utils](https://pypi.org/project/SQLAlchemy-Utils)
* [pymysql](https://pypi.org/project/PyMySQL)



Note that script [get_repeats_ensembl.sh](./get_repeats_ensembl.sh) has some more dependencies listed at its header.

## Argument lists

If you run 

    ./Red2Ensembl.py -h

you'll get the list of supported arguments and what they're for:

```
./Red2Ensembl.py -h
usage: Red2Ensembl.py [-h] [--exe EXE] [--cor COR] [--msk_file MSK_FILE]
                      [--bed_file BED_FILE] [--host HOST] [--user USER]
                      [--pw PW] [--port PORT] [--db DB]
                      [--logic_name LOGIC_NAME] [--description DESCRIPTION]
                      [--displaylabel DISPLAYLABEL]
                      fasta_file outdir

Script to run RepeatDetector (a fork of Red v2) to mask repeats,
and optionally feed results into an Ensembl core database.

positional arguments:
  fasta_file            path to FASTA file with top-level genomic sequences
  outdir                path to directory to store Red temp results

optional arguments:
  -h, --help            show this help message and exit
  --exe EXE             path to Red executable, default: ./../lib/Red/bin/Red
  --cor COR             number of cores for Red, default: 1
  --msk_file MSK_FILE   name of output FASTA file with soft-masked sequences
  --bed_file BED_FILE   name of output BED file with repeated ranges, uses
                        original sequence names
  --host HOST           name of the database host, required to store repeats
                        in Ensembl core
  --user USER           host user, required to store repeats in Ensembl core
  --pw PW               host password, required to store repeats in Ensembl
                        core
  --port PORT           host port, required to store repeats in Ensembl core
  --db DB               name of the core database, required to store repeats
                        in Ensembl core
  --logic_name LOGIC_NAME
                        logic name of Ensembl analysis, default:
                        repeatdetector
  --description DESCRIPTION
                        quoted string with Ensembl analysis description,
                        default: Repeats detected using <a href="https://bmcbi
                        oinformatics.biomedcentral.com/articles/10.1186/s12859
                        -015-0654-5">Red (REPeatDetector)</a>
  --displaylabel DISPLAYLABEL
                        string with Ensembl analysis display label, default:
                        Repeats:Red

Citation:
Contreras-Moreira et al (2021) https://doi.org/10.1002/tpg2.20143
Girgis HZ (2015) BMC Bioinformatics 16:227. doi: 10.1186/s12859-015-0654-5
```

Similarly, if you run

    ./AnnotRedRepeats.py -h

you'll get:

```
usage: AnnotRedRepeats.py [-h] [--exe EXE] [--cor COR] [--minlen MINLEN]
                          [--bed_file BED_FILE] [--host HOST] [--user USER]
                          [--pw PW] [--port PORT] [--db DB]
                          [--logic_name LOGIC_NAME]
                          [--description DESCRIPTION]
                          [--displaylabel DISPLAYLABEL]
                          repeat_fasta_file outdir

Script to annotate Red repeats and optionally
feed the new consensus_repeats into an Ensembl core database.

positional arguments:
  repeat_fasta_file     path to FASTA file with repeat sequences in RepBase
                        format
  outdir                path to directory with stored Red results

optional arguments:
  -h, --help            show this help message and exit
  --exe EXE             path to minimap2 executable, default:
                        ./../lib/minimap2/minimap2
  --cor COR             number of cores for minimap2, default: 1
  --minlen MINLEN       min length of repeats to be annotated, default: 90bp
  --bed_file BED_FILE   name of output BED file with annotated repeats
  --host HOST           name of the database host
  --user USER           host user
  --pw PW               host password
  --port PORT           host port
  --db DB               name of the core database
  --logic_name LOGIC_NAME
                        logic name of Ensembl analysis, default:
                        repeatdetector_annotated
  --description DESCRIPTION
                        quoted string with Ensembl analysis description,
                        default: Repeats detected using <a href="https://bmcbi
                        oinformatics.biomedcentral.com/articles/10.1186/s12859
                        -015-0654-5">Red (REPeatDetector)</a> and annotated by
                        alignment to a repeat library.
  --displaylabel DISPLAYLABEL
                        string with Ensembl analysis display label, default:
                        'Repeats:Red (annotated)'

Citation:
Contreras-Moreira et al (2021) https://doi.org/10.1002/tpg2.20143
Girgis HZ (2015) BMC Bioinformatics 16:227. doi: 10.1186/s12859-015-0654-5
Li H (2018) Bioinformatics 34(18):3094–3100. doi: 10.1093/bioinformatics/bty191
```

## Examples

Note that the input FASTA file can be GZIP/BZIP2 compressed.
The script *Red2Ensembl.py* will attempt to estimate the GB RAM needed for the input genome.

### i) Masking

```
## test run, saves results in folder 'test_Atha_chr4' 
./Red2Ensembl.py ../files/Arabidopsis_thaliana.fna.gz test_Atha_chr4 --msk_file Atha.sm.fna --bed_file Atha.bed

# parsing FASTA file
# genome length = 18585056 bp
...


## real example, with several chromosomes, taking 4 CPU cores 
./Red2Ensembl.py Brachypodium_distachyon_v3.0.dna.toplevel.fa Brachypodium_distachyon --cor 4 

## local run & loading repeats in core Ensembl db (will re-use previous Red results)
./Red2Ensembl.py Brachypodium_distachyon_v3.0.dna.toplevel.fa Brachypodium_distachyon \
	--host pl1 --user xyz --pw XYZ --port 123 --db brachypodium_distachyon_core_49_102
```

### ii) Annotating masked repeated sequences

The repeats called by Red can be optionally annotated by similarity to sequences in an external FASTA file, 
such as the library **nrTEplants**. The script does not load the resulting annotations in a core db just yet:
```
## test run, re-uses folder 'test_Atha_chr4'
./AnnotRedRepeats.py ../files/nrTEplantsJune2020.fna test_Atha_chr4 --bed_file test.nrTEplants.bed

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


## Error messages

+ ERROR: cannot run Red -9: This means the Red process was killed by the Operating system, usually for taking too much RAM. You will need more RAM to run this job.

