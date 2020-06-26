
# Programmatic analysis of plant genomes with Ensembl Plants

This folder contains code examples for interrogating Ensembl Plants from your own scripts.

<!-- [![Build Status](https://travis-ci.com/Ensembl/plant_tools.svg?branch=master)](https://travis-ci.com/Ensembl/plant_tools) -->

## Dependencies

Some of the scripts depend on additional software packages, see below to learn how to install them.

### FTP

The examples for bulk downloads from the FTP site require the software [wget](https://www.gnu.org/software/wget/), which is usually installed on most Linux distributions. For macOS it is available on [Homebrew](https://brew.sh). For Windows it ships with [MobaXterm](https://mobaxterm.mobatek.net).

### MySQL

The examples for SQL queries to Ensembl Genomes database servers require the [MySQL](https://www.mysql.com) client.

### Perl

As listed in [cpanfile](./cpanfile), three modules are required for the REST examples: [JSON](https://metacpan.org/pod/JSON), [JSON::XS](https://metacpan.org/pod/JSON::XS) and [HTTP::Tiny](https://metacpan.org/pod/HTTP::Tiny). You can install them as explained in [phylogenomics/README.md](../phylogenomics/README.md). 

For the recipes using the Ensembl Perl API please follow the installation [instructions](http://www.ensembl.org/info/docs/api/api_installation.html). If you use git please follow [these](http://www.ensembl.org/info/docs/api/api_git.html). There is also a debugging [guide](https://m.ensembl.org/info/docs/api/debug_installation_guide.html), which lists some extra dependencies that might not have, such as modules [DBI](https://metacpan.org/pod/DBI) and [DBD::mysql](https://metacpan.org/pod/DBD::mysql).

The software VEP has the following Perl dependencies: [DBI](https://metacpan.org/pod/DBI), [DBD::mysql](https://metacpan.org/pod/DBD::mysql) and [Archive::Zip](https://metacpan.org/pod/Archive::Zip). See full documentation at https://github.com/Ensembl/ensembl-vep

### Python

The REST recipes written in python require library [requests](https://pypi.org/project/requests).

### R

You will need BioConductor package [biomaRt](http://www.bioconductor.org/packages/release/bioc/html/biomaRt.html). See installation instructions [here](https://www.ensembl.org/info/data/biomart/biomart_r_package.html).


## List of recipes

A few scripts for phylogenomic studies are documented in folder [phylogenomics](../phylogenomics)

These are the recipes in this folder, obtained with:

```
grep -P "^## \w\d+" example*

exampleAPI.pl:## A1) Load the Registry object with details of genomes available
exampleAPI.pl:## A2) Check which analyses are available for a species
exampleAPI.pl:## A3) Get soft masked sequences from Arabidopsis thaliana
exampleAPI.pl:## A4) Get BED file with repeats in chr4
exampleAPI.pl:## A5) Find the DEAR3 gene
exampleAPI.pl:## A6) Get its canonical sequence
exampleAPI.pl:## A7) Find all orthologues among rosids
exampleAPI.pl:## A8) Get markers mapped on chr1D of bread wheat

exampleBiomart.R:## B1) Check plant marts and select dataset
exampleBiomart.R:## B2) Check available filters and attributes
exampleBiomart.R:## B3) Download GO terms associated to genes
exampleBiomart.R:## B4) Get Pfam domains annotated in genes
exampleBiomart.R:## B5) Get SNP consequences from a selected variation source

exampleFTP.sh:## F1) download peptide sequences in FASTA format
exampleFTP.sh:## F2) download CDS nucleotide sequences in FASTA format
exampleFTP.sh:## F3) download transcripts (cDNA) in FASTA format
exampleFTP.sh:## F4) download soft-masked genomic sequences
exampleFTP.sh:## F5) Upstream/downstream sequences
exampleFTP.sh:## F6) get mappings to UniProt proteins
exampleFTP.sh:## F7) get indexed, bgzipped VCF file with variants mapped
exampleFTP.sh:## F8) get precomputed VEP cache files
exampleFTP.sh:## F9) download all homologies in a single TSV file, several GBs
exampleFTP.sh:## F10) download UniProt report of Ensembl Plants, 
exampleFTP.sh:## F11) retrieve list of new species in current release
exampleFTP.sh:## F12) get current plant species tree (cladogram)

exampleMySQL.sh:## S1) Check currently supported Ensembl Genomes (EG) core schemas,
exampleMySQL.sh:## S2) Count protein-coding genes of a particular species
exampleMySQL.sh:## S3) Get canonical transcript stable_ids
exampleMySQL.sh:## S4) Get variants significantly associated to phenotypes
exampleMySQL.sh:## S5) get Triticum aestivum homeologous genes across A,B & D subgenomes

exampleREST.p[ly]:## R1) Create a HTTP client and a helper functions 
exampleREST.p[ly]:## R2) Get metadata for all plant species 
exampleREST.p[ly]:## R3) Find features overlapping genomic region
exampleREST.p[ly]:## R4) Fetch phenotypes overlapping genomic region
exampleREST.p[ly]:## R5) Find homologues of selected gene
exampleREST.p[ly]:## R6) Get annotation of orthologous genes/proteins
exampleREST.p[ly]:## R7) Fetch variant consequences for multiple variant ids
exampleREST.p[ly]:## R8) Check consequences of single SNP within CDS sequence

exampleVEP.sh:## V1) Download, install and update VEP
exampleVEP.sh:## V2) Unpack downloaded cache file & check SIFT support 
exampleVEP.sh:## V3) Call effect of variants 
exampleVEP.sh:## V4) Call effect of variants for species not in Ensembl
```
