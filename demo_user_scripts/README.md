
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

As listed in [cpanfile](./cpanfile), three modules are required: [JSON](https://metacpan.org/pod/JSON), [JSON::XS](https://metacpan.org/pod/JSON::XS) and [HTTP::Tiny](https://metacpan.org/pod/HTTP::Tiny). You can install them as explained in [phylogenomics/README.md](../phylogenomics/README.md)

### Python3



### R

You will need BioConductor package [biomaRt](http://www.bioconductor.org/packages/release/bioc/html/biomaRt.html). See installation instructions [here](https://www.ensembl.org/info/data/biomart/biomart_r_package.html).


## List of recipes

```
grep -P "^## \w\d+" example*

```

A few scripts for phylogenomic studies are documented in folder [phylogenomics](../phylogenomics)

exampleAPI.pl:## A1) Load the Registry object with details of genomes available
exampleAPI.pl:## A2) Find the DEAR3 gene from Arabidopsis thaliana
exampleAPI.pl:## A3) Find all orthologues among rosids
exampleAPI.pl:## A4) Get BED coordinates of all repeats in chr4 
exampleAPI.pl:## A5) Get markers mapped on chr1D of bread wheat
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
exampleFTP.sh:## F8) download all homologies in a single TSV file, several GBs
exampleFTP.sh:## F9) download UniProt report of Ensembl Plants, 
exampleFTP.sh:## F10) retrieve list of new species in current release
exampleFTP.sh:## F11) get current plant species tree (cladogram)
exampleMySQL.sh:## S1) check currently supported Ensembl Genomes (EG) core schemas,
exampleMySQL.sh:## S2) count protein-coding genes of a particular species
exampleMySQL.sh:## S3) get variants significantly associated to phenotypes
exampleMySQL.sh:## S4) get Triticum aestivum homeologous genes across A,B & D subgenomes
exampleREST.pl:## R1) Create an HTTP client and a helper function for invoking a
exampleREST.pl:## R2) Get metadata for all plant species 
exampleREST.pl:## R3) Find features overlapping genomic region
exampleREST.pl:## R4) Fetch phenotypes overlapping genomic region
exampleREST.pl:## R5) Find homologues of selected gene
exampleREST.pl:## R6) Get annotation of orthologous genes/proteins
exampleREST.pl:## R7) Fetch variant consequences for multiple variant ids
