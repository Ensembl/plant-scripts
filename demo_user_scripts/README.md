
# Programmatic analysis of plant genomes with Ensembl Plants

This folder contains code examples for interrogating Ensembl Plants from your own scripts.

[![Build Status](https://travis-ci.com/Ensembl/plant_tools.svg?branch=master)](https://travis-ci.com/Ensembl/plant_tools)

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



## Phylogenomics

A few scripts for phylogenomic studies are documented in folder [phylogenomics](../phylogenomics)

## List of recipes

``
grep -P "^## \w\d+" example*

```
