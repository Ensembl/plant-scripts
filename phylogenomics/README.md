
# Plant phylogenomics scripts

These scripts interrogate Ensembl Plants through [REST endpoints](https://rest.ensembl.org) 
and the FTP site to export data that might be useful for phylogenomic and pan-gene set studies.

These scripts were tested at the 
[CABANA workshop: Analysis of crop genomics data ](http://training.ensembl.org/events/2021/2021-03-01-CABANA).

## Documentation and examples

Run any of the scripts with argument -h to get instructions and examples.

## Dependencies

The following dependencies can be installed in the parent folder with:

    make install_REST

The scripts require the following non-core Perl modules:
* [HTTP::Tiny](https://metacpan.org/release/HTTP-Tiny)
* [JSON](https://metacpan.org/release/JSON)
* [DBI](https://metacpan.org/pod/DBI)
* [DBD::mysql](https://metacpan.org/pod/DBD::mysql)

which can be installed with: 
```
# install cpanminus installer, check more options at https://metacpan.org/pod/App::cpanminus
sudo cpan -i App::cpanminus  

# actually install modules
sudo apt-get install -y mysql-client libmysqlclient-dev
cpanm JSON JSON::XS HTTP::Tiny DBI DBD::mysql
```

In addition the scripts import module [PlantCompUtils.pm](./PlantCompUtils.pm), 
which is included in this folder.


### ens_single-copy_core_genes.pl

This script can be used to obtain single-copy core genes present within a clade.
Example calls include:

```
perl ens_single-copy_core_genes.pl -c Brassicaceae -f Brassicaceae
perl ens_single-copy_core_genes.pl -c Brassicaceae -f Brassicaceae -t cdna -o beta_vulgaris
perl ens_single-copy_core_genes.pl -f poaceae -c 4479 -r oryza_sativa -WGA 75
perl ens_single-copy_core_genes.pl -f all -c 33090 -m all -r physcomitrium_patens
```

Note option -f produces FASTA files of aligned peptide sequences, one per cluster. 
Such a task takes usually takes over an hour over the Ensembl REST API.


### ens_syntelogs.pl

This script is related to [ens_single-copy_core_genes.pl](ens_single-copy_core_genes.pl) but explicitely considers only orthogroups with Gene Order Conservation (GOC) score >= 75 by default. The output matrix contains also the genomic coordinates of genes of the reference genome:

```
perl ens_syntelogs.pl -c Brassicaceae -f Brassicaceae

```

A sample output matrix is available in [Brassicaceae.syntelogs.GOC75.tsv](./bench/Brassicaceae.syntelogs.GOC75.tsv). 
A benchmark is described in <https://github.com/Ensembl/plant_tools/tree/master/bench/synthelogs>.

Note option -f produces FASTA files of aligned peptide sequences, one per cluster. 
Such a task takes usually takes over an hour over the Ensembl REST API.

WARNING: not all species are included in the Compara gene-tree analysis. You can exclude them with -i.

### ens_sequences.pl

Produces a FASTA file with the canonical cds/pep sequences of species in a clade in Ensembl Plants:
```
perl ens_syntelogs.pl -c Brassicaceae -f Brassicaceae.fna

```


### ens_pangene_analysis.pl

This was a prototype which was eventually replaced by the scripts at 
[pangenes](https://github.com/Ensembl/plant-scripts/tree/master/pangenes).

