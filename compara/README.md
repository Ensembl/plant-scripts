
# Compara

These are scripts that interrogate Ensembl Genomes Compara through [REST endpoints](https://rest.ensembl.org).

The scripts in folder API/ use directly the [Perl API](https://www.ensembl.org/info/docs/api/index.html).

## Documentation and examples

Please run any of these scripts with argument -h to get instructions and examples.

For example, single-copy_core_genes.pl can be used to obtain single-copy core genes present within a clade. 
If you run ./single-copy_core_genes.pl you'll get:

```
usage: ./single-copy_core_genes.pl [options]

-c NCBI Taxonomy clade of interest      (required, example: -c Brassicaceae or -c 3700)
-l list supported species_names         (optional, example: -l)
-d Ensembl division                     (optional, default: -d Plants)
-r reference species_name               (optional, default: -r arabidopsis_thaliana)
-o outgroup species_name                (optional, example: -o brachypodium_distachyon)
-m multi-copy species_name(s)           (optional, example: -m brassica_napus -m ... -m all)
-i ignore species_name(s)               (optional, example: -i selaginella_moellendorffii -i ...)
-f folder to output FASTA files         (optional, example: -f myfolder)
-t sequence type [protein|cdna]         (optional, requires -f, default: -t protein)
-L allow low-confidence orthologues     (optional, by default these are skipped)
-v verbose                              (optional, example: -v

The following options are only available for some clades:

-G min Gene Order Conservation [0:100]  (optional, example: -G 75)
   see modules/Bio/EnsEMBL/Compara/PipeConfig/EBI/Plants/ProteinTrees_conf.pm
   at https://github.com/Ensembl/ensembl-compara

-W min Whole Genome Align score [0:100] (optional, example: -W 75)
   see ensembl-compara/scripts/pipeline/compara_plants.xml
   at https://github.com/Ensembl/ensembl-compara

Read about GOC and WGA at:
https://www.ensembl.org/info/genome/compara/Ortholog_qc_manual.html

Example calls:

 perl ./single-copy_core_genes.pl -c Brassicaceae -f Brassicaceae
 perl ./single-copy_core_genes.pl -c Brassicaceae -f Brassicaceae -t cdna -o theobroma_cacao
 perl ./single-copy_core_genes.pl -f poaceae -c 4479 -r oryza_sativa -WGA 75
 perl ./single-copy_core_genes.pl -f all -c 33090 -m all -r physcomitrella_patens

```

