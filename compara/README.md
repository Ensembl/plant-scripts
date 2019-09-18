
# Compara

These scripts interrogate Ensembl Genomes Compara through [REST endpoints](https://rest.ensembl.org) and 
the FTP site.

The scripts in folder API/ use directly the [Perl API](https://www.ensembl.org/info/docs/api/index.html).

## Documentation and examples

Run any of the scripts with argument -h to get instructions and examples.

### pangene_analysis.pl

Script pangene_analysis.pl can be used to analyse a clade-specific pan-gene set, with several options:

```
usage: pangene_analysis.pl [options]

-c NCBI Taxonomy clade of interest         (required, example: -c Brassicaceae or -c 3700)
-f output folder                           (required, example: -f myfolder)
-r reference species_name to name clusters (required, example: -r arabidopsis_thaliana)
-l list supported species_names            (optional, example: -l)
-d Ensembl division                        (optional, default: -d Plants)
-o outgroup species_name                   (optional, example: -o brachypodium_distachyon)
-i ignore species_name(s)                  (optional, example: -i selaginella_moellendorffii -i ...)
-L allow low-confidence orthologues        (optional, by default these are skipped)
-v verbose                                 (optional, example: -v

The following options are only available for some clades:

-G min Gene Order Conservation [0:100]  (optional, example: -G 75)
   see modules/Bio/EnsEMBL/Compara/PipeConfig/EBI/Plants/ProteinTrees_conf.pm
   at https://github.com/Ensembl/ensembl-compara

-W min Whole Genome Align score [0:100] (optional, example: -W 75)
   see ensembl-compara/scripts/pipeline/compara_plants.xml
   at https://github.com/Ensembl/ensembl-compara

Read about GOC and WGA at:
https://www.ensembl.org/info/genome/compara/Ortholog_qc_manual.html
```

These examples generate results in folders [Brassicaceae](./Brassicaceae) and [Oryza](./Oryza) and produce the following output:
```
perl pangene_analysis.pl -c Oryza -f Oryza -r oryza_sativa
perl pangene_analysis.pl -r arabidopsis_thaliana -c Brassicaceae -f Brassicaceae

...

# pangene_analysis.pl -d Plants -c Brassicaceae -r arabidopsis_thaliana -o  -f Brassicaceae -t protein -G 0 -W 0 -L 0

# supported species in NCBI taxon Brassicaceae : 6

# total selected species : 6

...

# arabidopsis_thaliana : sequences = 27628 clusters = 26478 (singletons = 1878)
# arabidopsis_lyrata : sequences = 32667 clusters = 29462 (singletons = 3787)
# brassica_rapa : sequences = 41018 clusters = 27880 (singletons = 1141)
# brassica_oleracea : sequences = 59220 clusters = 39457 (singletons = 8077)
# brassica_napus : sequences = 101040 clusters = 46335 (singletons = 11181)
# arabidopsis_halleri : sequences = 32158 clusters = 29172 (singletons = 3981)

# total sequences = 293731

# number_of_clusters = 71515 (core = 16609)

# cluster_list = Brassicaceae/arabidopsisthaliana_Brassicaceae_algEnsemblCompara.cluster_list
# cluster_directory = Brassicaceae/arabidopsisthaliana_Brassicaceae_algEnsemblCompara

# percent_conserved_proteins_file = Brassicaceae/POCP.matrix.tab

# pangenome_file = Brassicaceae/pangenome_matrix.tab tranposed = Brassicaceae/pangenome_matrix.tr.tab
# pangenome_genes = Brassicaceae/pangenome_matrix_genes.tab transposed = Brassicaceae/pangenome_matrix_genes.tr.tab
# pangenome_FASTA_file = Brassicaceae/pangenome_matrix.fasta

# genome composition report (samples=10,seed=12345)
## sample 0 (arabidopsis_thaliana | 0,1,2,3,4,5,)
## sample 1 (brassica_oleracea | 3,4,5,1,2,0,)
## sample 2 (arabidopsis_lyrata | 1,2,5,0,3,4,)
## sample 3 (brassica_oleracea | 3,2,4,5,1,0,)
## sample 4 (arabidopsis_thaliana | 0,2,1,5,3,4,)
## sample 5 (brassica_oleracea | 3,2,0,5,4,1,)
## sample 6 (arabidopsis_lyrata | 1,4,3,5,0,2,)
## sample 7 (arabidopsis_halleri | 5,3,1,0,2,4,)
## sample 8 (brassica_rapa | 2,4,5,3,1,0,)
## sample 9 (arabidopsis_thaliana | 0,2,4,3,5,1,)

# pan-gene (number of clusters) = Brassicaceae/pan_gene.tab
# core-gene (number of clusters) = Brassicaceae/core_gene.tab

# runtime: 72 wallclock secs (62.80 usr  2.80 sys + 27.61 cusr  1.35 csys = 94.56 CPU)
```
And those files can be used to obtain plots such as these with scripts from 
[GET-HOMOLOGUES](https://github.com/eead-csic-compbio/get_homologues): 




### single-copy_core_genes.pl

Instead, single-copy_core_genes.pl can be used to obtain single-copy core genes present within a clade. 
Example calls include:

```
 perl ./single-copy_core_genes.pl -c Brassicaceae -f Brassicaceae
 perl ./single-copy_core_genes.pl -c Brassicaceae -f Brassicaceae -t cdna -o theobroma_cacao
 perl ./single-copy_core_genes.pl -f poaceae -c 4479 -r oryza_sativa -WGA 75
 perl ./single-copy_core_genes.pl -f all -c 33090 -m all -r physcomitrella_patens

```

