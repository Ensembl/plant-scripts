# Benchmark

This section summarizes a benchmark experiment carried out to evaluate Ensembl Compara synthelogs.
These are obtained by i) computing [tree-based](https://www.ensembl.org/info/genome/compara/homology_method.html) orthologous genes 
and ii) filtering those with Gene Order Conservation (GOC) scores >= 75. GOC is computed on the four closest neighbours of a gene.

Our gold standard was a list of 22 genomic blocks building up the Ancestral Crucifer Karyotype, reported in https://www.ncbi.nlm.nih.gov/pubmed/26945766. Blocks are defined by *A. thaliana* intervals and the BAC clones that contain them. Those BAC clones were used to design chromosome painting probes used for comparative cytogenetics and define intervals that can be shorter than the gene-based intervals. Blocks are listed in file [Athaliana_blocks.uc.tsv](Athaliana_blocks.uc.tsv). 

The following commands were performed to complete the benchmark:

```
# retrieve synthelogs from Ensembl (GOC>=75), leave out Brassica napus for being allopolyploid
perl ../ens_synthelogs.pl -c Brassicaceae -r arabidopsis_thaliana -i brassica_napus \
	-a > Brassicaceae.synthelogs.tsv

# compare genes within pre-defined synthenic blocks with Ensembl synthelogs
perl _bench_blocks.pl Athaliana_blocks.uc.tsv ../downloads/Arabidopsis_thaliana.TAIR10.45.gtf.gz \
	Brassicaceae.synthelogs.tsv > Athaliana_blocks.report.tsv
```

The final report is in file [Athaliana_blocks.report.tsv](Athaliana_blocks.report.tsv).
