
# Pangene scripts

The scripts in this folder can be used to call collinear genes in genomes of your choice and produce clusters (files in FASTA format) of the corresponding encoded cDNA, CDS and peptide sequences. 

## Dependencies

* Perl
* https://github.com/lh3/minimap2 or alternatively https://github.com/ekg/wfmash
* https://bedtools.readthedocs.io/en/latest/
* https://github.com/gpertea/gffread

##  Example with 4 rice genomes

I'll demonstrate how to use them with a test dataset which includes the following assemblies & annotations for 4 rices. You should obtain similar files for your project:

    Oryza_sativa.IRGSP-1.0.dna.toplevel.fa
    all.chrs.con
    Oryza_indica.ASM465v1.dna.toplevel.fa
    Oryza_nivara.Oryza_nivara_v1.0.dna.toplevel.fa

	Oryza_sativa.IRGSP-1.0.51.gff3
    all.chr.gff3                   
	Oryza_indica.ASM465v1.51.gff3
	Oryza_nivara.Oryza_nivara_v1.0.51.gff3

These files correspond to the following custom named genomes, which can be binary or ternary:

    oryza_sativa
    oryza_sativa_MSU
    oryza_indica
    oriza_nivara

Note that the first are actually the same assembly with different annotated gene models.

### 1) Compute pairwise collinear genes with minimap2

The first step involves comparing all vs all genomes. Note that in the example all the results for a species are added to the same TSV file, as well as the logs. Also note that the call below assume a binary called minimap2 can be found in $PATH, but you can use -M to point to another location. In our benchmark we obtained best results with v2.17-r941:

```
# Osativa
perl get_collinear_genes.pl \
  -sp1 oryza_sativa -fa1 Oryza_sativa.IRGSP-1.0.dna.toplevel.fa -gf1 Oryza_sativa.IRGSP-1.0.51.gff3 \
  -sp2 oryza_nivara -fa2 Oryza_nivara.Oryza_nivara_v1.0.dna.toplevel.fa -gf2 Oryza_nivara.Oryza_nivara_v1.0.51.gff3 \
  -out Minimap2.homologies.oryza_sativa.overlap0.5.tsv &> log.Osativa

perl get_collinear_genes.pl \
  -sp1 oryza_sativa -fa1 Oryza_sativa.IRGSP-1.0.dna.toplevel.fa -gf1 Oryza_sativa.IRGSP-1.0.51.gff3 \
  -sp2 oryza_sativa_MSU -fa2 all.chrs.con -gf2 all.chr.gff3 -r -add \
  -out Minimap2.homologies.oryza_sativa.overlap0.5.tsv &>> log.Osativa

perl get_collinear_genes.pl \
  -sp1 oryza_sativa -fa1 Oryza_sativa.IRGSP-1.0.dna.toplevel.fa -gf1 Oryza_sativa.IRGSP-1.0.51.gff3 \
  -sp2 oryza_indica -fa2 Oryza_indica.ASM465v1.dna.toplevel.fa -gf2 Oryza_indica.ASM465v1.51.gff3 -r -add \
  -out Minimap2.homologies.oryza_sativa.overlap0.5.tsv &>> log.Osativa

# Onivara
perl get_collinear_genes.pl \
  -sp1 oryza_nivara -fa1 Oryza_nivara.Oryza_nivara_v1.0.dna.toplevel.fa -gf1 Oryza_nivara.Oryza_nivara_v1.0.51.gff3 \
  -sp2 oryza_sativa_MSU -fa2 all.chrs.con -gf2 all.chr.gff3 \
  -out Minimap2.homologies.oryza_nivara.overlap0.5.tsv &> log.Onivara

perl get_collinear_genes.pl \
  -sp1 oryza_nivara -fa1 Oryza_nivara.Oryza_nivara_v1.0.dna.toplevel.fa -gf1 Oryza_nivara.Oryza_nivara_v1.0.51.gff3 \
  -sp2 oryza_indica -fa2 Oryza_indica.ASM465v1.dna.toplevel.fa.gz -gf2 Oryza_indica.ASM465v1.51.gff3.gz -r -add \
  -out Minimap2.homologies.oryza_nivara.overlap0.5.tsv &>> log.Onivara

# Oindica
perl get_collinear_genes.pl \
  -sp1 oryza_indica -fa1 Oryza_indica.ASM465v1.dna.toplevel.fa.gz -gf1 Oryza_indica.ASM465v1.51.gff3.gz \
  -sp2 oryza_sativa_MSU -fa2 all.chrs.con -gf2 all.chr.gff3 \
  -out Minimap2.homologies.oryza_indica.overlap0.5.tsv &> log.Oindica
```

The TSV files produced contain data like this excerpt:

    gene_stable_id  protein_stable_id       species overlap homology_type   homology_gene_stable_id homology_protein_stable_id      homology_species        overlap dn      ds      goc_score       wga_coverage    is_high_confidence      coordinates
    gene-mag-r      gene-mat-r      oryza_sativa    525     ortholog_collinear      ONIVA02G03490   ONIVA02G03490   oryza_nivara    525     NULL    NULL    NULL    100.00  1       Mt:315811-317848;Mt:315625-316336
    gene-orf165     gene-orf165     oryza_sativa    498     ortholog_collinear      ONIVA02G03510   ONIVA02G03510   oryza_nivara    498     NULL    NULL    NULL    100.00  1       Mt:337076-337574;Mt:323628-338167
    gene-rps1       gene-rps1       oryza_sativa    519     ortholog_collinear      ONIVA02G03500   ONIVA02G03500   oryza_nivara    519     NULL    NULL    NULL    100.00  1       Mt:321367-321886;Mt:318304-322745
    Os01g0100100    Os01g0100100    oryza_sativa    7833    ortholog_collinear      ONIVA01G00100   ONIVA01G00100   oryza_nivara    7833    NULL    NULL    NULL    100.00  1       1:2982-10815;1:2929-12267

### 2) Extract cDNA, CDS and pep sequences

This script cuts all available sequences found for each gene. Note that not all cDNA sequences encode CDS/peptides. In the example the path to gffread is explicitely passed as an argument:

```
perl cut_sequences.pl \
  -sp oryza_sativa  -fa Oryza_sativa.IRGSP-1.0.dna.toplevel.fa -gf Oryza_sativa.IRGSP-1.0.51.gff3 \
  -nr -p ~/soft/gffread-0.12.7.Linux_x86_64/gffread

perl cut_sequences.pl \
  -sp oryza_sativa_MSU -fa all.chrs.con -gf all.chr.gff3 \
  -nr -p ~/soft/gffread-0.12.7.Linux_x86_64/gffread

perl cut_sequences.pl \
  -sp oryza_indica -fa Oryza_indica.ASM465v1.dna.toplevel.fa -gf Oryza_indica.ASM465v1.51.gff3 \
  -nr -p ~/soft/gffread-0.12.7.Linux_x86_64/gffread

perl cut_sequences.pl \
  -sp oryza_nivara -fa Oryza_nivara.Oryza_nivara_v1.0.dna.toplevel.fa -gf Oryza_nivara.Oryza_nivara_v1.0.51.gff3 \
  -nr -p ~/soft/gffread-0.12.7.Linux_x86_64/gffread
```

This step will produce three FASTA files per genome.

### 3) Define pangenes, clusters of collinear genes

Finally, this step uses the files computed in the previous steps to compile clusters of collinear genes. In this case the results are saved in folder results/ and oryza_sativa was selected as reference genome: 

```
perl pangene_analysis.pl -T Minimap2.homologies.oryza_sativa.overlap0.5.tsv \
  -T Minimap2.homologies.oryza_indica.overlap0.5.tsv \
  -T Minimap2.homologies.oryza_nivara.overlap0.5.tsv \
  -f results -r oryza_sativa
```

In this example, the clusters are stored in folder oryzasativa/ and a text file describing the clusters is also produced, named oryzasativa.cluster_list, which looks like this:

    cluster gene-mag-r size=1 taxa=2 cdnafile: gene-mag-r.cdna.fna cdsfile: gene-mag-r.cds.fna pepfile: gene-mag-r.cds.faa
    : oryza_nivara
    cluster gene-orf165 size=3 taxa=3 cdnafile: gene-orf165.cdna.fna cdsfile: gene-orf165.cds.fna pepfile: gene-orf165.cds.faa
    : oryza_nivara
    : oryza_sativa_MSU
    : oryza_sativa
    cluster gene-rps1 size=3 taxa=3 cdnafile: gene-rps1.cdna.fna cdsfile: gene-rps1.cds.fna pepfile: gene-rps1.cds.faa
    : oryza_nivara
    : oryza_indica
    : oryza_sativa

Each clusters is a FASTA file, which looks like this:

```
grep ">" oryzasativa/Os01g0100100.cdna.fna
>Os01t0100100-01 Os01g0100100 1:2983-10815 [oryza_sativa]
>Os01t0100200-01 Os01g0100200 1:11218-12435 [oryza_sativa]
>ONIVA01G00100.2 ONIVA01G00100 1:104921-115645 [oryza_nivara]
>ONIVA01G00100.3 ONIVA01G00100 1:104921-115645 [oryza_nivara]
>ONIVA01G00100.4 ONIVA01G00100 1:104921-115645 [oryza_nivara]
>ONIVA01G00100.1 ONIVA01G00100 1:104921-116326 [oryza_nivara]
>LOC_Os01g01010.1 LOC_Os01g01010 Chr1:2903-10817 [oryza_sativa_MSU]
>LOC_Os01g01010.2 LOC_Os01g01010 Chr1:2984-10562 [oryza_sativa_MSU]
>LOC_Os01g01019.1 LOC_Os01g01019 Chr1:11218-12435 [oryza_sativa_MSU]
>BGIOSGA002569-TA BGIOSGA002569 1:30220-36442 [oryza_indica]
>BGIOSGA002570-TA BGIOSGA002570 1:38569-39088 [oryza_indica]
```

Multiple alignments can be computed for each FASTA file to determine which is the most conserved gene structure.

The script also produces POCP and pangenomes matrices (described [here](../phylogenomics)). Note that currently clusters are not guaranteed to be sorted by chr position.


