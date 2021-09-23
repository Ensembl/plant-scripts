
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

The first step involves comparing all vs all genomes. Note that in the example all the results for a species are added to the same TSV file, as well as the logs:

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
