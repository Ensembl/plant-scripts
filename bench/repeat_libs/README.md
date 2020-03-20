

# nrTEplants

This Markdown document explains how to produce a library of non-redundant transposable elements (TE) found in plants and annotated in the following libraries, contained in FASTA format in folder [repeats/](./repeats/): 

|library|files downloaded|sequences|publication|
|-------|----------------|---------|-----------|
|[TREP](http://botserv2.uzh.ch/kelldata/trep-db)|http://botserv2.uzh.ch/kelldata/trep-db/downloads/trep-db_nr_Rel-19.fasta.gz|4162||
|[SINEbase](http://sines.eimb.ru)|http://sines.eimb.ru/banks/SINEs.bnk|60|https://www.ncbi.nlm.nih.gov/pubmed/23203982|
|[REdat](https://pgsb.helmholtz-muenchen.de/plant/recat)|ftp://ftpmips.helmholtz-muenchen.de/plants/REdat/mipsREdat_9.3p_ALL.fasta.gz|61730|https://www.ncbi.nlm.nih.gov/pubmed/23203886|
|[RepetDB](http://urgi.versailles.inra.fr/repetdb/begin.do)|Exported all Viridiplantae in i) FASTA and ii) TSV formats|33416|https://www.ncbi.nlm.nih.gov/pubmed/30719103|

## cDNA sequences

In order to gauge the overlap between TE libraries and coding sequences, transcripts from the best annotated plant species in Ensembl Plants were also downloaded with script [ens_sequences.pl](../../compara/ens_sequences.pl). The obtained nucleotide files were also put in folder [repeats/](./repeats/). Check the [README.txt](./repeats/README.txt) file there for details.

## Clustering sequences

All TE sequences and cDNA were clustered with [GET_HOMOLOGUES](https://github.com/eead-csic-compbio/get_homologues). This software runs BLASTN and the MCL algorithm, and computes coverage by combining local alignments. 

```
git clone https://github.com/eead-csic-compbio/get_homologues.git
```

Now, in file *get_homologues-est.pl* modify lines [L36-7](https://github.com/eead-csic-compbio/get_homologues/blob/0dce095527ba059c69fba3aa267162e17374f86d/get_homologues-est.pl#L36):
```
set my $MAXSEQLENGTH = 55000;
set my $MINSEQLENGTH = 90;
```

Sequence clustering can be done with the following commands:
```
# default sequence identity = 95%
# default alignment coverage = 75% of shortest sequence in pair
# OMCL algorithm
# skip redundant sequences within TE/cDNA lib overlap>100bp with no mismatches
perl get_homologues-est.pl -d repeats -m cluster -M -t 0 -i 100 &> log.M

# produce pangenome matrix in folder all_clusters
perl compare_clusters.pl -d repeats_est_homologues/SINEs_isoover100_0taxa_algOMCL_e0_ -o all_clusters -m -n

# intersection output directory: all_clusters
# intersection size = 602746 clusters

# pangenome_file = all_clusters/pangenome_matrix_t0.tab transposed = all_clusters/pangenome_matrix_t0.tr.tab
# pangenome_genes = all_clusters/pangenome_matrix_genes_t0.tab transposed = all_clusters/pangenome_matrix_genes_t0.tr.tab

# compute intersections among TE and cDNA libraries
perl parse_pangenome_matrix.pl -m all_clusters/pangenome_matrix_t0.tab -s -x

# intersection pangenome matrix: all_clusters/pangenome_matrix_t0__intersection.tab
# mean %cluster intersection: 5.70

# heatmap of cluster intersection
./plot_matrix_heatmap.sh -i pangenome_matrix_t0__intersection.tab -t "cDNAs from EG46 plants vs TE libraries" -o pdf

# count how many clusters do not contain TEs
perl parse_pangenome_matrix.pl -m all_clusters/pangenome_matrix_t0.tab -A repeats/cdna.list -B repeats/TE.list -a

# finding genes which are absent in B ...
# file with genes absent in B (532793): all_clusters/pangenome_matrix_t0__pangenes_list.txt
```

## Align TE clusters and annotate Pfam domains

We now concentrate on the subset of clusters containing TE sequences. Note that over two thousand clusters contain TE and cDNA sequences, and are thus called 'mixed clusters':

```
perl annot_TEs.pl all_clusters/pangenome_matrix_genes_t0.tr.tab &>log.annot

perl get_ambiguous_Pfam_domains.pl log.annot > Pfam.tsv

# TEclusters=69953
# mixedclusters=2109
```

The resulting file *Pfam.tsv* was then imported to Google Sheets, and columns 'frac_potgenes', 'TElibs_per_cluster' and 'notes'. Note that 'frac_potgenes' is rounded to two decimals and captures the fraction of clustered sequences that share a Pfam domain and have been annotated as 'PotentialHostGene' (PHG) in library RepetDB. 
The resulting TSV file is [Pfam_notes.tsv](./Pfam_notes.tsv).

## Removing ambiguous TE sequences

Two lists of Pfam domains were curated as controls. 
File [control_pos.list](./control_pos.list) contains a set of Pfam domains contained in coding sequences of bona fide TEs. 
Instead, file [control_neg_NLR.list](./control_neg_NLR.list) contains Pfam domains of NLR genes, curated by Carla Filippi. NLR genes are known to be often masked as a side-effect when masking repeated sequences in genomes (see https://www.nature.com/articles/s41477-018-0264-0).

Both files were used to compute the performance of removing TE sequences clustered with PHGs:

```
#TP = True positive = correctly identified
#FP = False positive = incorrectly identified
#TN = True negative = correctly rejected
#FN = False negative = incorrectly rejected

# check control +, column [7] is frac_potgenes
grep -f control_pos.list Pfam_notes.tsv | wc
22   
grep -f control_pos.list Pfam_notes.tsv | perl -F"\t" -lane 'print if($F[7]>0)' | wc
3     

#TP = 19
#FN = 3

# check control -
grep -f control_neg_NLR.list Pfam_notes.tsv | wc
25
grep -f control_neg_NLR.list Pfam_notes.tsv | perl -F"\t" -lane 'print if($F[7]>0)' | wc
21 

#TN = 21
#FP = 4

#sensitivity = TP / (TP + FN) = 19 / (19 + 3) = 0.86
#specificity = TN / (TN + FP) = 21 / 21 + 4) = 0.84
```

Based on this benchmark we decided to filter out TE sequences with Pfam associated to PHGs:
```
perl -F"\t" -lane 'print if ($F[7]>0)' Pfam_notes.tsv  | wc
332

perl -F"\t" -lane 'if($F[7]>0){ foreach $cl (split(/,/,$F[10])){print $cl}}' Pfam_notes.tsv | sort -u > clusters2remove.list
```

Finally, a non-redundant library of plant TEs was produced as follows:

```
./select_TE_clusters.pl log.annot clusters2remove.list nrTEplantsMar2020 &> log.select 
```

The resulting library contains 69,209 sequences from 4 TE collections and has the following contents:

```
# clusters=69209 sequences=69209
mipsREdat_9.3p_ALL      41848
repetDB.Mar2020 26495
trep-db_nr_Rel-19       819
SINEs.plants    47

LTR/Gypsy       19291
LTR     16445
LTR/Copia       8797
TIR     5372
MobileElement   4022
LINE    3229
Unclassified    2758
DNA     892
Helitron        743
SINE    704
DIRS    685
MITE    662
rRNA    568
DNA/En-Spm      548
Other   431
DNA/Mite        414
TRIM    397
Retroelement    393
LARD    298
TIR/hAT 297
DNA/MuDR        284
DNA/TcMar       242
Other/Simple    187
TIR/PIF-Harbinger       176
nonLTR  173
DNA/Stowaway    169
Satellite       160
TIR/Mutator     149
DNA/Harbinger   119
TIR/Mariner     89
TIR/Harbinger   72
DNA/hAT 60
TIR/CACTA       51
non-LTR(SINE)   43
Helitron/Helitron       32
DNA/Tourist     32
non-LTR(SINE)/I 32
Maverick        29
TIR/EnSpm/CACTA 25
LARD|TRIM       25
SINE|TRIM       10
non-LTR(SINE)/Jokey     10
LTR|TIR 10
DNA/hAT-Ac      9
non-LTR(SINE)/L1        8
non-LTR(SINE)/Pan       8
RC/Helitron     8
DIRS|TIR        7
Helitron|LARD   7
LTR/Echo        5
LTR/Halcyon     5
Other/Centromeric       5
PLE     5
Helitron|TRIM   3
LTR|DIRS        3
TIR/PiggyBac    2
TIR|Maverick    2
SINE|LARD       2
DNA/TcMar-Pogo  1
TIR/P   1
non-LTR(SINE)/R2        1
Crypton 1
non-LTR(SINE)/Chronos   1
```

## Clustering sequences from TE libraries with CD-HIT

In order to put the previous results in context a similar clustering experiment, including only TEs and no cDNAs, was carried out with [CD-HIT](http://weizhongli-lab.org/cd-hit):

```
cat mipsREdat_9.3p_ALL.fasta SINEs.bnk.fna trep-db_complete_Rel-19.fasta repetDB.Mar2020.fna >
all.fna

~/soft/cd-hit-v4.8.1-2019-0228/cd-hit-est -i all.fna -c 0.95 -T 3 -o TE.nr.fna

total seq: 99538
longest and shortest : 54107 and 25
Total letters: 632259937
Sequences have been sorted

# comparing sequences from          0  to         18
...
# comparing sequences from      98665  to      99538
..................---------- new table with      552 representatives

    99538  finished      84944  clusters

Approximated maximum memory consumption: 837M
Total CPU time 17598.12
```
