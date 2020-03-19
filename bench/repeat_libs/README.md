

# nrTEplants

This Markdown document explains how to produce a library of non-redundant transposable elements (TE) found in plants and annotated in the following libraries, contained in FASTA format in folder [repeats/](./repeats/): 

|library|URL|files downloaded|publication|
|-------|---|----------------|-----------|
|TREP|http://botserv2.uzh.ch/kelldata/trep-db|http://botserv2.uzh.ch/kelldata/trep-db/downloads/trep-db_nr_Rel-19.fasta.gz||
|SINEbase|http://sines.eimb.ru|http://sines.eimb.ru/banks/SINEs.bnk|https://www.ncbi.nlm.nih.gov/pubmed/23203982|
|REdat|https://pgsb.helmholtz-muenchen.de/plant/recat|ftp://ftpmips.helmholtz-muenchen.de/plants/REdat/mipsREdat_9.3p_ALL.fasta.gz|https://www.ncbi.nlm.nih.gov/pubmed/23203886|
|RepetDB|http://urgi.versailles.inra.fr/repetdb/begin.do|Exported all Viridiplantae in FASTA|https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6350395/|

## cDNA sequences

In order to gauge the overlap between TE libraries and coding sequences, transcripts from the best annotated plant species in Ensembl Plants were also downloaded with script [ens_sequences.pl](../../compara/ens_sequences.pl). The obtained nucleotide files were also put in folder [repeats/](./repeats/). Check the [README.txt](./repeats/README.txt) file there for details.

## Clustering sequences

All TE sequences and cDNA were clustered with software [GET_HOMOLOGUES](https://github.com/eead-csic-compbio/get_homologues), which can compute alignment coverage of split alignments.

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
perl get_homologues-est.pl -d repeats -m cluster -M -t 0 -i 100 &> log.M

perl compare_clusters.pl -d repeats_est_homologues/SINEs_isoover100_0taxa_algOMCL_e0_ -o all_clusters -m -n

# number of input cluster directories = 1

# parsing clusters in repeats_est_homologues/SINEs_isoover100_0taxa_algOMCL_e0_ ...
# cluster_list in place, will parse it (repeats_est_homologues/SINEs_isoover100_0taxa_algOMCL_e0_.cluster_list)
# number of clusters = 602746

# intersection output directory: all_clusters
# intersection size = 602746 clusters

# intersection list = all_clusters/intersection_t0.cluster_list

# pangenome_file = all_clusters/pangenome_matrix_t0.tab transposed = all_clusters/pangenome_matrix_t0.tr.tab
# pangenome_genes = all_clusters/pangenome_matrix_genes_t0.tab transposed = all_clusters/pangenome_matrix_genes_t0.tr.tab
# pangenome_phylip file = all_clusters/pangenome_matrix_t0.phylip
# pangenome_FASTA file = all_clusters/pangenome_matrix_t0.fasta
# pangenome CSV file (Scoary) = all_clusters/pangenome_matrix_t0.tr.csv

perl parse_pangenome_matrix.pl -m all_clusters/pangenome_matrix_t0.tab -s -x

# intersection pangenome matrix: all_clusters/pangenome_matrix_t0__intersection.tab
# mean %cluster intersection: 5.70

./plot_matrix_heatmap.sh -i pangenome_matrix_t0__intersection.tab -t "cDNAs from EG46 plants vs TE libraries" -o pdf

perl parse_pangenome_matrix.pl -m all_clusters/pangenome_matrix_t0.tab -A repeats/cdna.list -B repeats/TE.list -a

# matrix contains 602746 clusters and 18 taxa

# taxa included in group B = 4

# finding genes which are absent in B ...
# file with genes absent in B (532793): all_clusters/pangenome_matrix_t0__pangenes_list.txt
```

## Align and get TE sequences

perl annot_TEs.pl all_clusters/pangenome_matrix_genes_t0.tr.tab &>log.annot

# can be used to get Pfam table
#perl -lne 'if(/# Pfam domains: $/){ $tot{'noPfam'}++ } while(/(PF[\d\.]+)/g){ $tot{$1}++; } END{ foreach my $pf (keys(%tot)){ print "$pf\t$tot{$pf}"} }' log.annot | sort -rn -k2,2 

# Pfam domains found in clusters with both TE and transcripts
perl get_ambiguous_Pfam_domains.pl  log.annot > Pfam.tsv
# TEclusters=69953
# mixedclusters=2109

# import to https://docs.google.com/spreadsheets/d/1HzN7Hpr0ul9TxSwAFEaYEoR26cPcnIDogVyfzc_qOfo/edit#gid=1353841590,
# add columns frac_potgenes (round,2), TElibs_perl_cluster & notes, then download Pfam_notes.tsv

True positive = correctly identified
False positive = incorrectly identified
True negative = correctly rejected
False negative = incorrectly rejected

# check control + , column [7] is frac_potgenes
grep -f control_pos.list Pfam_notes.tsv | wc
22   
grep -f control_pos.list Pfam_notes.tsv | perl -F"\t" -lane 'print if($F[7]>0)' |wc
3     

TP = 19
FN = 3

# check control -
grep -f control_neg_NLR.list Pfam_notes.tsv | wc
25
grep -f control_neg_NLR.list Pfam_notes.tsv | perl -F"\t" -lane 'print if($F[7]>0)' | wc
21 

TN = 21
FP = 4

sensitivity = TP / (TP + FN) = 19 / (19 + 3) = 0.86
specificity = TN / (TN + FP) = 21 / 21 + 4) = 0.84

perl -F"\t" -lane 'print if ($F[7]>0)' Pfam_notes.tsv  | wc
332

# clusters to remove
perl -F"\t" -lane 'if($F[7]>0){ foreach $cl (split(/,/,$F[10])){print $cl}}' Pfam_notes.tsv | sort -u > clusters2remove.list

./select_TE_clusters.pl log.annot clusters2remove.list nrTEplantsMar2020 &> log.select 


