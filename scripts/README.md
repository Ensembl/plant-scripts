The script merge_genes.py merges gene from Nipponbare 3 annotation sources (Gramene, MSU and RAPDB) based on GFF (obtained after merging GFF from 3 annotations namely oryza_sativa_RAPDB.gff, oryza_sativa_MSU.gff and oryza_sativa_gramene.gff)

Input: sorted_all_gene.gff, this file contain genes from all three annotation sources and sorted based on chromosome number and position.
Output: merged_genes.tsv, contains merged genes in a tsv format.
 
Script options can be explored with help option.

`python merge_genes.py --help`

For example:
`python merge_genes.py --input_gff_file sorted_all_gene.gff --output_dir /homes/user/merged_genes.tsv`

After running merge_genes.py, the information of mRNA, CDS, exon etc. is added for each gene by script add_info.py.

Input: merged_genes.tsv, oryza_sativa_RAPDB.gff, oryza_sativa_MSU.gff and oryza_sativa_gramene.gff
Output: merged_genes.gff

Usage: `python add_info.py`
