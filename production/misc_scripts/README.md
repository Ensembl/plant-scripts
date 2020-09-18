#### add_gene_from_mRNA_gff.pl:
add genes to a GFF file containing only mRNA features

#### add_gene_from_mRNA_no_parent_gff.pl:
add genes to a GFF file containing only mRNA features with no explicit parent

#### add_transcript_from_CDS_gff.pl:
add transcripts to a GFF file containing only gene and CDS features (no exons)

#### add_transcript_from_gene_gff.pl:
add transcripts to a GFF file containing gene, exon and CDS features

#### BIGD2ENA.pl
takes a FASTA file with contigs from BIGD (non INSDC archive) and attempts to match them
 to sequences in a ENA sequence report asumming order is conserved. Uses sequence length as a control.

#### check_coding_translations.pl
this check is to mimic the dump process by Production
https://github.com/Ensembl/ensembl-production/blob/fd20cda83bf1808b3792ec227e1b29c00a705181/modules/Bio/EnsEMBL/Production/Pipeline/FASTA/DumpFile.pm

#### check_xref_objects.pl
gets the logic name for a list of xrefs objects with no ontology

#### compare_two_fasta_files.pl
Takes two FASTA nucleotide files and checks whether they contain the same headers & sequences

#### delete_genes.pl
Removes genes (or any user-selected logic names/biotypes) from a core db; adapted from
https://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/Core+database+production

#### extract_gff_description_xref.pl
Reads a GFF file and writes a 4-column TSV file with stable_id , gene name, external db and description

#### extract_RAPDB_description_xref.pl
reads a RAP-DB GFF file and writes a 4-column TSV file with stable_id , gene name, external db and description
 for hierarchical external dbs i) RAP-DB and ii) Oryzabase

#### fix_redundant_GFF_IDs.pl
script to fix GFF files 

#### gene_api_example.pl
example of how to get info about a gene from the API

#### get_assembly_provider.pl
get info about assemblies for a given division

#### get_exon_dna.pl
get exon data for a given gene using the API

#### get_slices.pl
get slices and sequence for a given list of seq_regions

#### get_species_with_analysis.pl
get list of species by logic name

#### get_translation_example.pl
translate gene using API

#### json2perlhash.pl
translate from json to hash

#### local_checkout_shell.sh
checkout local API

#### mask_stats.pl
repeat masking stats

#### merge_gff.pl
This script takes two GFF3 files, one with genes and one with transcripts and merges them

#### whosconnected.pl
adapted from https://ben.lobaugh.net/blog/202527/mysql-who-see-who-is-connected-to-your-mysql-server-with-this-script

