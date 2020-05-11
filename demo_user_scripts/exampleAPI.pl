# Example adapted from https://europepmc.org/article/MED/27987162

# 1. Install the Ensembl Perl API ( see Note 28 )

# full, instruction at https://ensembl.org/info/docs/api/api_git.html
# summary:
#git clone -b release-1-6-924 --depth 1 https://github.com/bioperl/bioperl-live.git
#git clone https://github.com/Ensembl/ensembl.git
#cd ensembl
#git checkout release/100
#cd ..

PERL5LIB=${PERL5LIB}:${HOME}/src/bioperl-1.6.924
PERL5LIB=${PERL5LIB}:${HOME}/src/ensembl/modules
export PERL5LIB

# 2. Load the Registry object with details of genomes available
# from the public Ensembl Genomes servers:
use warnings;
use strict;
use Bio::EnsEMBL::Registry;

Bio::EnsEMBL::Registry->load_registry_from_db(
	-USER => 'anonymous',
	-HOST => 'mysql-eg-publicsql.ebi.ac.uk',
	-PORT => '4157',
);

# 3. Find the DEAR3 gene from A. thaliana :

# gene to look for
my $gene_name = 'DEAR3';

# species to look for
my $species = 'arabidopsis_thaliana';

# get a gene adaptor to work with genes from
# the species
my $gene_adaptor = Bio::EnsEMBL::Registry->
	get_adaptor($species, 'core', 'gene');

# find the gene with the specifi ed name using
# the adaptor
my ($gene_obj) = @{$gene_adaptor->
	fetch_all_by_external_name($gene_name)};

# 4. Find all orthologues from Tracheophytes in the plant Compara:

# compara database to search in
my $division = 'plants';

# get an adaptor to work with genes from compara
my $gene_member_adaptor = Bio::EnsEMBL::Registry->
	get_adaptor($division, 'compara', 'GeneMember');

# find the corresponding gene in compara
my $gene_member = $gene_member_adaptor->
	fetch_by_source_stable_id('ENSEMBLGENE',
		$gene_obj->stable_id);

# get an adaptor to work with homologues in compara
my $homology_adaptor = Bio::EnsEMBL::Registry->
	get_adaptor($division, 'compara', 'Homology');

# find all homologues of the gene
my @homologies = @{$homology_adaptor->
	fetch_all_by_Member($gene_member)};

# filter out homologues based on taxonomy and type
@homologies = grep {
	$_->taxonomy_level eq 'Tracheophyta' &&
	$_->description =~ m/ortholog/
} @homologies;

# 5. Find each orthologous protein:
foreach my $homology (@homologies) {

	# get the protein from the target
	my $target = $homology->get_all_Members->[1];
	
	my $translation = $target->get_Translation;
	
	print $target->genome_db->name, ' orthologue ',
		$translation->stable_id, "\n";
	
# 6. For each translation, print information about GO annotation:

	# find all the GO terms for this translation
	foreach my $go_term (@{$translation->get_all_DBEntries('GO')}) {
	
	# print some information about each GO annotation
	print $go_term->primary_id, ' ', $go_term->description;
	
	# print the evidence for the GO annotation
	print ' Evidence: ', (join ', ', map {$_->[0]}
		@{$go_term->get_all_linkage_info}), "\n";
	}
}

