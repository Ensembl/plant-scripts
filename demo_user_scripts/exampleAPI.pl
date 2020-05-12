#!/usr/bin/env perl

# Examples of queries to Ensembl Plants using the native Perl API
# Check the tutorial fore more examples:
# https://m.ensembl.org/info/docs/api/core/core_tutorial.html
#
# Copyright [2017-2020] EMBL-European Bioinformatics Institute


# 1. Install the Ensembl Perl API and updated env as explained in
# http://www.ensembl.org/info/docs/api/api_installation.html
# http://www.ensembl.org/info/docs/api/api_git.html
# https://m.ensembl.org/info/docs/api/debug_installation_guide.html
# Note you might need some dependencies, such as libmysqlclient-dev
# or Perl modules DBI & DBD::mysql


# 2. Load the Registry object with details of genomes available
# from the public Ensembl Genomes servers:
use warnings;
use strict;
use Bio::EnsEMBL::Registry;

Bio::EnsEMBL::Registry->load_registry_from_db(
	-USER => 'anonymous',
	-HOST => 'mysql-eg-publicsql.ebi.ac.uk',
	-PORT => '4157',
	#-VERBOSE => 1 # uncomment to see dbs loaded
);


# 3. Find the DEAR3 gene from Arabidopsis thaliana

# gene of interest and species
my $gene_name = 'DEAR3';
my $species = 'arabidopsis_thaliana';

# get a gene adaptor to work with genes from
# the species
my $gene_adaptor = Bio::EnsEMBL::Registry->
	get_adaptor($species, 'core', 'gene');

# find the gene with the specified name using
# the adaptor
my ($gene_obj) = @{$gene_adaptor->
   fetch_all_by_external_name($gene_name)};

# 4. Find all orthologues among rosids

# compara database to search in
my $division = 'plants';

# get an adaptor to work with genes from compara
my $gene_member_adaptor = Bio::EnsEMBL::Registry->
	get_adaptor($division, 'compara', 'GeneMember');

# find the corresponding gene in compara
my $gene_member = $gene_member_adaptor->
	fetch_by_stable_id($gene_obj->stable_id());

# get an adaptor to work with homologues in compara
my $homology_adaptor = Bio::EnsEMBL::Registry->
	get_adaptor($division, 'compara', 'Homology');

# find all homologues of the gene
my @homologies = @{$homology_adaptor->
	fetch_all_by_Member($gene_member)};

# filter out homologues based on taxonomy and type
@homologies = grep {
	$_->taxonomy_level eq 'rosids' &&
	$_->description =~ m/ortholog/
} @homologies;

foreach my $homology (@homologies) {

	# get the protein from the target
	my $target = $homology->get_all_Members->[1];
	
	my $translation = $target->get_Translation;
	
	printf("%s\t%\t%s\t%s\n",
		$gene_obj->stable_id(), 
		$species, 
		$translation->stable_id(),
		$target->genome_db->name, 
		$translation->stable_id );
}

# 5. get BED coordinates of all repeats in chr1 
my $chrname = 'chr1';
my $slice_adaptor = Bio::EnsEMBL::Registry->
	get_adaptor( $species, 'Core', 'Slice' );
my $slice = $slice_adaptor->
	fetch_by_region( 'chromosome', $chrname );

my @repeats = @{ $slice->get_all_RepeatFeatures() };

foreach my $repeat (@repeats) {
	printf("%s\t%d\t%d\t%s\n",
		$chrname,
		$repeat->start(), 
		$repeat->end(), 
		$repeat->display_id() );
}

# 6. markers

