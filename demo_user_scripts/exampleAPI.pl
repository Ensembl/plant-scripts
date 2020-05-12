#!/usr/bin/env perl

# Examples of queries to Ensembl Plants using the native Perl API
# Check the tutorials for more examples:
# https://m.ensembl.org/info/docs/api/core/core_tutorial.html
# https://www.ensembl.org/info/docs/api/compara/compara_tutorial.html
#
# Copyright [2017-2020] EMBL-European Bioinformatics Institute


# Install the Ensembl Perl API and updated env as explained in
# http://www.ensembl.org/info/docs/api/api_installation.html
# http://www.ensembl.org/info/docs/api/api_git.html
# https://m.ensembl.org/info/docs/api/debug_installation_guide.html
# Note you might need some dependencies, such as libmysqlclient-dev
# or Perl modules DBI & DBD::mysql


## A1) Load the Registry object with details of genomes available
#      from the public Ensembl Genomes server:
use warnings;
use strict;
use Bio::EnsEMBL::Registry;

Bio::EnsEMBL::Registry->load_registry_from_db(
	-USER => 'anonymous',
	-HOST => 'mysql-eg-publicsql.ebi.ac.uk',
	-PORT => '4157',
	#-VERBOSE => 1 # uncomment to see dbs loaded
);


## A2) Find the DEAR3 gene from Arabidopsis thaliana

# gene of interest and species
my $gene_name = 'DEAR3';
my $species = 'arabidopsis_thaliana';
my $division = 'plants';

# get a gene adaptor to work with genes from
# the species
my $gene_adaptor = Bio::EnsEMBL::Registry->
	get_adaptor($species, 'core', 'gene');

# find the gene with the specified name using
# the adaptor
my ($gene_obj) = @{$gene_adaptor->
   fetch_all_by_external_name($gene_name)};

# A3) Find all orthologues among rosids

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
	
	printf("%s\t%s\t%s\t%s\n",
		$gene_obj->stable_id(), 
		$species, 
		$translation->stable_id(),
		$target->genome_db->name() );
}

# A4) Get BED coordinates of all repeats in chr4 

my $chrname = 'chr4';

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

# A5) Get markers mapped on chr1D of bread wheat

# Note: only a few plants have markers
# As of release EG47/100:
# triticum_aestivum, oryza_indica, brassica_rapa

$species = 'triticum_aestivum';
$chrname = '1D';

$slice_adaptor = Bio::EnsEMBL::Registry->
	get_adaptor( $species, 'Core', 'Slice' );

$slice = $slice_adaptor->
	fetch_by_region( 'chromosome', $chrname );

my $marker_features = $slice->get_all_MarkerFeatures();
while ( my $mf = shift @{$marker_features} ) {

	my $marker = $marker_feature->marker(); 

	printf("%s\t%d\t%d\t%s\t%s\t%s\t%d\n",
		$mf->seq_region_name(),
		$mf->start(),      
		$mf->end(), 
		$mf->display_id(),
		$marker->left_primer(),
		$marker->marker()->right_primer(),
		$marker->marker()->max_primer_dist() );
}
