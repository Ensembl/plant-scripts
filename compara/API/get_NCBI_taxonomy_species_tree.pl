#!/usr/bin/env perl

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::MetaData::DBSQL::GenomeInfoAdaptor;
use Bio::EnsEMBL::Compara::Utils::SpeciesTree;

# This script creates a NCBI Taxonomy species tree for a (plant) clade of 
# by querying the NCBITaxon facilities of the Compara database with production names
# This is based on 
# https://github.com/Ensembl/ensembl-compara/blob/release/97/scripts/examples/species_buildSpeciesTree.pl
# NOTE1: some production names are not valid such as oryza_indica or panicum_hallii_fil2
# NOTE2: requires the previous release of ensembl-compara; it use relase/97 to query Ensembl Plants 98
# Bruno Contreras Moreira 2019

my $NCBICLADE = 33090; # NCBI Taxonomy id, Viridiplantae
my $VERBOSE = 0;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
  -host=>'ensembldb.ensembl.org',
  -user=>'anonymous', 
);

my (@list_of_species);

## 1) check species in clade ####################################################################

# get a metadata adaptor
my $e_gdba = Bio::EnsEMBL::MetaData::DBSQL::GenomeInfoAdaptor->build_ensembl_genomes_adaptor();

# find and iterate over all genomes from Ensembl Plants
for my $genome (@{$e_gdba->fetch_all_by_taxonomy_branch($NCBICLADE)}) {
	push(@list_of_species, $genome->name());
	print $genome->name()."\n" if($VERBOSE);
}

printf("# number species in NCBICLADE %d : %d\n\n", $NCBICLADE, scalar(@list_of_species));

## 2) build NCBI Taxonomy species tree for this clade

my $taxonDBA = $reg->get_adaptor("Multi", "compara", "NCBITaxon");
my @taxon_ids = ();
foreach my $species_name (@list_of_species) {
  my $taxon = $taxonDBA->fetch_node_by_name($species_name);
  unless ($taxon) {
      warn "Cannot find '$species_name' in the NCBI taxonomy tables\n";
      next;
  }
  push @taxon_ids, $taxon->dbID;
}

my $root = Bio::EnsEMBL::Compara::Utils::SpeciesTree->create_species_tree(
    -COMPARA_DBA    => $reg->get_DBAdaptor("Multi", "compara"),
    -SPECIES_SET    => undef,
    -NO_PREVIOUS    => 1,
    -RETURN_NCBI_TREE       => 1,
    -EXTRATAXON_SEQUENCED   => \@taxon_ids,
);

$root->print_tree(5);

