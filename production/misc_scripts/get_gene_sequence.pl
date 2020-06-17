#!/usr/bin/env perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use feature qw /say/;

if(!$ARGV[2]){
	die "# usage: $0 <registry> <species> <gene stable id>\n";
}

my ($regfile,$species,$stable_id) = @ARGV;

# Load the registry
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($regfile);
    
# fetch a gene by its stable identifier
my $gene_adaptor = $registry->get_adaptor($species, "core", "gene");
my $gene = $gene_adaptor->fetch_by_stable_id($stable_id);
my $transcript = $gene->canonical_transcript();

printf("gene\t%s\n\n",$gene->seq());
printf("mrna\t%s\n\n",$transcript->seq());
