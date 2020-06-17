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
my $t_adaptor = $registry->get_adaptor($species, "core", "Transcript");
my $transcript = $t_adaptor->fetch_by_stable_id($stable_id);

printf("cDNA\t%s\n\n",$transcript->spliced_seq());

foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
	print  "exon: ", $exon->start(), " ", $exon->end(), "\n";
}

printf("CDS \t%s\n\n",$transcript->translateable_seq());


printf("pep \t%s\n\n",$transcript->translate()->seq());

