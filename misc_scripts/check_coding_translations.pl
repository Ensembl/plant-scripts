#!/usr/bin/env perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::SeqIO;

# this check is to mimic the dump process by Production
#https://github.com/Ensembl/ensembl-production/blob/fd20cda83bf1808b3792ec227e1b29c00a705181/modules/Bio/EnsEMBL/Production/Pipeline/FASTA/DumpFile.pm

my ($reg_conf, $species);
my $schema_type = 'core';

if(!$ARGV[1]){ die "# usage: $0 <reg_file> <species_name> [otherfeatures]\n" }
else{
	$reg_conf= $ARGV[0];
	$species = $ARGV[1];
}

if($ARGV[2]){
	$schema_type = $ARGV[1];
}

print "# $0 $reg_conf $species $schema_type\n\n";

#################################

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($reg_conf);

my $ga = $registry->get_adaptor( $species, $schema_type, "gene");

my $genes = $ga->fetch_all_by_biotype('protein_coding'); 

my $n_of_coding_genes = 0;
my ($tr,$seq, $gene);
foreach $gene (@$genes) {
	foreach $tr (@{ $gene->get_all_Transcripts }){
		$seq = $tr->translate();
		if(!defined($seq)) {
			printf("Gene %s [%s] Translation %s has no sequence\n", $gene->stable_id(),$gene->biotype, $tr->stable_id());
		} 
		#else { printf("Gene %s [%s] Translation %s has sequence\n", $gene->stable_id(),$gene->biotype, $tr->stable_id()) }
	}

	$n_of_coding_genes++;
}

