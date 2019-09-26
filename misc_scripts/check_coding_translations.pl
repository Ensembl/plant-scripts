#!/usr/bin/env perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::SeqIO;

# this check is to mimic the dump process by Production
#https://github.com/Ensembl/ensembl-production/blob/fd20cda83bf1808b3792ec227e1b29c00a705181/modules/Bio/EnsEMBL/Production/Pipeline/FASTA/DumpFile.pm

my $reg_conf = $ENV{'ENSAPIPATH'}."/Registries/p1pan.reg";
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($reg_conf);

###fetch a gene by its stable identifier
my $ga = $registry->get_adaptor("saccharum_spontaneum", "core", "gene");

my $genes = $ga->fetch_all_by_biotype('protein_coding'); 
#my $genes = $ga->fetch_all_by_logic_name('sugarcaneuil');

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

