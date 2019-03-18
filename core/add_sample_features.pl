#!/bin/env perl
use strict;
use warnings;
use Getopt::Std;
use Bio::EnsEMBL::Registry;

# This script adds sample features (location, gene, sample) to 1 or more core databases.
# Sample features are actually added as meta keys and values in table meta.
# 
# This might be useful after running load_GFF_hive.pl 
#
# Adapted from ensembl-genomeloader/modules/Bio/EnsEMBL/GenomeLoader/GenomeLoader.pm
# by Bruno Contreras Moreira EMBL-EBI 2019

my ($species_name,$reg_file,$registry);
my (%opts, @species);

getopts('hR:s:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)){
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-s species_name or file with 1name/line    (required, example: -s solanum_lycopersicum)\n";
  print "-R registry file of production db server   (required, example: -R \$p1panreg)\n\n";
  exit(0);
}

if($opts{'R'} && -e $opts{'R'}){ 
	$reg_file = $opts{'R'}; 

	# try a connection to production db
	$registry = 'Bio::EnsEMBL::Registry';	
	if($registry->load_all($reg_file) == 0){
		die "# EXIT: cannot connect to production db server\n";
	}
}
else{ die "# EXIT : need a valid -R file, such as -R \$p2panreg\n" }

if($opts{'s'}){ $species_name = $opts{'s'} }
else{ die "# EXIT : need a valid -s species_name, such as -d solanum_lycopersicum\n" }

if(-s $species_name){ # file with several
	open(LIST,'<',$species_name) || die "# ERROR: cannot read $species_name\n";
	while(<LIST>){
		my $this_sp_name = (split)[0];
		push(@species, $this_sp_name);
	}
	close(LIST);
} else{ # single species
	push(@species, $species_name);
}

## sample features one species at a time
foreach $species_name (@species){
	
	print "# updating $species_name ...\n";

	my $registry2 = 'Bio::EnsEMBL::Registry';
	$registry2->load_registry_from_db(
        -host => 'ensembldb.ensembl.org',
        -user => 'anonymous',
	);


	# get db adaptors from API
	my $meta_adaptor = $registry2->get_adaptor($species_name, "core", "MetaContainer");
	my $gene_adaptor = $registry2->get_adaptor($species_name, "core", "Gene");

	# get all protein-coding genes, if any
	my @genes = @{ $gene_adaptor->fetch_all_by_biotype( 'protein_coding' ) }; 
	if(scalar(@genes) == 0){
                print "# skip, no annotated genes\n";
                next;
        }
	 
	# choose a random gene, giving preference to named genes
  	my @named_genes = grep { defined $_->external_name() } @genes;
  	if ( scalar(@named_genes) > 0 ) { @genes = @named_genes }
	my $random_gene_index = int( rand(scalar(@genes)) ); #print "$random_gene_index\n";
	my $random_gene = $genes[$random_gene_index];

	# get the corresponding location
	my $sr_name     = $random_gene->seq_region_name();
	my $sr_start    = $random_gene->seq_region_start();
	my $sr_end      = $random_gene->seq_region_end();
	my $location    = "$sr_name:${sr_start}-${sr_end}";

	# get the corresponding transcript
	my $transcript = @{ $random_gene->get_all_Transcripts() }[0];

	printf("# Storing sample location, gene & transcript data: %s %s\n",
		$location, 
		$random_gene->external_name() || $random_gene->stable_id() );

	$meta_adaptor->delete_key( 'sample.location_param' );
	$meta_adaptor->store_key_value( 'sample.location_param', $location );
	$meta_adaptor->delete_key( 'sample.location_text' );
	$meta_adaptor->store_key_value( 'sample.location_text', $location );

	$meta_adaptor->delete_key( 'sample.gene_param' );
	$meta_adaptor->store_key_value( 'sample.gene_param', $random_gene->stable_id() );
	$meta_adaptor->delete_key( 'sample.gene_text' );
	$meta_adaptor->store_key_value( 'sample.gene_text', 
		$random_gene->external_name() || $random_gene->stable_id() );

	$meta_adaptor->delete_key( 'sample.transcript_param' );
	$meta_adaptor->store_key_value( 'sample.transcript_param', $transcript->stable_id() );
	$meta_adaptor->delete_key( 'sample.transcript_text' );
	$meta_adaptor->store_key_value( 'sample.transcript_text',
		$transcript->external_name() || $transcript->stable_id() );
}
