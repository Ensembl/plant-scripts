#!/usr/bin/env perl

# Example in Perl of queries to the Ensembl REST endpoints 
# Your usage of the data returned by the REST service is 
# subject to same conditions as laid out on the Ensembl website.
#
# The full set of endpoints is documented 
# at http://rest.ensembl.org
#
# See https://github.com/Ensembl/ensembl-rest/wiki for
# examples in other programming languages

# Copyright [2017-2020] EMBL-European Bioinformatics Institute

use strict;
use warnings;
use JSON;
use HTTP::Tiny;
use Data::Dumper;

## 1. Create an HTTP client and a helper function for invoking a
##  REST endpoint:

# create an HTTP client
my $http = HTTP::Tiny->new;
my $server = 'http://rest.ensembl.org';

# function for invoking endpoint, see other options at
# https://github.com/Ensembl/ensembl-rest/wiki
sub call_endpoint {
	my ($url, $verbose) = @_;
	
	print "Invoking $url\n" if($verbose);

	my $response = $http-> get($url, {headers =>	
		{'Content-type' => 'application/json'} });
	
	return decode_json($response->{content});
}

## 2. Get metadata for all plant species 
my $url =
    join('/', $server, '/info/genomes/division/EnsemblPlants').
	    "?content-type=application/json";

# call url endpoint and get an array back,
# note that actual data structure changes across endpoints
my $metadata = call_endpoint($url);

# parse the data from the response
foreach my $sp (@$metadata){
	printf("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n",
		$sp->{name},
		$sp->{strain} ||'NA', 
		$sp->{assembly_accession},
		$sp->{assembly_level},
		$sp->{has_peptide_compara},
		$sp->{has_variations},
		$sp->{has_genome_alignments},
		$sp->{has_synteny});
}

## 3. Find features overlapping genomic region, full list at
## https://rest.ensembl.org/documentation/info/overlap_region

my $species = 'triticum_aestivum';
my $region = '3D:379400000-379540000';

# genes
$url =
    join('/', $server, 'overlap/region', $species, $region).
	"?feature=gene;content-type=application/json";

print "# check it on the browser: $url\n\n";

my $overlap_data = call_endpoint($url);

foreach my $overlap_feat (@$overlap_data){
	printf("%s\t%s\t%s\n",$overlap_feat->{id},
		$overlap_feat->{start}, $overlap_feat->{end});
}

# now LTR repeats 
$url = join('/', $server, 'overlap/region', $species, $region).
		"?feature=repeat;content-type=application/json";

$overlap_data = call_endpoint($url);

foreach my $overlap_feat (@$overlap_data){

	next if($overlap_feat->{description} !~ 'LTR');

	printf("%s\t%s\t%s\n", $overlap_feat->{description},
		$overlap_feat->{start}, $overlap_feat->{end});
}

# finally get EMS variants
$url = join('/', $server, 'overlap/region', $species, $region).
	"?feature=variation;content-type=application/json";

$overlap_data = call_endpoint($url);

foreach my $overlap_feat (@$overlap_data){
	next if($overlap_feat->{source} ne 'EMS-induced mutation');
	print $overlap_feat->{id}."\n";
}

## 4. Find homologues of A. thaliana DEAR3 gene
$species = 'arabidopsis_thaliana';
my $division = 'plants';
my $gene = 'DEAR3';
my $homoltype = 'ortholog';

# optional define a target clade, such as 4479
# for Poaceae, see https://www.ncbi.nlm.nih.gov/taxonomy
my $target_clade = 4479; 

$url = 
	join('/', $server, 'homology/symbol', $species, $gene).
	"?content-type=application/json&compara=$division";

# restrict to homologues from this clade/taxon
if(defined($target_clade)){
	$url .= "&target_taxon=$target_clade";
}

my $homology_data = call_endpoint($url);

my @homologies = @{$homology_data->{data}[0]{homologies}};

# filter out homologues based on type
if(defined($homoltype)){
	@homologies = grep {
		$_->{type} =~ m/$homoltype/
	} @homologies;
}

## 5. Print some information about the orthologous protein
for my $homolog (@homologies) {

	my $target_species = $homolog->{target}{species};
	my $target_id = $homolog->{target}{protein_id};
	print "$target_species $homolog->{type} $target_id\n";

	# For each translation, print GO annotation
	# using the xrefs/id endpoint:
	$url = join('/', $server, "xrefs/id/$target_id").
		"?content-type=application/json;external_db=GO;all_levels=1";

	my $go_data = call_endpoint($url);
	for my $go (@{$go_data}) {
		print $go->{display_id}, ' ', $go->{description} || 'NA',
		' Evidence: ', join(', ', @{$go->{linkage_types}}),
		"\n";
	}
}
