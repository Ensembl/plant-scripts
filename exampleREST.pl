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

## R1) Create a HTTP client and helper functions 

# create a new HTTP client
my $http = HTTP::Tiny->new();
my $server = 'http://rest.ensembl.org';

# function for invoking endpoint, 
# see how to cope with failed request and retries at 
# https://github.com/Ensembl/ensembl-rest/wiki/Example-Perl-Client
sub call_endpoint {
	my ($http, $url, $verbose) = @_;
	
	print "Invoking $url (GET)\n" if($verbose);

	my $response = $http->get($url, {
		headers => 
			{'Content-type' => 'application/json'} });
	
	die "# Failed GET request $url\n" unless $response->{success};

	return decode_json($response->{content});
}

# function to post data to endpoint and 
# wait for response
sub post_endpoint {
	my ($http, $url, $data, $verbose) = @_;

	print "Invoking $url (POST)\n" if($verbose);

	my $response = $http->request('POST', $url, {	
		headers => { 
			'Content-type' => 'application/json',
			'Accept' => 'application/json'
		},
		content => $data
	});
	
	die "# Failed POST request $url\n" unless $response->{success};

	return decode_json($response->{content});
}

## R2) Get metadata for all plant species 

my $url =
    join('/', $server, '/info/genomes/division/EnsemblPlants').
	    "?content-type=application/json";

# call url endpoint and get an array back,
# note that actual data structure changes across endpoints
my $metadata = call_endpoint($http,$url);

# parse the data from the response
foreach my $sp (@$metadata){
	printf("%s\t%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\n",
		$sp->{name},
		$sp->{strain} ||'NA', 
		$sp->{assembly_accession},
		$sp->{base_count}, #golden path genome size
		$sp->{assembly_level},
		$sp->{has_peptide_compara},
		$sp->{has_variations},
		$sp->{has_genome_alignments},
		$sp->{has_synteny});
}

## R3) Find features overlapping genomic region

# full list at
# https://rest.ensembl.org/documentation/info/overlap_region

# Note 1-based inclusive coords are returned

my $species = 'triticum_aestivum';
my $region = '3D:379400000-379540000';

# genes
$url =
    join('/', $server, 'overlap/region', $species, $region).
	"?feature=gene;content-type=application/json";

print "# check it on the browser: $url\n\n";

my $overlap_data = call_endpoint($http,$url);

foreach my $overlap_feat (@$overlap_data){
	printf("%s\t%s\t%s\n",
		$overlap_feat->{id},
		$overlap_feat->{start}, 
		$overlap_feat->{end});
}

# now LTR repeats 
$url = join('/', $server, 'overlap/region', $species, $region).
		"?feature=repeat;content-type=application/json";

$overlap_data = call_endpoint($http,$url);

foreach my $overlap_feat (@$overlap_data){

	next if($overlap_feat->{description} !~ 'LTR');

	printf("%s\t%s\t%s\n", 
		$overlap_feat->{description},
		$overlap_feat->{start}, 
		$overlap_feat->{end});
}

# EMS variants

my $source_variation = 'EMS-induced mutation';
#$source_variation = 'CerealsDB';

$url = join('/', $server, 'overlap/region', $species, $region).
	"?feature=variation;content-type=application/json";

$overlap_data = call_endpoint($http,$url);

foreach my $overlap_feat (@$overlap_data){
	next if($overlap_feat->{source} ne $source_variation);
	printf("%s\t%s\n",
		$overlap_feat->{id},
		$overlap_feat->{source});
}

# protein-coding genes from additional annotation tracks,
# also called otherfeatures dbs
my $dbtype = 'otherfeatures';

$url = join('/', $server, 'overlap/region', $species, $region).
	"?feature=transcript;db_type=$dbtype;content-type=application/json";

$overlap_data = call_endpoint($http,$url);

foreach my $overlap_feat (@$overlap_data){
	next if($overlap_feat->{biotype} ne 'protein_coding');

	# get peptide sequences encoded by these genes
	$url = join('/', $server, 'sequence/id', $overlap_feat->{id}).
		"?db_type=$dbtype;type=protein;species=$species;object_type=transcript;".
		"content-type=application/json";
	
	my $result = call_endpoint($http,$url);
	printf(">%s\n%s\n",$result->{'id'},$result->{'seq'});
}


## R4) Fetch phenotypes overlapping genomic region

$species = 'arabidopsis_thaliana';
$region = '3:19431095-19434450';
my $p_cutoff = 0.0001; 

$url = join('/', $server, 'phenotype/region', $species, $region).
        "?content-type=application/json;feature_type=Variation";

my $pheno_data = call_endpoint($http,$url);

foreach my $feat (@$pheno_data){

	foreach my $assoc (@{ $feat->{phenotype_associations} }) {
		
		next if($assoc->{attributes}{p_value} > $p_cutoff);

		printf("%s\t%s\t%s\n",
		$feat->{id},
		$assoc->{location},
		$assoc->{description});
	}
}


## R5) Find homologues of selected gene

$species = 'triticum_aestivum'; 
my $division = 'plants';
my $gene = 'TraesCS1B02G195200'; 
my $homoltype = 'ortholog'; #paralog

# optionally define a target clade, such as 4479
# for Poaceae, see https://www.ncbi.nlm.nih.gov/taxonomy
my $target_clade=4479; # or $target_clade='Poaceae';

#see endpoint doc at 
#https://rest.ensembl.org/documentation/info/homology_symbol
$url = 
	join('/', $server, 'homology/symbol', $species, $gene).
	"?content-type=application/json&compara=$division";


# restrict to homologues from this clade/taxon
if(defined($target_clade)){
	$url .= "&target_taxon=$target_clade";
}

my $homology_data = call_endpoint($http,$url);

my @homologies = @{$homology_data->{data}[0]{homologies}};

# filter out homologues based on type
if(defined($homoltype)){
	@homologies = grep {
		$_->{type} =~ m/$homoltype/
	} @homologies;
}

## R6) Get annotation of orthologous genes/proteins

# using the xrefs/id endpoint
# https://rest.ensembl.org/documentation/info/xref_id

my $total_annots = 0;

for my $homolog (@homologies) {

	last if($total_annots++ > 5); # for brevity

	my $target_species = $homolog->{target}{species};
	my $target_id = $homolog->{target}{id};
	my $target_prot_id = $homolog->{target}{protein_id};
	
	print "$gene\t$species\t$target_id\t$target_species\n";

	# GO annotation (protein)
	$url = join('/', $server, "xrefs/id/$target_prot_id").
		"?content-type=application/json;external_db=GO;all_levels=1";

	my $go_data = call_endpoint($http,$url);
	for my $go (@{$go_data}) {
		print $go->{dbname},': ', 
			$go->{display_id}, ' ', 
			$go->{description} || 'NA',
			' Evidence: ', 
			join(', ', @{$go->{linkage_types}}),
			"\n";
	}

	# check KEGG Enzyme annotation (protein)
	$url = join('/', $server, "xrefs/id/$target_prot_id").
	        "?content-type=application/json;".
			"external_db=KEGG_Enzyme";

	my $KE_data = call_endpoint($http,$url);
	for my $ke (@{$KE_data}) {
		if(defined($ke->{description})){
			print $ke->{dbname}, ': ',
				$ke->{display_id}, ' ',
				$ke->{description}, ' ',
				'Evidence: ', $ke->{info_type},
				"\n";
		}
	}

	# now check Plant Reactome annotation (gene)
	$url = join('/', $server, "xrefs/id/$target_id").
		"?content-type=application/json;".
		"external_db=Plant_Reactome_Pathway";

	my $PR_data = call_endpoint($http,$url);
	for my $pr (@{$PR_data}) {
		print $pr->{dbname}, ': ', 
			$pr->{display_id}, ' ', 
			'Evidence: ', $pr->{info_type},
			"\n";
	}
}

## R7) Fetch variant consequences for multiple variant ids

# Note: unless previous examples, this is a POST REST request, 
# where user data is posted to the server and after some time
# a response in parsed. Read more at:
# https://github.com/Ensembl/ensembl-rest/wiki/POST-Requests

$species = 'oryza_sativa';

$url = join('/', $server, "/vep/$species/id");

# max one thousand ids
my $variants = '{ "ids" : [ "10522356134" ] }';

my $vep_data = post_endpoint($http,$url,$variants);

print Dumper $vep_data;


## R8) Check consequences of single SNP within CDS sequence

# Note: you need the relevant transcript id from species of interest
# This query involves 2 consecutive REST calls

$species = 'triticum_aestivum';
my $transcript_id = 'TraesCS4B02G042700.1';
my $SNPCDScoord = 812;
my $SNPbase = 'T'; 

# convert CDS coords to genomic coords
$url = join('/', $server, "map/cds/$transcript_id/$SNPCDScoord..$SNPCDScoord").
        "?content-type=application/json;species=$species";

my $map_cds = call_endpoint($http,$url);

if(defined($map_cds->{mappings}->[0]->{seq_region_name})){

	my $SNPgenome_coord = 
		$map_cds->{mappings}->[0]->{seq_region_name} . 
		':' .
		$map_cds->{mappings}->[0]->{start} .
		'-' .
		$map_cds->{mappings}->[0]->{end};
	
	# fetch VEP consequences for this region
	$url = join('/', $server, "vep/$species/region/$SNPgenome_coord/$SNPbase").
	        "?content-type=application/json";

	my $conseq = call_endpoint($http,$url);

	if(defined($conseq->[0]->{allele_string})){
		
		foreach my $tcons ($conseq->[0]->{transcript_consequences}){
			
			printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
				$transcript_id,
				$SNPCDScoord,
				$conseq->[0]->{allele_string},
				$tcons->[0]->{biotype},
				$tcons->[0]->{codons} || 'NA',
				$tcons->[0]->{amino_acids} || 'NA',
				$tcons->[0]->{protein_start} || 'NA',
				$tcons->[0]->{impact},
				# not all variants have SIFT scores precomputed
				$tcons->[0]->{sift_prediction} || 'NA',
				$tcons->[0]->{sift_score} || 'NA' 
			);
		}
	} 
}

## R9) Retrieve variation sources of a species

$url = join('/', $server, "info/variation/$species") .
		"?content-type=application/json";

my $var_source_data = call_endpoint($http,$url);

for my $pr (@{$var_source_data}) {
	printf("%s\t%s\t%s\n", 
		$species, 
		$pr->{name},
		$pr->{description} );
}
