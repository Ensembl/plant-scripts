#!/usr/bin/env Rscript

# Example in R of queries to the Ensembl REST endpoints 
# Your usage of the data returned by the REST service is 
# subject to same conditions as laid out on the Ensembl website.
#
# The full set of endpoints is documented 
# at http://rest.ensembl.org
#
# See https://github.com/Ensembl/ensembl-rest/wiki for
# examples in other programming languages

# Copyright [2017-2020] EMBL-European Bioinformatics Institute

library(httr)
library(jsonlite)

server = 'http://rest.ensembl.org'

## R1) helper functions 

# function for invoking endpoint
call_endpoint <- function(url, content_type) {

	result = GET(url, accept(content_type))
	stop_for_status(result)

	if(content_type == 'application/json'){
		return (fromJSON(content(result, "text", encoding="UTF-8")))
	} else {
		return (content(result, "text", encoding="UTF-8"))
	}
}

# function to post data to endpoint and 
# wait for response
#sub post_endpoint {
#	my ($http, $url, $data, $verbose) = @_;

#	print "Invoking $url (POST)\n" if($verbose);

#	my $response = $http->request('POST', $url, {	
#		headers => { 
#			'Content-type' => 'application/json',
#			'Accept' => 'application/json'
#		},
#		content => $data
#	});
	
#	die "# Failed POST request $url\n" unless $response->{success};

#	return decode_json($response->{content});
#}

## R2) Get metadata for all plant species 

url = paste(server, '/info/genomes/division/EnsemblPlants', sep="/")

# call url endpoint and get an array back,
# note that actual data structure changes across endpoints
metadata = call_endpoint(url,"application/json");

# parse the data from the response
for (row in 1:nrow(metadata)) {
	print(paste(
		metadata[ row, 'name'], 
		metadata[ row, 'strain'],
		metadata[ row, 'assembly_accession'],
		metadata[ row, 'base_count'],
		metadata[ row, 'assembly_level'],
		metadata[ row, 'has_peptide_compara'],
		metadata[ row, 'has_variations'],
		metadata[ row, 'has_genome_alignments'],
		metadata[ row, 'has_synteny'] , collapse="\t") );
}

