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
#see https://cran.r-project.org/web/packages/jsonlite/vignettes/json-apis.html
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

subset(metadata, select=c( 'name','strain','assembly_accession',
	'base_count','assembly_level','has_peptide_compara',
	'has_variations','has_genome_alignments','has_synteny'))

## R3) Find features overlapping genomic region

# full list at
# https://rest.ensembl.org/documentation/info/overlap_region

species = 'triticum_aestivum';
region = '3D:379400000-379540000';

# genes
url = paste( 
		paste(server, 'overlap/region', species, region, sep="/"),
		"feature=gene;content-type=application/json",
		sep="?")

overlap_data = call_endpoint(url,"application/json")
subset(overlap_data, select=c( 'id','start','end'))

# now LTR repeats
url = paste(
        paste(server, 'overlap/region', species, region, sep="/"),
		"feature=repeat;content-type=application/json",
		sep="?")

overlap_data = call_endpoint(url,"application/json")
subset(overlap_data, grepl("LTR", description), 
	select=c( 'description','start','end'))

# finally get EMS variants
url = paste(
        paste(server, 'overlap/region', species, region, sep="/"),
		"feature=variation;content-type=application/json",
		sep="?")

overlap_data = call_endpoint(url,"application/json")
subset(overlap_data, grepl("EMS-induced mutation", source),
    select=c( 'id','source'))

## R4) Fetch phenotypes overlapping genomic region

species = 'arabidopsis_thaliana'
region = '3:19431095-19434450'
p_cutoff = 0.0001

url = paste(
		paste(server, 'phenotype/region', species, region, sep="/"),
		"feature_type=Variation;content-type=application/json",
        sep="?")

pheno_data = call_endpoint(url,"application/json")
if(length(pheno_data) > 0) {
	for (row in 1:nrow(pheno_data)) {
		feat = pheno_data[[ row, 'phenotype_associations']]
		for (f in 1:nrow(feat)) {
			assoc = feat[ f, ]
			if(as.numeric(assoc$attributes$p_value) <= p_cutoff){
				print( paste(
					pheno_data[row,'id'],
					assoc$location,
					assoc$description,
					sep='\t') )
			}
		}
	}
}


## R5) Find homologues of selected gene

species = 'triticum_aestivum'
gene = 'TraesCS1B02G195200'
homoltype = 'ortholog' #paralog

# optional define a target clade, such as 4479
# for Poaceae, see https://www.ncbi.nlm.nih.gov/taxonomy
target_clade=4479

#see endpoint doc at
#https://rest.ensembl.org/documentation/info/homology_symbol

url = paste(
        paste(server, 'homology/symbol', species, gene, sep="/"),
	    "compara=plants;content-type=application/json",
		sep="?")

# restrict to homologues from this clade/taxon
if(!is.null(target_clade)){
	url = paste(url, "&target_taxon=", target_clade, sep="")
}

homology_data = call_endpoint(url, "application/json")

# extract homologies
homologies=homology_data$data$homologies[[1]]

# filter out homologues based on type
filt_homol = subset(homologies, grepl(homoltype, type))

## R6) Get annotation of orthologous genes/proteins

# using the xrefs/id endpoint
# https://rest.ensembl.org/documentation/info/xref_id

total_annots=0

for (hom in 1:nrow(filt_homol)) {

	total_annots = total_annots + 1
	if(total_annots > 10){
		break; # for brevity
	}

	target_species = filt_homol[ hom,'target']$species
	target_id = filt_homol[ hom,'target']$id
	target_prot_id = filt_homol[ hom,'target']$protein_id

    cat(paste(gene, species, target_id, target_species, "\n", sep="\t"))

	# GO annotation (protein)
	url = paste(
	        paste(server, 'xrefs/id', target_prot_id, sep="/"),
			"content-type=application/json;external_db=GO;all_levels=1",
			sep="?")

	go_data = call_endpoint(url, "application/json");
	if(length(go_data) > 0){
		print(subset(as.data.frame(go_data), select=c( 'dbname','display_id',
							'description','linkage_types')))
	}

	# check KEGG Enzyme annotation (protein)
	url = paste(
            paste(server, 'xrefs/id', target_prot_id, sep="/"),
            "content-type=application/json;external_db=KEGG_Enzyme",
            sep="?")

	KE_data = call_endpoint(url, "application/json");
	if(length(KE_data) > 0){
		print(subset(KE_data, select=c( 'dbname','display_id',
	                            'description','info_type')))
	}

	# now check Plant Reactome annotation (gene)
	url = paste(
			paste(server, 'xrefs/id', target_id, sep="/"),
			"content-type=application/json;external_db=Plant_Reactome_Pathway",
			sep="?")

	PR_data = call_endpoint(url, "application/json");
	if(length(PR_data) > 0){
		print(subset(PR_data, select=c( 'dbname','display_id', 'info_type')))
	}
}


