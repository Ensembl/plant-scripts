#!/usr/bin/env python

from __future__ import print_function ##Allow Python 3.X printing for Python 2
import requests, sys
import re

# Example of python queries to the Ensembl REST endpoints
# Your usage of the data returned by the REST service is
# subject to same conditions as laid out on the Ensembl website.
#
# The full set of endpoints is documented
# at http://rest.ensembl.org
#
# See https://github.com/Ensembl/ensembl-rest/wiki for
# examples in other programming languages

# Copyright [2020] EMBL-European Bioinformatics Institute

#=====================================================
# helper functions that fetch JSON from REST endpoints

def get_json(ext):
    '''simple get function'''	

    server = "https://rest.ensembl.org"
    url = server+ext
    r = requests.get(url, headers={ "Content-Type" : "application/json"})
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    decoded = r.json()
    return decoded

def get_json_post(ext,ids):
    '''simple post function'''

    server = "https://rest.ensembl.org"
    url = server+ext
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    r = requests.post(server+ext, headers=headers, data=ids)
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    decoded = r.json()
    return decoded


#===========================================
# example funtions that query REST endpoints 


def get_metadata():
    '''get metadata for plants division
    More examples: 
    https://rest.ensembl.org/documentation/info/info_genomes_division
    '''

    ext = "/info/genomes/division/EnsemblPlants?"
    decoded = get_json(ext)

    for d in decoded:
        
        # Printing relevant info, notice some items have been 
        # converted to str for easier printing
        meta = (d['name'],d['assembly_accession'],str(d['base_count']),
            d['assembly_level'],str(d['has_peptide_compara']),
            str(d['has_variations']),str(d['has_genome_alignments']),
            str(d['has_synteny']));
        separator = '\t'
        print(separator.join(meta))


def get_overlapping_features(species,region):
    '''get genes, repeats & variants overlapping a chosen region
    Note: produces 1-based inclusive coordinates
    More examples:
    https://rest.ensembl.org/documentation/info/overlap_region
	'''
    ## get the genes via the API
    overlap_url = ("/overlap/region/" + species + "/" + region)
    ext = (overlap_url + "?feature=gene;content-type=application/json")
    overlap_data = get_json(ext)
    for overlap_feat in overlap_data:
        print("%s\t%s\t%s" % 
		(overlap_feat['id'],overlap_feat['start'],overlap_feat['end']))

    ## now LTR repeats 
    ext = (overlap_url + "?feature=repeat;content-type=application/json");
    overlap_data = get_json(ext)

    for overlap_feat in overlap_data:
        ltr_match = re.search('LTR', overlap_feat['description'])
        if ltr_match:
            print("%s\t%s\t%s" % (overlap_feat['description'],overlap_feat['start'],
                overlap_feat['end']))


    ## search for a specific variation source
    source_variation = 'EMS-induced mutation';
    #source_variation = 'CerealsDB';
    
    ext = (overlap_url + "?feature=variation;content-type=application/json");
    overlap_data = get_json(ext)

    for overlap_feat in overlap_data:
        if overlap_feat['source'] == source_variation:
            print("%s\t%s" % (overlap_feat['id'],overlap_feat['source']))

    ## protein-coding genes from additional annotation tracks,
    ## also called otherfeatures dbs
    dbtype = 'otherfeatures'

    ext = (overlap_url + 
        "?feature=transcript;db_type=" + dbtype + 
        ";content-type=application/json")
    overlap_data = get_json(ext)    

    for overlap_feat in overlap_data:
        if overlap_feat['biotype'] == 'protein_coding':

            ext2 = ( "/sequence/id/" + overlap_feat['id'] +
                "?db_type=" + dbtype + ";type=protein;object_type=transcript;" +
                "species=" + species + ";content-type=application/json" )
           
            prot_data = get_json(ext2)
            print(">%s\n%s" % (prot_data['id'],prot_data['seq']))  


def get_phenotypes(species,region,p_cutoff):
    '''gets phenotypes associated to a genomic region under a P-value cutoff
    More examples: 
    https://rest.ensembl.org/documentation/info/phenotype_region
    '''

    ext = ('/phenotype/region/' + species + "/" + region + 
              "?feature_type=Variation;content-type=application/json")
    pheno_data = get_json(ext)
    
    for feat in pheno_data:
        for assoc in feat['phenotype_associations']:
            pval = float(assoc['attributes']['p_value']) #pval is a float
            if pval <= p_cutoff:
                print("%s\t%s\t%s" % 
					(feat['id'], assoc['location'], assoc['description']))


def find_homologues(species,division,gene,homoltype,target_clade=None):
    '''
    Creates a homologue list for a given gene and species
    More examples here: https://rest.ensembl.org/documentation/info/homology_symbol
    '''

    ext = ('/homology/symbol/' + species + "/" + gene + 
             "?content-type=application/json&compara="+division)

    if target_clade:
        ext = ext + "&target_taxon=" + target_clade
    
    homology_data = get_json(ext)
    homologies = homology_data['data'][0]['homologies']
    
    ##create a filtered list of homolgies according to homoltype
    filtered_homologies = []
    for hom in homologies:
        if re.search(homoltype,hom['type']):
            filtered_homologies.append(hom)
    return filtered_homologies


def parse_homologies(homologies):
    '''
    Parse homologs in a homology list
    '''
    for homolog in homologies:
        target_species = homolog['target']['species']
        target_id      = homolog['target']['id']
        target_prot_id = homolog['target']['protein_id']
        
        #GO annotation (protein)
        ext = ( "/xrefs/id/" + target_prot_id +
                "?content-type=application/json;external_db=GO;all_levels=1" )

        go_data = get_json(ext)
        for go in go_data:
            separator = ','
            linkage_types =  separator.join(go['linkage_types'])
            
            ##Default NA for description
            if 'description' not in go or not go['description']:
                go['description'] = 'NA'
                
            #print linkage_types
            go_vars = (go['dbname'],go['display_id'],go['description'],linkage_types)
            print("%s: %s,%s Evidence: %s" % go_vars)


        #check KEGG Enzyme annotation (protein)
        ext = ( "/xrefs/id/" + target_prot_id + 
                "?content-type=application/json;external_db=KEGG_Enzyme" )
        KE_data = get_json(ext)
        for ke in KE_data:
            if  'description' in ke and  ke['description']:
                ke_vars = (ke['dbname'],ke['display_id'],ke['description'],ke['info_type'])
                print("%s: %s,%s Evidence: %s" % ke_vars)
            
        #now check Plant Reactome annotation (gene)
        ext = ( "/xrefs/id/" + target_id + 
                "?content-type=application/json;external_db=Plant_Reactome_Pathway" )

        PR_data = get_json(ext)
        for pr in PR_data:
            print (pr['dbname'] + ': ' + pr['display_id'] + ' ' +
                        'Evidence: ' + pr['info_type'])

def get_vep_data(species,variants_ids):
    '''
    Note: unlike previous examples, this is a POST REST request, 
    where user data is posted to the server and after some time
    a response in parsed. Read more at:
    https://github.com/Ensembl/ensembl-rest/wiki/POST-Requests 
    '''
    ext = "/vep/" + species + "/id"

    vep_data = get_json_post(ext,variants_ids)

    print (vep_data)


def check_snp_consequences(species,transcript_id,SNPCDScoord,SNPbase):
    '''
    checks the consequences of a given SNP
    More examples: https://rest.ensembl.org/documentation/info/vep_region_get
    '''
    
    # convert CDS coords to genomic coords
    # using: https://rest.ensembl.org/documentation/info/assembly_cds
    ext = ("/map/cds/" + transcript_id + "/" + SNPCDScoord + ".." + SNPCDScoord 
           + "?content-type=application/json;species=" + species)
    map_cds = get_json(ext)
    
    #Check if we have all the data in the JSON output
    try:
        if map_cds['mappings'][0]['seq_region_name']:
            mapping = map_cds['mappings'][0]
    except:
        print("# ERROR: failed mapping CDS coords")
        return

    # fetch VEP consequences for this region
    SNPgenome_coord = ( mapping['seq_region_name'] + ':' + 
        str(mapping['start']) + '-' + str(mapping['end']) )
    ext = ("/vep/"+ species + "/region/" + SNPgenome_coord + "/" + 
        SNPbase + "?content-type=application/json")
    conseq = get_json(ext)
 
    ##Print all the relevant info for the given variant
    if conseq[0]['allele_string']:
        for tcons in conseq[0]['transcript_consequences']:
            
            ##Making sure NA is default in void parameters
            defaults = ['codons', 'amino_acids', 'protein_start', 
                'sift_prediction', 'sift_score']
            for d in defaults:
                if d not in tcons or not tcons[d]:
                    tcons[d] = 'NA'
            
            ##Values to print for the given variant
            values = (transcript_id,SNPCDScoord,conseq[0]['allele_string'],
                      tcons['biotype'],tcons['codons'],tcons['amino_acids'],
                      tcons['protein_start'],tcons['impact'],tcons['sift_prediction'],
                      tcons['sift_score'])
            for val in values:
                print (val, end="\t")
            print()

def get_variation_sources(species):
    '''retreive the variation sources of a species'''

    ext = ('/info/variation/' + species + 
        "?content-type=application/json")
    var_source_data = get_json(ext)

    for source in var_source_data:
        print("%s\t%s\t%s" % (species, 
            source['name'], source['description']))


#======================================== 
#Main

## R2) Get metadata for all plant species
get_metadata() #function call

## R3) Find features overlapping genomic region
species = 'triticum_aestivum';
region = '3D:379400000-379540000';
get_overlapping_features(species,region) #function call

## R4) Fetch phenotypes overlapping genomic region
species  = 'arabidopsis_thaliana'
region   = '3:19431095-19434450'
p_cutoff = 0.0001
get_phenotypes(species,region,p_cutoff) #function call

## R5) Find homologues of selected gene
#see endpoint doc at 
#https://rest.ensembl.org/documentation/info/homology_symbol
species   = 'triticum_aestivum' 
division  = 'plants'
gene      = 'TraesCS1B02G195200'
homoltype = 'ortholog' #can also be paralog

# optionally define a target clade, such as 4479
# for Poaceae, see https://www.ncbi.nlm.nih.gov/taxonomy
target_clade = '4479' # or target_clade = 'Poaceae' 
homologies = find_homologues(species,division,gene,homoltype,target_clade) #function call

## R6) Parse homologies and get annotation of orthologous genes/proteins
# using the xrefs/id endpoint
# https://rest.ensembl.org/documentation/info/xref_id

parse_homologies(homologies) #function call

## R7) Fetch variant consequences for multiple variant ids
species   = 'oryza_sativa' 
# no more than 1000 ids in one call, please
variants_json = '{ "ids" : [ "10522356134" ] }'
get_vep_data(species,variants_json) #function call

## R8) Check consequences of single SNP within CDS sequence
# Note: you need the relevant transcript id from species of interest
# This query involves 2 consecutive REST calls
species = 'triticum_aestivum'
transcript_id = 'TraesCS4B02G042700.1'
SNPCDScoord = '812'
SNPbase = 'T'
check_snp_consequences(species,transcript_id,SNPCDScoord,SNPbase) 

## R9) Retrieve variation sources of a species

get_variation_sources(species)

