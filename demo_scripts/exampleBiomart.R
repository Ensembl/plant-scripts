#!/usr/bin/env Rscript 

# Examples of R queries to Ensembl Plants Biomart
# Your usage of the data returned by the Biomart service is
# subject to same conditions as laid out on the Ensembl website.
#
# See documentation at 
# https://www.ensembl.org/info/data/biomart/how_to_use_biomart.html
# https://www.ensembl.org/info/data/biomart/biomart_r_package.html
#
# Copyright [2020] EMBL-European Bioinformatics Institute

# uncomment to install package biomaRt package
# if(!requireNamespace("BiocManager", quietly = TRUE))
# 	install.packages("BiocManager")
# BiocManager::install("biomaRt")

library("biomaRt")

## B1) Check plant marts and select dataset

listMarts( host="plants.ensembl.org" )

EPgenes = useEnsembl( biomart="plants_mart", 
			host="plants.ensembl.org")

dsets = listDatasets(EPgenes)

dsets[grep("Triticum aestivum", dsets$description),]
#              dataset                     description version
# 69 taestivum_eg_gene Triticum aestivum genes (IWGSC)   IWGSC

# take a note of the dataset name 'taestivum_eg_gene'


## B2) Check available filters and attributes

EPgenes = useMart( 
	biomart="plants_mart",
  	host="plants.ensembl.org",
	dataset="taestivum_eg_gene")

head( listFilters(EPgenes) )

head( listAttributes(EPgenes) )


## B3) Download GO terms associated to genes

# Note genes might appear in several rows

go = getBM( 
		attributes=c("ensembl_gene_id", "go_id"), 
		mart=EPgenes) 

head(go)

## B4) Get Pfam domains annotated in genes

EPgenes = useMart(  
	biomart="plants_mart",
	host="plants.ensembl.org",
	dataset="hannuus_eg_gene")

pfam = getBM(
		attributes=c("ensembl_gene_id", "pfam"),
		mart=EPgenes)

head(pfam)

## B5) Get SNP consequences from a selected variation source

# Note this requires connecting to a different mart (snp)
# Note this query takes a few minutes to run

EPvar = useMart( biomart="plants_variations",
        	host="plants.ensembl.org", 
			dataset="taestivum_eg_snp")

snp_source = c("EMS-induced mutation")

chrs = listFilterValues(mart=EPvar,
		filter="chr_name")

attribs = c(
	"refsnp_id", 
	"refsnp_source",
	"ensembl_gene_stable_id",
	"consequence_type_tv",
	"sift_prediction",
	"sift_score")

filts = c( 
	"variation_source", 
	"chr_name",
    "sift_prediction")

preds = c(
	"tolerated", # comment if unwanted
	"deleterious")

snps <- NULL
for(chr in chrs){
	print(chr) # show progress 

	for(pred in preds){
		print(pred) # show progress

		tmp_s <- getBM(
			attributes=attribs,
			filters=filts,
			values=list(
				variation_source=snp_source, 
				chr_name=chr, 
				sift_prediction=c(pred)),
			mart=EPvar)
		
		# append SNP batches to object snps
		if(is.null(snps)){
			snps<-tmp_s
		}else{
		    snps<-rbind(snps,tmp_s)
		}																			
	}
}

head(snps)
