#!/usr/bin/env Rscript 

# Examples of R queries to Ensembl Plants Biomart
# Your usage of the data returned by the REST service is
# subject to same conditions as laid out on the Ensembl website.
#
# See documentation at 
# https://www.ensembl.org/info/data/biomart/how_to_use_biomart.html
# https://www.ensembl.org/info/data/biomart/biomart_r_package.html
#
# Copyright [2020] EMBL-European Bioinformatics Institute

# uncomment to install package biomaRt
# if(!requireNamespace("BiocManager", quietly = TRUE))
# 	install.packages("BiocManager")
# BiocManager::install("biomaRt")

library("biomaRt")

## 1) download GO terms associated to genes
mart <- useMart(host="plants.ensembl.org", biomart="plants_mart", 
	dataset="taestivum_eg_gene")
go <- getBM(attributes=c("ensembl_gene_id", "go_id"), mart=mart) 
head(go_e46)


## 2)  
The other common task is to look at the effects of the SNPs, I was doing a global analysis, and it was timing out (I couldn't figure out how to do it in an SQL query). So my solution was like this:

				v_source <- c("EMS-induced mutation")
				snps <- NULL
				for(chr in listFilterValues(mart = varmart, filter = "chr_name")){
					    print(chr)
						    for(pred in c("tolerated", "deleterious")){
								        print(pred)
										        tmp_s <- getBM(attributes=c("refsnp_id", "refsnp_source", "ensembl_gene_stable_id", "consequence_type_tv","sift_prediction","sift_score"), 
												              filters=c( "variation_source", "chr_name","sift_prediction"),
															                values=list(variation_source=v_source, chr_name=chr, sift_prediction=c(pred)),
																			              mart=varmart)
												        if(is.null(snps)){
															            snps<-tmp_s
																		        }else{
																					        snps<-rbind(snps,tmp_s)
																							        }
																									    }
				}
				snps
