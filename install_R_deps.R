#!/usr/bin/env Rscript

# Installs missing R dependencies 

if(!requireNamespace("BiocManager", quietly=T))
	install.packages("BiocManager")
BiocManager::install("biomaRt")

#repository = 'https://cloud.r-project.org';

#required_packages = c("httr", "jsonlite")

#local_lib = "./lib/R/"

#.libPaths( c( .libPaths(), local_lib) )

#for (package in required_packages) {
#	if (!require(package, character.only=T, quietly=T)) {
#		sprintf("# cannot load %s, will get it from %s and install it in %s",
#			package,repository,local_lib)
#		install.packages(package, dependencies=TRUE, lib=local_lib, repos=repository)
#	}
#}

sessionInfo()
