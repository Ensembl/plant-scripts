#!/usr/bin/env Rscript

# Installs missing R dependencies 

local_lib = "./lib/R/"
.libPaths( c( .libPaths(), local_lib) )

if(!requireNamespace("BiocManager", quietly=T))
	install.packages("BiocManager", dependencies=T, lib=local_lib)

install.packages("dplyr", "stringi", "knitr", "httr", "jsonlite", dependencies=T, lib=local_lib)

BiocManager::install("biomaRt", lib=local_lib, dependencies=T)

sessionInfo()
