
#!/usr/bin/env Rscript 
# based on 
# https://github.com/eead-csic-compbio/get_homologues/blob/master/pfam_enrich.pl

args = commandArgs(trailingOnly=TRUE)

if(length(args)<3) {
  stop("# Usage: <exp> <control> <out> filenames", call.=FALSE)
}

# globals
direction="greater" #"two.sides" 
multitest="fdr"
verbose=F

##query_file="aegilops_tauschii.Red.genes.ovl50.tsv"
#control_file="aegilops_tauschii.tsv.count.tsv"
#out_file="kk"

query_file=args[1]
control_file=args[2]
out_file=args[3]

# parse data
que_data=read.csv(query_file, sep="\t", header=FALSE);
que_rows=nrow(que_data);
ref_data=read.csv(control_file, sep="\t", header=FALSE);
ref_rows=nrow(ref_data); 

# uses globals que_total, ref_total
enrich_test <- function(x){
  que_id=x[1];
  que_value=as.numeric(x[2]);
  ref_value=as.numeric(ref_data[ref_data$V1==que_id,2]);
  if (length(ref_value)==0){
    ref_value=0;
  }
  if(verbose==T){
    cat(paste(que_id, "\n"), file=stderr());
  }
  values=c(que_value, ref_value, que_total, ref_total);
  input_matrix=matrix(values, nrow = 2,
                      dimnames=list(c("exp", "control"), c("Pfam", "total")));
  if(verbose==T){
    cat(paste(input_matrix, "\n"), file=stderr());
  }
  fisher_htest=fisher.test(input_matrix, alternative=direction);
  if(verbose==T){
    cat(paste(fisher_htest, "\n"), file=stderr());
  }
  ret_value=c( que_id, fisher_htest$p.value );
  return(ret_value);
}

print_pvalues <- function(x) {
  cat(paste(x[1],"\t",x[2],"\t",x[3],"\n"),file=out_file,append=T);
}

que_total=sum(que_data[,2]);
ref_total=sum(ref_data[,2]);
pvalues=apply(que_data, 1, enrich_test);

## Multiple test adjustment
num_pvalues=as.numeric(pvalues[2,]);
adj_pvalues=p.adjust(num_pvalues, method=multitest);
result_pvalues=rbind(pvalues, adj_pvalues);
result_pvalues=t(result_pvalues); # rows to columns
apply(result_pvalues, 1, print_pvalues);
