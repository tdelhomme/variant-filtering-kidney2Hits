#! /usr/bin/env Rscript

args <- commandArgs(TRUE)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$vcf)) {stop("no input VCF file")} else {vcf = args$vcf}
if(is.null(args$min_q)) {min_q = 20} else {min_q = as.numeric(args$min_q)}
if(is.null(args$max_q)) {max_q = 50} else {max_q = as.numeric(args$max_q)}
if(is.null(args$chunk_size)) {chunk_size = 10000} else {chunk_size = as.numeric(args$chunk_size)}

output_vcf=paste(gsub(".vcf.gz","",vcf),"_q",min_q,"_q",max_q,".vcf",sep="")

library(VariantAnnotation)

#initiate the first chunk
vcf_con <- open(VcfFile(vcf,  yieldSize=chunk_size))
vcf_chunk = readVcf(vcf_con, "hg19")

#and continue
while(dim(vcf_chunk)[1] != 0) {
  # count the qvals
  QVAL_matrix = geno(vcf_chunk,"QVAL")
  QVAL_INV_matrix = geno(vcf_chunk,"QVAL_INV")
  
  n_qval=apply(QVAL_matrix,1,function(x){length(which(x>min_q & x<max_q))})
  n_qval_inv=apply(QVAL_INV_matrix,1,function(x){length(which(x>min_q & x<max_q))})
  
  #annotate the header of the chunk
  name_qval=paste("N_QVAL_",min_q,"_",max_q,sep="")
  info(header(vcf_chunk))[name_qval,]=list("1","Integer",paste("Number of samples with ",min_q,"<QVAL<",max_q,sep=""))
  name_qval_inv=paste("N_QVAL_INV_",min_q,"_",max_q,sep="")
  info(header(vcf_chunk))[name_qval_inv,]=list("1","Integer",paste("Number of samples with ",min_q,"<QVAL_INV<",max_q,sep=""))
  
  #annotate the chunk with computed values
  info(vcf_chunk)[,name_qval] = n_qval
  info(vcf_chunk)[,name_qval_inv] = n_qval_inv
  
  #write out the annotated VCF
  con = file(output_vcf, open = "a")
  writeVcf(vcf_chunk, con)
  vcf_chunk = readVcf(vcf_con, "hg19")
  close(con)
}
