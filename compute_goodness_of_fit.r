#! /usr/bin/env Rscript

library(VariantAnnotation)

args <- commandArgs(TRUE)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$vcf)) {stop("no input VCF file")} else {vcf = args$vcf}
if(is.null(args$chunk_size)) {chunk_size = 10000} else {chunk_size = as.numeric(args$chunk_size)}
if(is.null(args$max_QVAL)) {max_QVAL = 20} else {max_QVAL = as.numeric(args$max_QVAL)}
if(is.null(args$output_vcf)) {output_vcf=gsub(".vcf.gz","_GOF.vcf",vcf)} else {output_vcf=args$output_vcf}

#initiate the first chunk
vcf_con <- open(VcfFile(vcf,  yieldSize=chunk_size))
vcf_chunk = readVcf(vcf_con, "hg19")

#and continue
while(dim(vcf_chunk)[1] != 0) {
  # get AO, DP and QVAL
  AO_matrix = geno(vcf_chunk,"AO")
  DP_matrix = geno(vcf_chunk,"DP")
  QVAL_matrix = geno(vcf_chunk,"QVAL")
  QVAL_INV_matrix = geno(vcf_chunk,"QVAL_INV")
  
  # compute AO expected
  errors = info(vcf_chunk)$ERR
  AOexp_matrix = lapply(1:nrow(DP_matrix), function(i){ DP_matrix[i,] * errors[i] })
  AOexp_matrix = data.frame(matrix(unlist(AOexp_matrix), nrow=length(AOexp_matrix), byrow=T))
  
  # compute chi-squared statistic
  gof_stat = function(AOobs, AOexp){
    (AOobs - AOexp) ^2 / ((1 - AOexp)^2)
  }
  
  all_gof_stat = unlist(lapply(1:nrow(AOexp_matrix), function(i){
    if(is.na(as.character(geno(vcf_chunk[i,1])$QVAL))){ #here we have QVAL_INV
      ids = which(QVAL_INV_matrix[i,]>max_QVAL) } else { ids = which(QVAL_matrix[i,]<=max_QVAL) }
    all_AOexp = AOexp_matrix[i,][ids]
    all_AOobs = AO_matrix[i,][ids]
    sum(mapply(gof_stat, all_AOobs, all_AOexp)) / length(all_AOobs)
  }))
  
  #annotate the header of the chunk
  info(header(vcf_chunk))["GOF",]=list("1","Integer","Goodness of fit from the chi-squared statistic")

  #annotate the chunk with computed values
  info(vcf_chunk)[,"GOF"] = all_gof_stat

  #write out the annotated VCF
  con = file(output_vcf, open = "a")
  writeVcf(vcf_chunk, con)
  vcf_chunk = readVcf(vcf_con, "hg19")
  close(con)
}
