#! /usr/bin/env Rscript

library(VariantAnnotation)

args <- commandArgs(TRUE)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$vcf)) {stop("no input VCF file")} else {vcf = args$vcf}
if(is.null(args$chunk_size)) {chunk_size = 10000} else {chunk_size = as.numeric(args$chunk_size)}
if(is.null(args$max_QVAL_gof)) {max_QVAL_gof = 50} else {max_QVAL_gof = as.numeric(args$max_QVAL_gof)}
if(is.null(args$max_QVAL_count)) {max_QVAL_count = 50} else {max_QVAL_count = as.numeric(args$max_QVAL_count)}
if(is.null(args$min_QVAL_count)) {min_QVAL_count = 20} else {min_QVAL_count = as.numeric(args$min_QVAL_count)}
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
    (AOobs - AOexp) ^2 / ((1 + AOexp)^2) # one at the denominator to deal with value between 0 and 1
  }
  
  all_gof_stat = unlist(lapply(1:nrow(AOexp_matrix), function(i){
    if(is.na(as.character(geno(vcf_chunk[i,1])$QVAL))){ #here we have QVAL_INV
      ids = which(QVAL_INV_matrix[i,]>max_QVAL_gof) } else { ids = which(QVAL_matrix[i,]<=max_QVAL_gof) }
    all_AOexp = AOexp_matrix[i,][ids]
    all_AOobs = AO_matrix[i,][ids]
    sum(mapply(gof_stat, all_AOobs, all_AOexp)) / length(all_AOobs) #normalize by number of samples
  }))
  
  # compute number of low QVAL variants
  n_qval=apply(QVAL_matrix,1,function(x){length(which(x>min_QVAL_count & x<max_QVAL_count))})
  n_qval_inv=apply(QVAL_INV_matrix,1,function(x){length(which(x>min_QVAL_count & x<max_QVAL_count))})
  
  #annotate the header of the chunk
  info(header(vcf_chunk))["GOF",]=list("1","Float","Goodness of fit based on the chi-squared statistic")
  
  name_qval=paste("N_QVAL_",min_QVAL_count,"_",max_QVAL_count,sep="")
  info(header(vcf_chunk))[name_qval,]=list("1","Integer",paste("Number of samples with ",min_QVAL_count,"<QVAL<",max_QVAL_count,sep=""))
  name_qval_inv=paste("N_QVAL_INV_",min_QVAL_count,"_",max_QVAL_count,sep="")
  info(header(vcf_chunk))[name_qval_inv,]=list("1","Integer",paste("Number of samples with ",min_QVAL_count,"<QVAL_INV<",max_QVAL_count,sep=""))
  

  #annotate the chunk with computed values
  info(vcf_chunk)[,"GOF"] = all_gof_stat
  info(vcf_chunk)[,name_qval] = n_qval
  info(vcf_chunk)[,name_qval_inv] = n_qval_inv
  
  #write out the annotated VCF
  con = file(output_vcf, open = "a")
  writeVcf(vcf_chunk, con)
  vcf_chunk = readVcf(vcf_con, "hg19")
  close(con)
}
