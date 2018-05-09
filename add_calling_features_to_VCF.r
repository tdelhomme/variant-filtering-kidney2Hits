#! /usr/bin/env Rscript

library(VariantAnnotation)

args <- commandArgs(TRUE)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$vcf)) {stop("no input VCF file")} else {vcf = args$vcf}
if(is.null(args$chunk_size)) {chunk_size = 10000} else {chunk_size = as.numeric(args$chunk_size)}
if(is.null(args$output_vcf)) {output_vcf=gsub(".vcf.gz","_addVCFfeatures.vcf",vcf)} else {output_vcf=args$output_vcf}

#initiate the first chunk
vcf_con <- open(VcfFile(vcf,  yieldSize=chunk_size))
vcf_chunk = readVcf(vcf_con, "hg19")

#and continue
while(dim(vcf_chunk)[1] != 0) {
  # get DP
  DP_matrix = geno(vcf_chunk,"DP")

  # compute median DP
  median_DP = apply(DP_matrix, 1, function(x){ median(x) })
  
  #annotate the header of the chunk
  info(header(vcf_chunk))["medianDP",]=list("1","Float","Median DP over all samples")
  
  #annotate the chunk with computed values
  info(vcf_chunk)[,"medianDP"] = median_DP
  
  #write out the annotated VCF
  con = file(output_vcf, open = "a")
  writeVcf(vcf_chunk, con)
  vcf_chunk = readVcf(vcf_con, "hg19")
  close(con)
}