library(reprtree)
library(caret)
library(randomForest)
setwd("~/Documents/Models/variant-filtering")
set.seed(123) # here I fix the seed of the random generator to have the same random numbers if re-run 

test_table="downsampling/K2H_AllVariants_All_Lib_downsampling_noequal_GOF_addVCFfeatures_illuminaBED_WES_samples_annotated_with_coverage_INFO_GENO_status_supp_features.txt"
test_table = read.table(test_table, quote="\"", stringsAsFactors=F, sep="\t", header=T)
test_table[which(is.infinite(test_table$MIN_DIST)),"MIN_DIST"] = 1000000000
test_table[which(is.infinite(test_table$MaxRatioWin)),"MaxRatioWin"] = 1000000000
test_table$IoD = 1 + test_table$SIG * 10^(test_table$ERR)
test_table$IoD = 1 + test_table$SIG * 10^(test_table$ERR)

for(type in c("snv","indel")){
  load(paste("downsampling/RF_model_downsample_",type,".Rdata",sep=""))
  if(type == "snv") { ids=which(test_table$TYPE_INFO=="snv") ; rf = rf_snv }
  if(type == "indel") { ids=which(test_table$TYPE_INFO!="snv") ; rf = rf_indel }
  test_table[ids,"FP_prob"] = predict(rf, test_table[ids,], type="prob")[,"FP"]
}

