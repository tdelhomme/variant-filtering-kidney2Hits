#library(reprtree)
#library(caret)
library(randomForest)

args <- commandArgs(TRUE)
parseArgs <- function(x) {
  res = strsplit(sub("^--", "", x), "=")
  if(length(unlist(res))==1) res[[1]][2]=""
  return(res)
}
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$help)) {help=FALSE} else {help=TRUE}

if(is.null(args$model_snv) | is.null(args$model_snv) | is.null(args$target_table) | help) {
  cat("

      Method: add features to each variant of a table to be used in the classifier 

      Mandatory arguments:
      --model_snv=Robject         - trained random forest from snv (output of compute_model.r)
      --model_indel=Robject       - trained random forest from indel (output of compute_model.r)
      --target_table=file_name    - table to apply the model on (output of add_calling_features.r) 

      Optional arguments:
      --output_folder             - folder to save the models (default: .)
      --output_table              - file name for the table annotated with random forest predictions (default: RF_table_result.txt)
      --help                      - print this text

      example: add_calling_features.r --train_table=table_for_training.txt \n\n")
  q(save="no")
}

if(is.null(args$output_folder)) {output_folder="."} else {output_folder=args$output_folder}
if(is.null(args$output_table)) {output_table="RF_table_result.txt"} else {output_table=args$output_table}

model_snv = args$model_snv
model_indel = args$model_indel

set.seed(123) # here I fix the seed of the random generator to have the same random numbers if re-run 

target_table = read.table(args$target_table, quote="\"", stringsAsFactors=F, sep="\t", header=T)
target_table[which(is.infinite(target_table$MIN_DIST)),"MIN_DIST"] = 1000000000
target_table[which(is.infinite(target_table$MaxRatioWin)),"MaxRatioWin"] = 1000000000
target_table$IoD = 1 + target_table$SIG * 10^(target_table$ERR)
target_table$IoD = 1 + target_table$SIG * 10^(target_table$ERR)

for(type in c("snv","indel")){
  load(get(paste("model_",type,sep="")))
  if(type == "snv") { ids=which(target_table$TYPE_INFO=="snv") ; rf = rf_snv }
  if(type == "indel") { ids=which(target_table$TYPE_INFO!="snv") ; rf = rf_indel }
  target_table[ids,"FP_prob"] = predict(rf, target_table[ids,], type="prob")[,"FP"]
}


write.table(target_table, file=paste(output_folder,"/",output_table,sep=""), quote=F, sep="\t", row.names = F)
