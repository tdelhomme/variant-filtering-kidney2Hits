library(reprtree)
library(caret)
library(randomForest)
setwd("~/Documents/Models/variant-filtering")
set.seed(123) # here I fix the seed of the random generator to have the same random numbers if re-run 

train_table="tables/K2H_AllVariants_All_Lib_NoMinAF_noequal_GOF_addVCFfeatures_illuminaBED_WES_samples_annotated_with_coverage_INFO_GENO_status_supp_features.txt"
test_table="tables/K2H_AllVariants_All_Lib_NoMinAF_noequal_GOF_addVCFfeatures_illuminaBED_WES_samples_annotated_with_coverage_INFO_GENO_status_supp_features.txt"

# if we want to test the downsampling
#train_table="downsampling/K2H_AllVariants_All_Lib_downsampling_noequal_GOF_addVCFfeatures_illuminaBED_WES_samples_annotated_with_coverage_INFO_GENO_status_supp_features.txt"
#test_table="downsampling/K2H_AllVariants_All_Lib_downsampling_noequal_GOF_addVCFfeatures_illuminaBED_WES_samples_annotated_with_coverage_INFO_GENO_status_supp_features.txt"

train_table = read.table(train_table, quote="\"", stringsAsFactors=F, sep="\t", header=T)
test_table = read.table(test_table, quote="\"", stringsAsFactors=F, sep="\t", header=T)

train_table$IoD = 1 + train_table$SIG * 10^(train_table$ERR)
test_table$IoD = 1 + test_table$SIG * 10^(test_table$ERR)

type="snv"
if(type=="snv"){
  train_table = train_table[which(train_table$TYPE_INFO=="snv"),]  
}
if(type=="indel"){
  train_table = train_table[which(train_table$TYPE_INFO=="ins" | train_table$TYPE_INFO=="del"),]  
}

train_table[which(is.infinite(train_table$MIN_DIST)),"MIN_DIST"] = 1000000000
train_table[which(is.infinite(train_table$MaxRatioWin)),"MaxRatioWin"] = 1000000000
test_table[which(is.infinite(test_table$MIN_DIST)),"MIN_DIST"] = 1000000000
test_table[which(is.infinite(test_table$MaxRatioWin)),"MaxRatioWin"] = 1000000000

propTP = as.numeric(table(train_table$status)["TP"] / nrow(train_table))
propFP = as.numeric(table(train_table$status)["FP"] / nrow(train_table))

my_features=c("status","RVSB", "QVAL","AF","ERR_INFO","DP", "medianDP_INFO",
              "FS", "MIN_DIST", "AO", "QUAL", "MaxRatioWin", "NbVarWin", "IoD", "HpLength", "N_QVAL_20_50_INFO")

if(exists("dbsnp_status")){
  if(dbsnp){
    train_table$WES_status = train_table$status
    train_table$status = "FP"
    train_table[which(!is.na(train_table$avsnp150)),"status"] = "TP"
  }
}

rf = randomForest(as.factor(status) ~ .,
                  data = train_table[,my_features],
                  importance = TRUE, # to allow us to inspect variable importance
                  ntree = 500, sampsize = as.numeric(table(train_table$status)["TP"]) #classwt = c(propTP, propFP), # weighted FP by propTP and TP by propFP
                  # ,maxnodes=10, nodesize=20
                  ) 

test_table$prediction = predict(rf, test_table)
sens = sum(test_table$status == "TP" & test_table$prediction == "TP") / sum(test_table$status == "TP")
fdr = 1 - (sum(test_table$status == "TP" & test_table$prediction == "TP") / sum(test_table$prediction == "TP"))

# look at variable importance
varImpPlot(rf)

# contruct the representative tree and plot it
# definition: Representative trees are, in some sense, trees in the ensemble which are on average the "closest" to 
#             all the other trees in the ensemble.
tree = ReprTree(rf, train_table, metric='d2')
plot(tree)

# plot first tree from the forest : reprtree:::plot.getTree(rf, k = 1)

# use K-fold
kfold_spec = c()
kfold_sens = c()
kfold_TDR = c()

folds <- createFolds(test_table$status, 10)

for(i in 1:10){
  print(paste("fold: ",i,sep=""))
  test = test_table[folds[[i]],]
  train = test_table[-folds[[i]],]
  
  propTP = as.numeric(table(train$status)["TP"] / nrow(train))
  propFP = as.numeric(table(train$status)["FP"] / nrow(train))
  
  rf_fold = randomForest(as.factor(status) ~ .,
                         data = train[,my_features],
                         importance = TRUE, # to allow us to inspect variable importance
                         ntree = 500, sampsize = as.numeric(table(train$status)["TP"])) #classwt = c(propTP, propFP))
  test$prediction = predict(rf_fold, test)
  kfold_sens = c(kfold_sens, (sum(test$status == "TP" & test$prediction == "TP") / sum(test$status == "TP")) )
  kfold_spec = c(kfold_spec, (sum(test$status == "FP" & test$prediction == "FP") / sum(test$status == "FP")) )
  kfold_TDR = c(kfold_TDR, (sum(test$status == "TP" & test$prediction == "TP") / sum(test$prediction == "TP")) )
}

if(type=="snv") { boxplot(data.frame("TDR"=kfold_TDR, "sensitivity"=kfold_sens), col="lightgrey", ylim=c(0.85,1), outpch=20, 
        outcex=0.75, axes=F, staplelwd=2, main=type)
axis(side=3, at=c(1,2), labels=c("TDR","sensitivity"), col = NA, col.ticks = NA)
axis(side=2, at=c("0.85","0.90","0.95","1")) }

if(type=="indel") { boxplot(data.frame("TDR"=kfold_TDR, "sensitivity"=kfold_sens), col="lightgrey", ylim=c(0,1), outpch=20, 
                          outcex=0.75, axes=F, staplelwd=2, main=type)
  axis(side=3, at=c(1,2), labels=c("TDR","sensitivity"), col = NA, col.ticks = NA)
  axis(side=2, at=seq(0,1,by=0.1)) }

# load functions from MLP.r and run plot_ggviolin(kfold_TDR, kfold_sens, type)

# compare with home filters
train_table$prediction_home_filters = "FP"
train_table[which(train_table$AF>0.1 & train_table$RVSB<0.95 & train_table$DP>50 & train_table$ERR<(-2)),"prediction_home_filters"] = "TP"
sens_home_filter = sum(train_table$status=="TP" & train_table$prediction_home_filters=="TP") / sum(train_table$status=="TP")
TDR_home_filters = sum(train_table$status == "TP" & train_table$prediction_home_filters == "TP") / sum(train_table$prediction_home_filters == "TP")

points(1, TDR_home_filters, pch=8, col="darkred", lwd=2)
points(2, sens_home_filter, pch=8, col="darkred", lwd=2)

legend(x=2, y=0.88, c("AF>0.1","RVSB<0.95","DP>50","ERR<0.01"), col=c("black","white","white","white"), pch=4, bty="n", cex=0.75)

sens_rf = sum(train_table$status=="TP" & train_table$prediction=="TP") / sum(train_table$status=="TP")
TDR_rf = sum(train_table$status == "TP" & train_table$prediction == "TP") / sum(train_table$prediction == "TP")


