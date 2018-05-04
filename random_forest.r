library(reprtree)
library(caret)
setwd("~/Documents/Models/variant-filtering")
set.seed(123) # here I fix the seed of the random generator to have the same random numbers if re-run 


train_table="K2H_AllVariants_NoMinAF_WES_samples_illuminaBED_annotated_with_coverage_INFO_GENO_status_supp_features.txt"
train_table = read.table(train_table, quote="\"", stringsAsFactors=F, sep="\t", header=T)

train_table = train_table[which(train_table$TYPE=="snv"),]

train_table[which(is.infinite(train_table$MIN_DIST)),"MIN_DIST"] = 1000000000
train_table[which(is.infinite(train_table$MaxRatioWin)),"MaxRatioWin"] = 1000000000

propTP = as.numeric(table(train_table$status)["TP"] / nrow(train_table))
propFP = as.numeric(table(train_table$status)["FP"] / nrow(train_table))

train_table$IoD = 1 + train_table$SIG * 10^(train_table$ERR)

my_features=c("status","RVSB", "FS","AF","ERR","DP", "QVAL", "MIN_DIST", "AO", "QUAL", "MaxRatioWin", "NbVarWin", "IoD")

rf = randomForest(as.factor(status) ~ .,
                  data = train_table[,my_features],
                  importance = TRUE, # to allow us to inspect variable importance
                  ntree = 500, classwt = c(propTP, propFP)) # weighted FP by propTP and TP by propFP

train_table$prediction = predict(rf, train_table)

# look at variable importance
varImpPlot(rf)

# contruct the representative tree and plot it
# definition: Representative trees are, in some sense, trees in the ensemble which are on average the "closest" to 
#             all the other trees in the ensemble.
tree = ReprTree(rf, train_table, metric='d2')
plot(tree)

# plot first tree from the forest : reprtree:::plot.getTree(rf, k = 1)

# use K-fold
folds <- createFolds(train_table$status, 10)
kfold_spec = c()
kfold_sens = c()
kfold_TDR = c()

for(i in 1:10){
  print(paste("fold: ",i,sep=""))
  test = train_table[folds[[i]],]
  train = train_table[-folds[[i]],]
  
  propTP = as.numeric(table(train$status)["TP"] / nrow(train))
  propFP = as.numeric(table(train$status)["FP"] / nrow(train))
  
  rf_fold = randomForest(as.factor(status) ~ .,
                         data = train[,my_features],
                         importance = TRUE, # to allow us to inspect variable importance
                         ntree = 500, classwt = c(propTP, propFP), maxnodes=10)
  test$prediction = predict(rf_fold, test)
  kfold_sens = c(kfold_sens, (sum(test$status == "TP" & test$prediction == "TP") / sum(test$status == "TP")) )
  kfold_spec = c(kfold_spec, (sum(test$status == "FP" & test$prediction == "FP") / sum(test$status == "FP")) )
  kfold_TDR = c(kfold_TDR, (sum(test$status == "TP" & test$prediction == "TP") / sum(test$prediction == "TP")) )
}

boxplot(data.frame("TDR"=kfold_TDR,
                   "sensitivity"=kfold_sens,
                   "specificity"=kfold_spec), col="lightgrey", ylim=c(0.8,1))

# compare with home filters

train_table$prediction_home_filters = "FP"
train_table[which(train_table$AF>0.1 & train_table$RVSB<0.95 & train_table$DP>50 & train_table$ERR<(-2) & train_table$AO>3),"prediction_home_filters"] = "TP"
sens_home_filter = sum(train_table$status=="TP" & train_table$prediction_home_filters=="TP") / sum(train_table$status=="TP")
TDR_home_filters = sum(train_table$status == "TP" & train_table$prediction_home_filters == "TP") / sum(train_table$prediction_home_filters == "TP")

sens_rf = sum(train_table$status=="TP" & train_table$prediction=="TP") / sum(train_table$status=="TP")
TDR_rf = sum(train_table$status == "TP" & train_table$prediction == "TP") / sum(train_table$prediction == "TP")
