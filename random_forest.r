library(randomForest)
library(reprtree)
library(caret)
setwd("~/Documents/Models/variant-filtering")
set.seed(123) # here I fix the seed of the random generator to have the same random numbers if re-run 


train_table="K2H_AllVariants_NoMinAF_WES_samples_illuminaBED_annotated_with_coverage_INFO_GENO_status_supp_features.txt"
train_table = read.table(train_table, quote="\"", stringsAsFactors=F, sep="\t", header=T)

train_table = train_table[which(train_table$TYPE=="snv"),]

train_table[which(is.infinite(train_table$MIN)),"MIN_DIST"] = 1000000000

propTP = as.numeric(table(train_table$status)["TP"] / nrow(train_table))
propFP = as.numeric(table(train_table$status)["FP"] / nrow(train_table))

rf = randomForest(as.factor(status) ~ RVSB + AF + ERR + DP,
                  data = train_table,
                  importance = TRUE, # to allow us to inspect variable importance
                  ntree = 500, classwt = c(propTP, propFP), maxnodes=10) # weighted FP by propTP and TP by propFP

# look at variable importance
varImpPlot(rf)

# contruct the consensus tree and plot it
tree = ReprTree(rf, train_table, metric='d2')
plot(tree)

# plot first tree from the forest : reprtree:::plot.getTree(rf, k = 1)

# use K-fold
folds <- createFolds(train_table$status, 20)
kfold_spec = c()
kfold_sens = c()
kfold_FDR = c()
kfold_TDR = c()

for(i in 1:10){
  print(paste("fold: ",i,sep=""))
  test = train_table[folds[[i]],]
  train = train_table[-folds[[i]],]
  
  propTP = as.numeric(table(train$status)["TP"] / nrow(train))
  propFP = as.numeric(table(train$status)["FP"] / nrow(train))
  
  rf_fold = randomForest(as.factor(status) ~ RVSB + AF + DP + ERR + QVAL + AO,
                         data = train,
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


