library(foreach)
library(doParallel)

train_table="tables/K2H_AllVariants_All_Lib_NoMinAF_noequal_GOF_addVCFfeatures_illuminaBED_WES_samples_annotated_with_coverage_INFO_GENO_status_supp_features.txt"
train_table = read.table(train_table, quote="\"", stringsAsFactors=F, sep="\t", header=T)

train_table$IoD = 1 + train_table$SIG * 10^(train_table$ERR)
train_table = train_table[which(train_table$TYPE_INFO=="snv"),]
train_table[which(is.infinite(train_table$MIN_DIST)),"MIN_DIST"] = 1000000000
train_table[which(is.infinite(train_table$MaxRatioWin)),"MaxRatioWin"] = 1000000000

kfold_auc_rf = c()
kfold_auc_mlp = c()
kfold_auc_svm = c()

folds <- createFolds(train_table$status, 10)

cores=detectCores()
cl <- makeCluster(cores/2) #not to overload your computer
registerDoParallel(cl)

res = foreach(i=1:10, .combine=cbind, .packages = c("randomForest","monmlp","e1071","ROCR")) %dopar% {
  print(paste("fold: ",i,sep=""))
  test = train_table[folds[[i]],]
  train = train_table[-folds[[i]],]
  
  train$target = 0 ; train[which(train$status=="TP"),"target"] = 1
  test$target = 0 ; test[which(test$status=="TP"),"target"] = 1
  
  my_features=c("target","RVSB_INFO", "QVAL","AF","ERR_INFO","DP",
                "FS", "MIN_DIST", "AO", "QUAL", "MaxRatioWin", "NbVarWin", "IoD", "HpLength")
  
  # random forest
  print("PERFORMING RANDOM FORESTS")
  rf_fold = randomForest(target ~ ., data = train[,my_features], importance = TRUE,
                         ntree = 500, sampsize = as.numeric(table(train$status)["TP"])) 
  test$prediction = predict(rf_fold, test)
  auc_rf = performance( prediction(test$prediction, test$target), "auc" )@y.values[[1]]
  plot(performance(prediction(test$prediction, test$target), "prec","rec" ), colorize=T, xlim=c(1,0.95), colorize.palette="darkred")
  
  # support vector machine
  print("PERFORMING SVMs")
  svm_fold = svm(target ~ . , data=train[,my_features])  
  test$prediction = predict(svm_fold, test[,my_features])
  auc_svm = performance( prediction(test$prediction, test$target), "auc" )@y.values[[1]]
  plot(performance(prediction(test$prediction, test$target), "prec","rec" ), colorize=T, xlim=c(1,0.95), add=T, colorize.palette="darkgreen")

  # multilayers perceptron
  print("PERFORMING MLPs")
  x = data.matrix(train[,my_features[my_features!="target"]])
  y = data.matrix(train$target)
  model_mlp <- monmlp.fit(x, y, hidden1=5, n.ensemble=15, monotone=1, bag=TRUE)
  test$prediction = monmlp.predict(x = data.matrix(test[,my_features[my_features!="target"]]), weights = model_mlp)
  auc_mlp = performance( prediction(test$prediction, test$target), "auc" )@y.values[[1]]
  plot(performance(prediction(test$prediction, test$target), "prec","rec" ), colorize=T, xlim=c(1,0.95), add=T, colorize.palette="darkblue")
  
  c(auc_rf, auc_mlp, auc_svm)
}
stopCluster(cl)


boxplot(data.frame("RF"=res[1,], "SVM"=res[2,], "MLP"=res[3,]))


