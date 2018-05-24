library(e1071)

model_svm <- svm(as.factor(status) ~ . , data=train_table[,my_features])

pred = predict(model_svm, train_table[,my_features])

sens_svm = sum(train_table$status=="TP" & pred=="TP") / sum(train_table$status=="TP")
TDR_svm = sum(train_table$status == "TP" & pred == "TP") / sum(pred == "TP")

kfold_spec = c()
kfold_sens = c()
kfold_TDR = c()

folds <- createFolds(train_table$status, 10)

for(i in 1:10){
  print(paste("fold: ",i,sep=""))
  test = train_table[folds[[i]],]
  train = train_table[-folds[[i]],]
    
  rf_fold = randomForest(as.factor(status) ~ .,
                         data = train[,my_features],
                         importance = TRUE, # to allow us to inspect variable importance
                         ntree = 500, sampsize = as.numeric(table(train$status)["TP"])) #classwt = c(propTP, propFP))
  test$prediction = predict(rf_fold, test)
  kfold_sens = c(kfold_sens, (sum(test$status == "TP" & test$prediction == "TP") / sum(test$status == "TP")) )
  kfold_spec = c(kfold_spec, (sum(test$status == "FP" & test$prediction == "FP") / sum(test$status == "FP")) )
  kfold_TDR = c(kfold_TDR, (sum(test$status == "TP" & test$prediction == "TP") / sum(test$prediction == "TP")) )
}

boxplot(data.frame("TDR"=kfold_TDR, "sensitivity"=kfold_sens), col="lightgrey", ylim=c(0.85,1), outpch=20, outcex=0.75, axes=F, staplelwd=2)
axis(side=3, at=c(1,2), labels=c("TDR","sensitivity"), col = NA, col.ticks = NA)
axis(side=2, at=c("0.85","0.90","0.95","1"))
