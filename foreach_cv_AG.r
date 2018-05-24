if(method=="svm"){
  results <- foreach(fold = cv, .packages="e1071") %dopar% {
    testData <- data_to_classified[fold,] # Get the opposite of the test observations to train on
    trainData <- data_to_classified[-fold,]
    if(probs){
      svm_model=svm(trainData,cancer_to_classified[-fold],method = "C-classification", kernel=svm_method,class.weights = (100/summary(cancer_to_classified[-fold]))/sum(100/summary(cancer_to_classified[-fold])),probability=T,cost=cost)
      list(accuracy=data.frame(pred=predict(svm_model,testData,probability=T),cohort=cancer_to_classified[fold]),proba=attr(predict(svm_model,testData,probability=T),"probabilities"))
    }else{
      svm_model=svm(trainData,cancer_to_classified[-fold],method = "C-classification", kernel=svm_method,class.weights = (100/summary(cancer_to_classified[-fold]))/sum(100/summary(cancer_to_classified[-fold])),cost=cost)
      list(accuracy=data.frame(pred=predict(svm_model,testData),cohort=cancer_to_classified[fold]))
      
    }
  }
}

if(method=="rf"){
  results <- foreach(fold = cv, .packages="randomForest") %dopar% {
    testData <- data_to_classified[fold,] # Get the opposite of the test observations to train on
    trainData <- data_to_classified[-fold,]
    if(probs==F){
      model=randomForest(trainData,cancer_to_classified[-fold],importance=T,do.trace=F,classwt=(100/summary(cancer_to_classified[-fold]))/sum(100/summary(cancer_to_classified[-fold])))
      list(accuracy=data.frame(pred=predict(model,testData),cohort=cancer_to_classified[fold]))
    }else{
      model=randomForest(trainData,cancer_to_classified[-fold],importance=T,do.trace=F,classwt=(100/summary(cancer_to_classified[-fold]))/sum(100/summary(cancer_to_classified[-fold])),probs=T)
      list(accuracy=data.frame(pred=predict(model,testData),cohort=cancer_to_classified[fold]),proba=predict(model,testData,type="prob"))
    }
  }
}

merge_results=foreach(fold.result=results, fold.num=icount(), .combine=rbind) %do%{
  as.data.frame(fold.result$accuracy) 
}

stopCluster(cl)


probabilities=foreach(fold.result=results, fold.num=icount(), .combine=rbind) %do%{
  as.data.frame(fold.result$proba) 
}
#write("after probabilities","test.txt",append = T)
res=cbind(probabilities,merge_results)