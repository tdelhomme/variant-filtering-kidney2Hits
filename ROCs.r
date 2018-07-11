
### ROCs ###
x = data.matrix(train_table[,my_features[my_features!="status"]])

y = train_table$status
y[which(y=="TP")] = 1
y[which(y=="FP")] = 0
y=as.matrix(as.numeric(y))
model_mlp <- monmlp.fit(x, y, hidden1=3, n.ensemble=1, monotone=1, bag=TRUE)

pred = monmlp.predict(x = x, weights = model_mlp)

library(ROCR)
perf = performance( prediction( pred, y ), "prec","rec" )
plot(perf, colorize=T, xlim=c(1,0.9))
auc = performance( prediction( pred, y ), "auc" )@y.values[[1]]


# analyses #

coeffs_from_2_pts <- function(x1, y1, x2, y2){
  if(y1==y2) return(list(a=0, b=y1))
  if(x1==x2) return(list(a=NULL, b=NULL))
  a = (y2 - y1) / (x2 - x1)
  b = y1 - a * x1
  return(list(a=a, b=b))
}

sens_thr = 0.98

cf = coeffs_from_2_pts(x1 = unlist(perf@x.values)[rev(which(unlist(perf@x.values)<=sens_thr))[1]],
                       y1 = unlist(perf@y.values)[rev(which(unlist(perf@x.values)<=sens_thr))[1]],
                       x2 = unlist(perf@x.values)[which(unlist(perf@x.values)>=sens_thr)[1]],
                       y2 = unlist(perf@y.values)[which(unlist(perf@x.values)>=sens_thr)[1]])

fdr = cf[["a"]] * sens_thr + cf[["b"]]

points(sens_thr, fdr)

mlp_thr = 
  
  # compute FDR given sens #
  
  tdr_thr = 0.98
cf = coeffs_from_2_pts(x1 = unlist(perf@x.values)[rev(which(rev(unlist(perf@y.values))<=tdr_thr))[1]],
                       y1 = rev(unlist(perf@y.values))[rev(which(rev(unlist(perf@y.values))<=tdr_thr))[1]],
                       x2 = unlist(perf@x.values)[which(rev(unlist(perf@y.values))>=tdr_thr)[1]],
                       y2 = rev(unlist(perf@y.values))[which(rev(unlist(perf@y.values))>=tdr_thr)[1]])

sens = cf[["a"]] * tdr_thr + cf[["b"]]
if(length(sens)==0) sens = unlist(perf@x.values)[rev(which(rev(unlist(perf@y.values))<=tdr_thr))[1]]

points(sens_thr, fdr)


#### using K-fold ####
#### run the beginning of MLP.r script to have the table status ####

folds <- createFolds(train_table$status, 10)

#all_pred = c()
#all_status = c()

all_pred_rf = c()
all_status_rf = c()

for(i in 1:10){
  print(paste("fold: ",i,sep=""))
  test = train_table[folds[[i]],]
  train = train_table[-folds[[i]],]
  
  #xtest = data.matrix(test[,my_features[my_features!="status"]])
  #ytest = test$status
  #ytest[which(ytest=="TP")] = 1
  #ytest[which(ytest=="FP")] = 0
  #ytest=as.matrix(as.numeric(ytest))
  
  #xtrain = data.matrix(train[,my_features[my_features!="status"]])
  #ytrain = train$status
  #ytrain[which(ytrain=="TP")] = 1
  #ytrain[which(ytrain=="FP")] = 0
  #ytrain=as.matrix(as.numeric(ytrain))
  
  #model_mlp <- monmlp.fit(xtrain, ytrain, hidden1=3, n.ensemble=1, monotone=1, bag=TRUE)
  
  #test$prediction = monmlp.predict(x = xtest, weights = model_mlp)
  #all_pred=c(all_pred, test$prediction)
  #all_status = c(all_status, ytest)
  
  train$old_status=train$status; test$old_status=test$status
  train$status = 0 ; train[which(train$old_status=="TP"),"status"] = 1
  test$status = 0 ; test[which(test$old_status=="TP"),"status"] = 1
  rf_fold = randomForest(as.factor(status) ~ ., data = train[,my_features], importance = TRUE,
                         ntree = 500, sampsize = as.numeric(table(train$status)["1"])) 
  prediction = predict(rf_fold, test, type="prob")[,2]
  all_pred_rf = c(all_pred_rf, prediction)
  all_status_rf = c(all_status_rf, test$status)
}

library(ROCR)
#perf = performance( prediction( all_pred, all_status ), "prec","rec" )
#plot(perf, colorize=T, xlim=c(1,0), lwd=3)
#auc = performance( prediction( all_pred, all_status ), "auc" )@y.values[[1]]

perf = performance( prediction( all_pred_rf, all_status_rf ), "rec" ,"ppv" ) #sens in y and tdr in x
if(type=="snv"){
  plot(perf, colorize=T, lwd=3, xlab="1-FDR", ylab="sensitivity", xaxt='n', main=type, 
       xlim=c(1,0.5), ylim=c(0.5,1))
} else {
  plot(perf, colorize=T, xlim=c(1,0), lwd=3, xlab="1-FDR", ylab="sensitivity", xaxt='n', main=type, colorkey=F)
}

auc = performance( prediction( all_pred_rf, all_status_rf ), "auc" )@y.values[[1]]

# compare with home filters
train_table$prediction_home_filters = "FP"
rvsbthr = ifelse(type=="snv", 0.95, 0.85)
train_table[which(train_table$AF>0.1 & train_table$RVSB<rvsbthr & train_table$DP>50 & train_table$ERR<(-2)),"prediction_home_filters"] = "TP"
sens_home_filter = sum(train_table$status=="TP" & train_table$prediction_home_filters=="TP") / sum(train_table$status=="TP")
TDR_home_filters = sum(train_table$status == "TP" & train_table$prediction_home_filters == "TP") / sum(train_table$prediction_home_filters == "TP")

points(TDR_home_filters, sens_home_filter, pch=4, lwd=2)
legend(x=0.6, y=0.7, c("AF>0.1","RVSB<0.95","DP>50","ERR<0.01"), 
       col=c("black","white","white","white"), pch=4, bty="n", cex=0.75)

id = which(unlist(perf@x.values) <= TDR_home_filters)[1]
points(unlist(perf@x.values)[id], unlist(perf@y.values)[id])



