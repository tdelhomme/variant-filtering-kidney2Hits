library(monmlp)

setwd("~/Documents/Models/variant-filtering")
set.seed(123) # here I fix the seed of the random generator to have the same random numbers if re-run 

train_table="tables/K2H_AllVariants_All_Lib_NoMinAF_noequal_GOF_addVCFfeatures_illuminaBED_WES_samples_annotated_with_coverage_INFO_GENO_status_supp_features.txt"
train_table = read.table(train_table, quote="\"", stringsAsFactors=F, sep="\t", header=T)

train_table$IoD = 1 + train_table$SIG * 10^(train_table$ERR)
train_table = train_table[which(train_table$TYPE_INFO=="snv"),]
#train_table = train_table[which(train_table$TYPE_INFO=="ins" | train_table$TYPE_INFO=="del"),]
train_table[which(is.infinite(train_table$MIN_DIST)),"MIN_DIST"] = 1000000000
train_table[which(is.infinite(train_table$MaxRatioWin)),"MaxRatioWin"] = 1000000000

propTP = as.numeric(table(train_table$status)["TP"] / nrow(train_table))
propFP = as.numeric(table(train_table$status)["FP"] / nrow(train_table))

my_features=c("status","RVSB", "QVAL","AF","ERR_INFO","DP", "medianDP_INFO",
              "FS", "MIN_DIST", "AO", "QUAL", "MaxRatioWin", "NbVarWin", "IoD", "HpLength")


# use K-fold
kfold_spec = c()
kfold_sens = c()
kfold_TDR = c()

folds <- createFolds(train_table$status, 10)

for(i in 1:10){
  print(paste("fold: ",i,sep=""))
  test = train_table[folds[[i]],]
  train = train_table[-folds[[i]],]
  
  xtest = data.matrix(test[,my_features[my_features!="status"]])
  ytest = test$status
  ytest[which(ytest=="TP")] = 1
  ytest[which(ytest=="FP")] = 0
  ytest=as.matrix(as.numeric(ytest))

  xtrain = data.matrix(train[,my_features[my_features!="status"]])
  ytrain = train$status
  ytrain[which(ytrain=="TP")] = 1
  ytrain[which(ytrain=="FP")] = 0
  ytrain=as.matrix(as.numeric(ytrain))
  
  model_mlp <- monmlp.fit(xtrain, ytrain, hidden1=3, n.ensemble=1, monotone=1, bag=TRUE)
  
  test$prediction = monmlp.predict(x = xtest, weights = model_mlp)
  thr = 0.5
  kfold_sens = c(kfold_sens, (sum(test$status == "TP" & test$prediction >= thr) / sum(test$status == "TP")) )
  kfold_spec = c(kfold_spec, (sum(test$status == "FP" & test$prediction <= thr) / sum(test$status == "FP")) )
  kfold_TDR = c(kfold_TDR, (sum(test$status == "TP" & test$prediction >= thr) / sum(test$prediction >= thr)) )
}

boxplot(data.frame("TDR"=kfold_TDR, "sensitivity"=kfold_sens), col="lightgrey", ylim=c(0.85,1), outpch=20, outcex=0.75, axes=F, staplelwd=2)
axis(side=3, at=c(1,2), labels=c("TDR","sensitivity"), col = NA, col.ticks = NA)
axis(side=2, at=c("0.85","0.90","0.95","1"))



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