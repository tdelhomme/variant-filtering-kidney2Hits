library(monmlp)

x = data.matrix(train_table[,my_features[my_features!="status"]])

y = train_table$status
y[which(y=="TP")] = 1
y[which(y=="FP")] = 0
y=as.matrix(as.numeric(y))

model_mlp <- monmlp.fit(x, y, hidden1=3, n.ensemble=15, monotone=1, bag=TRUE)

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


# compute FDR given sens #

tdr_thr = 0.98
cf = coeffs_from_2_pts(x1 = unlist(perf@x.values)[rev(which(rev(unlist(perf@y.values))<=tdr_thr))[1]],
                       y1 = rev(unlist(perf@y.values))[rev(which(rev(unlist(perf@y.values))<=tdr_thr))[1]],
                       x2 = unlist(perf@x.values)[which(rev(unlist(perf@y.values))>=tdr_thr)[1]],
                       y2 = rev(unlist(perf@y.values))[which(rev(unlist(perf@y.values))>=tdr_thr)[1]])

sens = cf[["a"]] * tdr_thr + cf[["b"]]
if(length(sens)==0) sens = unlist(perf@x.values)[rev(which(rev(unlist(perf@y.values))<=tdr_thr))[1]]

points(sens_thr, fdr)