library(monmlp)

x = as.matrix(train_table[,c("RVSB","QVAL")])

y = train_table$status
y[which(y=="TP")] = 1
y[which(y=="FP")] = 0
y=as.matrix(as.numeric(y))

model_mlp <- monmlp.fit(x, y, hidden1=3, n.ensemble=15, monotone=1, bag=TRUE)

pred = monmlp.predict(x = x, weights = model_mlp)

library(ROCR)
plot( performance( prediction( pred, y ), "tpr","fpr" ),colorize=T )
performance( prediction( pred, y ), "auc" )@y.values[[1]]
