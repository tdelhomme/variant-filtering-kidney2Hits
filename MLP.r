library(monmlp)

x = data.matrix(train_table[,my_features[my_features!="status"]])

y = train_table$status
y[which(y=="TP")] = 1
y[which(y=="FP")] = 0
y=as.matrix(as.numeric(y))

model_mlp <- monmlp.fit(x, y, hidden1=3, n.ensemble=15, monotone=1, bag=TRUE)

pred = monmlp.predict(x = x, weights = model_mlp)

library(ROCR)
plot( performance( prediction( pred, y ), "prec","rec" ), colorize=T, xlim=c(1,0.9) )
performance( prediction( pred, y ), "auc" )@y.values[[1]]
