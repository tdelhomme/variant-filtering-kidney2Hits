library(monmlp)
library(caret)

setwd("~/Documents/Models/variant-filtering")
set.seed(123) # here I fix the seed of the random generator to have the same random numbers if re-run 

###### functions #####
data_summary <- function(x) {
  m <- median(x)
  ymin <- as.numeric(quantile(x)[2])
  ymax <- as.numeric(quantile(x)[4])
  return(c(y=m,ymin=ymin,ymax=ymax))
}

plot_ggviolin <- function(kfold_TDR, kfold_sens, type){
  dat = data.frame(var=c(rep("TDR", length(kfold_TDR)),rep("sensitivity", length(kfold_sens))),
                   val=c(kfold_TDR, kfold_sens))
  ggplot(dat, aes(x=var, y=val, fill=var)) + theme_bw() + theme(legend.position="none",
                                                                plot.title = element_text(size=18,hjust = 0.5),
                                                                axis.ticks.x = element_blank(),
                                                                panel.border=element_blank(),
                                                                axis.text=element_text(size=12),
                                                                axis.title=element_text(size=14)) + ggtitle(type) +
    geom_violin() + stat_summary(fun.data=data_summary, color="black") +
    scale_fill_manual(values=c("darkorchid1","chocolate1")) + labs(x="", y="") +
    scale_y_continuous(limits = c(0.5, 1))
}  
##############################

for(type in c("snv", "indel")){
  train_table="tables/K2H_AllVariants_All_Lib_NoMinAF_noequal_GOF_addVCFfeatures_illuminaBED_WES_samples_annotated_with_coverage_INFO_GENO_status_supp_features.txt"
  train_table = read.table(train_table, quote="\"", stringsAsFactors=F, sep="\t", header=T)
  train_table$IoD = 1 + train_table$SIG * 10^(train_table$ERR)
  if(type=="snv") train_table = train_table[which(train_table$TYPE_INFO=="snv"),]
  if(type=="indel") train_table = train_table[which(train_table$TYPE_INFO=="ins" | train_table$TYPE_INFO=="del"),]
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
  
  #boxplot(data.frame("TDR"=kfold_TDR, "sensitivity"=kfold_sens), col="lightgrey", ylim=c(0.5,1), outpch=20, outcex=0.75, axes=F, staplelwd=2) 
  #title(type, line = 2.5)
  #axis(side=3, at=c(1,2), labels=c("TDR","sensitivity"), col = NA, col.ticks = NA)
  #axis(side=2, at=seq(0,1,by=0.1))
  plot_ggviolin(kfold_TDR, kfold_sens, type)
}



cols=c(rev(rainbow(400,start=0/6, end=4/6))[1:125],
       rep(rev(rainbow(400,start=0/6, end=4/6))[200],1000),
       rev(rainbow(400,start=0/6, end=4/6))[250:400])

cols=rev(rainbow(256,start=0/6, end=4/6))
  
plot(perf, colorize=T, xlim=c(1,0), lwd=3,  colorize.palette=cols)
       
       
