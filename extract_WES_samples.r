setwd("~/Documents/Models/variant-filtering/")
library(xlsx)

data = read.xlsx("K2H_SamplesInfo_03Apr2018.xlsx", sheetIndex = 1, startRow = 2)

sample_to_exclude = c("K2150075","K2110533","K2110519")

WES_samples = paste(as.character(data[which(data$Sent_to_CCG==1),"Barcode"]),
                    as.character(data[which(data$Sent_to_CCG==1),"ID"]),sep="-")

Lib = as.character(data[which(data$Sent_to_CCG==1),"Library.ID"])

WES_samples_name = unlist(lapply(WES_samples, function(s) unlist(strsplit(s,"-"))[2]))

WES_samples = WES_samples[which(! WES_samples_name %in% sample_to_exclude)]
Lib = Lib[which(! WES_samples_name %in% sample_to_exclude)]
Lib = as.numeric(substr(Lib,8,8))

#length(WES_samples) # 55 + 1 for the QC duplicate

dat = data.frame("Lib"=Lib, "sample"=WES_samples)

for(lib in unique(Lib)){
  datlib = dat[which(dat$Lib == lib),"sample"]
  write.table(datlib, file=paste("K2H_WES_sample_ID_Lib",lib,".txt",sep=""), quote = F, row.names = F, col.names = F)
}
