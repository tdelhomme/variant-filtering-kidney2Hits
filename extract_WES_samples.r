setwd("~/Documents/Models/variant-filtering/")
library(xlsx)

data = read.xlsx("K2H_SamplesInfo_03Apr2018.xlsx", sheetIndex = 1, startRow = 2)

sample_to_exclude = c("K2150075","K2110533","K2110519")

WES_samples = paste(as.character(data[which(data$Sent_to_CCG==1),"Barcode"]),
                    as.character(data[which(data$Sent_to_CCG==1),"ID"]),sep="-")

WES_samples_name = unlist(lapply(WES_samples, function(s) unlist(strsplit(s,"-"))[2]))

WES_samples = WES_samples[which(! WES_samples_name %in% sample_to_exclude)]

length(WES_samples) # 55 + 1 for the QC duplicate

write.table(WES_samples, file="K2H_WES_sample_ID.txt", quote = F, row.names = F, col.names = F)
