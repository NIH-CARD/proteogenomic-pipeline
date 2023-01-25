#read gtf files
#read in required packages
require(readxl)
require(tidyverse)


if (file.exists("principal_tx.csv")) {
  #Delete file if it exists
  file.remove("principal_tx.csv")
}

#create a list of the files from your target directory
file_list <- list.files(path="iPSC_gtfs",pattern = "\\.csv$")
ddf <- data.frame(TxID=character(),GeneID=character(),Gene_Name=character(), cov=numeric(),FPKM=numeric(),TPM=numeric())
all_ddf <- data.frame(TxID=character(),GeneID=character(),Gene_Name=character(), cov=numeric(),FPKM=numeric(),TPM=numeric())
for (i in 1:length(file_list))
{
  print(paste0('reading file: ',file_list[i]))
  rec <- as.data.frame(read.csv(paste0("iPSC_gtfs/",file_list[i]),sep="\t", header=FALSE))
  rec=subset (rec, select = -V7)
  colnames(rec) = c("TxID","GeneID","Gene_Name", "cov","FPKM","TPM")

  all_ddf=rbind(all_ddf,rec)
}
#remove spaces
all_ddf <- as.data.frame(apply(all_ddf, 2, str_remove_all, " "))
#and sort for fast retrieval
all_ddf<-all_ddf[order(all_ddf$Gene_Name),]
#get unique gene names
unique_genes=unique(all_ddf[c("Gene_Name")])
#
Tx_ddf <- data.frame(TxID=character(),GeneID=character(),Gene_Name=character())
#select row with max cov for each gene
#for (i in 1:dim(all_ddf)[1])
print(paste0('Now Generating Txs Table, Will take a while !!!!!!!'))
for (i in 1:length(unique_genes[[1]]))
{
  TxSubset <- all_ddf[grep(unique_genes[[1]][i], all_ddf$Gene_Name), ]
  tx_gene=TxSubset[which.max(TxSubset$cov),][1:3]
  Tx_ddf=rbind(Tx_ddf,tx_gene)
  #if((i %% 5000)==0)
  #{print(paste0("reached:",i))}
}
#write in file
#write.table(Tx_ddf,file=paste0("principal_tx",".csv"), quote=F, sep=",", row.names=F, col.names=F)
#also remove trailing version numbers
Tx_ddf$GeneID = unlist(lapply(Tx_ddf$GeneID, function (x) strsplit(as.character(x), ".", fixed=TRUE)[[1]][1]))
#write.table(Tx_ddf,file=paste0("principal_tx1",".csv"), quote=F, sep=",", row.names=F, col.names=F)

Tx_ddf$TxID = unlist(lapply(Tx_ddf$TxID, function (x) strsplit(as.character(x), ".", fixed=TRUE)[[1]][1]))
write.table(Tx_ddf,file=paste0("principal_txs",".csv"), quote=F, sep=",", row.names=F, col.names=F)

print(paste0('Done With Txs Table: principal_txs.csv'))