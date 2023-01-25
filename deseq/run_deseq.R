#!/usr/bin/env Rscript


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
# generate dds object by loading count matrix geenrated from stringtie2 (after merging and re-estimation)
countData <- read.csv("gene_count_matrix.csv")
#now load metadata file
metaData <- read.csv("new_samples.csv", sep="\t", row.names=1)
# run dds functions

dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~condition, tidy = TRUE)

dds <- DESeq(dds)

# if want to save the normalized count file
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, file="normalized_counts.csv")

#View(normalized_counts)

# run differential gene expression - here condition column in the metaData contains two conditions, Control and TDP43
#note: for all RBP samples, add all samples in col1 of new_samples.csv and condition names in col2
#res <- results(dds, contrast = c("condition", "TDP43", "Control"))

#write.csv(res, "TDPVsControl_test.csv")

#vextract results
conditions <- unique(metaData$condition)

#repeat for multiple conditions

for (i in 1:length(conditions)){
	if (conditions[i] != "ctrl"){
	col = conditions[i]
	print(paste("Processing ", col))
	res <- results(dds, contrast=c("condition", col, "ctrl"))
	write.csv(x=res, file=paste("control_vs_", col, ".csv", sep=""))
	# dfList[[i]] <- res
	} else {
	NULL
	}
}

