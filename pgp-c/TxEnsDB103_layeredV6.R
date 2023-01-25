###!/usr/bin/env Rscript

#using ENSG (commented out this - junction start and end range) to get transcripts
#Starting from previous version (TxEnsDB103_layeredV1.R), it has 41 events for which Tx selected does not encapsulate the event
#An example is event: chr12	22515151	22517985	AC053513.1
#So coordinate based search is back on to see if it makes difference
#THIS SCRIPT HAS BEEN MODIFIED - SO REPLACE IT IN ALL VERSIONS
#04/20/2022 - removing 5'utr
rm(list = ls())
library(GenomicRanges)
#USE EnsDb V86
#library(EnsDb.Hsapiens.v86)
#USE EnsDB V99 from Annotation hub
library(AnnotationHub)
library(GenomicFeatures)
library(Repitools)
library(ensembldb)
ah <- AnnotationHub()
#also load gtf file fron V86
#edb <- EnsDb.Hsapiens.v86
#get object of EnsDBV99
#edb=query(ah, c("EnsDb", "Hsapiens", "99"))[[1]]
edb=query(ah, c("EnsDb", "Hsapiens", "103"))[[1]]

#Also get Tx lengths to remove 5'utr
tx_lens=transcriptLengths(edb,with.utr5_len = TRUE,with.utr3_len = TRUE)
#read bed file from MAJIQ



args = commandArgs(trailingOnly=TRUE)



GeneIDField=6

#{
  #set current working directory
  #setwd("/Volumes/SYEDSHAH/MichaelLab/proteogenomic_pipeline/github_final_pgp/phase-b/ce_incl/v3")
  
  #Read Peaks File
  
  SpliceData <- read.csv(args[1], header=FALSE) #this gave poblem and I do not know why
  #Also read Tx list
  Tx_list <- read.csv(args[2], header=FALSE)
  
  
  #SpliceData <- read.csv("sorted_remaining_events.csv", header=FALSE)
  #Tx_list <- unique(read.csv("principal_txs.csv", header=FALSE)) #get unique data
  
  #also read appris annoation data to get principal 1 isoform
  appr_anno1 <- as.data.frame(read.csv("GRCh38_appris_data.principal.txt",sep="\t", header=FALSE))
  #V1 - hugo_symbol, V2 - ENSG_ID, V3 - TX_ID, V4 - , V2 - PRINCIPAL/ALTERNATE
  colnames(appr_anno1)=c('V1','V2','V3','V4','V5')
  #select only PRINCIPAL.1 ISOFORMS
  appr_anno=appr_anno1[with(appr_anno1,grepl("PRINCIPAL:1", V5)),] #APPRIS has multiple P1 Tx for a gene- e.g RALYL
  
  if (file.exists("all_tx_events.csv")) {
    #Delete file if it exists
    file.remove("all_tx_events.csv")
  }
  
  if (file.exists("all_events_bed_sashimi.tab")) {
    #Delete file if it exists
    file.remove("all_events_bed_sashimi.tab")
  }

  if (file.exists("events_to_tx_mapping_valid.csv")) {
    #Delete file if it exists
    file.remove("events_to_tx_mapping_valid.csv")
  }
  if (file.exists("events_to_tx_mapping_invalid.csv")) {
    #Delete file if it exists
    file.remove("events_to_tx_mapping_invalid.csv")
  }
  
  
#}


write(paste0('                                 '),file=args[3],append=TRUE)
write(paste0('Starting From TxEnsDB103_layeredV6.R: --------------- Processing file: ',args[1],' with: ',dim(SpliceData)[1],' events to generate each event .bed files in event_bedfiles/ folder: '),file=args[3],append=TRUE)

print(paste0('Started Generating BED files for Splicing Events in folder event_bedfiles/ from File: ',args[1]))

trackj = 1;
temp_gene="";
current_gene="";
tx_lengths <- c()
df_notfound <- data.frame(seqnames=character(), start=numeric(),end=numeric(),strand=character(),genename=character(),junc_type=character())
df_zeroutr <- data.frame(seqnames=character(), start=numeric(),end=numeric(),strand=character(),genename=character(),junc_type=character())
repeated_entries=0
iPSC_events=0
appris_events=0
principalTx_events=0
events_xTx=0
Tx_str=0 #0 for iPSC, 1 for APPRIS and 2 for EnsDB
Tx_valid=0
Total_Events = dim(SpliceData)[1]
probable_noise_events=0
probable_noncoding_events=0
utr5_events=0
for (i in 1:dim(SpliceData)[1])
{
  #i=1
  #print(paste('processing gene: ',SpliceData[i,5]))
  #print(paste('processing gene: ',SpliceData[i,6]))
  #Step 0 - get transcripts for each gene (MAJIQ only reports for some)
  #Deal with multiple events for the same gene
  rm(Tx_name) #clear Tx_name
  if (i==1)
  {
    temp_gene=SpliceData[i,5]
    trackj = 1
  }
  else if (temp_gene==SpliceData[i,5])
  {
    trackj = trackj+1
    
  }
  else
  {
    temp_gene=SpliceData[i,5]
    trackj=1
  }
    #get gene name and gene_id (from granges filter) using event coordinates
    grf <- GRangesFilter(GRanges(substr(SpliceData[i,1],4,nchar(as.character(SpliceData[i,1]))), ranges = IRanges(SpliceData[i,2], SpliceData[i,3]),strand = SpliceData[i,4]), type = "any")
    gn <- genes(edb, filter = grf)
    #First check if gene_name exactly matches upto gene_version and then get all its Txs
    genes_data=annoGR2DF(gn)
    #Txs <- transcripts(edb, filter = GeneNameFilter(gn$gene_name))
    #if(length(gn)>1)
    #  events_xTx=events_xTx+1
    #print(paste0('gene_names',gn))
    gn_id=strsplit(SpliceData[i,6],"[.]")[[1]][1]
    #print(paste0('gene ID: ',gn_id))
    #reset Tx flag
    Tx_flg=0
    #print(paste('num of gene_ids:',dim(genes_data)[1]))
    #if(length(gn_id)>0)
    if(dim(genes_data)[1]>0)
    {
      #First try iPSC Tx
      #tmp_tx=Tx_list[grep(gn_id, Tx_list$V2), 1]
      #tmp_tx=Tx_list[grep(genes_data$gene_id, Tx_list$V2), 1]
      #if(length(strsplit(Tx_list[grep(gn$gene_id, Tx_list$V2), 1],"[.]")[[1]][1])>0)
      if(dim(genes_data)[1]>1)
      {
        #some events (coordinates) are mapped to multiple gene_id's, so select one that has same gene_id as majiq and spans
        #the genomic range or the one which spans the genomic range
        flg_ex=0
        for(ti in 1:dim(genes_data)[1])
        {
          if(SpliceData[i,2]>=genes_data[ti,]$start && SpliceData[i,3]<=genes_data[ti,]$end && length(Tx_list[grep(genes_data[ti,]$gene_id, Tx_list$V2), 1])>0)
          {
            Tx_name = Tx_list[grep(genes_data[ti,]$gene_id, Tx_list$V2), 1]
            flg_ex=1
            break
          }
          if(flg_ex==1)
            break
        }
      }
      else 
      {
        if(SpliceData[i,2]>=genes_data$start && SpliceData[i,3]<=genes_data$end && length(Tx_list[grep(genes_data$gene_id, Tx_list$V2), 1])>0)
        {
          Tx_name = Tx_list[grep(genes_data$gene_id, Tx_list$V2), 1]
        }
      }
        
      if(exists("Tx_name") && length(Tx_name)>0)
      {
        #Tx_name = Tx_list[grep(gn$gene_id, Tx_list$V2), 1]
        #Tx_name = strsplit(Tx_list[grep(gn_id, Tx_list$V2), 1],"[.]")[[1]][1]
        tl1=transcriptLengths(edb, filter = TxIdFilter(Tx_name))
        #make sure that we have transcript in EnsDB 
        if (dim(tl1)[1]>0) 
        {
          #and event coordinates lies within the Tx
          dat=data.frame(exonsBy(edb, by = "tx", filter = TxIdFilter(Tx_name))@unlistData)
          if (SpliceData[i,2]>=min(dat[,2]) && SpliceData[i,3]<=max(dat[,3])) #make sure that event lies within the transcript
          {
            iPSC_events=iPSC_events+1
            Tx_str = 'iPSC'
            Tx_flg=1
          }
        }
      

      }
      
      #if ((length(tmp_tx)==0 || length(Tx_name) == 0 || dim(tl1)[1]==0) && length(appr_anno[which(gn_id==appr_anno["V2"]),"V3"])>0) #now APPRIS  
      #if ((Tx_flg==0) && length(appr_anno[which(gn_id==appr_anno["V2"]),"V3"])>0) #now APPRIS    
      if ((Tx_flg==0) && length(appr_anno[which(genes_data$gene_id==appr_anno["V2"]),"V3"])>0) #now APPRIS    
      #else if(length(appr_anno[which(gn$gene_id==appr_anno["V2"]),"V3"])>0) #now APPRIS
      {
        #Tx_name=appr_anno[which(gn$gene_id==appr_anno["V2"]),"V3"]
        #Tx_name=appr_anno[which(gn_id==appr_anno["V2"]),"V3"]
        Tx_name=appr_anno[which(genes_data$gene_id==appr_anno["V2"]),"V3"]
        tltt=transcriptLengths(edb, filter = TxIdFilter(Tx_name))
        #As APPRIS has multiple P1 Txs for some genes (like RALYL), make sure that the one with maximum exons and highest Tx length (in bp) is selected
        
        if (dim(tltt)[1]>0)
        {
          #if multiple P1 Txs, then select the one with highet num of exons and largest size (in bp) and encapsulates event
          tflg=0
          for (ti in 1:dim(tltt)[1])
          {
            dat=data.frame(exonsBy(edb, by = "tx", filter = TxIdFilter(tltt[ti,]$tx_id[1]))@unlistData)
            if (SpliceData[i,2]>=min(dat[,2]) && SpliceData[i,3]<=max(dat[,3])) #make sure that event lies within the transcript
            {
              
              if(tflg==0)
              {
                tlt=tltt[ti,]
              }
              else
              {
                tlt=rbind(tlt,tltt[ti,])
              }
              tflg=1
            }
          }
          ########
          #now make sure that longest Tx is selected
          if (tflg==1)
          {
          if(dim(tlt)[1]>0)
          {
            jj1=which(tlt$nexon==max(tlt$nexon)) #for nexon
            tl1t=tlt[jj1,]
            jj=which(tl1t$tx_len==max(tl1t$tx_len)) #which(tl$nexon==max(tl$nexon)&&tl$tx_len==max(tl$tx_len)) #for nexon
            tl1=tl1t[jj,]
            Tx_name=tl1[1,2]
            if (dim(tl1)[1]>0)
            {
              
          ########
                #and event coordinates lies within the Tx
                dat=data.frame(exonsBy(edb, by = "tx", filter = TxIdFilter(Tx_name))@unlistData)
                if (SpliceData[i,2]>=min(dat[,2]) && SpliceData[i,3]<=max(dat[,3])) #make sure that event lies within the transcript
                {
                  appris_events=appris_events+1
                  Tx_str = 'APPRIS'
                  Tx_flg=1
                }
            }
          }
        }
        }

        
      }
      #if ((length(tmp_tx)==0 || length(Tx_name) == 0 || dim(tl1)[1]==0) && dim(transcriptLengths(edb, filter = GeneNameFilter(gn$gene_name)))[1]>0) #and Finally longest Tx
      if ((Tx_flg==0) && dim(transcriptLengths(edb, filter = GeneNameFilter(gn$gene_name)))[1]>0) #and Finally longest Tx
        #else if (dim(transcriptLengths(edb, filter = GeneNameFilter(gn$gene_name)))[1]>0) #and Finally longest Tx
      {
        if(dim(transcriptLengths(edb, filter = GeneNameFilter(gn$gene_name)))[1]>1)
          events_xTx=events_xTx+1
        
        tltt=transcriptLengths(edb, filter = GeneNameFilter(gn$gene_name))
        #first get list of Tx's which encapsulates the event
        tflg=0
        for (ti in 1:dim(tltt)[1])
        {
          #make sure that transcipt starting with ENST is used - was case with gene FH, where one of the Tx_id is LRG_504t1
          if(grepl('ENST', tltt[ti,]$tx_id[1], fixed = TRUE))
          {
            dat=data.frame(exonsBy(edb, by = "tx", filter = TxIdFilter(tltt[ti,]$tx_id[1]))@unlistData)
            if (SpliceData[i,2]>=min(dat[,2]) && SpliceData[i,3]<=max(dat[,3])) #make sure that event lies within the transcript
            {
              if(tflg==0)
              {
                tlt=tltt[ti,]
              }
              else
              {
                tlt=rbind(tlt,tltt[ti,])
              }
              tflg=1
            }
          }
        }
        #now make sure that longest Tx is selected
        if(dim(tlt)[1]>0)
        {
          jj1=which(tlt$nexon==max(tlt$nexon)) #for nexon
          tl1t=tlt[jj1,]
          jj=which(tl1t$tx_len==max(tl1t$tx_len)) #which(tl$nexon==max(tl$nexon)&&tl$tx_len==max(tl$tx_len)) #for nexon
          tl1=tl1t[jj,]
          Tx_name=tl1[1,2]
          if (dim(tl1)[1]>0)
          {
            #and event coordinates lies within the Tx
            dat=data.frame(exonsBy(edb, by = "tx", filter = TxIdFilter(Tx_name))@unlistData)
            if (SpliceData[i,2]>=min(dat[,2]) && SpliceData[i,3]<=max(dat[,3])) #make sure that event lies within the transcript
            {
              principalTx_events=principalTx_events+1
              Tx_str = 'EnsDB'
              Tx_flg=1
            }
          }
        }
      }
      
      
          ################3
      #Step 1 - for each transcript, get all exons for it
      #ddf <- data.frame(seqnames=str(''), start=numeric(),end=numeric(),width=numeric(),exon_rank=numeric(),strand=str(''))
      ddf <- data.frame(seqnames=character(), start=numeric(),end=numeric(),width=numeric(),exon_rank=numeric(),strand=character())
      #includes Tx id to use with sashimi plots
      #ddf1 <- data.frame(seqnames=str(''), start=numeric(),end=numeric(),width=numeric(),exon_rank=numeric(),strand=str(''),TxID=str(''))
      ddf1 <- data.frame(seqnames=character(), start=numeric(),end=numeric(),width=numeric(),exon_rank=numeric(),strand=character(),TxID=character())
      #check if Tx is available from either of the three resources
      if ((Tx_flg==1) && dim(tl1)[1]>0 && length(Tx_name)>0 && length(exonsBy(edb, by = "tx", filter = TxIdFilter(Tx_name)))>0)
      #if(length(Tx_name)>0 && length(exonsBy(edb, by = "tx", filter = TxIdFilter(Tx_name)))>0)
      {
            #Now use Tx to 
            #EnsGenes1 <- exonsBy(edb, by = "tx", filter = TxNameFilter(Tx_name))
            EnsGenes1 <- exonsBy(edb, by = "tx", filter = TxIdFilter(Tx_name))
            EnsGenes_cds = exonsBy(edb, by = "tx", filter = TxIdFilter(Tx_name),columns=c("tx_seq_start","tx_cds_seq_start","tx_seq_end","tx_cds_seq_end","tx_biotype"))
            #get 5utr length
            utr5l=tx_lens[Tx_name,]$utr5_len
            utr3l=tx_lens[Tx_name,]$utr3_len
            #get exons
            #EnsGenes1 <- exonsBy(edb, by = "tx", filter = TxIdFilter(tl[j,2]))
            #now get data frame and subtract 5utr and 3utr
            g_dat=data.frame(EnsGenes1@unlistData)
            g_dat1=data.frame(EnsGenes_cds@unlistData)
            #print(paste0('processing gene: ',temp_gene, ' @ line: ',i))
            if(g_dat1$tx_biotype[1]=='protein_coding')
            {
              if (g_dat$strand[1]=="+") 
              {
                #First check if event lies within the cds
                if(SpliceData[i,]$V2 >= g_dat1$tx_cds_seq_start[1])
                {
                   if (g_dat$width[1] > utr5l) #AARS1 has trl5l > exon 1 width
                   {
                      g_dat$start[1] = g_dat$start[1]+utr5l-1
                      g_dat$width[1] = g_dat$width[1]-utr5l
                      #utr5_events=utr5_events+1
                      #write(paste0('5 UTR events#: ',utr5_events,' ** gene: ', SpliceData[i,]$V5, ' @ line: ',i,' lies outside of 5 UTR'),file=args[3],append=TRUE)
                      
                   }
                  #else
                  #{
                  #  print(paste0('=== gene:', SpliceData[i,]$V5, ' @ line: ',i,' has 5 exon width: ',g_dat[1,]$width, ' and 5 utr length: ' ,utr5l))
                  #}
                }
                else
                {
                  probable_noncoding_events=probable_noncoding_events+1
                  print(paste0('*** gene: ', SpliceData[i,]$V5, ' @ line: ',i,' lies outside of 5 UTR'))
                  write(paste0('Probable non-cds events#: ',probable_noncoding_events,' ** gene: ', SpliceData[i,]$V5, ' @ line: ',i,' lies outside of 5 UTR'),file=args[3],append=TRUE)
                }
              }
              else
              {
                #First check if event lies within the cds
                if(SpliceData[i,]$V3 <= g_dat1$tx_cds_seq_end[1]) #for -ve strand, ensdb return tx_cds_end as 5' utr, why??
                {
                  
                  if (g_dat$width[1] > utr5l)
                  {
                    g_dat$end[1] = g_dat$end[1]-utr5l
                    g_dat$width[1]=g_dat$width[1] - utr5l
                    #utr5_events=utr5_events+1
                    #write(paste0('5 UTR events#: ',utr5_events,' ** gene: ', SpliceData[i,]$V5, ' @ line: ',i,' lies outside of 5 UTR'),file=args[3],append=TRUE)
                    
                  }
                  #else
                  #{
                  #  print(paste0('=== gene:', SpliceData[i,]$V5, ' @ line: ',i,' has 5 exon width: ',g_dat[1,]$width, ' and 5 utr length: ' ,utr5l))
                  #}
                  
                }
                else
                {
                  probable_noncoding_events=probable_noncoding_events+1
                  print(paste0('** gene: ', SpliceData[i,]$V5, ' @ line: ',i,' lies outside of 5 UTR'))
                  write(paste0('Probable non-cds events#: ',probable_noncoding_events,' ** gene: ', SpliceData[i,]$V5, ' @ line: ',i,' lies outside of 5 UTR'),file=args[3],append=TRUE)
                }
                
              }
            }
            if(g_dat1$tx_biotype[1]!='protein_coding')
            {
              probable_noise_events=probable_noise_events+1
              print(paste0('Noise Event#: ',probable_noise_events,'--- Transcript: ',Tx_name,' gene: ', SpliceData[i,]$V5, ' @ line: ',i,' is: ',g_dat1$tx_biotype[1]))
              write(paste0('Noise Event#: ',probable_noise_events,'--- Transcript: ',Tx_name,' gene: ', SpliceData[i,]$V5, ' @ line: ',i,' is: ',g_dat1$tx_biotype[1]),file=args[3],append=TRUE)
            }
            #df=data.frame(EnsGenes1@unlistData)
            df=g_dat
            df$seqnames = paste0('chr', df$seqnames)
            #cols <- c(1:3,5)
            #include junction width as column 4 and exon rank as column 5
            
            cols <- c(1:4,8,5)
            ddf = rbind(ddf,df[,cols]) #concatenate all
            cols1 <- c(1:4,8,5,7)
            ddf1 = rbind(ddf1,df[,cols1]) #   rbind(ddf1,,df[,cols])
            if (NROW(ddf)>0)
            {
              #save bedfile for current junction
              SJdat <- data.frame(seqnames=SpliceData[i,1],
                                  start=SpliceData[i,2],
                                  end=SpliceData[i,3],
                                  name=1,
                                  score=0,
                                  strand=SpliceData[i,4])
      
              SJdat1 <- data.frame(seqnames=SpliceData[i,1],
                                  start=SpliceData[i,2],
                                  end=SpliceData[i,3],
                                  name=1,
                                  score=0,
                                  strand=SpliceData[i,4],
                                  strand=SpliceData[i,5],
                                  strand=Tx_name)
              if (SpliceData[i,2]>=min(ddf[,2]) && SpliceData[i,3]<=max(ddf[,3])) #make sure that event lies within the transcript
              {
                Tx_valid=Tx_valid+1
                SJdat2 <- data.frame(seqnames=SpliceData[i,1],
                                     start=SpliceData[i,2],
                                     end=SpliceData[i,3],
                                     gene_name=SpliceData[i,5],
                                     gene_id=SpliceData[i,6],
                                     Tx=Tx_name,
                                     Tx_start=min(ddf[,2]),
                                     Tx_end=max(ddf[,3]),
                                     Tx_type=Tx_str,
                                     valid=1)
              }
              else
              {
                SJdat2 <- data.frame(seqnames=SpliceData[i,1],
                                     start=SpliceData[i,2],
                                     end=SpliceData[i,3],
                                     gene_name=SpliceData[i,5],
                                     gene_id=SpliceData[i,6],
                                     Tx=Tx_name,
                                     Tx_start=min(ddf[,2]),
                                     Tx_end=max(ddf[,3]),
                                     Tx_type=Tx_str,
                                     valid=0)
                
              }
              
              write.table(SJdat, file=paste("event_bedfiles/temp_",SpliceData[i,5],"-",trackj,".bed",sep=""), quote=F, sep="\t", row.names=F, col.names=F)
              #also save all_events bed file with Tx selected for sashimi plots
              write.table(SJdat1, file='all_events_bed_sashimi.tab', append=TRUE, quote=F, sep="\t", row.names=F, col.names=F)
              write.table(ddf, file=paste("event_bedfiles/",SpliceData[i,5],"-",trackj,".bed",sep=""), quote=F, sep="\t", row.names=F, col.names=F)
              #TxID is used with sashimi plots
              #write.table(unique(ddf1), file=paste("event_bedfiles/TxID",SpliceData[i,5],"-",trackj,".bed",sep=""), quote=F, sep="\t", row.names=F, col.names=F)
              write.table(ddf1, file=paste("event_bedfiles/TxID",SpliceData[i,5],"-",trackj,".bed",sep=""), quote=F, sep="\t", row.names=F, col.names=F)
              
              #save for ce script
              write.table(SpliceData[i,],file='all_tx_events.csv',append=TRUE, quote=F, sep=",", row.names=F, col.names=F)
              #save for Tx mapping
              #if (SpliceData[i,2]>=min(ddf[,2]) && SpliceData[i,3]<=max(ddf[,3])) #make sure that event lies within the transcript
              {
                write.table(SJdat2,file='events_to_tx_mapping_valid.csv',append=TRUE, quote=F, sep=",", row.names=F, col.names=F)
              }
              #else
              #{
               # write.table(SJdat2,file='events_to_tx_mapping_invalid.csv',append=TRUE, quote=F, sep=",", row.names=F, col.names=F)
              #}
            }
            else
            {
              print(paste('got null at record: ',i,' for gene: ',SpliceData[i,5]))
              df_notfound = rbind(df_notfound,SpliceData[i,]) #concatenate all
            }
      }
      else
      {
        print(paste('got null at record: ',i,' for gene: ',SpliceData[i,5]))
        df_notfound = rbind(df_notfound,SpliceData[i,]) #concatenate all
        #save for ce script
        #write.table(SpliceData[i,],file='all_non_es.csv',append=TRUE, quote=F, sep=",", row.names=F, col.names=F)
      }
    }
    else
    {
      print(paste('got null at record: ',i,' for gene: ',SpliceData[i,5]))
      df_notfound = rbind(df_notfound,SpliceData[i,]) #concatenate all
      #save for ce script
      #write.table(SpliceData[i,],file='all_non_es.csv',append=TRUE, quote=F, sep=",", row.names=F, col.names=F)
    }
}
#finally write it
if(dim(df_notfound)[1]>0)
{
  write.table(df_notfound,file=paste0("EnsDB_tx_not_found",".csv"), quote=F, sep=",", row.names=F, col.names=F)
  #write.table(df_notfound,file=paste0("res_intronic_range/EnsDB_tx_not_found",".csv"), quote=F, sep=",", row.names=F, col.names=F)
}
#write.table(df_notfound,file=paste0("res_ce_boundary/still_all_not_found",".csv"), quote=F, sep=",", row.names=F, col.names=F)
#write.table(df_zeroutr,file=paste0("zero_utr",".csv"), quote=F, sep=",", row.names=F, col.names=F)
#print(paste('repeated entries are: ',repeated_entries))
print(paste('Out of a total of: ',Total_Events,' Events: Transcripts for: ',dim(df_notfound)[1],' Events are not_found are'))
print(paste('Total : ',iPSC_events, ' iPSC Princiapl Tx are selected'))
print(paste('Total : ',appris_events,' APPRIS Principal Tx are selected'))
print(paste('Total : ',principalTx_events, ' EnsDB Princiapl Tx (maximum number of Exons and largest size (in bp)) are selected'))
print(paste('Out of Total : ',Total_Events,' A total of: ',Tx_valid,' Events have Valid Transcripts (event lies between Transcript ends) and ',(Total_Events-Tx_valid),' Events have Transcripts that does not Encapsulate events'))


write(paste0('Finally total Noise Events detected are: ',probable_noise_events),file=args[3],append=TRUE)
write(paste0('Finally total Probable non-cds events#: ',probable_noncoding_events),file=args[3],append=TRUE)
#print(paste('genes with zero UTR are: ',dim(df_zeroutr)[1]))
#}
