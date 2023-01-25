library(data.table)
#library(readr)
#read bed file from MAJIQ
args = commandArgs(trailingOnly=TRUE)
library("stringr")                       # Load stringr package
#if (length(args)==0)
#{
#  stop("Please provide Path to ce_all_scan_range_junctions.bed and gene_name_TDP.cov.bed files, quitting", call.=FALSE)
#}else
#{
#  print(paste0("",args[1]))

#setwd("/Volumes/SYEDSHAH/MichaelLab/proteogenomic_pipeline/github_final_pgp/phase-a/new_bed_file_version")

#folder='res_ce_0.6'
#csv_record <- as.data.frame(read.csv(paste0(folder,"/IGV_ce_inclusion.csv"), sep="\t", header=FALSE,colClasses = "character"))  
#ce_coord <- as.data.frame(read.csv(paste0(folder,"/ce_all_scan_range_junctions.bed"), sep="\t", header=FALSE,colClasses = "character"))  
#intron_coord <- as.data.frame(read.csv(paste0(folder,"/ce_all_scan_intron.bed"), sep="\t", header=FALSE,colClasses = "character"))  
#pcnt=-.6 #percent drop in coverage - 40% of 10 = 4 = (10,4) -> %drop = -6/10=-.6


#print(paste0('got: ',args[1]))
csv_record=as.data.frame(read.csv(args[1], sep="\t", header=FALSE,colClasses = "character"))  
#print(paste0('got: ',args[2]))
ce_coord =as.data.frame(read.csv(args[2], sep="\t", header=FALSE,colClasses = "character"))  
#print(paste0('got: ',args[3]))
intron_coord <- as.data.frame(read.csv(args[4], sep="\t", header=FALSE,colClasses = "character"))  
#print(paste0('got: ',args[4]))
folder=str_split_fixed(args[1],"/",2)[1,1]
pcnt = -1*as.numeric(args[3])


pct <- function(x) diff(x)/x[-(length(x))]
#print(paste0(pct))
pct1 <- function(x) diff(x)/x[-(1)] #for - strands scanning reverse
#Check its existence
if (file.exists(paste0(folder,"/ce_inclusion_coord.bed"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/ce_inclusion_coord.bed"))
}
if (file.exists(paste0(folder,"/ce_extension_coord.bed"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/ce_extension_coord.bed"))
}
if (file.exists(paste0(folder,"/coverage_file_not_found.txt"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/coverage_file_not_found.txt"))
}

if (file.exists(paste0(folder,"/coverage_file_not_found.bed"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/coverage_file_not_found.bed"))
}

if (file.exists(paste0(folder,"/skipped_ce.csv"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/skipped_ce.csv"))
}

if (file.exists(paste0(folder,"/IGV_R_returned_ce_inclusion.csv"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/IGV_R_returned_ce_inclusion.csv"))
}

if (file.exists(paste0(folder,"/IGV_coverage_file_notfound_ce.csv"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/IGV_coverage_file_notfound_ce.csv"))
}

if (file.exists(paste0(folder,"/IGV_R_returned_ce_extension.csv"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/IGV_R_returned_ce_extension.csv"))
}

if (file.exists(paste0(folder,"/IGV_R_ce1_ce2.csv"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/IGV_R_ce1_ce2.csv"))
}
if (file.exists(paste0(folder,"/ce_inclusion_coord_sashimi.bed"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/ce_inclusion_coord_sashimi.bed"))
}
if (file.exists(paste0(folder,"/ce_extension_coord_sashimi.bed"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/ce_extension_coord_sashimi.bed"))
}

if (file.exists(paste0(folder,"/IGV_R_returned_IR.csv"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/IGV_R_returned_IR.csv"))
}

if (file.exists(paste0(folder,"/IR_coord_sashimi.bed"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/IR_coord_sashimi.bed"))
}

if (file.exists(paste0(folder,"/IR_coord.bed"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/IR_coord.bed"))
}

if (file.exists(paste0(folder,"/skipped_ce_sashimi.csv"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/skipped_ce_sashimi.csv"))
}




#pcnt=-.6 #percent drop in coverage - 40% of 10 = 4 = (10,4) -> %drop = -6/10=-.6

last_genename=""
k=0
IRFlag=0
scan_flag=0 #0 for start<end, 1 for start>end
start_cov=0
end_cov=0
ce_coord$V2 <- as.numeric(ce_coord$V2)
ce_coord$V3 <- as.numeric(ce_coord$V3)
#read every 3rd entry starting from 2
cnt_ce=0
cnt_ir=0
cnt_ext=0
cnt_f=0
cnt_cf=0
ii=0 #for csv_record
#print(paste0('managed here'))
for (i in seq(2,dim(ce_coord)[1],3)) #every 2nd entry (out of group of 3 - up_ex-ce-dn_ex) is probably ce to decide
{
  ii=ii+1
  #print(paste('i: ',i,', ii: ',ii))
  gene_name=ce_coord[i,]$V7
  chr_n=ce_coord[i,]$V1
  start=ce_coord[i,]$V2
  end=ce_coord[i,]$V3
  strnd=ce_coord[i,]$V6
  
  star_int=intron_coord[ii,]$V2
  end_int=intron_coord[ii,]$V3
  
  
  
  #also get intron start and end
  if(strnd=="+")
  {
    cstart=ce_coord[i-1,]$V3
    cend=ce_coord[i+1,]$V2
  }
  else
  {
    cstart=ce_coord[i+1,]$V3
    cend=ce_coord[i-1,]$V2
  }

  #print(paste('start and end is: ',start,end))
  if(start<end)
  {
    #ce1=start
    scan_flag=0
  }else
  {
    #ce1=end
    scan_flag=1
  }
  #make sure events with same gene_names are counted properly
#  if(last_genename == gene_name)
#  {
#    k=k+1
#    fname=paste0("coverages/",gene_name,"_",k,"_TDP.cov.bed")
#    last_genename = gene_name
#  }else
#  {
    k=1
    #fname=paste0("coverages/",gene_name,"_",k,"_TDP.cov.bed")
    #NOW COVERAGE FILE NAME FOR EACH EVENT IS CHANGED TO chr_start_end_gene_name to make each cov file unique
    fname=paste0("coverages/",chr_n,"_",star_int,"_",end_int,"_",gene_name,".cov.bed")
    
    #last_genename = gene_name
#  }

  #now read coverage file
  #if (file.exists(fname)&&fflg==1)
  if (file.exists(fname))    
  {
  CovData <- as.data.frame(read.csv(fname, sep="\t", header=FALSE))
  #now get ce coverage data for the two ends
  if (scan_flag==0)
  {
    data_val= CovData[,9]
    #start_cov=data_val[1]
    #change start_cov position as we are scanning whole intron now
    begt=as.numeric(CovData[,2])
    beg=begt[1]
    offv=start-beg
    start_cov=data_val[offv]
    data_val1=data_val[- (1:offv-1)]
    #distance=which(data_val <= floor(pcnt*start_cov))
    #end_cov_index= min(which(data_val <= floor(pcnt*start_cov)))
    #print(paste0('scanflg=0, pct is:',pct(data_val)))
    #end_cov_index=as.numeric(min(which(pct(data_val) <= -pcnt))) #new criterion
    
    end_cov_index=as.numeric(min(which(pct(data_val1) <= pcnt))) #same as befor - just used -sign with pcnt
    
    
    #print(paste0('0 end_cov_ind: ',end_cov_index+offv-1))

    ce1=start
    if(end_cov_index==Inf)
    {
      
      #print(paste0(gene_name,' is probably IR event, Reporting half range: '))
      #NOW CHECK FOR WHOLE INTRON COVERAGE
      #FIRST CHECK FOR THE COVERAGES ON THE RIGHT TILL NEXT EXON
      data_val1=data_val[(offv:1)] #get data on the right
      #end_cov_index=as.numeric(min(which(pct(data_val1) <= pcnt))) #same as befor - just used -sign with pcnt
      end_cov_index1=as.numeric(min(which(pct(data_val1) <= pcnt))) #same as befor - just used -sign with pcnt
      
      #print(paste0('0 end_cov_ind: ',end_cov_index+offv-1))
      if(end_cov_index1==Inf) #THIS IS IR EVENT
      {
        if(strnd=="+")
        {
        ce1=ce_coord[i-1,]$V3+1
        ce2=ce_coord[i+1,]$V2-1
        #print(paste0('IR EVENT DETECTED for GENE_NAME: ',gene_name,' ce1: ',ce1,' ce2: ',ce2))
        }else
        {
          ce1=ce_coord[i-1,]$V2-1
          ce2=ce_coord[i+1,]$V3+1
          #print(paste0('IR EVENT DETECTED for GENE_NAME: ',gene_name,' ce1: ',ce1,' ce2: ',ce2))
          
        }
        IRFlag=2
      }
      else #CE_EXTENSION
      {
        #if(end_cov_index>1)
        #{
        if(strnd=="+")
        {
          ce2=ce_coord[i+1,]$V2-1
        #ce2=ce1+end_cov_index #rev_data_frame[end_cov_index,]$V8
        
        }
        else
        {
          ce2=start
          ce1=ce_coord[i-1,]$V2-1
        }
        IRFlag=1
        #}
        #else
        #{
         # ce2=end #rev_data_frame[end_cov_index,]$V8
          #IRFlag=1
          
        #}
      }
    }
    else
    {
      IRFlag=0
      #ce2=ce1+CovData[end_cov_index,]$V8
      #changed due to whole intron
      ce2=ce1+end_cov_index-1 #CovData[end_cov_index+offv-1,]$V8
    }
    if (ce1>ce2) #avoid start>end
    {
      tempt=ce1
      ce1=ce2
      ce2=tempt
      
    }
    
  }
  else
  {
    #here we first reverse the dataframe w.r.to rows
    # transpose of dataframe 
    transpose <- t(CovData)
    
    # converting the result to dataframe
    transpose <- as.data.frame(transpose)
    
    # calculating reverse of dataframe
    rev_data_frame <- rev(transpose)
    
    # transpose of reverse dataframe 
    rev_data_frame <- t(rev_data_frame)
    
    # converting the result to dataframe
    rev_data_frame <- as.data.frame(rev_data_frame)
    rev_data_frame<-transform(rev_data_frame, V8 = as.numeric(V8), 
                              V9 = as.numeric(V9))
    
    
    data_val= rev_data_frame[,9]
    #start_cov=data_val[1]
    #change start_cov position as we are scanning whole intron now
    endt=as.numeric(rev_data_frame[,3]) #get end position
    et=endt[1]
    offv=et-start
    start_cov=data_val[offv]
    data_val1=data_val[- (1:offv-1)]
    #data_val1=rev(data_val[1:offv-1]) #reverse
    #data_val1=data_val[offv-1:1]
    #distance=which(data_val <= floor(pcnt*start_cov))
    #end_cov_index= min(which(data_val <= floor(pcnt*start_cov))) #which(data_val < floor(.3*start_cov) )
    #print(paste0('scanflg = 1, pct is:',pct(data_val)))
    end_cov_index=as.numeric(min(which(pct(data_val1) <= pcnt))) #new criterion
    #end_cov_index=as.numeric(min(which(pct1(data_val1) <= pcnt))) #new criterion
    #print(paste0('1 end_cov_ind: ',end_cov_index+offv-1))
    ce2=start
    if(end_cov_index==Inf)
    {
      #print(paste0(gene_name,' is probably IR event, Reporting half range: '))
      #NOW CHECK FOR WHOLE INTRON COVERAGE
      #end_cov_index=as.numeric(min(which(pct(data_val) <= pcnt))) #same as befor - just used -sign with pcnt
      #NOTE HERE WE ARE SCANNING IN THE FORWARD DIRECTION, SO USING pct instead of pct1
      data_val1=data_val[- (1:offv-1)]
      end_cov_index1=as.numeric(min(which(pct(data_val1) <= pcnt))) #same as befor - just used -sign with pcnt
      #print(paste0('0 end_cov_ind: ',end_cov_index1+offv-1))
      if(end_cov_index1==Inf) #THIS IS IR EVENT
      {
        #print(paste0('IR EVENT DETECTED for GENE_NAME: ',gene_name))
        IRFlag=2
        
        if(strnd=="+")
        {
          ce1=ce_coord[i-1,]$V3+1
          ce2=ce_coord[i+1,]$V2-1
        }
        else
        {
          ce1=ce_coord[i+1,]$V3+1
          ce2=ce_coord[i-1,]$V2-1
          
        }
      }
      else #CE_EXTENSION
      {
        #if (end_cov_index>1)
        #{
          if(strnd=="+")
          {
            ce2=start
            ce1=ce_coord[i-1,]$V3+1
            #ce1=ce2-end_cov_index #end #rev_data_frame[end_cov_index,]$V8
          }else
          {
            ce1=ce_coord[i+1,]$V3+1
            ce2=start
            #ce1=ce2-end_cov_index #end #rev_data_frame[end_cov_index,]$V8
          }
          IRFlag=1
        #}
        #else
        #{
         # ce1=end #rev_data_frame[end_cov_index,]$V8
          #IRFlag=1
          
        #}
      }
      
    }
    else
    {
	    #print(paste0('end of coverage is: ',end_cov_index))
      IRFlag=0
      if(strnd=="-")
      {
      	#ce1=ce2-end_cov_index +1 #rev_data_frame[end_cov_index,]$V8
        #chnaged
      	ce1=ce2-end_cov_index+2
      }
      else
      {
	      #ce1=ce2-end_cov_index+3 #3#do not know why i have add 3 here
	      ce1=ce2-end_cov_index+4
      }
      #save bedfile for current junction
      #write.table(SJdat, file="ce_inclusion_coord.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    }
    if (ce1>ce2) #avoid star>end
    {
      tempt=ce1
      ce1=ce2
      ce2=tempt
      
    }
    
    
  }
  #print(paste0("ce range for ",gene_name,'_',k,': ',ce1,'-',ce2))
  
  #save bedfile for current junction
  if(IRFlag==0 && ce1!=ce2) #avoid ce1=ce2
  {
    cnt_ce=cnt_ce+1
    SJdat <- data.frame(seqnames=ce_coord[i,]$V1,
                        #start=SpliceData[i,2]-1,
                        start=ce1-1, #no idea why -1
                        end=ce2,
                        name=1, #reverted back ce_coord[i,]$V4, #1, #changed to exon size
                        score=0,
                        strand=ce_coord[i,]$V6,
                        gene_name=ce_coord[i,]$V7)
    
    #write.table(SJdat, file="ce_inclusion_coord.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    #also new file with 5',ce and 3'
    recrd_1=ce_coord[i-1,]
    recrd_1$V4=1 #change exon length to 1 - causes trouble
    write.table(recrd_1, file=paste0(folder,"/ce_inclusion_coord.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    write.table(SJdat, file=paste0(folder,"/ce_inclusion_coord.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    recrd_2=ce_coord[i+1,]
    recrd_2$V4=1
    #write.table(ce_coord[i+1,], file="folder/ce_inclusion_coord.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    write.table(recrd_2, file=paste0(folder,"/ce_inclusion_coord.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    #also write bed file for sashimi plots including TxID
    c1=ce_coord[i-1,]
    c1$V8=strsplit(csv_record[ii,],',')[[1]][9]
    SJdat$TxID=strsplit(csv_record[ii,],',')[[1]][9]
    c2=ce_coord[i+1,]
    c2$V8=strsplit(csv_record[ii,],',')[[1]][9]
    write.table(c1, file=paste0(folder,"/ce_inclusion_coord_sashimi.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    write.table(SJdat, file=paste0(folder,"/ce_inclusion_coord_sashimi.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    write.table(c2, file=paste0(folder,"/ce_inclusion_coord_sashimi.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
        #also write csv file
    write.table(csv_record[ii,],file=paste0(folder,"/IGV_R_returned_ce_inclusion.csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)    
  }
  else if (IRFlag==1 && ce1!=ce2) #extension event
  {
    cnt_ext=cnt_ext+1
    SJdat <- data.frame(seqnames=ce_coord[i,]$V1,
                        #start=SpliceData[i,2]-1,
                        start=ce1,
                        end=ce2-1, #to match with DNM1 - for exton extension events
                        name=1, #ce_coord[i,]$V4, #1, 
                        score=0,
                        strand=ce_coord[i,]$V6,
                        gene_name=ce_coord[i,]$V7)
    
    #write.table(SJdat, file="ce_extension_coord.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    
    #also new file with 5',ce and 3'
    recrd_1=ce_coord[i-1,]
    recrd_1$V4=1
    write.table(recrd_1, file=paste0(folder,"/ce_extension_coord.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    write.table(SJdat, file=paste0(folder,"/ce_extension_coord.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    recrd_2=ce_coord[i+1,]
    recrd_2$V4=1
    
    write.table(recrd_2,file=paste0(folder,"/ce_extension_coord.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    #also write bed file for sashimi plots including TxID
    c1=ce_coord[i-1,]
    c1$V8=strsplit(csv_record[ii,],',')[[1]][9]
    SJdat$TxID=strsplit(csv_record[ii,],',')[[1]][9]
    c2=ce_coord[i+1,]
    c2$V8=strsplit(csv_record[ii,],',')[[1]][9]
    write.table(c1, file=paste0(folder,"/ce_extension_coord_sashimi.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    write.table(SJdat, file=paste0(folder,"/ce_extension_coord_sashimi.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    write.table(c2, file=paste0(folder,"/ce_extension_coord_sashimi.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    
    
    #also write csv file
    write.table(csv_record[ii,],file=paste0(folder,"/IGV_R_returned_ce_extension.csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)    
    
  }
  else if (IRFlag==2 && ce1!=ce2) #IR events
  {
    cnt_ir=cnt_ir+1
    SJdat <- data.frame(seqnames=ce_coord[i,]$V1,
                        #start=SpliceData[i,2]-1,
                        start=ce1,
                        end=ce2-1, #to match with DNM1 - for exton extension events
                        name=1, #ce_coord[i,]$V4, #1, 
                        score=0,
                        strand=ce_coord[i,]$V6,
                        gene_name=ce_coord[i,]$V7)
    
    #write.table(SJdat, file="ce_extension_coord.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    
    #also new file with 5',ce and 3'
    #recrd_1=ce_coord[i-1,]
    #recrd_1$V4=1
    #write.table(recrd_1, file=paste0(folder,"/ce_extension_coord.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    write.table(SJdat, file=paste0(folder,"/IR_coord.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    #recrd_2=ce_coord[i+1,]
    #recrd_2$V4=1
    
    #write.table(recrd_2,file="folder/ce_extension_coord.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    #also write bed file for sashimi plots including TxID
    #c1=ce_coord[i-1,]
    #c1$V8=strsplit(csv_record[ii,],',')[[1]][9]
    #SJdat$TxID=strsplit(csv_record[ii,],',')[[1]][9]
    #c2=ce_coord[i+1,]
    #c2$V8=strsplit(csv_record[ii,],',')[[1]][9]
    #write.table(c1, file="folder/ce_extension_coord_sashimi.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    SJdat$TxID=strsplit(csv_record[ii,],',')[[1]][9]
    write.table(SJdat, file=paste0(folder,"/IR_coord_sashimi.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    #write.table(c2, file="folder/ce_extension_coord_sashimi.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    
    
    #also write csv file
    write.table(csv_record[ii,],file=paste0(folder,"/IGV_R_returned_IR.csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)    
    
  }
  else
  {
    cnt_f=cnt_f+1
    print(paste('ce1==ce2 ',ce1,ce2,' for gene: ',gene_name, ' so skipping for ce and events are in folder/skipped_ce.csv'))
    SJdat <- data.frame(seqnames=ce_coord[i,]$V1,
                        #start=SpliceData[i,2]-1,
                        start=ce1,
                        end=ce2,
                        name=1, #ce_coord[i,]$V4, #1,
                        score=0,
                        strand=ce_coord[i,]$V6,
                        gene_name=ce_coord[i,]$V7)
    
    #write.table(SJdat, file="ce_extension_coord.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    
    #also new file with 5',ce and 3'
    recrd_1=ce_coord[i-1,]
    recrd_1$V4=1
    
    write.table(recrd_1, file=paste0(folder,"/skipped_ce.csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
    write.table(SJdat, file=paste0(folder,"/skipped_ce.csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
    recrd_2=ce_coord[i+1,]
    recrd_2$V4=1
    
    write.table(recrd_2,file=paste0(folder,"/skipped_ce.csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
    #also write bed file for sashimi plots including TxID
    c1=ce_coord[i-1,]
    c1$V8=strsplit(csv_record[ii,],',')[[1]][9]
    SJdat$TxID=strsplit(csv_record[ii,],',')[[1]][9]
    c2=ce_coord[i+1,]
    c2$V8=strsplit(csv_record[ii,],',')[[1]][9]
    write.table(c1, file=paste0(folder,"/skipped_ce_sashimi.csv"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    write.table(SJdat, file=paste0(folder,"/skipped_ce_sashimi.csv"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    write.table(c2, file=paste0(folder,"/skipped_ce_sashimi.csv"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    
    
    
    write.table(csv_record[ii,],file=paste0(folder,"/IGV_R_ce1_ce2.csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)    
  
  }
  }
  else
  {
    cnt_cf=cnt_cf+1
    print(paste0("coverage file: ",fname," does not exists, so skipping"))
    write.table(paste0("coverage file: ",fname," does not exists, so skipping"), file=paste0(folder,"/coverage_file_not_found.txt"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    #also write bed file to generate coverage files for missing junctions
    #also new file with 5',ce and 3'
    #write.table(ce_coord[i-1,], file="ce_extension_coord.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    write.table(ce_coord[i,], file=paste0(folder,"/coverage_file_not_found.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    #write.table(ce_coord[i+1,],file="ce_extension_coord.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    #also write csv file
    write.table(csv_record[ii,],file=paste0(folder,"/IGV_coverage_file_notfound_ce.csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)    
  }
}
print(paste('total events read are: ',dim(ce_coord)[1]/3))

#cat(paste('Total events unsuccessful (due to ce1==ce2) are: ',cnt_f, 'please see folder/skipped_ce.csv'),file="folder/Summary_stats.txt",sep="\n",append=TRUE)
cat(paste('Finished R session for CE_boundary calculations of a total of: ',dim(ce_coord)[1]/3,' successful events, please see folder/IGV_ce_inclusion.csv'),file=paste0(folder,"/Summary_stats.txt"),sep="\n",append=TRUE)
print(paste('Total events for which ce_coverage files not found are : ',cnt_cf, ' - please see folder/coverage_file_not_found.bed and folder/IGV_coverage_file_notfound_ce.csv'))
cat(paste('Total events for which ce_coverage files not found are : ',cnt_cf,' (please see folder/coverage_file_not_found.bed)'),file=paste0(folder,"/Summary_stats.txt"),sep="\n",append=TRUE)
print(paste('Total events unsuccessful (due to ce1==ce2) are: ',cnt_f, ' - Please see file folder/skipped_ce.csv file'))
cat(paste('Total events unsuccessful (due to ce1==ce2)   are: ',cnt_f,' Please see file folder/skipped_ce.csv file. It contains us, ce and ds exon coordinates for each event'),file=paste0(folder,"/Summary_stats.txt"),sep="\n",append=TRUE)
print(paste('Total events successfully processed are: ',cnt_ce+cnt_ir+cnt_ext))
#print(paste('Total events successfully processed are: ',cnt_ce+cnt_ir))
cat(paste('Total events successfully processed are: ',cnt_ce+cnt_ir+cnt_ext),file=paste0(folder,"/Summary_stats.txt"),sep="\n",append=TRUE)
print(paste('Out of total ',cnt_ce+cnt_ir+cnt_ext,' successfull events processed, ce_inclusion events are: ',cnt_ce, 'ce_extension events are: ',cnt_ext, ' and IR events are: ',cnt_ir))
cat(paste('Out of total: ',cnt_ce+cnt_ir+cnt_ext,' successfull events processed'),file=paste0(folder,"/Summary_stats.txt"),sep="\n",append=TRUE)
cat(paste('            ce_inclusion events are: ',cnt_ce),file=paste0(folder,"/Summary_stats.txt"),sep="\n",append=TRUE)
cat(paste('            ce_extension events are: ',cnt_ext),file=paste0(folder,"/Summary_stats.txt"),sep="\n",append=TRUE)
cat(paste('            IR events are: ',cnt_ir),file=paste0(folder,"/Summary_stats.txt"),sep="\n",append=TRUE)
#}
