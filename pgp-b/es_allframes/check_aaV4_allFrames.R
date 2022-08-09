library(stringr)
library(sjmisc)
library(tidyverse)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

#Read cds_aa fasta file
cds_aa <- readLines(args[1])
#Read all frames FASTA file
all_frames_aa <- readLines(args[2])
folder=str_split_fixed(args[1],"/",2)[1,1]
#also read ce_range as well
if (length(args)>3)
{
all_ce_range <-read.csv(args[4],sep="\t", header=FALSE)

all_ce_events <-read.csv(args[5],sep="\t", header=FALSE)
}
#print(paste0(all_ce_range))

#setwd("/Volumes/SYEDSHAH/MichaelLab/proteogenomic_pipeline/github_final_pgp/phase-b/ce_incl/v3")
#Read cds_aa fasta file
#cds_aa <- readLines("res_ce_0.6/aa.fasta")
#Read all frames FASTA file
#all_frames_aa <- readLines("res_ce_0.6/cds_CE_INCLUSION_FUSED_AA.fasta")

if (file.exists(paste0(folder,"/final_aa.fasta"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/final_aa.fasta"))
}

if (file.exists(paste0(folder,"/all_frames_final_aa.fasta"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/all_frames_final_aa.fasta"))
}


if (file.exists(paste0(folder,"/cds_PEAKS_coord_only.txt"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/cds_PEAKS_coord_only.txt"))
}

if (file.exists(paste0(folder,"/cds_PEAKS_IGV_events.csv"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/cds_PEAKS_IGV_events.csv"))
}


write(paste0('                                 '),file=args[3],append=TRUE)
write(paste0('Starting From check_aaV1.R: --------------- Processing file: ',args[2],' against reference AA (for protein coding genes only) in file: ',args[1]),file=args[3],append=TRUE)
#index in cds_aa file
cds_indx=1
not_found=0
eventn=0
all_frames_flag=1
for (i in seq(1,length(all_frames_aa),6))
{
  #aa_seq1=substr(all_frames_aa[i+1],2,15)
  #if (grepl(substr(all_frames_aa[i+1],2,15), cds_aa[cds_indx+1],fixed=T) )
  #NOW CHECKING ALL FRAME RETURNED BY cds AA
  eventn=eventn+1
  #first check upstream exon length
  if (length(args)>3)
  {
  strnd=str_split(all_frames_aa[i], "_")[[1]][10]
  if(strnd=="plus")
  {
    length1=as.numeric(str_split(all_frames_aa[i], "_")[[1]][4]) - as.numeric(str_split(all_frames_aa[i], "_")[[1]][3])
    length=length1%/%3
    if(length>15 && length <=20)length=16
  }else
  {
    length1=as.numeric(str_split(all_frames_aa[i], "_")[[1]][8]) - as.numeric(str_split(all_frames_aa[i], "_")[[1]][7])
    length=length1%/%3
    if(length>15 && length <=20)length=16
  }
  
  }else
  {
    strnd=str_split(all_frames_aa[i], "_")[[1]][8]
  
  if(strnd=="plus")
  {
    length1=as.numeric(str_split(all_frames_aa[i], "_")[[1]][4]) - as.numeric(str_split(all_frames_aa[i], "_")[[1]][3])
    length=length1%/%3
    if(length>15 && length <=20)length=16
  }else
  {
    length1=as.numeric(str_split(all_frames_aa[i], "_")[[1]][6]) - as.numeric(str_split(all_frames_aa[i], "_")[[1]][5])
    length=length1%/%3
    #if(length>15 && length <=20)length=16
    if(length>15)length=16
  }
  }
  if (all_frames_flag == 0)
  {
    #if (grep(substr(all_frames_aa[i+1],1,15), cds_aa[cds_indx+1]))
    if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+1],2,length-1))||str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+1],nchar(all_frames_aa[i+1])-19,nchar(all_frames_aa[i+1])-1)))
    #if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+1],2,15)) || str_contains(cds_aa[cds_indx+3],substr(all_frames_aa[i+1],2,15)) || str_contains(cds_aa[cds_indx+5],substr(all_frames_aa[i+1],2,15)) )
    {
      write.table(all_frames_aa[i],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      write.table(all_frames_aa[i+1],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      
    }  
    #else if (grep(substr(all_frames_aa[i+3],1,15), cds_aa[cds_indx+1]))
    else if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+3],2,length-1))||str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+3],nchar(all_frames_aa[i+3])-19,nchar(all_frames_aa[i+3])-1)))
    #else if (grepl(substr(all_frames_aa[i+3],2,15), cds_aa[cds_indx+1])||grepl(substr(all_frames_aa[i+3],2,15), cds_aa[cds_indx+3])||grepl(substr(all_frames_aa[i+3],2,15), cds_aa[cds_indx+5]))
    #else if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+3],2,15)) || str_contains(cds_aa[cds_indx+3],substr(all_frames_aa[i+3],2,15)) || str_contains(cds_aa[cds_indx+5],substr(all_frames_aa[i+3],2,15)) )
    {
      write.table(all_frames_aa[i+2],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      write.table(all_frames_aa[i+3],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      
    }
    #else if (grep(substr(all_frames_aa[i+5],1,15), cds_aa[cds_indx+1]))
    else if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+5],2,length-1))||str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+5],nchar(all_frames_aa[i+5])-19,nchar(all_frames_aa[i+5])-1)))
    #else if (grepl(substr(all_frames_aa[i+5],2,15), cds_aa[cds_indx+1])||grepl(substr(all_frames_aa[i+5],2,15), cds_aa[cds_indx+3])||grepl(substr(all_frames_aa[i+5],2,15), cds_aa[cds_indx+5]))
    #else if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+5],2,15)) || str_contains(cds_aa[cds_indx+3],substr(all_frames_aa[i+5],2,15)) || str_contains(cds_aa[cds_indx+5],substr(all_frames_aa[i+5],2,15)) )
    {
      write.table(all_frames_aa[i+4],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      write.table(all_frames_aa[i+5],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      
    }
    else
    {
      not_found=not_found+1
      print(paste0(not_found,': no aa frame found for gene: ',cds_aa[cds_indx], ' and event num: ',eventn))
      write(paste0('From check_aaV1.R: ',not_found,': no aa frame found for gene: ',cds_aa[cds_indx], ' and event num: ',eventn, ' in file: ',args[2]),file=args[3],append=TRUE)
    }
    cds_indx=cds_indx+2
  }
  else
  {
    #if (grep(substr(all_frames_aa[i+1],1,15), cds_aa[cds_indx+1]))
    #if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+1],2,15))||str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+1],nchar(all_frames_aa[i+1])-19,nchar(all_frames_aa[i+1])-1)))
    #if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+1],2,length-1)) || str_contains(cds_aa[cds_indx+3],substr(all_frames_aa[i+1],2,length-1)) || str_contains(cds_aa[cds_indx+5],substr(all_frames_aa[i+1],2,length-1)) )
    if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+1],1,length)) || str_contains(cds_aa[cds_indx+3],substr(all_frames_aa[i+1],1,length)) || str_contains(cds_aa[cds_indx+5],substr(all_frames_aa[i+1],1,length)))      
    {
      
      #write.table(all_frames_aa[i],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      #if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+1],2,length-1)))
      if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+1],1,length))) 
      {
        #find position in the string
        #start_end=str_locate(cds_aa[cds_indx+1],substr(all_frames_aa[i+1],2,length-1))  #position of AA seq in fasta file
        start_end=str_locate(cds_aa[cds_indx+1],substr(all_frames_aa[i+1],1,length))  #position of AA seq in fasta file
        
        #also modify fastaID as well. Append aa length from cds exon1 upto us_exon (excluding us_exon length)
        tem_fastaID=paste0(all_frames_aa[i],'_',nchar(substr(cds_aa[cds_indx+1],1,start_end[1]-1)))
        print(paste0('new ID: ',tem_fastaID, ' Old ID: ',all_frames_aa[i]))
        write.table(tem_fastaID,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)

        #now get the string from cds_aa upto start_end[2] and concatenate with all_frames_aa upto end
        #new_aa=paste0(substr(cds_aa[cds_indx+1],1,start_end[2]),substr(all_frames_aa[i+1],16,nchar(all_frames_aa[i+1])))
        #new_aa=paste0(substr(cds_aa[cds_indx+1],1,start_end[1]-2),all_frames_aa[i+1])
        print(paste0('exon1 length is: ',nchar(substr(cds_aa[cds_indx+1],1,start_end[1]-1))))
        new_aa=paste0(substr(cds_aa[cds_indx+1],1,start_end[1]-1),all_frames_aa[i+1])
        #write.table(cds_aa[cds_indx+1],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        write.table(new_aa,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      }
      #else if (str_contains(cds_aa[cds_indx+3],substr(all_frames_aa[i+1],2,length-1)) )
      else if (str_contains(cds_aa[cds_indx+3],substr(all_frames_aa[i+1],1,length)) )
      {
        #find position in the string
        #start_end=str_locate(cds_aa[cds_indx+3],substr(all_frames_aa[i+1],2,length-1))  #position of AA seq in fasta file
        start_end=str_locate(cds_aa[cds_indx+3],substr(all_frames_aa[i+1],1,length))  #position of AA seq in fasta file
        #now get the string from cds_aa upto start_end[2] and concatenate with all_frames_aa upto end
        #new_aa=paste0(substr(cds_aa[cds_indx+3],1,start_end[2]),substr(all_frames_aa[i+1],16,nchar(all_frames_aa[i+1])))
        #new_aa=paste0(substr(cds_aa[cds_indx+3],1,start_end[1]-2),all_frames_aa[i+1])
        
        #also modify fastaID as well. Append aa length from cds exon1 upto us_exon (excluding us_exon length)
        tem_fastaID=paste0(all_frames_aa[i],'_',nchar(substr(cds_aa[cds_indx+3],1,start_end[1]-1)))
        print(paste0('new ID: ',tem_fastaID, ' Old ID: ',all_frames_aa[i]))
        write.table(tem_fastaID,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        
        print(paste0('exon1 length is: ',nchar(substr(cds_aa[cds_indx+3],1,start_end[1]-1))))
        new_aa=paste0(substr(cds_aa[cds_indx+3],1,start_end[1]-1),all_frames_aa[i+1])
        #write.table(cds_aa[cds_indx+3],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        write.table(new_aa,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      }
      #else if (str_contains(cds_aa[cds_indx+5],substr(all_frames_aa[i+1],2,length-1)) )
      else if (str_contains(cds_aa[cds_indx+5],substr(all_frames_aa[i+1],1,length)) )
      {
        #find position in the string
        #start_end=str_locate(cds_aa[cds_indx+5],substr(all_frames_aa[i+1],2,length-1))  #position of AA seq in fasta file
        start_end=str_locate(cds_aa[cds_indx+5],substr(all_frames_aa[i+1],1,length))  #position of AA seq in fasta file
        #now get the string from cds_aa upto start_end[2] and concatenate with all_frames_aa upto end
        #new_aa=paste0(substr(cds_aa[cds_indx+5],1,start_end[2]),substr(all_frames_aa[i+1],16,nchar(all_frames_aa[i+1])))
        #new_aa=paste0(substr(cds_aa[cds_indx+5],1,start_end[1]-2),all_frames_aa[i+1])
        
        #also modify fastaID as well. Append aa length from cds exon1 upto us_exon (excluding us_exon length)
        tem_fastaID=paste0(all_frames_aa[i],'_',nchar(substr(cds_aa[cds_indx+5],1,start_end[1]-1)))
        print(paste0('new ID: ',tem_fastaID, ' Old ID: ',all_frames_aa[i]))
        write.table(tem_fastaID,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        
        print(paste0('exon1 length is: ',nchar(substr(cds_aa[cds_indx+5],1,start_end[1]-1))))
        new_aa=paste0(substr(cds_aa[cds_indx+5],1,start_end[1]-1),all_frames_aa[i+1])
        #write.table(cds_aa[cds_indx+5],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        write.table(new_aa,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      }
      
      #write.table(all_frames_aa[i+1],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      if (length(args)>3)
      {
      write.table(all_ce_range[eventn,],file=paste0(folder,"/cds_PEAKS_coord_only.txt"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
      
      write.table(all_ce_events[eventn,],file=paste0(folder,"/cds_PEAKS_IGV_events.csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
      }
      #print(paste0('got: ',all_ce_range[eventn,]))
      
      
      #ALSO SAVE ALL FRAMES FOR NEW REQUIREMENT

      write.table(tem_fastaID,file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      write.table(new_aa,file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      
      #also modify fastaID as well. Append aa length from cds exon1 upto us_exon (excluding us_exon length)
      zr=0
      tem_fastaID1=paste0(all_frames_aa[i+2],'_',zr)
      #print(paste0('new ID: ',tem_fastaID, ' Old ID: ',all_frames_aa[i]))
      
      write.table(tem_fastaID1,file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      write.table(all_frames_aa[i+3],file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      
      tem_fastaID1=paste0(all_frames_aa[i+4],'_',zr)
      
      write.table(tem_fastaID1,file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      write.table(all_frames_aa[i+5],file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      
      
      
      
    }  
    #else if (grep(substr(all_frames_aa[i+3],1,15), cds_aa[cds_indx+1]))
    #else if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+3],2,15))||str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+3],nchar(all_frames_aa[i+3])-19,nchar(all_frames_aa[i+3])-1)))
    #else if (grepl(substr(all_frames_aa[i+3],2,15), cds_aa[cds_indx+1])||grepl(substr(all_frames_aa[i+3],2,15), cds_aa[cds_indx+3])||grepl(substr(all_frames_aa[i+3],2,15), cds_aa[cds_indx+5]))
    #else if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+3],2,length-1)) || str_contains(cds_aa[cds_indx+3],substr(all_frames_aa[i+3],2,length-1)) || str_contains(cds_aa[cds_indx+5],substr(all_frames_aa[i+3],2,length-1)) )
    else if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+3],1,length)) || str_contains(cds_aa[cds_indx+3],substr(all_frames_aa[i+3],1,length)) || str_contains(cds_aa[cds_indx+5],substr(all_frames_aa[i+3],1,length)) )
    {
      #write.table(all_frames_aa[i+2],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      #if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+3],2,length-1)))
      if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+3],1,length)))
      {
        #find position in the string
        #start_end=str_locate(cds_aa[cds_indx+1],substr(all_frames_aa[i+3],2,length-1))  #position of AA seq in fasta file
        start_end=str_locate(cds_aa[cds_indx+1],substr(all_frames_aa[i+3],1,length))  #position of AA seq in fasta file
        #now get the string from cds_aa upto start_end[2] and concatenate with all_frames_aa upto end
        #new_aa=paste0(substr(cds_aa[cds_indx+1],1,start_end[2]),substr(all_frames_aa[i+3],16,nchar(all_frames_aa[i+3])))
        #new_aa=paste0(substr(cds_aa[cds_indx+1],1,start_end[1]-2),all_frames_aa[i+3])
        
        #also modify fastaID as well. Append aa length from cds exon1 upto us_exon (excluding us_exon length)
        tem_fastaID=paste0(all_frames_aa[i+2],'_',nchar(substr(cds_aa[cds_indx+1],1,start_end[1]-1)))
        print(paste0('new ID: ',tem_fastaID, ' Old ID: ',all_frames_aa[i+2]))
        write.table(tem_fastaID,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        
        print(paste0('exon1 length is: ',nchar(substr(cds_aa[cds_indx+1],1,start_end[1]-1))))
        new_aa=paste0(substr(cds_aa[cds_indx+1],1,start_end[1]-1),all_frames_aa[i+3])
        #write.table(cds_aa[cds_indx+1],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        write.table(new_aa,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      }
      #else if (str_contains(cds_aa[cds_indx+3],substr(all_frames_aa[i+3],2,length-1)) )
      else if (str_contains(cds_aa[cds_indx+3],substr(all_frames_aa[i+3],1,length)) )
      {
        #find position in the string
        start_end=str_locate(cds_aa[cds_indx+3],substr(all_frames_aa[i+3],1,length))  #position of AA seq in fasta file
        #now get the string from cds_aa upto start_end[2] and concatenate with all_frames_aa upto end
        #new_aa=paste0(substr(cds_aa[cds_indx+3],1,start_end[2]),substr(all_frames_aa[i+3],16,nchar(all_frames_aa[i+3])))
        
        #also modify fastaID as well. Append aa length from cds exon1 upto us_exon (excluding us_exon length)
        tem_fastaID=paste0(all_frames_aa[i+2],'_',nchar(substr(cds_aa[cds_indx+3],1,start_end[1]-1)))
        print(paste0('new ID: ',tem_fastaID, ' Old ID: ',all_frames_aa[i+2]))
        write.table(tem_fastaID,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        
        print(paste0('exon1 length is: ',nchar(substr(cds_aa[cds_indx+3],1,start_end[1]-1))))
        new_aa=paste0(substr(cds_aa[cds_indx+3],1,start_end[1]-1),all_frames_aa[i+3])
        #write.table(cds_aa[cds_indx+3],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        write.table(new_aa,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      }
      else if (str_contains(cds_aa[cds_indx+5],substr(all_frames_aa[i+3],1,length)) )
      {
        #find position in the string
        start_end=str_locate(cds_aa[cds_indx+5],substr(all_frames_aa[i+3],1,length))  #position of AA seq in fasta file
        #now get the string from cds_aa upto start_end[2] and concatenate with all_frames_aa upto end
        #new_aa=paste0(substr(cds_aa[cds_indx+5],1,start_end[2]),substr(all_frames_aa[i+3],16,nchar(all_frames_aa[i+3])))
        
        #also modify fastaID as well. Append aa length from cds exon1 upto us_exon (excluding us_exon length)
        tem_fastaID=paste0(all_frames_aa[i+2],'_',nchar(substr(cds_aa[cds_indx+5],1,start_end[1]-1)))
        print(paste0('new ID: ',tem_fastaID, ' Old ID: ',all_frames_aa[i+2]))
        write.table(tem_fastaID,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        
        print(paste0('exon1 length is: ',nchar(substr(cds_aa[cds_indx+5],1,start_end[1]-1))))
        new_aa=paste0(substr(cds_aa[cds_indx+5],1,start_end[1]-1),all_frames_aa[i+3])
        #write.table(cds_aa[cds_indx+5],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        write.table(new_aa,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      }
      
      #write.table(all_frames_aa[i+3],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      if (length(args)>3)
      {
        
      write.table(all_ce_range[eventn,],file=paste0(folder,"/cds_PEAKS_coord_only.txt"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
      write.table(all_ce_events[eventn,],file=paste0(folder,"/cds_PEAKS_IGV_events.csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
      }

      #ALSO SAVE ALL FRAMES FOR NEW REQUIREMENT
      
      zr=0
      tem_fastaID1=paste0(all_frames_aa[i],'_',zr)
      
      write.table(tem_fastaID1,file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      write.table(all_frames_aa[i+1],file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)

      write.table(tem_fastaID,file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      write.table(new_aa,file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      
      tem_fastaID1=paste0(all_frames_aa[i+4],'_',zr)      
      write.table(tem_fastaID1,file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      write.table(all_frames_aa[i+5],file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      
      
      
    }
    #else if (grep(substr(all_frames_aa[i+5],1,15), cds_aa[cds_indx+1]))
    #else if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+5],2,15))||str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+5],nchar(all_frames_aa[i+5])-19,nchar(all_frames_aa[i+5])-1)))
      #else if (grepl(substr(all_frames_aa[i+5],2,15), cds_aa[cds_indx+1])||grepl(substr(all_frames_aa[i+5],2,15), cds_aa[cds_indx+3])||grepl(substr(all_frames_aa[i+5],2,15), cds_aa[cds_indx+5]))
    else if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+5],1,length)) || str_contains(cds_aa[cds_indx+3],substr(all_frames_aa[i+5],1,length)) || str_contains(cds_aa[cds_indx+5],substr(all_frames_aa[i+5],1,length)) )
    {
      #write.table(all_frames_aa[i+4],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      if (str_contains(cds_aa[cds_indx+1],substr(all_frames_aa[i+5],1,length)))
      {
        #find position in the string
        start_end=str_locate(cds_aa[cds_indx+1],substr(all_frames_aa[i+5],1,length))  #position of AA seq in fasta file
        #now get the string from cds_aa upto start_end[2] and concatenate with all_frames_aa upto end
        #new_aa=paste0(substr(cds_aa[cds_indx+1],1,start_end[2]),substr(all_frames_aa[i+5],16,nchar(all_frames_aa[i+5])))
        
        #also modify fastaID as well. Append aa length from cds exon1 upto us_exon (excluding us_exon length)
        tem_fastaID=paste0(all_frames_aa[i+4],'_',nchar(substr(cds_aa[cds_indx+1],1,start_end[1]-1)))
        print(paste0('new ID: ',tem_fastaID, ' Old ID: ',all_frames_aa[i+4]))
        write.table(tem_fastaID,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        
        print(paste0('exon1 length is: ',nchar(substr(cds_aa[cds_indx+1],1,start_end[1]-1))))
        new_aa=paste0(substr(cds_aa[cds_indx+1],1,start_end[1]-1),all_frames_aa[i+5])
        #write.table(cds_aa[cds_indx+1],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        write.table(new_aa,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      }
      else if (str_contains(cds_aa[cds_indx+3],substr(all_frames_aa[i+5],1,length)))
      {
        #find position in the string
        start_end=str_locate(cds_aa[cds_indx+3],substr(all_frames_aa[i+5],1,length))  #position of AA seq in fasta file
        #now get the string from cds_aa upto start_end[2] and concatenate with all_frames_aa upto end
        #new_aa=paste0(substr(cds_aa[cds_indx+3],1,start_end[2]),substr(all_frames_aa[i+5],16,nchar(all_frames_aa[i+5])))
        
        #also modify fastaID as well. Append aa length from cds exon1 upto us_exon (excluding us_exon length)
        tem_fastaID=paste0(all_frames_aa[i+4],'_',nchar(substr(cds_aa[cds_indx+3],1,start_end[1]-1)))
        print(paste0('new ID: ',tem_fastaID, ' Old ID: ',all_frames_aa[i+4]))
        write.table(tem_fastaID,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        
        print(paste0('exon1 length is: ',nchar(substr(cds_aa[cds_indx+3],1,start_end[1]-1))))
        new_aa=paste0(substr(cds_aa[cds_indx+3],1,start_end[1]-1),all_frames_aa[i+5])
        #write.table(cds_aa[cds_indx+3],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        write.table(new_aa,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      }
      else if (str_contains(cds_aa[cds_indx+5],substr(all_frames_aa[i+5],1,length)) )
      {
        #find position in the string
        start_end=str_locate(cds_aa[cds_indx+5],substr(all_frames_aa[i+5],1,length))  #position of AA seq in fasta file
        #now get the string from cds_aa upto start_end[2] and concatenate with all_frames_aa upto end
        #new_aa=paste0(substr(cds_aa[cds_indx+5],1,start_end[2]),substr(all_frames_aa[i+5],16,nchar(all_frames_aa[i+5])))
        
        #also modify fastaID as well. Append aa length from cds exon1 upto us_exon (excluding us_exon length)
        tem_fastaID=paste0(all_frames_aa[i+4],'_',nchar(substr(cds_aa[cds_indx+5],1,start_end[1]-1)))
        print(paste0('new ID: ',tem_fastaID, ' Old ID: ',all_frames_aa[i+4]))
        write.table(tem_fastaID,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        
        print(paste0('exon1 length is: ',nchar(substr(cds_aa[cds_indx+5],1,start_end[1]-1))))
        new_aa=paste0(substr(cds_aa[cds_indx+5],1,start_end[1]-1),all_frames_aa[i+5])
        #write.table(cds_aa[cds_indx+5],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        write.table(new_aa,file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
        
        
        
      }
      
      #write.table(all_frames_aa[i+5],file=paste0(folder,"/final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      if (length(args)>3)
      {
        
      write.table(all_ce_range[eventn,],file=paste0(folder,"/cds_PEAKS_coord_only.txt"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
      write.table(all_ce_events[eventn,],file=paste0(folder,"/cds_PEAKS_IGV_events.csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
      }
      
      #ALSO SAVE ALL FRAMES FOR NEW REQUIREMENT
      zr=0
      tem_fastaID1=paste0(all_frames_aa[i],'_',zr)
      
      write.table(tem_fastaID1,file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      write.table(all_frames_aa[i+1],file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      
      
      tem_fastaID1=paste0(all_frames_aa[i+2],'_',zr)
      
      write.table(tem_fastaID1,file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      write.table(all_frames_aa[i+3],file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      
      
      write.table(tem_fastaID,file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      write.table(new_aa,file=paste0(folder,"/all_frames_final_aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
      
      
    }
    else
    {
      not_found=not_found+1
      print(paste0(not_found,': no aa frame found for gene: ',cds_aa[cds_indx], ' and event num: ',eventn))
      write(paste0('From check_aaV1.R: ',not_found,': no aa frame found for gene: ',cds_aa[cds_indx], ' and event num: ',eventn, ' in file: ',args[2]),file=args[3],append=TRUE)
    }
    
    cds_indx=cds_indx+6
  }
}

write(paste0('Done From check_aaV1.R:---------------- Processing file: ',args[2], ' Total events for which no canonical frame found are: ',not_found),file=args[3],append=TRUE)
write(paste0('!!!!!!!!*~@#$%^ PLEAS NOTE THAT Removed: ',not_found,' events are still in cds_IGV_unique_skiptics_translated.csv file*~@#$%^!!!!!!!'),file=args[3],append=TRUE)
write(paste0('                                 '),file=args[3],append=TRUE)
