#IMPORTANT
#ignoring whether string matches to 1,2 or 3 frame translated sequence
#TO RUN: Type
# Rscript aa_to_event_type.R PEAKS_***.csv FUSED_AA_***.fasta
# Example: Rscript aa_to_event_type.R PEAKS_results_for_Syed.csv CE_FUSED_AA.fasta

#THIS VERSION:
#1. SAVES bed file for the events to get AA fasta file
#2. SAVES .csv file for SASHIMI's (PLEASE NOTE THAT WE CONCATENATE CE_INCL and CE_EXT csv files here)
#3. As of 02/13/2022 removing /3 for range sizes and catering for up_exon_ce and ce_ds_exon mismatch when exon lengths are not multiple of 3

library(stringr)
header=c('peptide','gene_name','peptide_coord','IGV_coord','fasta_id','event_type')
command_line_flg=1
#ENABLE THIS IF READING FILE FROM terminal
if (command_line_flg==1)
{
  args = commandArgs(trailingOnly=TRUE)
  if (length(args)==0)
  {
    stop("Please provide a valid input csv file (output of PEAKS searc) and FUSED AA fasta file, quitting", call.=FALSE)
  } 
  #Read Peaks File
  peaks_data <- as.data.frame(read.csv(args[1], header=TRUE))
  #Read ce FASTA file
  ce_fasta_strings <- readLines(args[2])
  #Also Read majiq csv file of the selected events (Here we concatenate ce_inclusion and extension csv files)
  ce_incl_ext_csv <- readLines(args[3])
  #Also Read nt seq file ce events
  ce_incl_ext_nt_seq <- readLines(args[4])
  
  #Read skiptics FASTA file
  es_fasta_strings <- readLines(args[5])
  #Also Read majiq csv file of the selected events (Here we concatenate ce_inclusion and extension csv files)
  es_csv <- readLines(args[6])

  #Also Read nt seq file skiptics events
  es_nt_seq <- readLines(args[7])
  
  #Read IR FASTA file
  IR_fasta_strings <- readLines(args[8])
  #Also Read majiq csv file of IR events
  IR_csv <- readLines(args[9])
  
  #Also Read nt seq file for IR events
  IR_nt_seq <- readLines(args[10])
  
} else
{
  #set current working directory
  setwd("/Volumes/SYEDSHAH/MichaelLab/proteogenomic_pipeline/github_final_pgp/phase-c/test_modifiedFastaID")
  #Read Peaks File
  peaks_data <- as.data.frame(read.csv("inputs/6cryptics.csv", header=TRUE))
  #Read corresponding FASTA file
  #ce_fasta_strings <- readLines("small_new_CE_FUSED_AA.fasta")
  ce_fasta_strings <- readLines("inputs/PEAKS_AA.fasta")
  
  #Also Read majiq csv file of the selected events (Here we concatenate ce_inclusion and extension csv files)
  ce_incl_ext_csv <- readLines("inputs/cds_merged_igv_ce.csv")
  
  #Also Read nt seq file ce events
  ce_incl_ext_nt_seq <- readLines("inputs/cds_merged_nt_ce.fasta")
  
  #Read skiptics FASTA file
  es_fasta_strings <- readLines("inputs/cds_PEAKS_SKIPTICS_FUSED_AA.fasta")
  #Also Read majiq csv file of the selected events (Here we concatenate ce_inclusion and extension csv files)
  es_csv <- readLines("inputs/cds_IGV_unique_skiptics_translated.csv")

  #Also Read nt seq file ce events
  es_nt_seq <- readLines("inputs/cds_skiptics_fused_transeq_in.fasta")
  
  #Read IR FASTA file
  IR_fasta_strings <- readLines("inputs/FINAL_IR_AA.fasta")
  #Also Read majiq csv file of IR events
  IR_csv <- readLines("inputs/IGV_unique_IR.csv")
  
  #Also Read nt seq file for IR events
  IR_nt_seq <- readLines("inputs/IR_coord_uniq_nt.transeq_in.fasta")
  
  
  
}
#first get rows with empty accession
empty_accession=peaks_data[which(peaks_data$Accession=="",),]
#save 
if (file.exists("empty_accession.csv")) {
  #Delete file if it exists
  file.remove("empty_accession.csv")
}

write.table(empty_accession,file=paste0("empty_accession",".csv"), quote=F, sep=",")#, row.names=F, col.names=F)
#get clean peaks data
clean_peaks_data=peaks_data[which(!peaks_data$Accession=="",),]
#now get each element
#df_notfound <- data.frame(seqnames=str(''),igv=str(''),chrom=str(''), start=numeric(),end=numeric(),genename=str(''))
#df_peptide_category_all <- data.frame(seqnames=str(''),category=str(''),genename=str(''))


if (file.exists("PEAKS_STATS.txt")) {
  #Delete file if it exists
  file.remove("PEAKS_STATS.txt")
}



if (file.exists("undecided_category.csv")) {
  #Delete file if it exists
  file.remove("undecided_category.csv")
}

if (file.exists("undecided_skiptics.csv")) {
  #Delete file if it exists
  file.remove("undecided_skiptics.csv")
}

if (file.exists("undecided_cryptics.csv")) {
  #Delete file if it exists
  file.remove("undecided_cryptics.csv")
}


if (file.exists("peptide_category.csv")) {
  #Delete file if it exists
  file.remove("peptide_category.csv")
}
if (file.exists("ce_res/us_exon_ce.csv")) {
  #Delete file if it exists
  file.remove("ce_res/us_exon_ce.csv")
}

if (file.exists("ce_res/sashimi_us_exon_ce.csv")) {
  #Delete file if it exists
  file.remove("ce_res/sashimi_us_exon_ce.csv")
}

if (file.exists("ce_res/us_exon_ce.bed")) {
  #Delete file if it exists
  file.remove("ce_res/us_exon_ce.bed")
}


if (file.exists("ce_res/ce.csv")) {
  #Delete file if it exists
  file.remove("ce_res/ce.csv")
}

if (file.exists("ce_res/sashimi_ce.csv")) {
  #Delete file if it exists
  file.remove("ce_res/sashimi_ce.csv")
}

if (file.exists("ce_res/ce.bed")) {
  #Delete file if it exists
  file.remove("ce_res/ce.bed")
}

if (file.exists("ce_res/ce_ds_exon.csv")) {
  #Delete file if it exists
  file.remove("ce_res/ce_ds_exon.csv")
}

if (file.exists("ce_res/sashimi_ce_ds_exon.csv")) {
  #Delete file if it exists
  file.remove("ce_res/sashimi_ce_ds_exon.csv")
}

if (file.exists("ce_res/ce_ds_exon.bed")) {
  #Delete file if it exists
  file.remove("ce_res/ce_ds_exon.bed")
}



if (file.exists("ce_res/ds_exon.csv")) {
  #Delete file if it exists
  file.remove("ce_res/ds_exon.csv")
}

if (file.exists("ce_res/sashimi_ds_exon.csv")) {
  #Delete file if it exists
  file.remove("ce_res/sashimi_ds_exon.csv")
}
if (file.exists("ce_res/ds_exon.bed")) {
  #Delete file if it exists
  file.remove("ce_res/ds_exon.bed")
}


if (file.exists("ce_res/us_exon.csv")) {
  #Delete file if it exists
  file.remove("ce_res/us_exon.csv")
}

if (file.exists("ce_res/us_exon.bed")) {
  #Delete file if it exists
  file.remove("ce_res/us_exon.bed")
}

if (file.exists("undecided.csv")) {
  #Delete file if it exists
  file.remove("undecided.csv")
}

if (file.exists("ce_res/up_exon_annotated.csv")) {
  #Delete file if it exists
  file.remove("ce_res/up_exon_annotated.csv")
}

if (file.exists("es_res/annotated_es.csv")) {
  #Delete file if it exists
  file.remove("es_res/annotated_es.csv")
}


if (file.exists("ce_res/sashimi_us_exon.csv")) {
  #Delete file if it exists
  file.remove("ce_res/sashimi_us_exon.csv")
}

if (file.exists("ce_res/peaks_sashimi_us_exon.csv")) {
  #Delete file if it exists
  file.remove("ce_res/peaks_sashimi_us_exon.csv")
}

if (file.exists("ce_res/peaks_sashimi_us_exon_ce.csv")) {
  #Delete file if it exists
  file.remove("ce_res/peaks_sashimi_us_exon_ce.csv")
}

if (file.exists("ce_res/peaks_sashimi_ce.csv")) {
  #Delete file if it exists
  file.remove("ce_res/peaks_sashimi_ce.csv")
}

if (file.exists("ce_res/peaks_sashimi_ce_ds_exon.csv")) {
  #Delete file if it exists
  file.remove("ce_res/peaks_sashimi_ce_ds_exon.csv")
}

if (file.exists("ce_res/peaks_sashimi_ds_exon.csv")) {
  #Delete file if it exists
  file.remove("ce_res/peaks_sashimi_ds_exon.csv")
}

if (file.exists("es_res/es_us_exon.csv")) {
  #Delete file if it exists
  file.remove("es_res/es_us_exon.csv")
}

if (file.exists("es_res/es_sashimi_us_exon.csv")) {
  #Delete file if it exists
  file.remove("es_res/es_sashimi_us_exon.csv")
}

if (file.exists("es_res/es_us_exon.bed")) {
  #Delete file if it exists
  file.remove("es_res/es_us_exon.bed")
}

if (file.exists("es_res/es_peaks_sashimi_us_exon.csv")) {
  #Delete file if it exists
  file.remove("es_res/es_peaks_sashimi_us_exon.csv")
}

if (file.exists("es_res/es_ds_exon.csv")) {
  #Delete file if it exists
  file.remove("es_res/es_ds_exon.csv")
}

if (file.exists("es_res/es_sashimi_ds_exon.csv")) {
  #Delete file if it exists
  file.remove("es_res/es_sashimi_ds_exon.csv")
}

if (file.exists("es_res/es_ds_exon.bed")) {
  #Delete file if it exists
  file.remove("es_res/es_ds_exon.bed")
}

if (file.exists("es_res/es_peaks_sashimi_ds_exon.csv")) {
  #Delete file if it exists
  file.remove("es_res/es_peaks_sashimi_ds_exon.csv")
}

if (file.exists("es_res/es_us_exon_ds_exon.csv")) {
  #Delete file if it exists
  file.remove("es_res/es_us_exon_ds_exon.csv")
}

if (file.exists("es_res/es_us_exon_ds_exon.bed")) {
  #Delete file if it exists
  file.remove("es_res/es_us_exon_ds_exon.bed")
}

if (file.exists("es_res/es_sashimi_us_exon_ds_exon.csv")) {
  #Delete file if it exists
  file.remove("es_res/es_sashimi_us_exon_ds_exon.csv")
}

if (file.exists("es_res/es_peaks_sashimi_us_exon_ds_exonn.csv")) {
  #Delete file if it exists
  file.remove("es_res/es_peaks_sashimi_us_exon_ds_exonn.csv")
}

if (file.exists("IR_res/IR_events.csv")) {
  #Delete file if it exists
  file.remove("IR_res/IR_events.csv")
}

if (file.exists("IR_res/IR_sashimi.csv")) {
  #Delete file if it exists
  file.remove("IR_res/IR_sashimi.csv")
}

if (file.exists("IR_res/IR_events.bed")) {
  #Delete file if it exists
  file.remove("IR_res/IR_events.bed")
}

if (file.exists("IR_res/IR_peaks_sashimi.csv")) {
  #Delete file if it exists
  file.remove("IR_res/IR_peaks_sashimi.csv")
}




tot=0
tot_un=0
tupex=0
tupexce=0
tce=0
tcednex=0
tdnex=0
skip_uex=0
skip_uex_dex=0
skip_dex=0
ir_cnt=0
unknow_events_cnt=0
ce_unknown_category=0
skiptic_unknown_category=0
ir_unknown_category=0

Total_Events_Read = dim(clean_peaks_data)[1]
#print(paste0('read: ', clean_peaks_data))
upex_annotated_ce_incl=0
annotated_es=0
for (i in 1:dim(clean_peaks_data)[1])
{
  #print(paste0('record: ',i, dim(clean_peaks_data)[1]))
  aa_seq=clean_peaks_data[i,]$Peptide
  print(paste0(aa_seq))
  fasta_id_all_temp=strsplit(clean_peaks_data[i,]$Accession,':')[[1]][1] #get only first fasta_id if an aa appear in more than one events
  print(paste0(fasta_id_all_temp))
  #FIRST CHECK IF 
  #FIRST CHECK WHICH CATEGORY IT BELOGS TO (CE, SKIPTICS or IR)
  #if (length(str_split(fasta_id_all, "_")[[1]])==11) #CRYPTICS (inclusion and extension)
  #if (length(str_split(fasta_id_all, "_")[[1]])==12) #CRYPTICS (inclusion and extension) #now 12 as an extra event# flag is added
  if (length(str_split(fasta_id_all_temp, "_")[[1]])==13) #CRYPTICS (inclusion and extension) #now 13 as an extra flag is added for length of aa from exon1 to upstream exon
  {
    print(paste0('CE Accession ID: ',clean_peaks_data[i,]$Accession, ' @ Line: ',i)) #we use only first fastaID from the list
    
    ############# This section retrieves aa length from exon 1 to upstream exon
    AALengthEx1_upexon=as.numeric(str_split(fasta_id_all_temp, "_")[[1]][12])
    
    ntLengthEx1_upexon =0# 3*AALengthEx1_upexon
    fasta_id_all=paste(str_split(fasta_id_all_temp, "_")[[1]][1],str_split(fasta_id_all_temp, "_")[[1]][2],str_split(fasta_id_all_temp, "_")[[1]][3],str_split(fasta_id_all_temp, "_")[[1]][4],str_split(fasta_id_all_temp, "_")[[1]][5],
                   str_split(fasta_id_all_temp, "_")[[1]][6],str_split(fasta_id_all_temp, "_")[[1]][7],str_split(fasta_id_all_temp, "_")[[1]][8],str_split(fasta_id_all_temp, "_")[[1]][9],str_split(fasta_id_all_temp, "_")[[1]][10],str_split(fasta_id_all_temp, "_")[[1]][11],str_split(fasta_id_all_temp, "_")[[1]][13],sep="_",collapse = NULL)
    
    print(paste0('New Fasta ID: ',fasta_id_all, ' and AA exo1_upexon length is: ',AALengthEx1_upexon)) #we use only first fastaID from the list
    ################### End section retrieves aa length from exon 1 to upstream exon
    
    
    #First get fasta ID
    #get event ID
    eventID=as.numeric(str_split(fasta_id_all, "_")[[1]][9])
    #get strand
    strand=paste0(str_split(fasta_id_all, "_")[[1]][10],'_') 
    
    #Verify that we have valid strand type
    if(strand !="plus_" && strand !="minus_")
    {
      print(paste0('UNKNOWN STRAND TYPE: ',strand, ' @ Line Number: ',i, ' SO ABANDONING IDENTIFICATION FOR THIS ACCESSION'))
      next
    }
    if(strand =="plus_") strnd='+'
    if(strand =="minus_") strnd='-'
    
    #frame=as.numeric(str_split(fasta_id_all, "_")[[1]][10])
    frame=as.numeric(str_split(fasta_id_all, "_")[[1]][11])
    #Also first remove _eventID part of the string
    #fasta_id_all=paste(str_split(fasta_id_all, "_")[[1]][1],str_split(fasta_id_all, "_")[[1]][2],str_split(fasta_id_all, "_")[[1]][3],str_split(fasta_id_all, "_")[[1]][4],str_split(fasta_id_all, "_")[[1]][5],
    #                    str_split(fasta_id_all, "_")[[1]][6],str_split(fasta_id_all, "_")[[1]][7],str_split(fasta_id_all, "_")[[1]][8],str_split(fasta_id_all, "_")[[1]][10],str_split(fasta_id_all, "_")[[1]][11],str_split(fasta_id_all, "_")[[1]][12],sep="_",collapse = NULL)
    
    #now retain upto frame number and remove last character as it is coming from splitting final aa seq in lengths >8 and removing *
    #fasta_id = substr(unlist(fasta_id_all)[1],1,nchar(unlist(fasta_id_all)[1])-2)
    fasta_id=paste(str_split(fasta_id_all, "_")[[1]][1],str_split(fasta_id_all, "_")[[1]][2],str_split(fasta_id_all, "_")[[1]][3],str_split(fasta_id_all, "_")[[1]][4],str_split(fasta_id_all, "_")[[1]][5],
                   str_split(fasta_id_all, "_")[[1]][6],str_split(fasta_id_all, "_")[[1]][7],str_split(fasta_id_all, "_")[[1]][8],str_split(fasta_id_all, "_")[[1]][9],str_split(fasta_id_all, "_")[[1]][10],str_split(fasta_id_all, "_")[[1]][11],sep="_",collapse = NULL)
    
    #Also get fasta ID for 
    fasta_id_nt=paste(str_split(fasta_id, "_")[[1]][1],str_split(fasta_id, "_")[[1]][2],str_split(fasta_id, "_")[[1]][3],str_split(fasta_id, "_")[[1]][4],str_split(fasta_id, "_")[[1]][5],
                      str_split(fasta_id, "_")[[1]][6],str_split(fasta_id, "_")[[1]][7],str_split(fasta_id, "_")[[1]][8],str_split(fasta_id, "_")[[1]][9],str_split(fasta_id, "_")[[1]][10],sep="_",collapse = NULL)
    
    #print(paste0('CE, Finally I have Accession ID: ',fasta_id, ' @ Line: ',i))
    #ignoring whether string matches to 1,2 or 3 frame translated sequence
    
    flag=0
    flag_event=0
    #go through fasta file to see if this aa_seq is in the fasta_file
    #k=1 #to track lines in nt file
    for(j in seq(2,length(ce_fasta_strings),2))
    {
      
      #search if aa_seq exists in line j and fasta_id exists in line i-1 if the fasta file
      #if (grepl(aa_seq, ce_fasta_strings[j],fixed=T) && grepl(fasta_id, ce_fasta_strings[j-1],fixed=T))
      if (grepl(aa_seq, ce_fasta_strings[j],fixed=T) && grepl(fasta_id_all_temp, ce_fasta_strings[j-1],fixed=T))
      {
        category='failed'
        flag=1
        #Also get the nt sequence (as we have same number of lines in both nt and AA seq fasta files)
        #if (grepl(fasta_id_nt, ce_incl_ext_nt_seq[k],fixed=T))
        #Also search for nt seq
        nt_flg=0
        for(knt in seq(2,length(ce_incl_ext_nt_seq),2))
        {
          if (grepl(fasta_id_nt, ce_incl_ext_nt_seq[knt-1],fixed=T))
          {
            nt_seq= ce_incl_ext_nt_seq[knt]
            print(paste0('CE Fatsta ID: ',fasta_id_nt, ' and nt_seq: ',nt_seq))
            #print(paste0('CE event nt seq: ',nt_seq, ' for Line: ',i, ' in AA fasta and line ',knt,' in nt fasta')) #we use only first fastaID from the list
            nt_flg=1
            break
          }
        }
        if(nt_flg==0)
        {
          print(paste0('CE Fatsta ID:',fasta_id_nt))
          print(paste0('For CE event @ line : ',i, ' fastaID: ',clean_peaks_data[i,]$Accession, ' in AA and nt fasta files do not match, so skipping')) #we use only first fastaID from the list
          break
        }
        
        #fasta_id=''
        start_end=str_locate(ce_fasta_strings[j],aa_seq)  #position of AA seq in fasta file
        
        #also decide which category it belongs to
        coord_list=strsplit(fasta_id,"_")
        gene_name1=sapply(coord_list,"[",1)
        gene_name=substr(gene_name1,4,nchar(gene_name1[1]))
        chr=sapply(coord_list,"[",2)
        
        # mul position by three to match to genomic coordinates and subtract 3 
        #start_end1 = 3*start_end-3 #position starts from 0 #+ 3 #to avoid 3 frame mismatch - should check
        start_end=start_end-1
        #start_end1=3*start_end
        start_end1=3*start_end-3*AALengthEx1_upexon
        start_end1[1] = start_end1[1]+1
        start_end1[2] = start_end1[2]+frame-1
        
        
        
        #get up/dn and ce start and end positions
        uex_st=as.numeric(sapply(coord_list,"[",3))
        uex_end=as.numeric(sapply(coord_list,"[",4))
        ce_st=as.numeric(sapply(coord_list,"[",5))
        ce_end=as.numeric(sapply(coord_list,"[",6))
        dex_st=as.numeric(sapply(coord_list,"[",7))
        dex_end=as.numeric(sapply(coord_list,"[",8))
        
        
        #fasta_id_save=paste0(chr,':',uex_st,'-',uex_end,'-',chr,':',ce_st,'-',ce_end,'-',chr,':',dex_st,'-',dex_end)
        
        #check upstream exon
        start_endu=start_end1+uex_st
        uex_rng = (uex_end-uex_st)/3
        upex_len=uex_end-uex_st
        
        start_endce=start_end1+ce_st
        ce_rng = (ce_end-ce_st)/3
        boundary2=uex_rng+ce_rng
        ce_len=ce_end-ce_st
        
        start_endd=start_end1+dex_st
        dex_rng = (dex_end-dex_st)/3
        boundary3=boundary2+dex_rng
        dex_len=dex_end-dex_st
        
        
        #also adjust for different length up and down exons for up_exon_ce and ce_ds_exon events
        if ((upex_len+ce_len)%%3==0 )
        {
          if (frame==1)
          {
            frame_offset=0
          }else if (frame==2)
          {
            frame_offset=2
          }else if (frame==3)
          {
            frame_offset=1
          }
        }else if ((upex_len+ce_len)%%3==1)
        {
          if (frame==1)
          {
            frame_offset=1
          }else if (frame==2)
          {
            frame_offset=0
          }else if (frame==3)
          {
            frame_offset=2
          }
          
        }else if ((upex_len+ce_len)%%3==2)
        {
          if (frame==1)
          {
            frame_offset=2
          }else if (frame==2)
          {
            frame_offset=1
          }else if (frame==3)
          {
            frame_offset=0
          }
          
        }
        
        ##########THIS CHECKS FOR ANNOTATED EVENTS
        if(start_end[2] <=AALengthEx1_upexon)
        {
          upex_annotated_ce_incl=upex_annotated_ce_incl+1
          print(paste0('Event: ',i,' lies in upstream exon1 annotated ce', ' and total upex_annotated_ce_incl are: ',upex_annotated_ce_incl))
          write.table(clean_peaks_data[i,],file=paste0("ce_res/up_exon_annotated",".csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        }
        #check which range this event lies in
        else if (start_end[1] >=AALengthEx1_upexon && start_end[2]<=uex_rng+AALengthEx1_upexon)  #completely in uex
        {
          flag_event=1
          tot=tot+1
          tupex=tupex+1
          category = 'us_exon'
          if(strand=="plus_")
          {
            #start_coord=start_end1[1]+uex_st#-1#-frame+1
            #end_coord=start_end1[2]+uex_st+3#+frame-1
            
            start_coord=start_end1[1]+uex_st-ntLengthEx1_upexon#-1#-frame+1
            end_coord=start_end1[2]+uex_st+3-ntLengthEx1_upexon#+frame-1
            
          }else
          {
            start_coord=uex_end-start_end1[2]-2#-frame+1 #2
            end_coord=uex_end-start_end1[1]+1
          }
          #also get nt seq
          nt1=substr(nt_seq,start_end1[1],start_end1[2]+3)
          nt1=paste0(nt1,'_',frame)
          # print(paste0('up_exon got nt: ',nt1,' and length is: ', nchar(nt1)-2))
          # print(paste0('up_exon got AA seq: ',aa_seq,' and length is: ',nchar(aa_seq)))
          
          
          
        }
        #else if(start_end[1]<=uex_rng+AALengthEx1_upexon && start_end[2]>uex_rng+AALengthEx1_upexon && start_end[2]<=boundary2+AALengthEx1_upexon) #overlaps between up_ex and ce
        else if(start_end[1]<uex_rng+AALengthEx1_upexon && start_end[2]>uex_rng+AALengthEx1_upexon && start_end[2]<=boundary2+AALengthEx1_upexon) #overlaps between up_ex and ce          
        {
          flag_event=1
          tot=tot+1
          tupexce=tupexce+1
          
          category = 'us_exon_ce'
          #WE SAVE 2 RANGES FOR us_exon_ce - 1 for us_ex and other for ce
          if(strand=="plus_")
          {
            #how much event is inside ce
            ce_in= start_end1[2] - 3*uex_rng
            
            #if (frame!=1)
            #{
            #now find start and end for complete range
            #range in us_ex - start and end
            #start_coord1=start_end1[1]+uex_st
            #end_coord1=uex_end#+frame-1##frame -1 is to compensate for frame translation

            start_coord1=start_end1[1]+uex_st-ntLengthEx1_upexon
            end_coord1=uex_end-ntLengthEx1_upexon#+frame-1##frame -1 is to compensate for frame translation
            
            
            
            #range in ce - start and end
            #start_coord2=ce_st+1#+frame##frame is to compensate for frame translation
            #end_coord2=ce_st+ce_in+3 #+frame-1 #start_end1[2]
            
            start_coord2=ce_st+1-ntLengthEx1_upexon#+frame##frame is to compensate for frame translation
            end_coord2=ce_st+ce_in+3-ntLengthEx1_upexon #+frame-1 #start_end1[2]
            
                        
            
            #also write compelet range - 
            start_coord = start_coord1
            end_coord = end_coord2
            
            #also get nt seq
            #nt1=substr(nt_seq,start_end1[1],start_end1[2]+3) 
            #nt1=paste0(nt1,'_',frame)
            #print(paste0('up_exon_ce got nt: ',nt1,' and length is: ', nchar(nt1)-2))
            #print(paste0('up_exon_ce got AA seq: ',aa_seq,' and length is: ',nchar(aa_seq)))
            
            #also get up_exon and ce nt
            #nt2=substr(nt_seq,start_end1[1],3*uex_rng+frame-1)#frame -1 is to compensate for frame translation
            nt2=substr(nt_seq,start_end1[1],3*uex_rng)#+frame-1)#frame -1 is to compensate for frame translation
            nt2=paste0(nt2,'_',frame)
            #nt3=substr(nt_seq,3*uex_rng+1+frame-1,start_end1[2]+3)
            nt3=substr(nt_seq,3*uex_rng+1,start_end1[2]+3)
            nt3=paste0(nt3,'_',frame)
            #also get aa seq
            #aa_1=substr(aa_seq,1,uex_rng-start_end[1])
            #aa_2=substr(aa_seq,uex_rng-start_end[1]+1,nchar(aa_seq))
            
            aa_1=substr(aa_seq,1,uex_rng-(start_end[1]-AALengthEx1_upexon))
            aa_2=substr(aa_seq,uex_rng-start_end[1]+1+AALengthEx1_upexon,nchar(aa_seq))
            
          }else
          {
            
            ###############
            start_coord1=uex_st+1 #-start_end1[2]-2#-frame+1 #2
            end_coord1=uex_end-start_end1[1]+1
            #get dn_ex range
            start_coord2=ce_end-(start_end1[2]-upex_len)-2#-frame+1
            end_coord2=ce_end#-(start_end1[1]-upex_len)+1
            #also write compelet range - 
            start_coord = start_coord2
            end_coord = end_coord1
            
            #also get up_exon and ce nt
            nt2=substr(nt_seq,start_end1[1],upex_len)
            nt2=paste0(nt2,'_',frame)
            nt3=substr(nt_seq,upex_len+1,start_end1[2]+3)
            nt3=paste0(nt3,'_',frame)
            #also get aa seq
            #aa_1=substr(aa_seq,1,uex_rng-start_end[1])
            #aa_2=substr(aa_seq,uex_rng-start_end[1]+1,nchar(aa_seq))

            aa_1=substr(aa_seq,1,uex_rng-(start_end[1]-AALengthEx1_upexon))
            aa_2=substr(aa_seq,uex_rng-start_end[1]+1+AALengthEx1_upexon,nchar(aa_seq))
            
            ###############
            
          }
          # #also get nt seq
          #also get nt seq
          nt1=substr(nt_seq,start_end1[1],start_end1[2]+3)
          nt1=paste0(nt1,'_',frame)
          #print(paste0('up_exon_ce got nt: ',nt1,' and length is: ', nchar(nt1)-2))
          #print(paste0('up_exon_ce got AA seq: ',aa_seq,' and length is: ',nchar(aa_seq)))
          
          #also print ntup and ntce
          
          #print(paste0('up_exon_ce got AA seq for up_exon: ',aa_1,' and length is: ',nchar(aa_1)))
          #print(paste0('up_exon_ce got nt seq for up_exon: ',nt2,' and length is: ',nchar(nt2)-2))
          #print(paste0('up_exon_ce got AA seq for ce: ',aa_2,' and length is: ',nchar(aa_2)))
          #print(paste0('up_exon_ce got nt seq for ce: ',nt3,' and length is: ',nchar(nt3)-2))
          
        }
        else if (start_end[1]>=uex_rng+AALengthEx1_upexon && start_end[2]<=boundary2+AALengthEx1_upexon) #ce region
        {
          flag_event=1
          tot=tot+1
          tce=tce+1
          
          category = 'ce'
          if(strand=="plus_")
          {
            #start_coord=(start_end1[1]-3*uex_rng)+ce_st#-1#+frame-2#2
            #end_coord=(start_end1[2]-3*uex_rng)+ce_st+3#+frame-1#-1
            

            start_coord=(start_end1[1]-3*uex_rng)+ce_st-ntLengthEx1_upexon#-1#+frame-2#2
            end_coord=(start_end1[2]-3*uex_rng)+ce_st+3-ntLengthEx1_upexon#+frame-1#-1
            
            
          }else
          {
            
            start_coord=ce_end-(start_end1[2]-3*uex_rng)-2#-frame+1 #+3-frame
            end_coord=ce_end-(start_end1[1]-3*uex_rng)+1#-frame+1 #+3-frame
          }
          #also get nt seq
          nt1=substr(nt_seq,start_end1[1],start_end1[2]+3)
          nt1=paste0(nt1,'_',frame)
          
          # print(paste0('ce got nt: ',nt1,' and length is: ', nchar(nt1)-2))
          # print(paste0('ce got AA seq: ',aa_seq,' and length is: ',nchar(aa_seq)))
          
          
        }
        else if(start_end[1]<=boundary2+AALengthEx1_upexon && start_end[2]>boundary2+AALengthEx1_upexon && start_end[2]<=boundary3+AALengthEx1_upexon) #overlaps between ce and ds_ex
        {
          flag_event=1
          tot=tot+1
          tcednex=tcednex+1
          
          category = 'ce_ds_exon'
          if(strand=="plus_")
          {
            # start_coord1=(start_end1[1]-3*uex_rng)+ce_st#+frame-1
            # end_coord1=ce_end#+frame
            # #dex_st = dex_end-60
            # start_coord2=dex_st+1#+frame#-1
            # end_coord2=dex_st+(start_end1[2]-3*uex_rng-3*ce_rng)+3#+frame
            # #also write compelet range - 
            # start_coord = start_coord1
            # end_coord = end_coord2
            # 
            # #also get up_exon and ce nt
            # nt2=substr(nt_seq,start_end1[1],3*(uex_rng+ce_rng))
            # nt2=paste0(nt2,'_',frame)
            # nt3=substr(nt_seq,3*(uex_rng+ce_rng)+1,start_end1[2]+3)
            # nt3=paste0(nt3,'_',frame)
            # #also get aa seq
            # aa_1=substr(aa_seq,1,uex_rng+ce_rng-start_end[1])
            # aa_2=substr(aa_seq,uex_rng+ce_rng-start_end[1]+1,nchar(aa_seq))
            #Also changing to lengths instead of range
            #start_coord1=(start_end1[1]-upex_len)+ce_st#+frame-1
            #end_coord1=ce_end#+frame
            
            
            start_coord1=(start_end1[1]-upex_len)+ce_st-ntLengthEx1_upexon#+frame-1
            end_coord1=ce_end-ntLengthEx1_upexon#+frame
            
            
            
            #dex_st = dex_end-60
            #start_coord2=dex_st+1#+frame#-1
            #end_coord2=dex_st+(start_end1[2]-upex_len-ce_len)+3#+frame
            
            start_coord2=dex_st+1-ntLengthEx1_upexon#+frame#-1
            end_coord2=dex_st+(start_end1[2]-upex_len-ce_len)+3-ntLengthEx1_upexon#+frame
            
            
            #also write compelet range - 
            start_coord = start_coord1
            end_coord = end_coord2
            
            #also get up_exon and ce nt
            nt2=substr(nt_seq,start_end1[1],upex_len+ce_len)
            nt2=paste0(nt2,'_',frame)
            nt3=substr(nt_seq,upex_len+ce_len+1,start_end1[2]+3)
            nt3=paste0(nt3,'_',frame)
            #also get aa seq
            aa_1=substr(aa_seq,1,uex_rng+ce_rng-(start_end[1]-AALengthEx1_upexon))
            aa_2=substr(aa_seq,uex_rng+ce_rng-start_end[1]+1+AALengthEx1_upexon,nchar(aa_seq))
            
            
            
            
          }else
          {
            start_coord1=ce_st+1
            end_coord1=ce_end-(start_end1[1]-upex_len)+1
            
            end_coord2=dex_end
            start_coord2=dex_end-(start_end1[2]-upex_len-ce_len)-2 #-frame-1
            #also write compelet range - 
            start_coord = start_coord2
            end_coord = end_coord1
            
            #also get up_exon and ce nt
            nt2=substr(nt_seq,start_end1[1],upex_len+ce_len)#+frame_offset-1)
            nt2=paste0(nt2,'_',frame)
            nt3=substr(nt_seq,upex_len+ce_len+1,start_end1[2]+3)
            nt3=paste0(nt3,'_',frame)
            
            # nt2=substr(nt_seq,start_end1[1],3*(uex_rng+ce_rng)-1)
            # nt2=paste0(nt2,'_',frame)
            # nt3=substr(nt_seq,3*(uex_rng+ce_rng),start_end1[2]+3)
            # nt3=paste0(nt3,'_',frame)
            #also get aa seq
            aa_1=substr(aa_seq,1,uex_rng+ce_rng-(start_end[1]-AALengthEx1_upexon))
            aa_2=substr(aa_seq,uex_rng+ce_rng-start_end[1]+1+AALengthEx1_upexon,nchar(aa_seq))
            
          }
          #also get nt seq
          nt1=substr(nt_seq,start_end1[1],start_end1[2]+3)
          nt1=paste0(nt1,'_',frame)
          
          #also print ntup and ntce
          #print(paste0('ce_ds_exon got nt: ',nt1,' and length is: ', nchar(nt1)-2))
          #print(paste0('ce_ds_exon got AA seq: ',aa_seq,' and length is: ',nchar(aa_seq)))
          
          #print(paste0('ce_ds_exon got AA seq for ce: ',aa_1,' and length is: ',nchar(aa_1)))
          #print(paste0('ce_ds_exon got nt seq for ce: ',nt2,' and length is: ',nchar(nt2)-2))
          #print(paste0('ce_ds_exon got AA seq for ds_exon: ',aa_2,' and length is: ',nchar(aa_2)))
          #print(paste0('ce_ds_exon got nt seq for ds_exon: ',nt3,' and length is: ',nchar(nt3)-2))
          
        }
        else if (start_end[1]>=boundary2+AALengthEx1_upexon && start_end[2]<=boundary3+AALengthEx1_upexon) #ds_ex
        {
          flag_event=1
          tot=tot+1
          tdnex=tdnex+1
          
          category ='ds_exon'
          #start_coord=start_end1[1]+dex_st
          #end_coord=start_end1[2]+dex_end
          
          if(strand=="plus_")
          {
            #if (frame!=1)
            #{
            #start_coord=start_end1[1]-3*uex_rng-3*ce_rng+dex_st
            #end_coord=start_end1[2]-3*uex_rng-3*ce_rng + dex_st+3
            
            start_coord=start_end1[1]-3*uex_rng-3*ce_rng+dex_st-ntLengthEx1_upexon
            end_coord=start_end1[2]-3*uex_rng-3*ce_rng + dex_st+3-ntLengthEx1_upexon
            
            
          }else
          {
            start_coord=dex_end-(start_end1[2]-3*uex_rng-3*ce_rng)-2#+3-frame
            end_coord=dex_end-(start_end1[1]-3*uex_rng-3*ce_rng)+1 #+3-frame
          }
          #also get nt seq
          nt1=substr(nt_seq,start_end1[1],start_end1[2]+3)
          nt1=paste0(nt1,'_',frame)
          #print(paste0('ds_exon got nt: ',nt1,' and length is: ', nchar(nt1)-2))
          #print(paste0('ds_exon got AA seq: ',aa_seq,' and length is: ',nchar(aa_seq)))
          
          
        }
        else #undecided
        {
          
          ce_unknown_category=ce_unknown_category+1
          
          print(paste0('range for this event is: ',start_end[1],'-',start_end[2], ' and bounds are: ',boundary3))
          
          print(paste('SEQ: ', aa_seq,' @ line: ',i,' has undecided category'))
          write.table(clean_peaks_data[i,],file=paste0("undecided_category",".csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
          
          write.table(clean_peaks_data[i,],file=paste0("undecided_cryptics",".csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
          
          
        }
      }
      if(flag==1)
      {break}
    }
    if(flag==1 && flag_event==1)
    {
      event_fastaID=paste0(chr,':',start_coord,'-',end_coord)
      SJdat <- data.frame(aa_seq,
                          gene_name,
                          event_fastaID,
                          fasta_id,
                          fasta_id_all,
                          category,
                          nt1)
      
      if (category == 'us_exon_ce' || category == 'ce_ds_exon')
      {
        event_fastaID=paste0(chr,':',start_coord,'-',end_coord)
        SJdat_mixed <- data.frame(aa_seq,
                                  gene_name,
                                  event_fastaID,
                                  fasta_id,
                                  fasta_id_all,
                                  category,
                                  nt1,
                                  nt2,
                                  nt3,
                                  aa_1,
                                  aa_2,
                                  chr,
                                  start_coord1,
                                  end_coord1,
                                  start_coord2,
                                  end_coord2)
        
      }
      
      SashimiDat <- data.frame(strsplit(ce_incl_ext_csv[eventID],',')[[1]][2],
                               strsplit(ce_incl_ext_csv[eventID],',')[[1]][3],
                               strsplit(ce_incl_ext_csv[eventID],',')[[1]][4],
                               strsplit(ce_incl_ext_csv[eventID],',')[[1]][7],
                               strsplit(ce_incl_ext_csv[eventID],',')[[1]][8],
                               strsplit(ce_incl_ext_csv[eventID],',')[[1]][9])
      
      
      if (category=='us_exon' || category=='ce' || category=='ds_exon')
      {
        SJdat1 <- data.frame(chr,
                             start_coord,
                             end_coord,
                             1,
                             0,
                             strnd,
                             gene_name)
      }
      if (category == 'us_exon_ce' || category == 'ce_ds_exon')
      {
        SJdat1 <- data.frame(chr,
                             start_coord1,
                             end_coord1,
                             1,
                             0,
                             strnd,
                             gene_name)
        SJdat2 <- data.frame(chr,
                             start_coord2,
                             end_coord2,
                             1,
                             0,
                             strnd,
                             gene_name)
        
        
        
        
      }
      
      #save us_ex
      if (category=='us_exon')
      {
        write.table(SJdat,file="ce_res/us_exon.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        write.table(SJdat1,file="ce_res/sashimi_us_exon.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        write.table(SJdat1,file="ce_res/us_exon.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
        #also save for sashimi plots
        write.table(SashimiDat,file="ce_res/peaks_sashimi_us_exon.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        
      }
      if (category == 'us_exon_ce')
      {
        #write.table(SJdat,file=paste0("us_exon_ce",".csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        write.table(SJdat_mixed,file="ce_res/us_exon_ce.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        #write.table(SJdat1,file=paste0("sashimi_us_exon_ce",".csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        write.table(SJdat2,file="ce_res/sashimi_us_exon_ce.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        write.table(SJdat1,file="ce_res/us_exon_ce.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
        write.table(SJdat2,file="ce_res/us_exon_ce.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
        
        #also save for sashimi plots
        write.table(SashimiDat,file="ce_res/peaks_sashimi_us_exon_ce.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        
      }
      
      if(category=='ce')
      {
        write.table(SJdat,file="ce_res/ce.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        write.table(SJdat1,file="ce_res/sashimi_ce.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        write.table(SJdat1,file="ce_res/ce.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
        #also save for sashimi plots
        write.table(SashimiDat,file="ce_res/peaks_sashimi_ce.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        
      }
      
      if(category=='ce_ds_exon')
      {
        
        #write.table(SJdat,file=paste0("ce_ds_exon",".csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        write.table(SJdat_mixed,file="ce_res/ce_ds_exon.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        
        #write.table(SJdat1,file=paste0("sashimi_ce_ds_exon",".csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        write.table(SJdat2,file="ce_res/sashimi_ce_ds_exon.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        write.table(SJdat1,file="ce_res/ce_ds_exon.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
        write.table(SJdat2,file="ce_res/ce_ds_exon.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
        #also save for sashimi plots
        write.table(SashimiDat,file="ce_res/peaks_sashimi_ce_ds_exon.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        
      }
      if(category=='ds_exon')
      {
        
        write.table(SJdat,file="ce_res/ds_exon.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        write.table(SJdat1,file="ce_res/sashimi_ds_exon.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        write.table(SJdat1,file="ce_res/ds_exon.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
        #also save for sashimi plots
        write.table(SashimiDat,file="ce_res/peaks_sashimi_ds_exon.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        
      }
      
      #save all
      write.table(SJdat,file="peptide_category.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
    }
    
    
  }
  #else if (length(str_split(fasta_id_all, "_")[[1]])==10) #SKIPTCS - was 9
  else if (length(str_split(fasta_id_all_temp, "_")[[1]])==11) #SKIPTCS - was 10 now 11 extra flag for aa length from exon 1 to upstream exon
  {
    print(paste0('SKIPTICS Accession ID: ',clean_peaks_data[i,]$Accession, ' @ Line: ',i))
    
    
    ############# This section retrieves aa length from exon 1 to upstream exon
    AALengthEx1_upexon=as.numeric(str_split(fasta_id_all_temp, "_")[[1]][10])
    fasta_id_all=paste(str_split(fasta_id_all_temp, "_")[[1]][1],str_split(fasta_id_all_temp, "_")[[1]][2],str_split(fasta_id_all_temp, "_")[[1]][3],str_split(fasta_id_all_temp, "_")[[1]][4],str_split(fasta_id_all_temp, "_")[[1]][5],
                       str_split(fasta_id_all_temp, "_")[[1]][6],str_split(fasta_id_all_temp, "_")[[1]][7],str_split(fasta_id_all_temp, "_")[[1]][8],str_split(fasta_id_all_temp, "_")[[1]][9],str_split(fasta_id_all_temp, "_")[[1]][11],sep="_",collapse = NULL)
    
    print(paste0('New Fasta ID: ',fasta_id_all, ' and AA exo1_upexon length is: ',AALengthEx1_upexon)) #we use only first fastaID from the list
    ################### End section retrieves aa length from exon 1 to upstream exon
    
    
    #First get fasta ID
    
    #get event ID
    eventID=as.numeric(str_split(fasta_id_all, "_")[[1]][7])
    #get strand
    strand=paste0(str_split(fasta_id_all, "_")[[1]][8],'_') 
    
    #Verify that we have valid strand type
    if(strand !="plus_" && strand !="minus_")
    {
      print(paste0('UNKNOWN STRAND TYPE: ',strand, ' @ Line Number: ',i, ' SO ABANDONING IDENTIFICATION FOR THIS ACCESSION'))
      next
    }
    if(strand =="plus_") strnd='+'
    if(strand =="minus_") strnd='-'
    
    frame=as.numeric(str_split(fasta_id_all, "_")[[1]][9])
    #Also first remove _eventID part of the string
    fasta_id_all=paste(str_split(fasta_id_all, "_")[[1]][1],str_split(fasta_id_all, "_")[[1]][2],str_split(fasta_id_all, "_")[[1]][3],str_split(fasta_id_all, "_")[[1]][4],str_split(fasta_id_all, "_")[[1]][5],
                       str_split(fasta_id_all, "_")[[1]][6],str_split(fasta_id_all, "_")[[1]][7],str_split(fasta_id_all, "_")[[1]][8],str_split(fasta_id_all, "_")[[1]][9],str_split(fasta_id_all, "_")[[1]][10],sep="_",collapse = NULL)
    
    #now retain upto frame number and remove last character as it is coming from splitting final aa seq in lengths >8 and removing *
    #fasta_id = substr(unlist(fasta_id_all)[1],1,nchar(unlist(fasta_id_all)[1])-2)
    fasta_id=paste(str_split(fasta_id_all, "_")[[1]][1],str_split(fasta_id_all, "_")[[1]][2],str_split(fasta_id_all, "_")[[1]][3],str_split(fasta_id_all, "_")[[1]][4],str_split(fasta_id_all, "_")[[1]][5],
                   str_split(fasta_id_all, "_")[[1]][6],str_split(fasta_id_all, "_")[[1]][7],str_split(fasta_id_all, "_")[[1]][8],sep="_",collapse = NULL)

    #Also get fasta ID for 
    fasta_id_nt=paste(str_split(fasta_id, "_")[[1]][1],str_split(fasta_id, "_")[[1]][2],str_split(fasta_id, "_")[[1]][3],str_split(fasta_id, "_")[[1]][4],str_split(fasta_id, "_")[[1]][5],
                      str_split(fasta_id, "_")[[1]][6],str_split(fasta_id, "_")[[1]][7],str_split(fasta_id, "_")[[1]][8],sep="_",collapse = NULL)
    
        
    #ignoring whether string matches to 1,2 or 3 frame translated sequence
    
    flag=0
    flag_event=0
    print(paste0('SKIPTICS, Finally I have Accession ID: ',fasta_id, ' @ Line: ',i, ' and fasta_id_nt: ',fasta_id_nt))
    
    #go through fasta file to see if this aa_seq is in the fasta_file
    for(j in seq(2,length(es_fasta_strings),2))
    {
      
      #search if aa_seq exists in line j and fasta_id exists in line i-1 if the fasta file
      if (grepl(aa_seq, es_fasta_strings[j],fixed=T) && grepl(fasta_id_all_temp, es_fasta_strings[j-1],fixed=T))
      {
        category='failed'
        flag=1
        
        #Also get the nt sequence (as we have same number of lines in both nt and AA seq fasta files)
        #if (grepl(fasta_id_nt, ce_incl_ext_nt_seq[k],fixed=T))
        #Also search for nt seq
        nt_flg=0
        for(knt in seq(2,length(es_nt_seq),2))
        {
          if (grepl(fasta_id_nt, es_nt_seq[knt-1],fixed=T))
          {
            nt_seq= es_nt_seq[knt]
            #print(paste0('ES Fatsta ID:',fasta_id_nt))
            #print(paste0('ES event nt seq: ',nt_seq, ' for Line: ',i, ' in AA fasta and line ',knt,' in nt fasta')) #we use only first fastaID from the list
            nt_flg=1
            break
          }
        }
        if(nt_flg==0)
        {
          #print(paste0('ES Fatsta ID:',fasta_id_nt))
          print(paste0('For ES event @ line : ',i, ' fastaID: ',clean_peaks_data[i,]$Accession, ' in AA and nt fasta files do not match, so skipping')) #we use only first fastaID from the list
          break
        }
        #fasta_id=''
        start_end=str_locate(es_fasta_strings[j],aa_seq)  #position of AA seq in fasta file
        
        #also decide which category it belongs to
        coord_list=strsplit(fasta_id,"_")
        gene_name1=sapply(coord_list,"[",1)
        gene_name=substr(gene_name1,4,nchar(gene_name1[1]))
        chr=sapply(coord_list,"[",2)
        
        # mul position by three to match to genomic coordinates and subtract 3 
        
        start_end=start_end-1
        start_end1=3*start_end -3*AALengthEx1_upexon
        start_end1[1] = start_end1[1]+1
        start_end1[2] = start_end1[2]+frame-1
        
        
        #start_end1 = 3*start_end-3 #position starts from 0 #+ 3 #to avoid 3 frame mismatch - should check
        #start_end[1]=start_end[1]-1
        #start_end1=3*start_end
        
        #get up/dn and ce start and end positions
        uex_st=as.numeric(sapply(coord_list,"[",3))
        uex_end=as.numeric(sapply(coord_list,"[",4))
        #ce_st=as.numeric(sapply(coord_list,"[",5))
        #ce_end=as.numeric(sapply(coord_list,"[",6))
        dex_st=as.numeric(sapply(coord_list,"[",5))
        dex_end=as.numeric(sapply(coord_list,"[",6))
        
        #fasta_id_save=paste0(chr,':',uex_st,'-',uex_end,'-',chr,':',dex_st,'-',dex_end)
        
        #check upstream exon
        start_endu=start_end1+uex_st
        uex_rng = (uex_end-uex_st)/3
        upex_len=uex_end-uex_st
        
        
        
        start_endd=start_end1+dex_st
        dex_rng = (dex_end-dex_st)/3
        dex_len=dex_end-dex_st
        
        boundary2=upex_len+dex_len #uex_rng+dex_rng
        
        #also adjust for different length up and down exons for up_exon_ce and ce_ds_exon events
        if (upex_len%%3==0 )
        {
          if (frame==1)
          {
            es_frame_offset=0
          }else if (frame==2)
          {
            es_frame_offset=2
          }else if (frame==3)
          {
            fes_rame_offset=1
          }
        }else if (upex_len%%3==1)
        {
          if (frame==1)
          {
            es_frame_offset=1
          }else if (frame==2)
          {
            es_frame_offset=0
          }else if (frame==3)
          {
            es_frame_offset=2
          }
          
        }else if (upex_len%%3==2)
        {
          if (frame==1)
          {
            es_frame_offset=2
          }else if (frame==2)
          {
            fes_rame_offset=1
          }else if (frame==3)
          {
            es_frame_offset=0
          }
          
        }
        
        ##########THIS CHECKS FOR ANNOTATED EVENTS
        if(start_end[2] <=AALengthEx1_upexon)
        {
          annotated_es=annotated_es+1
          print(paste0('Event: ',i,' lies in annotated es', ' and total annotated_es are: ',annotated_es))
          write.table(clean_peaks_data[i,],file=paste0("es_res/annotated_es",".csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        }
        #check which range this event lies in
        
        #check which range this event lies in
        else if (start_end[1] >=AALengthEx1_upexon && start_end[2]<=uex_rng+AALengthEx1_upexon)  #completely in uex
        {
          flag_event=1
          tot=tot+1
          skip_uex=skip_uex+1
          
          category = 'es_us_exon'
          if(strand=="plus_")
          {
            start_coord=start_end1[1]+uex_st#+frame-1
            end_coord=start_end1[2]+uex_st+3#+frame-1
          }else
          {
            start_coord=uex_end-start_end1[2]-2#-frame+1 #2
            end_coord=uex_end-start_end1[1]+1
          }
          #also get nt seq
          nt1=substr(nt_seq,start_end1[1],start_end1[2]+3)
          nt1=paste0(nt1,'_',frame)
          
        }
        else if (start_end[1] <uex_rng+AALengthEx1_upexon && start_end[2]>uex_rng+AALengthEx1_upexon && start_end[2]<=boundary2+AALengthEx1_upexon)  # both up_ex-ds_ex
        {
          flag_event=1
          tot=tot+1
          skip_uex_dex=skip_uex_dex+1
          
          
          category = 'es_us_exon_ds_exon'
          if(strand=="plus_")
          {
            #range in us_ex
            start_coord1=start_end1[1]+uex_st#+frame-1
            end_coord1=uex_end#+es_frame_offset-1#frame-3
            #range in ds_ex
            start_coord2=dex_st+1#+es_frame_offset#frame-2#+frame
            end_coord2=dex_st+start_end1[2]-upex_len+3#+frame-1
            #also write compelet range - 
            start_coord = start_coord1
            end_coord = end_coord2
            
            #also get up_exon and ce nt
            
            nt2=substr(nt_seq,start_end1[1],upex_len)#frame-3)#+frame-1)#frame -1 is to compensate for frame translation
            nt2=paste0(nt2,'_',frame)
            #nt3=substr(nt_seq,3*uex_rng+1+frame-1,start_end1[2]+3)
            #nt3=substr(nt_seq,upex_len+frame-2,start_end1[2]+3)
            nt3=substr(nt_seq,upex_len+1,start_end1[2]+3)
            nt3=paste0(nt3,'_',frame)
            #also get aa seq
            aa_1=substr(aa_seq,1,uex_rng-(start_end[1]-AALengthEx1_upexon))
            aa_2=substr(aa_seq,uex_rng-start_end[1]+1+AALengthEx1_upexon,nchar(aa_seq))
            
            

          }else
          {
            #get up_ex range
            #start_coord1=uex_end-start_end1[2]-2#-frame+1 #2
            start_coord1=uex_st+1 #-start_end1[2]-2#-frame+1 #2
            end_coord1=uex_end-start_end1[1]+1
            #get dn_ex range
            start_coord2=dex_end-(start_end1[2]-upex_len)-2#-frame+1
            end_coord2=dex_end#-(start_end1[1]-upex_len)+1
            #also write compelet range - 
            start_coord = start_coord2
            end_coord = end_coord1
            
            #also get up_exon and ce nt
            nt2=substr(nt_seq,start_end1[1],upex_len)
            nt2=paste0(nt2,'_',frame)
            nt3=substr(nt_seq,upex_len+1,start_end1[2]+3)
            nt3=paste0(nt3,'_',frame)
            #also get aa seq
            aa_1=substr(aa_seq,1,uex_rng-(start_end[1]-AALengthEx1_upexon))
            aa_2=substr(aa_seq,uex_rng-start_end[1]+1+AALengthEx1_upexon,nchar(aa_seq))
            
          }
          #also get nt seq
          nt1=substr(nt_seq,start_end1[1],start_end1[2]+3)
          nt1=paste0(nt1,'_',frame)
          
        }
        
        else if (start_end[1] >uex_rng+AALengthEx1_upexon && start_end[2]<=boundary2+AALengthEx1_upexon)  #completely in dex
        {
          flag_event=1
          tot=tot+1
          skip_dex=skip_dex+1
          category = 'es_ds_exon'
          if(strand=="plus_")
          {
            # start_coord=start_end1[1]-3*uex_rng+dex_st+frame-1
            # end_coord=start_end1[2]-3*uex_rng+dex_st+frame-1
            start_coord=start_end1[1]-upex_len+dex_st#+frame-1
            end_coord=start_end1[2]-upex_len+dex_st+3#+frame-1
            
          }else
          {
            start_coord=dex_end-(start_end1[2]-upex_len)-2#-frame+1
            end_coord=dex_end-(start_end1[1]-upex_len)+1
          }
          #also get nt seq
          nt1=substr(nt_seq,start_end1[1],start_end1[2]+3)
          nt1=paste0(nt1,'_',frame)
          
        }
        else #undecided
        {
          
          
          print(paste0('range for this event is: ',start_end[1],'-',start_end[2], ' and bounds are: ',boundary2))
          skiptic_unknown_category=skiptic_unknown_category+1
          print(paste('Skiptic SEQ: ', aa_seq,' @ line: ',i,' has undecided category'))
          write.table(clean_peaks_data[i,],file=paste0("undecided_category",".csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
          write.table(clean_peaks_data[i,],file=paste0("undecided_skiptics",".csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        }
        
        
      }#IF FOUND IN FASTA FILE
      if(flag==1)
      {break}
      
    }#FATSA FILE LOOP
    if(flag==1 && flag_event==1)
    {
      event_fastaID=paste0(chr,':',start_coord,'-',end_coord)
      SJdat <- data.frame(aa_seq,
                          gene_name,
                          event_fastaID,
                          fasta_id,
                          fasta_id_all,
                          category,
                          nt1)
      SashimiDat <- data.frame(strsplit(es_csv[eventID],',')[[1]][2],
                               strsplit(es_csv[eventID],',')[[1]][3],
                               strsplit(es_csv[eventID],',')[[1]][4],
                               strsplit(es_csv[eventID],',')[[1]][7],
                               strsplit(es_csv[eventID],',')[[1]][8],
                               strsplit(es_csv[eventID],',')[[1]][9])
      
      

      if (category=='es_us_exon' || category=='es_ds_exon')
      {
        SJdat1 <- data.frame(chr,
                             start_coord,
                             end_coord,
                             1,
                             0,
                             strnd,
                             gene_name)
      }
      if (category == 'es_us_exon_ds_exon' )
      {
        SJdat1 <- data.frame(chr,
                             start_coord1,
                             end_coord1,
                             1,
                             0,
                             strnd,
                             gene_name)
        SJdat2 <- data.frame(chr,
                             start_coord2,
                             end_coord2,
                             1,
                             0,
                             strnd,
                             gene_name)
        
        
        event_fastaID=paste0(chr,':',start_coord,'-',end_coord)
        SJdat_mixed <- data.frame(aa_seq,
                                  gene_name,
                                  event_fastaID,
                                  fasta_id,
                                  fasta_id_all,
                                  category,
                                  nt1,
                                  nt2,
                                  nt3,
                                  aa_1,
                                  aa_2,
                                  chr,
                                  start_coord1,
                                  end_coord1,
                                  start_coord2,
                                  end_coord2)
        
        
      }
    #}
    #save us_ex
    if (category=='es_us_exon')
    {
      write.table(SJdat,file="es_res/es_us_exon.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
      write.table(SJdat1,file="es_res/es_sashimi_us_exon.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
      write.table(SJdat1,file="es_res/es_us_exon.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
      #also save for sashimi plots
      write.table(SashimiDat,file="es_res/es_peaks_sashimi_us_exon.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
      
      
    }
    if(category=='es_ds_exon')
    {
      
      write.table(SJdat,file="es_res/es_ds_exon.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
      write.table(SJdat1,file="es_res/es_sashimi_ds_exon.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
      write.table(SJdat1,file="es_res/es_ds_exon.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
      #also save for sashimi plots
      write.table(SashimiDat,file="es_res/es_peaks_sashimi_ds_exon.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
      
    }
    if(category=='es_us_exon_ds_exon')
    {
      write.table(SJdat_mixed,file="es_res/es_us_exon_ds_exon.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
      
      
      write.table(SJdat2,file="es_res/es_sashimi_us_exon_ds_exon.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
      write.table(SJdat1,file="es_res/es_us_exon_ds_exon.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
      write.table(SJdat2,file="es_res/es_us_exon_ds_exon.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
      #also save for sashimi plots
      write.table(SashimiDat,file="es_res/es_peaks_sashimi_us_exon_ds_exonn.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
      
    }
    #save all
    write.table(SJdat,file=paste0("peptide_category",".csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
    }
    
  }
  else if (length(str_split(fasta_id_all_temp, "_")[[1]])==8) #IR was 7
  {
    print(paste0('IR Accession ID: ',clean_peaks_data[i,]$Accession, ' @ Line: ',i))
    fasta_id_all = fasta_id_all_temp
    eventID=as.numeric(str_split(fasta_id_all, "_")[[1]][5])
    #get strand
    strand=paste0(str_split(fasta_id_all, "_")[[1]][6],'_') 
    
    #Verify that we have valid strand type
    if(strand !="plus_" && strand !="minus_")
    {
      print(paste0('UNKNOWN STRAND TYPE: ',strand, ' @ Line Number: ',i, ' IN SKIPTICS, SO ABANDONING IDENTIFICATION FOR THIS ACCESSION'))
      next
    }
    if(strand =="plus_") strnd='+'
    if(strand =="minus_") strnd='-'
    
    frame=as.numeric(str_split(fasta_id_all, "_")[[1]][7])
    #Also first remove _eventID part of the string
    fasta_id_all=paste(str_split(fasta_id_all, "_")[[1]][1],str_split(fasta_id_all, "_")[[1]][2],str_split(fasta_id_all, "_")[[1]][3],str_split(fasta_id_all, "_")[[1]][4],str_split(fasta_id_all, "_")[[1]][5],
                       str_split(fasta_id_all, "_")[[1]][6],str_split(fasta_id_all, "_")[[1]][7],str_split(fasta_id_all, "_")[[1]][8],sep="_",collapse = NULL)
    
    #now retain upto frame number and remove last character as it is coming from splitting final aa seq in lengths >8 and removing *
    #fasta_id = substr(unlist(fasta_id_all)[1],1,nchar(unlist(fasta_id_all)[1])-2)
    fasta_id=paste(str_split(fasta_id_all, "_")[[1]][1],str_split(fasta_id_all, "_")[[1]][2],str_split(fasta_id_all, "_")[[1]][3],str_split(fasta_id_all, "_")[[1]][4],str_split(fasta_id_all, "_")[[1]][5],
                   str_split(fasta_id_all, "_")[[1]][6],str_split(fasta_id_all, "_")[[1]][7],sep="_",collapse = NULL)
    
    #Also get fasta ID for 
    fasta_id_nt=paste(str_split(fasta_id, "_")[[1]][1],str_split(fasta_id, "_")[[1]][2],str_split(fasta_id, "_")[[1]][3],str_split(fasta_id, "_")[[1]][4],str_split(fasta_id, "_")[[1]][5],
                      str_split(fasta_id, "_")[[1]][6],sep="_",collapse = NULL)
    
    
    #ignoring whether string matches to 1,2 or 3 frame translated sequence
    
    flag=0
    flag_event=0
    print(paste0('IR, Finally I have Accession ID: ',fasta_id, ' @ Line: ',i))
    
    #go through fasta file to see if this aa_seq is in the fasta_file
    for(j in seq(2,length(IR_fasta_strings),2))
    {
      
      #search if aa_seq exists in line j and fasta_id exists in line i-1 if the fasta file
      if (grepl(aa_seq, IR_fasta_strings[j],fixed=T) && grepl(fasta_id, IR_fasta_strings[j-1],fixed=T))
      {
        category='failed'
        flag=1
        
        #Also get the nt sequence (as we have same number of lines in both nt and AA seq fasta files)
        #if (grepl(fasta_id_nt, ce_incl_ext_nt_seq[k],fixed=T))
        #Also search for nt seq
        nt_flg=0
        for(knt in seq(2,length(IR_nt_seq),2))
        {
          if (grepl(fasta_id_nt, IR_nt_seq[knt-1],fixed=T))
          {
            nt_seq= IR_nt_seq[knt]
            print(paste0('IR Fatsta ID:',fasta_id_nt))
            #print(paste0('IR event nt seq: ',nt_seq, ' for Line: ',i, ' in AA fasta and line ',knt,' in nt fasta')) #we use only first fastaID from the list
            nt_flg=1
            break
          }
        }
        if(nt_flg==0)
        {
          print(paste0('IR Fatsta ID:',fasta_id_nt))
          print(paste0('For IR event @ line : ',i, ' fastaID: ',clean_peaks_data[i,]$Accession, ' in AA and nt fasta files do not match, so skipping')) #we use only first fastaID from the list
          break
        }
        #fasta_id=''
        start_end=str_locate(IR_fasta_strings[j],aa_seq)  #position of AA seq in fasta file
        
        #also decide which category it belongs to
        coord_list=strsplit(fasta_id,"_")
        gene_name1=sapply(coord_list,"[",1)
        gene_name=substr(gene_name1,4,nchar(gene_name1[1]))
        chr=sapply(coord_list,"[",2)
        
        # mul position by three to match to genomic coordinates and subtract 3 
        
        start_end=start_end-1
        start_end1=3*start_end
        start_end1[1] = start_end1[1]+1
        start_end1[2] = start_end1[2]+frame-1
        
        #get up/dn and ce start and end positions
        uex_st=as.numeric(sapply(coord_list,"[",3))
        uex_end=as.numeric(sapply(coord_list,"[",4))
        upex_len=(uex_end-uex_st)/3
        #print(paste0('I have '))
        if (start_end[1] >=0 && start_end[2]<=upex_len)  #completely in uex
        {
          flag_event=1
          tot=tot+1
          ir_cnt=ir_cnt+1
          category='IR'
          
          if(strand=="plus_")
          {
            start_coord=start_end1[1]+uex_st#+frame-1
            end_coord=start_end1[2]+uex_st+3#+frame-1
          }else
          {
            start_coord=uex_end-start_end1[2]-2#-frame+1 #2
            end_coord=uex_end-start_end1[1]+1
          }
          #also get nt seq
          nt1=substr(nt_seq,start_end1[1],start_end1[2]+3)
          nt1=paste0(nt1,'_',frame)
        }
        else #undecided
        {
          
          ir_unknown_category=ir_unknown_category+1
          print(paste('IR SEQ: ', aa_seq,' @ line: ',i,' has undecided category'))
          write.table(clean_peaks_data[i,],file=paste0("undecided_category",".csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
        }
      }
      if(flag==1)
      {break}
    }
    if(flag==1 && flag_event==1)
    {
      event_fastaID=paste0(chr,':',start_coord,'-',end_coord)
      SJdat <- data.frame(aa_seq,
                          gene_name,
                          event_fastaID,
                          fasta_id,
                          fasta_id_all,
                          category,
                          nt1)
      SashimiDat <- data.frame(#strsplit(IR_csv[eventID],',')[[1]][1],
                               strsplit(IR_csv[eventID],',')[[1]][2],
                               strsplit(IR_csv[eventID],',')[[1]][3],
                               strsplit(IR_csv[eventID],',')[[1]][4],
                               strsplit(IR_csv[eventID],',')[[1]][5],
                               strsplit(IR_csv[eventID],',')[[1]][6],
                               strsplit(IR_csv[eventID],',')[[1]][7])
      SJdat1 <- data.frame(chr,
                           start_coord,
                           end_coord,
                           1,
                           0,
                           strnd,
                           gene_name)
      
    
      write.table(SJdat,file="IR_res/IR_events.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
      write.table(SJdat1,file="IR_res/IR_sashimi.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
      write.table(SJdat1,file="IR_res/IR_events.bed",append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
      #also save for sashimi plots
      write.table(SashimiDat,file="IR_res/IR_peaks_sashimi.csv",append = TRUE, quote=F, sep=",", row.names=F, col.names=F)

      #save all
      write.table(SJdat,file=paste0("peptide_category",".csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
      
    }
    
  }
  else  #UNKNOWN CATEGORY
  {
    unknow_events_cnt=unknow_events_cnt+1
    print(paste0('Invalid Accession ID: ',clean_peaks_data[i,]$Accession, ' @ Line: ',i))
  }
  
  
}
#Finally sort all file on gene_name
if (file.exists("ce_res/us_exon.csv")) {
  dat1 <- read.table("ce_res/us_exon.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V2) , ]
  write.table(dat2,file="ce_res/us_exon.csv", quote=F, sep=",", row.names=F, col.names=F)
  dat1 <- read.table("ce_res/peaks_sashimi_us_exon.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="ce_res/peaks_sashimi_us_exon.csv", quote=F, sep=",", row.names=F, col.names=F)

  dat1 <- read.table("ce_res/sashimi_us_exon.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="ce_res/sashimi_us_exon.csv", quote=F, sep=",", row.names=F, col.names=F)
  
  
  dat1 <- read.table("ce_res/us_exon.bed", header=F, sep="\t")
  dat2 <- dat1[order(dat1$V7) , ]
  write.table(dat2,file="ce_res/us_exon.bed", quote=F, sep="\t", row.names=F, col.names=F)
}
if (file.exists("ce_res/us_exon_ce.csv")) {
  dat1 <- read.table("ce_res/us_exon_ce.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V2) , ]
  write.table(dat2,file="ce_res/us_exon_ce.csv", quote=F, sep=",", row.names=F, col.names=F)
  dat1 <- read.table("ce_res/peaks_sashimi_us_exon_ce.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="ce_res/peaks_sashimi_us_exon_ce.csv", quote=F, sep=",", row.names=F, col.names=F)

  dat1 <- read.table("ce_res/sashimi_us_exon_ce.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="ce_res/sashimi_us_exon_ce.csv", quote=F, sep=",", row.names=F, col.names=F)
  
  
    dat1 <- read.table("ce_res/us_exon_ce.bed", header=F, sep="\t")
  dat2 <- dat1[order(dat1$V7) , ]
  write.table(dat2,file="ce_res/us_exon_ce.bed", quote=F, sep="\t", row.names=F, col.names=F)
}
if (file.exists("ce_res/ce.csv")) {
  dat1 <- read.table("ce_res/ce.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V2) , ]
  write.table(dat2,file="ce_res/ce.csv", quote=F, sep=",", row.names=F, col.names=F)
  dat1 <- read.table("ce_res/peaks_sashimi_ce.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="ce_res/peaks_sashimi_ce.csv", quote=F, sep=",", row.names=F, col.names=F)

  dat1 <- read.table("ce_res/sashimi_ce.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="ce_res/sashimi_ce.csv", quote=F, sep=",", row.names=F, col.names=F)
  
  
    dat1 <- read.table("ce_res/ce.bed", header=F, sep="\t")
  dat2 <- dat1[order(dat1$V7) , ]
  write.table(dat2,file="ce_res/ce.bed", quote=F, sep="\t", row.names=F, col.names=F)
}
if (file.exists("ce_res/ce_ds_exon.csv")) {
  dat1 <- read.table("ce_res/ce_ds_exon.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V2) , ]
  write.table(dat2,file="ce_res/ce_ds_exon.csv", quote=F, sep=",", row.names=F, col.names=F)
  dat1 <- read.table("ce_res/peaks_sashimi_ce_ds_exon.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="ce_res/peaks_sashimi_ce_ds_exon.csv", quote=F, sep=",", row.names=F, col.names=F)

  dat1 <- read.table("ce_res/sashimi_ce_ds_exon.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="ce_res/sashimi_ce_ds_exon.csv", quote=F, sep=",", row.names=F, col.names=F)
  
  
    dat1 <- read.table("ce_res/ce_ds_exon.bed", header=F, sep="\t")
  dat2 <- dat1[order(dat1$V7) , ]
  write.table(dat2,file="ce_res/ce_ds_exon.bed", quote=F, sep="\t", row.names=F, col.names=F)
}
if (file.exists("ce_res/ds_exon.csv")) {
  dat1 <- read.table("ce_res/ds_exon.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V2) , ]
  write.table(dat2,file="ce_res/ds_exon.csv", quote=F, sep=",", row.names=F, col.names=F)
  dat1 <- read.table("ce_res/peaks_sashimi_ds_exon.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="ce_res/peaks_sashimi_ds_exon.csv", quote=F, sep=",", row.names=F, col.names=F)

  dat1 <- read.table("ce_res/sashimi_ds_exon.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="ce_res/sashimi_ds_exon.csv", quote=F, sep=",", row.names=F, col.names=F)
  
  
  
    dat1 <- read.table("ce_res/ds_exon.bed", header=F, sep="\t")
  dat2 <- dat1[order(dat1$V7) , ]
  write.table(dat2,file="ce_res/ds_exon.bed", quote=F, sep="\t", row.names=F, col.names=F)
}

if (file.exists("es_res/es_us_exon.csv")) {
  dat1 <- read.table("es_res/es_us_exon.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V2) , ]
  write.table(dat2,file="es_res/es_us_exon.csv", quote=F, sep=",", row.names=F, col.names=F)
  dat1 <- read.table("es_res/es_peaks_sashimi_us_exon.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="es_res/es_peaks_sashimi_us_exon.csv", quote=F, sep=",", row.names=F, col.names=F)

  dat1 <- read.table("es_res/es_sashimi_us_exon.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="es_res/es_sashimi_us_exon.csv", quote=F, sep=",", row.names=F, col.names=F)
  
  
  dat1 <- read.table("es_res/es_us_exon.bed", header=F, sep="\t")
  dat2 <- dat1[order(dat1$V7) , ]
  write.table(dat2,file="es_res/es_us_exon.bed", quote=F, sep="\t", row.names=F, col.names=F)
}

if (file.exists("es_res/es_ds_exon.csv")) {
  dat1 <- read.table("es_res/es_ds_exon.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V2) , ]
  write.table(dat2,file="es_res/es_ds_exon.csv", quote=F, sep=",", row.names=F, col.names=F)
  dat1 <- read.table("es_res/es_peaks_sashimi_ds_exon.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="es_res/es_peaks_sashimi_ds_exon.csv", quote=F, sep=",", row.names=F, col.names=F)
  
  dat1 <- read.table("es_res/es_sashimi_ds_exon.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="es_res/es_sashimi_ds_exon.csv", quote=F, sep=",", row.names=F, col.names=F)
  
  
  dat1 <- read.table("es_res/es_ds_exon.bed", header=F, sep="\t")
  dat2 <- dat1[order(dat1$V7) , ]
  write.table(dat2,file="es_res/es_ds_exon.bed", quote=F, sep="\t", row.names=F, col.names=F)
}

if (file.exists("es_res/es_us_exon_ds_exon.csv")) {
  dat1 <- read.table("es_res/es_us_exon_ds_exon.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V2) , ]
  write.table(dat2,file="es_res/es_us_exon_ds_exon.csv", quote=F, sep=",", row.names=F, col.names=F)
  dat1 <- read.table("es_res/es_peaks_sashimi_us_exon_ds_exonn.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="es_res/es_peaks_sashimi_us_exon_ds_exonn.csv", quote=F, sep=",", row.names=F, col.names=F)
  
  dat1 <- read.table("es_res/es_sashimi_us_exon_ds_exon.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="es_res/es_sashimi_us_exon_ds_exon.csv", quote=F, sep=",", row.names=F, col.names=F)
  
  
  dat1 <- read.table("es_res/es_us_exon_ds_exon.bed", header=F, sep="\t")
  dat2 <- dat1[order(dat1$V7) , ]
  write.table(dat2,file="es_res/es_us_exon_ds_exon.bed", quote=F, sep="\t", row.names=F, col.names=F)
}

if (file.exists("IR_res/IR_events.csv")) {
  dat1 <- read.table("IR_res/IR_events.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V2) , ]
  write.table(dat2,file="IR_res/IR_events.csv", quote=F, sep=",", row.names=F, col.names=F)
  dat1 <- read.table("IR_res/IR_peaks_sashimi.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="IR_res/IR_peaks_sashimi.csv", quote=F, sep=",", row.names=F, col.names=F)
  
  dat1 <- read.table("IR_res/IR_sashimi.csv", header=F, sep=",")
  dat2 <- dat1[order(dat1$V5) , ]
  write.table(dat2,file="IR_res/IR_sashimi.csv", quote=F, sep=",", row.names=F, col.names=F)
  
  
  dat1 <- read.table("IR_res/IR_events.bed", header=F, sep="\t")
  dat2 <- dat1[order(dat1$V7) , ]
  write.table(dat2,file="IR_res/IR_events.bed", quote=F, sep="\t", row.names=F, col.names=F)
}




print(paste0('Total Input Events: ',Total_Events_Read))
write(paste0('Total Input Events: ',Total_Events_Read),file="PEAKS_STATS.txt",append=TRUE)
print(paste0('Total Unknown Events: ',unknow_events_cnt))
write(paste0('Total Unknown Events: ',unknow_events_cnt),file="PEAKS_STATS.txt",append=TRUE)

print(paste0('All Undecided Events Are: ',ce_unknown_category+skiptic_unknown_category+ir_unknown_category))
write(paste0('All Undecided Events Are: ',ce_unknown_category+skiptic_unknown_category+ir_unknown_category, ' Please see undecided_category.csv file'),file="PEAKS_STATS.txt",append=TRUE)


print(paste0('total upex_annotated_ce_incl: ',upex_annotated_ce_incl))

print(paste0('Total Undecided CE Events: ',ce_unknown_category, ' Please see undecided_cryptics'))
write(paste0('Total Undecided CE Events: ',ce_unknown_category),file="PEAKS_STATS.txt",append=TRUE)

print(paste0('Total Undecided SKIPTICS Events: ',skiptic_unknown_category, ' Please see undecided_skiptics.csv'))
write(paste0('Total Undecided SKIPTICS Events: ',skiptic_unknown_category),file="PEAKS_STATS.txt",append=TRUE)

print(paste0('Total Undecided IR Events: ',ir_unknown_category))
write(paste0('Total Undecided IR Events: ',ir_unknown_category),file="PEAKS_STATS.txt",append=TRUE)




print(paste0('Finally Total Processed Events are: ',tot))
write(paste0('Finally Total Processed Events are: ',tot),file="PEAKS_STATS.txt",append=TRUE)
print(paste0('Finally Total UNProcessed Events are: ',tot_un))
write(paste0('Finally Total UNProcessed Events are: ',tot_un),file="PEAKS_STATS.txt",append=TRUE)

print(paste0('Total CRYPTICS PROCESSED are: ',tupex+tupexce+tce+tcednex+tdnex))
write(paste0('Total CRYPTICS PROCESSED are: ',tupex+tupexce+tce+tcednex+tdnex),file="PEAKS_STATS.txt",append=TRUE)

print(paste0('Total up_exon (in ce) Events are: ',tupex))
write(paste0('Total up_exon (in ce) Events are: ',tupex),file="PEAKS_STATS.txt",append=TRUE)

print(paste0('Total up_exon_ce Events are: ',tupexce))
write(paste0('Total up_exon_ce Events are: ',tupexce),file="PEAKS_STATS.txt",append=TRUE)

print(paste0('Total ce Events are: ',tce))
write(paste0('Total ce Events are: ',tce),file="PEAKS_STATS.txt",append=TRUE)
print(paste0('Total ce_dn_exon Events are: ',tcednex))
write(paste0('Total ce_dn_exon Events are: ',tcednex),file="PEAKS_STATS.txt",append=TRUE)

print(paste0('Total dn_exon (in ce) Events are: ',tdnex))
write(paste0('Total dn_exon (in ce) Events are: ',tdnex),file="PEAKS_STATS.txt",append=TRUE)

print(paste0('Total SKIPTICS PROCESSED are: ',skip_uex+skip_uex_dex+skip_dex))
write(paste0('Total SKIPTICS PROCESSED are: ',skip_uex+skip_uex_dex+skip_dex),file="PEAKS_STATS.txt",append=TRUE)

print(paste0('Total up_exon (in skiptics) Events are: ',skip_uex))
write(paste0('Total up_exon (in skiptics) Events are: ',skip_uex),file="PEAKS_STATS.txt",append=TRUE)

print(paste0('Total up_exon_dn_exon (in skiptics) Events are: ',skip_uex_dex))
write(paste0('Total up_exon_dn_exon (in skiptics) Events are: ',skip_uex_dex),file="PEAKS_STATS.txt",append=TRUE)

print(paste0('Total dn_exon (in skiptics) Events are: ',skip_dex))
write(paste0('Total dn_exon (in skiptics) Events are: ',skip_dex),file="PEAKS_STATS.txt",append=TRUE)

print(paste0('Total IR Events are: ',ir_cnt))
write(paste0('Total IR Events are: ',ir_cnt),file="PEAKS_STATS.txt",append=TRUE)


print(paste0('ALL CRYPTICS+SKIPTICS+IR Events are: ',tupex+tupexce+tce+tcednex+tdnex+skip_uex+skip_uex_dex+skip_dex+ir_cnt))
write(paste0('ALL CRYPTICS+SKIPTICS+IR Events are: ',tupex+tupexce+tce+tcednex+tdnex+skip_uex+skip_uex_dex+skip_dex+ir_cnt),file="PEAKS_STATS.txt",append=TRUE)
