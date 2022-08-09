library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(data.table)
library(GenomicFeatures)
library(plyranges)
library(msa)
library(annotables)
library(GeneStructureTools)
## 'cds_by_tx' must be a GRangesList object obtained with cdsBy(txdb, by="tx")
addCdsPhase <- function(cds_by_tx)
{
  cds_phase <- pc(rep(IntegerList(0), length(cds_by_tx)),
                  heads((3L - (cumsum(width(cds_by_tx)) %% 3L)) %% 3L, n=-1L))
  unlisted_cds_by_tx <- unlist(cds_by_tx, use.names=FALSE)
  mcols(unlisted_cds_by_tx)$cds_phase <- unlist(cds_phase, use.names=FALSE)
  relist(unlisted_cds_by_tx, cds_by_tx)
}


my_orf = function (transcripts, BSgenome = NULL, returnLongestOnly = TRUE, 
                   allFrames = FALSE,  
                   uORFs = FALSE) 
{
  if (allFrames == TRUE) {
    returnLongestOnly = FALSE
    longest = 1
  }
  
  transcripts$exon_number <- as.numeric(transcripts$exon_number)
  order <- order(transcripts$transcript_id, transcripts$exon_number)
  transcripts <- transcripts[order]
  transcripts$seq <- as.character(Biostrings::getSeq(BSgenome, 
                                                     transcripts))
  seqCat <- aggregate(seq ~ transcript_id, mcols(transcripts), 
                      function(x) (paste(x, collapse = "")))
  ids <- as.character(seqCat$transcript_id)
  
  seqCat <- seqCat$seq
  rm <- which(grepl("N", seqCat))
  if (length(rm) > 0) {
    seqCat <- seqCat[-rm]
    removeId <- ids[rm]
    ids <- ids[-rm]
    transcripts <- transcripts[-which(transcripts$transcript_id %in% 
                                        removeId)]
  }
  seqCat <- c(seqCat, stringr::str_sub(seqCat, 2), stringr::str_sub(seqCat, 
                                                                    3))
  frames <- rep(c(1, 2, 3), each = length(ids))
  ids <- c(ids, ids, ids)
  orf <- suppressWarnings(unlist(lapply(seqCat, function(x) as.character(Biostrings::translate(Biostrings::DNAString(x))))))
  orfDF <- data.frame(id = ids, aa_sequence = orf, frame = frames, 
                      stringsAsFactors = FALSE)
  orfDF$seq_length <- nchar(orfDF$aa_sequence)
  orfDF$seq_length_nt <- nchar(seqCat) + orfDF$frame - 1
  
  
  return(orfDF)
}




txdb = transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
#setwd("/Volumes/SYEDSHAH/MichaelLab/proteogenomic_pipeline/github_final_pgp/phase-b/ce_incl/v3")
#cryptic_regions = fread("res_ce_0.15/ce_inclusion_coord_only.txt",sep = "\t")
#folder='res_ce_all'
#folder="res_ce_all"
args = commandArgs(trailingOnly=TRUE)
print(paste0('get file: ',args[1]))
cryptic_regions = fread(args[1],sep = "\t")
folder=str_split_fixed(args[1],"/",2)[1,1]

  
  
#write(paste0('                                 '),file="res_skiptics/FINAL_STATS_SKIPTICS.txt",append=TRUE)
write(paste0('                                 '),file=args[2],append=TRUE)
write(paste0('From get_orf_cds.R: -------------------- Starting processing file: ',args[1]),file=args[2],append=TRUE)

cryptic_regions_protein_coding = cryptic_regions %>% 
  mutate(V2 = ifelse(V6 == "-",V2 + 1,V2)) %>% 
  mutate(V2 = ifelse(V6 == "+",V2 + 1,V2)) %>% 
  left_join(annotables::grch38_tx2gene, 
            by = c("V8" = "enstxp")) %>% 
  left_join(annotables::grch38 %>% dplyr::select(ensgene,biotype)) %>% 
  left_join(annotables::grch38 %>% dplyr::select(symbol,biotype),
            by = c("V7" = "symbol")) %>% 
  mutate(biotype_full = case_when(is.na(biotype.y) ~ biotype.x,
                                  is.na(biotype.x) ~ biotype.y,
                                  is.na(biotype.x) & is.na(biotype.y) ~ NA_character_,
                                  TRUE ~ biotype.y)) %>% 
  dplyr::select(-biotype.x,-biotype.y) %>% 
  dplyr::filter(biotype_full == 'protein_coding')

cryptic_regions_protein_coding = cryptic_regions_protein_coding %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                           seqnames.field = 'V1',
                           start.field = 'V2',
                           end.field = 'V3',
                           strand.field = 'V6')

cryptic_regions_protein_coding$transcript_id = cryptic_regions_protein_coding$V8

cds_regions = cdsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "tx",use.names = TRUE)
cds_regions = unlist(cds_regions)

cds_regions$transcript_id = gsub("\\..*", "", names(cds_regions))

#ALSO TRYING EnsDB
#Also get from EnsDB 103
#library(AnnotationHub)
#ah <- AnnotationHub()
#edb=query(ah, c("EnsDb", "Hsapiens", "103"))[[1]]

#Txs <- as_tibble(as.list(transcripts(edb,
#                                     columns = c(listColumns(edb , "tx"), "gene_name"),
#                                     return.type = "DataFrame")))

#cryptic_regions_protein_coding1 = cryptic_regions %>% 
#  mutate(V2 = ifelse(V6 == "-",V2 + 1,V2)) %>% 
#  mutate(V2 = ifelse(V6 == "+",V2 + 1,V2)) %>% 
#  left_join(Txs, 
#            by = c("V8" = "tx_id"))%>% 
#  left_join(Txs %>% dplyr::select(gene_id,tx_biotype)) %>% 
#  left_join(Txs %>% dplyr::select(gene_name,tx_biotype),
#            by = c("V7" = "gene_name"))%>% 
#  mutate(biotype_full = case_when(is.na(tx_biotype.y) ~ tx_biotype.x,
#                                  is.na(tx_biotype.x) ~ tx_biotype.y,
#                                  is.na(tx_biotype.x) & is.na(tx_biotype.y) ~ NA_character_,
#                                  TRUE ~ tx_biotype.y)) %>% 
#  dplyr::select(-tx_biotype.x,-tx_biotype.y) %>% 
#  dplyr::filter(biotype_full == 'protein_coding')

#cryptic_regions_protein_coding1$transcript_id = cryptic_regions_protein_coding1$V8



#cds_regions1 = cdsBy(edb, "tx",use.names = TRUE)
#cds_regions1 = unlist(cds_regions1)

#cds_regions1$transcript_id = gsub("\\..*", "", names(cds_regions1))
#convert from NCBI to UCSC
#cds_regions1=renameSeqlevels(cds_regions1, na.omit(mapSeqlevels(seqlevels(cds_regions1),"UCSC")))

#END EnsDB

if (file.exists(paste0(folder,"/cds_successful_frames_list.csv"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/cds_successful_frames_list.csv"))
}

if (file.exists(paste0(folder,"/cds_unsuccessful_frames_list.csv"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/cds_unsuccessful_frames_list.csv"))
}

if (file.exists(paste0(folder,"/protein_coding.bed"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/protein_coding.bed"))
}
#this is for selecting protein coding genes in main pipeline
if (file.exists(paste0(folder,"/protein_coding_a.bed"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/protein_coding_a.bed"))
}

if (file.exists(paste0(folder,"/aa.fasta"))) {
  #Delete file if it exists
  file.remove(paste0(folder,"/aa.fasta"))
}


cryptic_regions_protein_coding$PTC = NA_character_
#get unique list only
cryptic_regions_protein_coding=unique(cryptic_regions_protein_coding)
zero_cds=0
all_frames_flag=1
for(i in 1:length(cryptic_regions_protein_coding)){
  
  #gene_name = cryptic_regions_protein_coding[i]$V7
  strand=as.character(strand(cryptic_regions_protein_coding[i]))
  gene_name = cryptic_regions_protein_coding[i]$V7
  parent_transcript = cryptic_regions_protein_coding[i]$transcript_id
  
  cds_parent = cds_regions %>% filter(transcript_id == parent_transcript)
  #ALSO MAKE SURE THAT cds covers EVENT RNAGE. WAS CASE FOR ATP6V1B2
  min_event=min(start(range(cryptic_regions_protein_coding[i])),end(range(cryptic_regions_protein_coding[i])))
  max_event=min(start(range(cryptic_regions_protein_coding[i])),end(range(cryptic_regions_protein_coding[i])))
  #if(length(cds_parent) == 0){
  if((length(cds_parent) == 0)||(min(start(ranges(cds_parent))))>=min_event && max(end(ranges(cds_parent))) <= max_event){
    print(paste0('got zero length cds for gene: ',gene_name,' @ line: ',i))
    SJdat1 <- data.frame(as.character(seqnames(cryptic_regions_protein_coding[i])),
                         start(range(cryptic_regions_protein_coding[i]))-1,
                         end(range(cryptic_regions_protein_coding[i])),
                         cryptic_regions_protein_coding[i]$V4,
                         0,
                         as.character(strand(cryptic_regions_protein_coding[i])),
                         cryptic_regions_protein_coding[i]$V7,
                         cryptic_regions_protein_coding[i]$V8)
    write.table(SJdat1,file=paste0(folder,"/cds_unsuccessful_frames_list.csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
    zero_cds=zero_cds+1
    write(paste0('From get_orf_cds.R: got zero length cds for gene: ',gene_name,' see file: ',args[1]),file=args[2],append=TRUE)
    
    next()
  }
  else
  {
    print(paste0('processing cds for gene: ',gene_name,' @ line: ',i))
    SJdat1 <- data.frame(as.character(seqnames(cryptic_regions_protein_coding[i])),
                         start(range(cryptic_regions_protein_coding[i]))-1,
                         end(range(cryptic_regions_protein_coding[i])),
                         cryptic_regions_protein_coding[i]$V4,
                         0,
                         as.character(strand(cryptic_regions_protein_coding[i])),
                         cryptic_regions_protein_coding[i]$V7,
                         cryptic_regions_protein_coding[i]$V8)
    
    SJdat2 <- data.frame(as.character(seqnames(cryptic_regions_protein_coding[i])),
                         start(range(cryptic_regions_protein_coding[i]))-1,
                         end(range(cryptic_regions_protein_coding[i])),
                         1,
                         0,
                         as.character(strand(cryptic_regions_protein_coding[i])),
                         cryptic_regions_protein_coding[i]$V7)
    
    
    write.table(SJdat1,file=paste0(folder,"/protein_coding.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    write.table(SJdat2,file=paste0(folder,"/protein_coding_a.bed"),append = TRUE, quote=F, sep="\t", row.names=F, col.names=F)
    
    #print(paste0("processing gene: ",gene_name, ' line #: ',i))
  }
  
  #######FOR ORF
  
  cds_parent1 = GeneStructureTools::reorderExonNumbers(sort(cds_parent))
  
  if(all_frames_flag == 0)
  {
    aaseq=getOrfs(cds_parent1,BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, returnLongestOnly = TRUE, allFrames = FALSE,longest = 1)
    write.table(gene_name,file=paste0(folder,"/aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
    write.table(aaseq$orf_sequence[1],file=paste0(folder,"/aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
    
  }
  else
  {
    #now writing all orf as we ran into problem for GRIP2
    aaseq=getOrfs(cds_parent1,BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, returnLongestOnly = TRUE, allFrames = TRUE)#,longest = 1)

    write.table(gene_name,file=paste0(folder,"/aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
    write.table(aaseq$orf_sequence[1],file=paste0(folder,"/aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
    write.table(gene_name,file=paste0(folder,"/aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
    write.table(aaseq$orf_sequence[2],file=paste0(folder,"/aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
    write.table(gene_name,file=paste0(folder,"/aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
    write.table(aaseq$orf_sequence[3],file=paste0(folder,"/aa.fasta"),append = TRUE, quote=F, row.names=F, col.names=F)
    
  }
  
  #########END
  
#  exon_one = cryptic_regions_protein_coding[i]
  
#  new_model = sort(c(exon_one,cds_parent))
  
#  new_model = GeneStructureTools::reorderExonNumbers(new_model)
  
 # names(new_model) = NULL
#  new_model = unlist(addCdsPhase(GRangesList(new_model)))
  
  #get_frame=new_model %>%filter(str_detect(V7, exon_one$V7))

#  cryptic_number = new_model %>% plyranges::filter(is.na(cds_id)) %>% 
#   as.data.frame() %>% pull(exon_number)
  
#  cryptic_up_down = new_model %>% filter(exon_number %in% c(cryptic_number - 1, 
 #                                                          cryptic_number,cryptic_number + 1))
  #print(paste0('got cryptic_up_down length: ',length(cryptic_up_down), ' for gene: ',gene_name, ' starnd: ',strand,' @ line: ',i))
  
#  if(strand == '+')
#  {
#    get_frame=cryptic_up_down[1]
#  }
#  else
#  {
#    if (length(cryptic_up_down)==2)
#      get_frame=cryptic_up_down[2]
#    else get_frame=cryptic_up_down[3]
#  }
  #for exon skip, we always provide up_stream exon
#  get_frame=cryptic_up_down[which(cryptic_up_down$V7==gene_name)]
  
#  cds_parent = GeneStructureTools::reorderExonNumbers(cds_parent)
#  normal = my_orf(cds_parent,BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38) %>% 
#    dplyr::slice(1) %>% pull(aa_sequence)
  
#  tmp = GeneStructureTools::getOrfs(new_model,BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
 #                                   returnLongestOnly = FALSE, 
  #                               allFrames = TRUE,
   #                               uORFs = TRUE)
  SJdat <- data.frame(as.character(seqnames(cryptic_regions_protein_coding[i])),
                     start(range(cryptic_regions_protein_coding[i]))-1,
                    end(range(cryptic_regions_protein_coding[i])),
                   cryptic_regions_protein_coding[i]$V4,
                   0,
                 as.character(strand(cryptic_regions_protein_coding[i])),
                cryptic_regions_protein_coding[i]$V7,
               cryptic_regions_protein_coding[i]$V8)
              #which.max(tmp$orf_length),
              #get_frame$cds_phase)
  
  #chr=as.character(seqnames(cryptic_regions_protein_coding[i]))
  #start=start(range(cryptic_regions_protein_coding[i]))
  #end=end(range(cryptic_regions_protein_coding[i]))
  #size=cryptic_regions_protein_coding[i]$V4
  #strand=as.character(strand(cryptic_regions_protein_coding[i]))
  #gen=cryptic_regions_protein_coding[i]$V7
  #tx=cryptic_regions_protein_coding[i]$V8
  #cdsphase=new_model[which(new_model$biotype_full=="protein_coding")]$cds_phase
  #frame=which.max(tmp$orf_length)
  
  #SJdat <- data.frame(chr,
   #                   start,
    #                  end,
     #                 size,
      #                0,
       #               strand,
        #              gen,
         #             tx,
          #            frame)
  
  
  #nr=cbind(cryptic_regions[k],which.max(tmp$orf_length))
  write.table(SJdat,file=paste0(folder,"/cds_successful_frames_list.csv"),append = TRUE, quote=F, sep=",", row.names=F, col.names=F)
}
write(paste0('From get_orf_cds.R: ---------------- Done processing file: ',args[1], ' zero_cds events are: ',zero_cds),file=args[2],append=TRUE)
write(paste0('                                 '),file=args[2],append=TRUE)
