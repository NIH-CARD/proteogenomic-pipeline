
#INLINE README FOR SASHIMI PLOTS FOR ALL EVENTS
#FIRST CALLS peaks_mappingV3.R to get csv file for each case to get bed files

#Rscript pgp-c_mapping.R synth_small_peptide_Final_CE_2021127.csv synth_6tfi_CE_INCLUSION_FUSED_AA.fasta IGV_unique_ce_inclusion.csv synth_ce_inclusion_fused.transeq_in.fasta synth_6tfi_SKIPTICS_FUSED_AA.fasta IGV_unique_skiptics_translated.csv synth_6tfi_skiptics_fused_transeq_in.fasta IR_AA.fasta IR_list.csv IR_nt.transeq_in.fasta #synth_6tfi_FINAL_CE_INCLUSION_AA.fasta IGV_unique_ce_inclusion.csv #small_peptide_Final_CE_2021127.csv new_CE_FUSED_AA.fasta
#mkdir -p peaks_beds
mkdir -p ce_res
mkdir -p es_res
mkdir -p IR_res

#Rscript pgp-c_mappingV2.R inputs/all25events_synth.csv inputs/PEAKS_AA.fasta inputs/cds_merged_igv_ce.csv inputs/cds_merged_nt_ce.fasta inputs/cds_PEAKS_SKIPTICS_FUSED_AA.fasta inputs/cds_IGV_unique_skiptics_translated.csv inputs/cds_skiptics_fused_transeq_in.fasta inputs/FINAL_IR_AA.fasta inputs/IGV_unique_IR.csv inputs/IR_coord_uniq_nt.transeq_in.fasta #synth_6tfi_FINAL_CE_INCLUSION_AA.fasta IGV_unique_ce_inclusion.csv #small_peptide_Final_CE_2021127.csv new_CE_FUSED_AA.fasta

Rscript pgp-c_mappingV3.R inputs/clean_no_humans.csv inputs/Allframes_PEAKS_AA_truncated.fasta inputs/cds_merged_igv_ce.csv inputs/cds_merged_nt_ce.fasta inputs/AllFrames_PEAKS_SKIPTICS_FUSED_AA.fasta inputs/cds_IGV_unique_skiptics_translated.csv inputs/cds_skiptics_fused_transeq_in.fasta inputs/FINAL_IR_AA.fasta inputs/IGV_unique_IR.csv inputs/IR_coord_uniq_nt.transeq_in.fasta #synth_6tfi_FINAL_CE_INCLUSION_AA.fasta IGV_unique_ce_inclusion.csv #small_peptide_Final_CE_2021127.csv new_CE_FUSED_AA.fasta
#NOW CALL AA SCRIPT FOR EACH EVENT TYPE

#SKIPTICS
if [ -f es_res/es_us_exon.bed ]
then
  echo came to es_res/es_us_exon.bed
  #source pgp-c_gc_aa.sh es_res/es_us_exon.bed
  bash run_sashimiV1.sh es_res/es_peaks_sashimi_us_exon.csv es_res/es_us_exon.bed es_res/es_us_exon.csv 1
fi
if [ -f es_res/es_us_exon_ds_exon.bed ]
then
  echo came to es_res/es_us_exon_ds_exon.bed
  #bash pgp-c_gc_aa.sh es_res/es_us_exon_ds_exon.bed
  bash run_sashimiV1.sh es_res/es_peaks_sashimi_us_exon_ds_exonn.csv es_res/es_us_exon_ds_exon.bed es_res/es_us_exon_ds_exon.csv 2 3
fi
if [ -f es_res/es_ds_exon.bed ]
then
  echo came to es_res/es_ds_exon.bed
  #bash pgp-c_gc_aa.sh es_res/es_ds_exon.bed
  bash run_sashimiV1.sh es_res/es_peaks_sashimi_ds_exon.csv es_res/es_ds_exon.bed es_res/es_ds_exon.csv 1
fi

#CRYPTICS
if [ -f ce_res/us_exon.bed ]
then
  echo came to ce_res/us_exon.bed
  #bash pgp-c_gc_aa.sh ce_res/us_exon.bed
  bash run_sashimiV1.sh ce_res/peaks_sashimi_us_exon.csv ce_res/us_exon.bed ce_res/us_exon.csv 1
fi
if [ -f ce_res/us_exon_ce.bed ]
then
  echo came to ce_res/us_exon_ce.bed
  #bash pgp-c_gc_aa.sh ce_res/us_exon_ce.bed
  bash run_sashimiV1.sh ce_res/peaks_sashimi_us_exon_ce.csv ce_res/us_exon_ce.bed ce_res/us_exon_ce.csv 2 0
fi



if [ -f ce_res/ce.bed ]
then
  echo came to ce_res/ce.bed
  #bash pgp-c_gc_aa.sh ce_res/ce.bed
  bash run_sashimiV1.sh ce_res/peaks_sashimi_ce.csv ce_res/ce.bed ce_res/ce.csv 1
fi
if [ -f ce_res/ce_ds_exon.bed ]
then
  echo came to ce_res/ce_ds_exon.bed
  #bash pgp-c_gc_aa.sh ce_res/ce_ds_exon.bed
  bash run_sashimiV1.sh ce_res/peaks_sashimi_ce_ds_exon.csv ce_res/ce_ds_exon.bed ce_res/ce_ds_exon.csv 2 1
fi
if [ -f ce_res/ds_exon.bed ]
then
  echo came to ce_res/ds_exon.bed
  #bash pgp-c_gc_aa.sh ce_res/ds_exon.bed
  bash run_sashimiV1.sh ce_res/peaks_sashimi_ds_exon.csv ce_res/ds_exon.bed ce_res/ds_exon.csv 1
fi
#IR EVENTS
if [ -f IR_res/IR_events.bed ]
then
  echo came to IR_res/IR_events.bed
  #bash pgp-c_gc_aa.sh IR_res/IR_events.bed
  bash run_sashimiV1.sh IR_res/IR_peaks_sashimi.csv IR_res/IR_events.bed IR_res/IR_events.csv 1
fi
#NOW CALL get_UpDnEx_sashimi.sh script to get up and down stream exon coordinates for sashimi plots
#INPUT: events csv file from peaks script, peaks bed file, peaks csv file and type flag. flag = 1 (un_exon, ce and ds_exon), 2 (us_exon_ce, ce_ds_exon),
#for mixed event cases, us_exon_ce.csv, flag=0 - up_Exon followed by ce and for ce_ds_exon.csv, flag=1 means, ce followed by ds_exon

#CE
#bash run_sashimiV1.sh peaks_sashimi_us_exon.csv us_exon.bed us_exon.csv 1
#bash run_sashimiV1.sh peaks_sashimi_us_exon_ce.csv us_exon_ce.bed us_exon_ce.csv 2 0
#bash run_sashimiV1.sh peaks_sashimi_ce.csv ce.bed ce.csv 1
#bash run_sashimiV1.sh peaks_sashimi_ce_ds_exon.csv ce_ds_exon.bed ce_ds_exon.csv 2 1
#bash run_sashimiV1.sh peaks_sashimi_ds_exon.csv ds_exon.bed ds_exon.csv 1

#SKIPTICS

#bash run_sashimiV1.sh es_peaks_sashimi_us_exon.csv es_us_exon.bed es_us_exon.csv 1
#bash run_sashimiV1.sh es_peaks_sashimi_us_exon_ds_exonn.csv es_us_exon_ds_exon.bed es_us_exon_ds_exon.csv 2 3
#bash run_sashimiV1.sh es_peaks_sashimi_ds_exon.csv es_ds_exon.bed es_ds_exon.csv 1

#IR EVENTS

#bash run_sashimiV1.sh IR_peaks_sashimi.csv IR_events.bed IR_events.csv 1
