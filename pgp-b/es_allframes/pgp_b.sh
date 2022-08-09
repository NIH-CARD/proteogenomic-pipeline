#!/bin/bash
#HOW TO RUN:
#1. skiptics: bash pgp_b.sh 3 events.csv
######## Proteogenomic Pipeline for Biomarker Discovery in iPSC neurons
######### Inline Readme
# 1. Get List of most abundant Transcripts from KD iPSC samples
# This step assumes that all bam files are in current folder and follow have
#
if [ $1 -eq 0 ] || [ $1 -eq 2 ] || [ $1 -eq 8 ] #8 is for start to finish
then
  echo STARTED SringTie Calculations - Will take a while depending on number of samples
mkdir -p iPSC_gtfs
#[ -e samples_nt_tdp.txt ] && rm samples_nt_tdp.txt
samples=$(ls *.bam | cut -d'.' -f1 | awk '{ print length, $0 }' | sort -n | cut -d" " -f2-)
#check if we got any samples
if [ -n "$samples" ]
then
for sample in $samples
do
    echo STARTED PROCESSING SAMPLE $sample
    #sam=$(echo $sample | cut -d. -f1)
    #echo $sam
    #stringtie $sample.pass2Aligned.sortedByCoord.out.bam -l gene_abundance.$sample.tab -C cov_$sample.gtf -b cov_tx_$sample.gtf -p 8 -A $sample -G /home/ext_syed_datatecnica_com/mountfolder/ref_genome/gencode.v38.annotation.gtf -o iPSC_gtfs/$sample.gtf
    stringtie $sample.pass2Aligned.sortedByCoord.out.bam  -p 8 -G gencode.v38.annotation.gtf -o iPSC_gtfs/$sample.gtf
    #echo assembly/"$sam".gtf >> samples_nt_tdp.txt
done
echo DONE WITH SringTie Calculations - NOW GENERATING LIST of ABUNDANT Txs
#Now get unique transcripts in csv files for each sample
for i in iPSC_gtfs/*.gtf; do
    [ -f "$i" ] || break
	#echo $i
	samp=$(echo $i|cut -d/ -f2|cut -d. -f1)
	echo PROCESSING SAMPLE $i
	#Step 1. Get all Txs having reference_id (ENST), ref_gene_id (ENSG) and ref_gene_name (HUGO SYMBOL), these are 24 column lines in stringtie's gtf file
	cat $i | awk 'NF==24{print $14,$16,$18,$20,$22,$24}' | awk '{gsub(";","\t");print}' > iPSC_gtfs/"$samp".csv
	#cat $i | awk 'NF==24{print}{}' > "$samp"_24.tab
#Step 2. Now select lines with Highest coverage (col 20) for each ref_gene_name (col 18)
#awk '$20>max[$18]{max[$18]=$20; row[$18]=$0} END{for (i in row) print row[i]}' Cont-A_S1_24.gtf > Cont-A_S1_Tx.gtf
done
#now select most abundant transcripts from all samples
Rscript abundant_tx.R
echo DONE WITH ABUNDANT Txs Generation for input BAM SAMPLES
else #ABORT if we got no samples
	echo No BAM SAMPLES FOUND, ABORTING!!!! Please check BAM samples path and Rerun again
	exit 1
fi
fi
#echo got flag "$1"
#exit 1
#NEXT SORT FILE ON COLUMN 5 - THIS IS IMPORTANT AS PIPELINE DEPENDS ON IT
if [ $1 -eq 1 ] || [ $1 -eq 2 ] || [ $1 -eq 8 ] #8 is for start to finish
then

if [ -f "$2" ]
then
	splicing_events_file=$(echo "$2"|cut -d. -f1)
	sort -t$',' -k5 $2 > sorted_$2
else
	#echo line 54 NO SPLICING_EVENTS csv file is provided, ABORTING!!!!
	echo PLEASE PROVIDE SPLICING_EVENTS csv file  "'"bash pgp_0.sh flag splicing_events.csv"'" and RERUN
	exit 1
fi

echo CAME IN FOR SASHIMI PLOTS FOR ALL EVENTS
#ALSO get bed file for sashimi plots for all events.
#FORMAT IS: chr, start of up_ex,end of ds_ex, 1, 0, strand, gene_name, TxID
generate_all_events_sashimi_beds=1

if [ $generate_all_events_sashimi_beds -eq 1 ]
then
	mkdir -p temp_all_events_sashimi
	echo "################################ GENERATING BED FILE FOR SASHIMI PLOTS FOR ALL EVENTS"
	echo "################################ GENERATING BED FILE FOR SASHIMI PLOTS FOR ALL EVENTS" >> ALL_STATS.txt
	#IMPORTANT - THIS CODE RELIES ON THE OUTPUT OF BEDTOOLS CLOSEST FUNCTION
	#GIVEN A MAJIQ EVENT (as a bed file) and BED FILE FOR THE TRANSCRIPT IT BELIEVES TO BE PART OF
	#BEDTOOLS CLOSEST FUNCTION USES INPUT RANGE (AMJIQ JUNTION) and finds closest exons (and their ranges) from the TRANSCRIPT BED FILE
	#BEDTOOLS CLOSEST RETURNS A BED FILE WITH FOLLOWING OUTPUT
	#INPUT: chr,start,end,1,0,+ (majiq event with strand)
	#OUTPUT: chr,start,end,1,0,+,chr,start,end,size,exon_num,+,distance
	#here first 6 entries are the original majiq input and next 7 entries are the resulting closest exon coordinates, its size (in bp),exon_num,strand and distance from reference
	#SHOULD ADD CHECKS ON EXON_NUM READING FROM FILE

		[ -e temp_all_events_sashimi/"$splicing_events_file"_all_sashimi.bed ] && rm temp_all_events_sashimi/"$splicing_events_file"_all_sashimi.bed #for abs(ex1-ex2)>1
		[ -e temp_all_events_sashimi/"$splicing_events_file"_all_sashimi.csv ] && rm temp_all_events_sashimi/"$splicing_events_file"_all_sashimi.csv
		[ -e temp_all_events_sashimi/"$splicing_events_file"_all_sashimi2.bed ] && rm temp_all_events_sashimi/"$splicing_events_file"_all_sashimi2.bed #for abs(ex1-ex2)>1
		[ -e temp_all_events_sashimi/"$splicing_events_file"_all_sashimi2.csv ] && rm temp_all_events_sashimi/"$splicing_events_file"_all_sashimi2.csv
		[ -e temp_all_events_sashimi/"$splicing_events_file"_all_sashimi01.bed ] && rm temp_all_events_sashimi/"$splicing_events_file"_all_sashimi01.bed #+ strand for abs(ex1-ex2)=0
		[ -e temp_all_events_sashimi/"$splicing_events_file"_all_sashimi01.csv ] && rm temp_all_events_sashimi/"$splicing_events_file"_all_sashimi01.csv
		[ -e temp_all_events_sashimi/"$splicing_events_file"_all_sashimi02.bed ] && rm temp_all_events_sashimi/"$splicing_events_file"_all_sashimi02.bed #- strand for abs(ex1-ex2)=0
		[ -e temp_all_events_sashimi/"$splicing_events_file"_all_sashimi02.csv ] && rm temp_all_events_sashimi/"$splicing_events_file"_all_sashimi02.csv
		[ -e temp_all_events_sashimi/"$splicing_events_file"_progress2.txt ] && rm temp_all_events_sashimi/"$splicing_events_file"_progress2.txt
		[ -e temp_all_events_sashimi/"$splicing_events_file"_progress01.txt ] && rm temp_all_events_sashimi/"$splicing_events_file"_progress01.txt
		[ -e temp_all_events_sashimi/"$splicing_events_file"_progress02.txt ] && rm temp_all_events_sashimi/"$splicing_events_file"_progress02.txt
		[ -e temp_all_events_sashimi/"$splicing_events_file"_progress_all.txt ] && rm temp_all_events_sashimi/"$splicing_events_file"_progress_all.txt
		#[ -e temp_all_events_sashimi/"$splicing_events_file"_progress1.txt ] && rm temp_all_events_sashimi/"$splicing_events_file"_progress1.txt

		generate_event_beds=1
		if [ $generate_event_beds -eq 1 ]
		then
				mkdir -p event_bedfiles
				[ "$(ls -A event_bedfiles/)" ] && rm event_bedfiles/*.*

				Rscript TxEnsDB103_layeredV5.R sorted_$2 principal_txs.csv #$tx_lst
		fi

		readarray -t csv_data < all_tx_events.csv
		csvi=0

		samples=$(ls event_bedfiles/temp_*.bed)

		for sample in $samples
		do
			#read csv entry
			csv_ln=${csv_data[$csvi]}
			((csvi=csvi+1))

			#echo processing "$sample"
			allexons=$(echo "$sample" | cut -d'/' -f2 |cut -d'_' -f2)
		  gene_name1=$(echo "$allexons" | cut -d'.' -f1)
		  gene_name=$(echo "$gene_name1" | cut -d'-' -f1)

		  #first sort the bed
			sortBed -i event_bedfiles/$allexons > event_bedfiles/t$allexons

			#Also read Tx Files to retrieve selected Tx - should find better ways
			TxID=$(head -1 event_bedfiles/TxID$allexons | awk '{print $7}')
			#echo got TxID $TxID

		  strnd=$(cat "$sample" | awk '{print $6}')
			#get distance to downstream exon (for ties, report first) from current reference and pick start, end and d
			ds=$(bedtools closest -a $sample -b event_bedfiles/t$allexons -s -D a -iu -d -t first )
		  #also get distance to upstream exon from current reference and pick start, end and d
			us=$(bedtools closest -a $sample -b event_bedfiles/t$allexons -s  -D a -id -d -t last)

			#get up and down stream exon numbers
			upexon=$(echo "$us" | awk '{print $11}')
			dnexon=$(echo "$ds" | awk '{print $11}')
			#events star and end
			event_st=$(echo "$us" | awk '{print $2}')
			event_end=$(echo "$us" | awk '{print $3}')

			diff_exon=$(($upexon-$dnexon))
			#take absolute value
			diff_exon_abs=${diff_exon#-}
			if [ "$diff_exon_abs" -ge 1 ] #ALL EVENTS THAT SPANS 2 OR MORE EXONS
			then
					if [ "$strnd" == "+" ]
					then
						start=$(echo "$us" | awk '{print $8}')
						end=$(echo "$ds" | awk '{print $9}')
					else
						start=$(echo "$ds" | awk '{print $8}')
						end=$(echo "$us" | awk '{print $9}')
					fi
					#also save
					#FORMAT IS: chr, start of up_ex,end of ds_ex, 1, 0, strand, gene_name, TxID
					#First check if event lies between selected exons
					if [ "$start" -le "$event_st" ] && [ "$end" -ge "$event_end" ]
					then
						echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi2.bed
						echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi.bed
						echo $csv_ln >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi2.csv
						echo $csv_ln >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi.csv
					else
						#echo ds $ds >> "$splicing_events_file"_progress2.txt
						#echo us $us >> "$splicing_events_file"_progress2.txt
						if [ "$strnd" == "+" ]
						then
							start=$(echo "$us" | awk '{print $8}')
							end=$(echo "$ds" | awk '{print $9}') #both are same
							#first check if star > event_start, then select upstream exon
							if [ "$start" -ge "$event_st" ]
							then
								#get exon (line number in bed file) to read
								exon=$(echo "$ds" | awk '{print $11}')
								exon=$(($exon-2))
								readarray -t bed_data < event_bedfiles/$allexons
								bed_ln=${bed_data[$exon]}
								#update start
								start=$(echo "$bed_ln" | awk '{print $2}')
							fi
							#now check if end <event_end
							if [ "$end" -le "$event_end" ]
							then
								#get exon (line number in bed file) to read
								exon=$(echo "$us" | awk '{print $11}')
								#exon=$(($exon+1)) #reading line for readarray starts from 0
								readarray -t bed_data < event_bedfiles/$allexons
								echo bed_data $bed_data
								bed_ln=${bed_data[$exon]}
								echo exon $exon bed_ln $bed_ln
								#update start
								end=$(echo "$bed_ln" | awk '{print $3}')
							fi
							#now one more time check if event lies between selected exons
							if [ "$start" -le "$event_st" ] && [ "$end" -ge "$event_end" ]
							then
								echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi2.bed
								echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi.bed
								echo $csv_ln >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi2.csv
								echo $csv_ln >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi.csv
							else
								echo ds 1 $ds >> temp_all_events_sashimi/"$splicing_events_file"_progress2.txt
								echo us 1 $us >> temp_all_events_sashimi/"$splicing_events_file"_progress2.txt

								echo diff_exon_abs is $diff_exon_abs selected event $sample has event_st $event_st selected start $start event end $event_end selected end $end - please check
							fi

						else #THIS IS FOR NEGATIVE STRAND
							start=$(echo "$ds" | awk '{print $8}')
							end=$(echo "$us" | awk '{print $9}')
							#echo ds $ds >> majiq_events_progress0.txt
							#echo us $us >> majiq_events_progress0.txt

							#first check if star > event_start, then select upstream exon
							if [ "$start" -ge "$event_st" ]
							then
								#get exon (line number in bed file) to read
								exon=$(echo "$ds" | awk '{print $11}')
								exon=$(($exon)) #Tx on -ve strand has exons listed from bottom to top in increasing order
								readarray -t bed_data < event_bedfiles/$allexons
								bed_ln=${bed_data[$exon]}
								#update start
								start=$(echo "$bed_ln" | awk '{print $2}')
							fi
							#now check if end <event_end
							if [ "$end" -le "$event_end" ]
							then
								#get exon (line number in bed file) to read
								exon=$(echo "$us" | awk '{print $11}')
								exon=$(($exon-2)) #Tx on -ve strand has exons listed from bottom to top in increasing order
								readarray -t bed_data < event_bedfiles/$allexons
								bed_ln=${bed_data[$exon]}
								#update start
								end=$(echo "$bed_ln" | awk '{print $3}')
							fi
							#now one more time check if event lies between selected exons
							if [ "$start" -le "$event_st" ] && [ "$end" -ge "$event_end" ]
							then
								echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi2.bed
								echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi.bed
								echo $csv_ln >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi2.csv
								echo $csv_ln >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi.csv
							else
								echo ds 2 $ds >> temp_all_events_sashimi/"$splicing_events_file"_progress2.txt
								echo us 2 $us >> temp_all_events_sashimi/"$splicing_events_file"_progress2.txt

								echo diff_exon_abs is $diff_exon_abs selected event $sample has event_st $event_st selected start $start event end $event_end selected end $end - please check
							fi
						fi
						#echo diff_exon_abs is $diff_exon_abs selected event $sample has event_st $event_st selected start $start event end $event_end selected end $end - please check
					fi
				elif [ "$diff_exon_abs" -eq 0 ]
				then
					if [ "$strnd" == "+" ]
					then
						start=$(echo "$us" | awk '{print $8}')
						end=$(echo "$ds" | awk '{print $9}') #both are same
						#first check if star > event_start, then select upstream exon
						if [ "$start" -ge "$event_st" ]
						then
							#get exon (line number in bed file) to read
							exon=$(echo "$us" | awk '{print $11}')
							exon=$(($exon-2))
							readarray -t bed_data < event_bedfiles/$allexons
							bed_ln=${bed_data[$exon]}
							#update start
							start=$(echo "$bed_ln" | awk '{print $2}')
							#now go on the other sied
							if [ "$start" -ge "$event_st" ] #get the other exon
							then
								#get exon (line number in bed file) to read
								exon=$(echo "$us" | awk '{print $11}')
								#exon=$(($exon+1))
								readarray -t bed_data < event_bedfiles/$allexons
								bed_ln=${bed_data[$exon]}
								#update start
								start=$(echo "$bed_ln" | awk '{print $2}')
							fi

						fi
						#now check if end <event_end
						if [ "$end" -le "$event_end" ]
						then
							#get exon (line number in bed file) to read
							exon=$(echo "$us" | awk '{print $11}')
							#exon=$(($exon+1)) #reading line for readarray starts from 0
							readarray -t bed_data < event_bedfiles/$allexons
							echo bed_data $bed_data
							bed_ln=${bed_data[$exon]}
							echo exon $exon bed_ln $bed_ln
							#update start
							end=$(echo "$bed_ln" | awk '{print $3}')
							if [ "$end" -le "$event_end" ]
							then
								#get exon (line number in bed file) to read
								exon=$(echo "$us" | awk '{print $11}')
								exon=$(($exon-2)) #reading line for readarray starts from 0
								readarray -t bed_data < event_bedfiles/$allexons
								bed_ln=${bed_data[$exon]}
								#update start
								end=$(echo "$bed_ln" | awk '{print $3}')
							fi

						fi
						#now one more time check if event lies between selected exons
						if [ "$start" -le "$event_st" ] && [ "$end" -ge "$event_end" ]
						then
							echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi01.bed
							echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi.bed

							echo $csv_ln >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi01.csv
							echo $csv_ln >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi.csv
						else
							echo ds $ds >> temp_all_events_sashimi/"$splicing_events_file"_progress01.txt
							echo us $us >> temp_all_events_sashimi/"$splicing_events_file"_progress01.txt

							echo diff_exon_abs is $diff_exon_abs selected event $sample has event_st $event_st selected start $start event end $event_end selected end $end - please check
						fi

					else #THIS IS FOR NEGATIVE STRAND
						start=$(echo "$ds" | awk '{print $8}')
						end=$(echo "$us" | awk '{print $9}')
						#echo ds $ds >> majiq_events_progress0.txt
						#echo us $us >> majiq_events_progress0.txt

						#first check if star > event_start, then select upstream exon
						if [ "$start" -ge "$event_st" ]
						then
							#get exon (line number in bed file) to read
							exon=$(echo "$us" | awk '{print $11}')
							#exon=$(($exon+1)) #Tx on -ve strand has exons listed from bottom to top in increasing order
							readarray -t bed_data < event_bedfiles/$allexons
							bed_ln=${bed_data[$exon]}
							#update start
							start=$(echo "$bed_ln" | awk '{print $2}')
						fi
						#now check if end <event_end
						if [ "$end" -le "$event_end" ]
						then
							#get exon (line number in bed file) to read
							exon=$(echo "$us" | awk '{print $11}')
							exon=$(($exon-2)) #Tx on -ve strand has exons listed from bottom to top in increasing order
							readarray -t bed_data < event_bedfiles/$allexons
							bed_ln=${bed_data[$exon]}
							#update start
							end=$(echo "$bed_ln" | awk '{print $3}')
						fi
						#now one more time check if event lies between selected exons
						if [ "$start" -le "$event_st" ] && [ "$end" -ge "$event_end" ]
						then
							echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi02.bed
							echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi.bed
							echo $csv_ln >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi02.csv
							echo $csv_ln >> temp_all_events_sashimi/"$splicing_events_file"_all_sashimi.csv
						else
							echo ds $ds >> temp_all_events_sashimi/"$splicing_events_file"_progress02.txt
							echo us $us >> temp_all_events_sashimi/"$splicing_events_file"_progress02.txt

							echo diff_exon_abs is $diff_exon_abs selected event $sample has event_st $event_st selected start $start event end $event_end selected end $end - please check
						fi
					fi
					#first check if star > event_start, then select upstream exon
				else
					echo ds $ds >>  temp_all_events_sashimi/"$splicing_events_file"_progress_all.txt
					echo us $us >>  temp_all_events_sashimi/"$splicing_events_file"_progress_all.txt
			fi


			#echo $lines1 | awk -v tx=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,tx}' >> res_skiptics/skiptics_uniq_sashimi.bed
			#echo $lines2 | awk -v tx=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,tx}' >> res_skiptics/skiptics_uniq_sashimi.bed

		done


		echo "################################ DONE "
		echo "################################ DONE " >> ALL_STATS.txt
		#if [ $generate_all_events_sashimi_plots -eq 1 ]
    if [ $1 -eq 1 ] || [ $1 -eq 2 ] || [ $1 -eq 8 ] #8 is for start to finish
		then
			echo "################################ STARTED CREATING SASHIMI PLOTS FOR ALL EVENTS - MAY TAKE MANY HOURS " >> ALL_STATS.txt
			echo "################################ STARTED CREATING SASHIMI PLOTS FOR ALL EVENTS - MAY TAKE MANY HOURS "

			#CALL ggsashimi for SASHIMI PLOTS
			source run_sashimi.sh 2 temp_all_events_sashimi/"$splicing_events_file"

			echo "################################ DONE - RESULTANT pdf is sashimi_plots/all_events_sashimi.pdf " >> ALL_STATS.txt
		else
			echo "################################ SKIPPING SASHIMI PLOTS " >> ALL_STATS.txt
		fi


fi

fi

[ -e ALL_STATS.txt ] && rm ALL_STATS.txt
if [ $1 -eq 3 ] || [ $1 -eq 4 ] || [ $1 -eq 8 ] #8 is for start to finish
then

  if [ -f "$2" ]
  then
  	splicing_events_file=$(echo "$2"|cut -d. -f1)
  	sort -t$',' -k5 $2 > sorted_"$splicing_events_file".csv
  else
  	#echo line 402 NO SPLICING_EVENTS csv file is provided, ABORTING!!!!
  	echo PLEASE PROVIDE SPLICING_EVENTS csv file  "'"bash pgp_b.sh flag splicing_events.csv"'" and RERUN
  	exit 1
  fi


#EVENTS CSV FILE CLEANING STARTS HERE
step00_flg=0
if [ $step00_flg -eq 1 ]
then
		echo STEP 0 - CLEANING INPUT FILE FOR DUPLICATES "("EVENTS WITH SAME CHR"#", START AND END")" AND STAR ">" END OR START==END EVENTS- WILL TAKE AROUND 30 MINUTES
		[ -e file_clean_report.txt ] && rm file_clean_report.txt

		#get total records

		t_recrds=$(cat sorted_"$splicing_events_file".csv | wc -l)
		echo original file has $t_recrds records >> file_clean_report.txt
		#delete if already exist
		[ -e clean_"$splicing_events_file".csv ] && rm clean_"$splicing_events_file".csv
		readarray -t all_splicing_events < sorted_$2 #for saving gene_ids as well to generate sashimi compatible csv file
		i=0

		while [ $i -lt $t_recrds ]
		do
		          line=${all_splicing_events[$i]}
							chr=$(echo $line | awk 'BEGIN {FS=","} {print $1}')
							start=$(echo $line | awk 'BEGIN {FS=","} {print $2}')
							end=$(echo $line | awk 'BEGIN {FS=","} {print $3}')
							gene_id=$(echo $line | awk 'BEGIN {FS=","} {print $6}')
							lastline=$i
							((i=i+1))
							#now go through each record
							j=$i
				      #now go through rest of the data
				      flg=0

				         while [ $j -lt $t_recrds ]
				         do
									  nextline=$(($j+1))
									 	nline=${all_splicing_events[$j]}
										nchr=$(echo $nline | awk 'BEGIN {FS=","} {print $1}')
				 						nstart=$(echo $nline | awk 'BEGIN {FS=","} {print $2}')
				 						nend=$(echo $nline | awk 'BEGIN {FS=","} {print $3}')
										ngene_id=$(echo $nline | awk 'BEGIN {FS=","} {print $6}')
				           	if [ $chr == $nchr ] && [ $start -eq $nstart ] && [ $end -eq $nend ] && [ $gene_id == $ngene_id ]
				           	then
											echo Line $i - $line - has same start: $start and end: $end as line $nextline - $nline -, so ignoring >> file_clean_report.txt
				            	flg=1
				           	fi
				           	((j=j+1))
				         done

							if [ $flg -eq 0 ] && [ $start -lt $end ]
							then
								#ALSO CHECK IF START < END
								if [ $start -lt $end ]
								then
									echo $line >> clean_"$splicing_events_file".csv
								else
									echo Line $i has problem with start: $start and end: $end >> file_clean_report.txt
								fi
							fi

		done
		c_records=$(cat clean_"$splicing_events_file".csv | wc -l)
		echo
		echo
		echo DONE WITH CLEANING EVENTS FILE
		echo CLEANED FILE clean_"$splicing_events_file".csv HAS $c_records RECORDS >> file_clean_report.txt
		echo PLEASE CHECK file_clean_report.txt FILE FOR CLEANING PROCESSES AND clean_"$splicing_events_file".csv FOR RESULTING CLEAN FILE
else
  cat sorted_"$splicing_events_file".csv > clean_"$splicing_events_file".csv
fi

  #statements
#EVENTS CSV FILE CLEANING ENDS HERE


#STEP 1. Call esV4_layered.sh for exon_skip events
echo "################################ STARTED SKIPTICS IDENTIFICATION "
echo "################################ STARTED SKIPTICS IDENTIFICATION " >> ALL_STATS.txt
[ -e EnsDB_tx_not_found.csv ] && rm EnsDB_tx_not_found.csv
skiptics_flg=1
#skiptics_sashimi_flg=0
if [ $skiptics_flg -eq 1 ]
then

    #if [ -f "$2" ]
    if [ -f clean_"$splicing_events_file".csv ]
    then
    	splicing_events_file=$(echo "$2"|cut -d. -f1)
    	sort -t$',' -k5 $2 > sorted_$2
    else
    	echo line 494 NO SPLICING_EVENTS csv file is provided, ABORTING!!!!
    	echo PLEASE PROVIDE SPLICING_EVENTS csv file  "'"bash pgp_0.sh flag splicing_events.csv"'" and RERUN
    	exit 1
    fi

		echo '#############################'
    echo INVOKING SKIPTICS Script esV5_layered_CDSV1.sh FOR FILE clean_"$splicing_events_file".csv from pgp_b.sh
		echo INVOKING SKIPTICS Script esV5_layered_CDSV1.sh FOR FILE clean_"$splicing_events_file".csv from pgp_b.sh >> ALL_STATS.txt
		echo #############################
		#source esV5_layered.sh sorted_$1 principal_tx2_merged_samples.csv #$2 #clean_"$splicing_events_file".csv #$1 #$2
    source esV5_layered_CDSV3.sh clean_"$splicing_events_file".csv principal_txs.csv #principal_tx2_merged_samples.csv #$2 #clean_majiq_events.csv #$1 #$2

		echo '#############################'
		echo BACK FROM SKIPTICS Script
		echo '#############################'
		echo "################################ DONE SKIPTICS IDENTIFICATION ", PLEASE SEE res_skiptics/FINAL_STATS_SKIPTICS.txt for skiptics details >> ALL_STATS.txt
		#CALL ggsashimi for SASHIMI PLOTS
		if [ $1 -eq 4 ] || [ $1 -eq 8 ] #8 is for start to finish
		then
			echo "################################ STARTED CREATING SASHIMI PLOTS FOR SKIPTICS EVENTS - MAY TAKE MANY HOURS " >> ALL_STATS.txt
			echo "################################ STARTED CREATING SASHIMI PLOTS FOR SKIPTICS EVENTS - MAY TAKE MANY HOURS "
			source run_sashimi.sh 1
			echo "################################ DONE - RESULTANT pdf is sashimi_plots/all_skiptics_sashimi.pdf " >> ALL_STATS.txt
		fi
fi
  echo PLEASE CHECK res_skiptics/SKIPTICS_FINAL_STATS.txt file all statistics
fi
#if [ $1 -eq 5 ] || [ $1 -eq 6] || [ $1 -eq 5]
#then
#Step 2 - STARTING ce calculations
if [ $1 -eq 5 ] || [ $1 -eq 8 ] #8 is for start to finish #create bed files for all bam samples
then
echo "############### NOW GENERATING BED FILES FOR ALL BAM SAMPLES, WILL TAKE A WHILE #########" >> ALL_STATS.txt
echo "############### NOW GENERATING BED FILES FOR ALL BAM SAMPLES, WILL TAKE A WHILE #########"

#Step 1. Calculate TDP_SAMPLE_XXX.bam.bed files for each TDP43KD BAM file
#var to disable following
#bam_bed_flg=1
#if [ $1 -eq 5 ]
#then
  mkdir -p bam_beds
  #for file in bam_files/*.bam
  for file in *.bam
  do
  	echo processing "$file"
  	#sample_bam=$(echo "$file" | cut -d'.' -f3 | cut -d'/' -f3 | awk '{print $1}')
  	#sample_bam=$(echo "$file" | cut -d'/' -f2 | cut -d'.' -f1)
    sample_bam=$(echo "$file" | cut -d'.' -f1)
  	echo generating bam_beds/"$sample_bam".bam.bed

  	bedtools bamtobed -split -i "$file" > bam_beds/"$sample_bam".bam.bed
  	#also create genome file containing chromosomes start and end for this bam file
  	samtools idxstats "$file" | cut -f 1-2 > bam_beds/"$sample_bam".chromosomes1.txt
    #I have to do this to delete last line which is always 0, do not know why?
    sed '$d' bam_beds/"$sample_bam".chromosomes1.txt > bam_beds/"$sample_bam".chromosomes.txt
    rm bam_beds/"$sample_bam".chromosomes1.txt
  	#Now sort the BED file according to this genome file
  	bedtools sort -faidx bam_beds/"$sample_bam".chromosomes.txt -i bam_beds/"$sample_bam".bam.bed > bam_beds/"$sample_bam"-sorted.bam.bed
  done
  if [ $1 -ne 8 ]
  then
    echo DONE CREATING BED FIELS FOR BAM SAMPLES - NOW EXITING
    exit 1
  else
    echo "############### DONE CREATING BED FIELS FOR BAM SAMPLES #########" >> ALL_STATS.txt
    echo "############### DONE CREATING BED FIELS FOR BAM SAMPLES #########"

  fi
fi

#Step 2. Generate intronic_range and ce__all_scan_range_junctions.bed and ce_all_scan_range.bed files
if [ $1 -eq 6 ] || [ $1 -eq 7 ] || [ $1 -eq 8 ] #8 is for start to finish # 6 is to also create cov bed files and 7 is to create sashimi plots with out re-generating coverage files
then
if [ $1 -eq 8 ]
then
  if [ -f all_non_skiptics.csv ]
  then
	   inpfile=all_non_skiptics.csv
   else
     echo FOUND EMPTY all_non_skiptics.csv SPLICING EVENTS FILE, PLEAES FIX THE PROBLEM AND RERUN
     echo EXITING NOW
   fi
else
	inpfile=$2
#	inpfile=all_non_skiptics.csv
fi
#tx_lst=$2
#flag to run this section
cryptics_flg=1
if [ $cryptics_flg -eq 1 ]
then
  echo "############### NOW STARTING CE IDENTIFICATION #########" >> ALL_STATS.txt

  echo NOW STARTING CE IDENTIFICATION

#mkdir -p res_intronic_range
mkdir -p res_ce_all
mkdir -p event_bedfiles


[ -e events_to_tx_mapping_valid.csv ] && rm events_to_tx_mapping_valid.csv
[ -e events_tx_mapping_invalid.csv ] && rm events_tx_mapping_invalid.csv

########IR SECTION

[ -e res_ce_all/IR_coord.bed ] && rm res_ce_all/IR_coord.bed
[ -e res_ce_all/IR_coord_sashimi.bed ] && rm res_ce_all/IR_coord_sashimi.bed
[ -e res_ce_all/IGV_R_returned_IR.csv ] && rm res_ce_all/IGV_R_returned_IR.csv
[ -e res_ce_all/IR_coord_uniq.bed ] && rm res_ce_all/IR_coord_uniq.bed
[ -e res_ce_all/IR_coord_uniq_sashimi.bed ] && rm res_ce_all/IR_coord_uniq_sashimi.bed
[ -e res_ce_all/IR_coord_only.bed ] && rm res_ce_all/IR_coord_only.bed
[ -e res_ce_all/IGV_unique_IR.csv ] && rm res_ce_all/IGV_unique_IR.csv



######END IR SECTION

#[ "$(ls -A event_bedfiles/)" ] && rm event_bedfiles/*.*
echo CLEARING UP res_ce_all folder from previous run if any. SASHIMI FOLDERS WILL ONLY BE DELETED BEFORE CREATING NEW SASHIMI PLOTS.
[ -e res_ce_all/ce_all_scan_intron.bed ] && rm res_ce_all/ce_all_scan_intron.bed
[ -e res_ce_all/ce_inclusion_coord_sashimi.bed ] && rm res_ce_all/ce_inclusion_coord_sashimi.bed
[ -e res_ce_all/ce_inclusion_coord_only.bed ] && rm res_ce_all/ce_inclusion_coord_only.bed
[ -e res_ce_all/ce_extension_coord_only.bed ] && rm res_ce_all/ce_extension_coord_only.bed
[ -e res_ce_all/IGV_R_returned_ce_extension.csv ] && rm res_ce_all/IGV_R_returned_ce_extension.csv
[ -e res_ce_all/IGV_R_returned_ce_inclusion.csv ] && rm res_ce_all/IGV_R_returned_ce_inclusion.csv
[ -e res_ce_all/IGV_skiptics.csv ] && rm res_ce_all/IGV_skiptics.csv
[ -e res_ce_all/IGV_unique_ce_extension.csv ] && rm res_ce_all/IGV_unique_ce_extension.csv
[ -e res_ce_all/IGV_unique_ce_inclusion.csv ] && rm res_ce_all/IGV_unique_ce_inclusion.csv
[ -e res_ce_all/non_ce_events.csv ] && rm res_ce_all/non_ce_events.csv

[ -e res_ce_all/ce_extension_coord_sashimi.bed ] && rm res_ce_all/ce_extension_coord_sashimi.bed
[ -e res_ce_all/ce_extension_coord_uniq_sashimi.bed ] && rm res_ce_all/ce_extension_coord_uniq_sashimi.bed
[ -e res_ce_all/ce_extension_coord.bed ] && rm res_ce_all/ce_extension_coord.bed
[ -e res_ce_all/ce_inclusion_coord.bed ] && rm res_ce_all/ce_inclusion_coord.bed

[ -e res_ce_all/ce_all_scan_unique_range.bed ] && rm res_ce_all/ce_all_scan_unique_range.bed
[ -e res_ce_all/ce_inclusion_coord_uniq_sashimi.bed ] && rm res_ce_all/ce_inclusion_coord_uniq_sashimi.bed
[ -e res_ce_all/ce_extension_coord_uniq_sashimi.bed ] && rm res_ce_all/ce_extension_coord_uniq_sashimi.bed
[ -e res_ce_all/IGV_skiptics.csv ] && rm res_ce_all/IGV_skiptics.csv
[ -e res_ce_all/IGV_unique_ce_inclusion.csv ] && rm res_ce_all/IGV_unique_ce_inclusion.csv
[ -e res_ce_all/IGV_unique_ce_extension.csv ] && rm res_ce_all/IGV_unique_ce_extension.csv
[ -e res_ce_all/ex1_ex2_ce.txt ] && rm res_ce_all/ex1_ex2_ce.txt

[ -e res_ce_all/ce_all_scan_range_junctions.bed ] && rm res_ce_all/ce_all_scan_range_junctions.bed
[ -e res_ce_all/ce_all_scan_range.bed ] && rm res_ce_all/ce_all_scan_range.bed


[ -e res_ce_all/ce_extension_coord_uniq.bed ] && rm res_ce_all/ce_extension_coord_uniq.bed
[ -e res_ce_all/ce_extension_coord_repeated.bed ] && rm res_ce_all/ce_extension_coord_repeated.bed


[ -e res_ce_all/CE_all_FINAL_STATS.txt ] && rm res_ce_all/CE_all_FINAL_STATS.txt

[ -e res_ce_all/IGV_ce_inclusion.csv ] && rm res_ce_all/IGV_ce_inclusion.csv


[ -e res_ce_all/ce_inclusion_coord_uniq.bed ] && rm res_ce_all/ce_inclusion_coord_uniq.bed
[ -e res_ce_all/ce_inclusion_coord_repeated.bed ] && rm res_ce_all/ce_inclusion_coord_repeated.bed

[ -e res_ce_all/non_ce_events.txt ] && rm res_ce_all/non_ce_events.txt
[ -e res_ce_all/non_ce_events.csv ] && rm res_ce_all/non_ce_events.csv

[ -e res_ce_all/IGV_problematic_junctions.csv ] && rm res_ce_all/IGV_problematic_junctions.csv
[ -e res_ce_all/problematic_junctions.txt ] && rm res_ce_all/problematic_junctions.txt


events_bed_create_flg=1
if [ $events_bed_create_flg -eq 1 ]
then
	[ "$(ls -A event_bedfiles/)" ] && rm event_bedfiles/*.*
	[ -e res_ce_all/EnsDB_tx_not_found.csv ] && rm res_ce_all/EnsDB_tx_not_found.csv
	Rscript TxEnsDB103_layeredV5.R $inpfile principal_txs.csv #$tx_lst
fi

readarray -t all_csv_data < $inpfile #for saving gene_ids as well to generate sashimi compatible csv file
event_i=0

samples=$(ls event_bedfiles/temp_*.bed)

overlap_allowed=5 #overlap allowed (intron/exon) for ce events

for sample in $samples
do
	#read line from csv file
	line_csv=${all_csv_data[$event_i]}
	((event_i=event_i+1))
	echo processing "$sample"
	#echo processing "$sample" >> res_intronic_range/all_processed.txt
	#get  all exons bed
	#allexons=$(echo "$sample" | cut -d'_' -f2)
	allexons=$(echo "$sample" | cut -d'/' -f2 |cut -d'_' -f2)
	#also get 5UTR and 3UTR bed files
	#utr=$(echo "$allexons" | cut -d'.' -f1)
	#echo $allexons
	#echo $utr
  #gene_name=$(echo "$allexons" | cut -d'.' -f1-2)
  gene_name1=$(echo "$allexons" | cut -d'.' -f1)
  gene_name=$(echo "$gene_name1" | cut -d'-' -f1)

  #echo $gene_name
  #exit 1
  #first sort the bed
	sortBed -i event_bedfiles/$allexons > event_bedfiles/t$allexons
	#Also read Tx Files to retrieve selected Tx - should find better ways
	TxID=$(head -1 event_bedfiles/TxID$allexons | awk '{print $7}')

  #echo got TxID $TxID

  strnd=$(cat "$sample" | awk '{print $6}')
  #echo strand "$strnd"
  cat "$sample"

	#get distance to downstream exon (for ties, report first) from current reference and pick start, end and d
	#ds=$(bedtools closest -a $sample -b t$allexons -t last -D ref -iu -d -s)
	ds=$(bedtools closest -a $sample -b event_bedfiles/t$allexons -s -D a -iu -d -t first )
  #also get distance to upstream exon from current reference and pick start, end and d
  #us=$(bedtools closest -a $sample -b t$allexons -t first -D ref -id -d -s)
	us=$(bedtools closest -a $sample -b event_bedfiles/t$allexons -s  -D a -id -d -t first)

  #echo ds $ds
	#echo us $us

  #get up/dn exon lengths
	upexonl=$(echo "$us" | awk '{print $10}')
	dnexonl=$(echo "$ds" | awk '{print $10}')

  #get up and down stream exon numbers
	upexon=$(echo "$us" | awk '{print $11}')
	dnexon=$(echo "$ds" | awk '{print $11}')

	#get overlap with up/dn exon - 0 means complete overlap which is assumed as exon skip event
	dsovlp1=$(echo "$ds" | awk '{print $13}')
	usovlp1=$(echo "$us" | awk '{print $13}')
	#take absolute values
	dsovlp=${dsovlp1#-}
	usovlp=${usovlp1#-}


  diff_exon=$(($upexon-$dnexon))
	#take absolute value
	diff_exon_abs=${diff_exon#-}
  if [[ "$upexon" == "$dnexon" ]] && [ $upexonl -ne $dnexonl ] # condition [ $upexonl -ne $dnexonl ] REMOVES EXON_SKIP events
	#if ([ "$upexon" == "$dnexon" ] || [ "$diff_exon" == 1 ]) && ([ "$diff21" -eq 1 ] || [ "$diff12" -eq 1 ] || [ "$ds12" -eq 1 -] || [ "$us12" -eq -1 ])
	then
		echo $gene_name is a CE event >> res_ce_all/ex1_ex2_ce.txt  #>> res_intronic_range/intronic_range_events.txt
		echo ds $ds >> res_ce_all/ex1_ex2_ce.txt
		echo us $us >> res_ce_all/ex1_ex2_ce.txt

		#get distance to downstream exon (for ties, report first) from current reference and pick start, end and d
		#ds=$(bedtools closest -a $sample -b t$allexons -t first -io -D ref -iu -d)
		dsn=$(bedtools closest -a $sample -b event_bedfiles/t$allexons -t first -io -D a -iu -d -s)
		#also get ds overlap
		dso=$(bedtools closest -a $sample -b event_bedfiles/t$allexons -t first -D a -iu -d -s)
		#get start and end of dso
		dso1=$(echo $dso | awk '{print $2}')
		dso2=$(echo $dso | awk '{print $9}')
		dsodiff=$(($dso2-$dso1))
		echo got dsodiff $dsodiff
	  #also get distance to upstream exon from current reference and pick start, end and d
	  #us=$(bedtools closest -a $sample -b t$allexons -t last -io -D ref -id -d)
		usn=$(bedtools closest -a $sample -b event_bedfiles/t$allexons -t last -io -D a -id -d -s)
		#also get us overlap
		uso=$(bedtools closest -a $sample -b event_bedfiles/t$allexons -t last -D a -id -d -s)
		#get up/dn exon lengths
		usnexonl=$(echo "$usn" | awk '{print $10}')
		dsnexonl=$(echo "$dsn" | awk '{print $10}')

		#echo dso $dso and dso1 $dso1 and dso2 $dso2 and dsodiff is $dsodiff
		#echo uso $uso
	  echo now dsn $dsn >> res_ce_all/ex1_ex2_ce.txt
		echo now usn $usn >> res_ce_all/ex1_ex2_ce.txt
    echo dsodiff $dsodiff >> res_ce_all/ex1_ex2_ce.txt
		echo up exon length $usnexonl and dsnexonl $dsnexonl

		if [ $usnexonl == "." ] || [ $dsnexonl == "." ]
		then
			#echo $gene_name is ce with up exon length $usnexonl and dsnexonl $dsnexonl , LEAVING it for now PLEASE MAKE SURE IT IS PROPERLY IDENTIFIED >> res_intronic_range/problematic_junctions.txt
			echo $gene_name is ce with up exon length $usnexonl and dsnexonl $dsnexonl , LEAVING it for now PLEASE MAKE SURE IT IS PROPERLY IDENTIFIED >> res_ce_all/problematic_junctions.txt
			#cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_intronic_range/IGV_problematic_junctions.csv
			cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_ce_all/IGV_problematic_junctions.csv

		else
    if [[ "$strnd" == "-" ]]
    then
      if [[ $dsodiff == 1 ]]
      then
				if [[ $upexonl -lt 60 ]]
				then
        	upex=$(echo "$us" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8+1,$9,$10,$5,$6,gen}')
				else
					upex=$(echo "$us" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8+1,$8+60+1,$10,$5,$6,gen}')
				fi
        #echo upex $upex
        #ce=$(echo "$us" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$2+1,$8+1,$4,$5,$6,gen}')
        ce=$(echo "$us" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$2+2,$8,$10,$5,$6,gen}')
        #echo ce $ce
				if [[ $dnexonl -lt 60 ]]
				then
					echo need to change this part of the code as i do not have second end of this down stream exon - strand 1
        	dsex=$(echo "$ds" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$2-60+1,$2+1,$10,$5,$6,gen}')
				else
					dsex=$(echo "$ds" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$2-60+1,$2+1,$10,$5,$6,gen}')
				fi
        #echo dsex $dsex
        #also get data for start and stop codons - only take -/+ 60 not the whole exon
        #upexc=$(echo "$us" | awk  'BEGIN {OFS="\t"} {print $1,$8+1,$8+60+1,$4,$5,$6}')
        #dsexc=$(echo "$ds" | awk 'BEGIN {OFS="\t"} {print $1,$2-60+1,$2+1,$4,$5,$6}')
      else
				if [[ $upexonl -lt 60 ]]
				then
					echo need to change this part of the code as i do not have first end of this up stream exon - strand 2
        	upex=$(echo "$us" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$3,$3+60,$10,$5,$6,gen}')
				else
					upex=$(echo "$us" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$3,$3+60,$10,$5,$6,gen}')
				fi

        #echo upex $upex
        ce=$(echo "$ds" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$9+1,$3-1,$10,$5,$6,gen}')
        echo ce $ce
				if [[ $dnexonl -lt 60 ]]
				then
        	dsex=$(echo "$ds" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8,$9,$10,$5,$6,gen}')
				else
					dsex=$(echo "$ds" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$9-60,$9,$10,$5,$6,gen}')
				fi
        #echo dsex $dsex
        #also get data for start and stop codons - only take -/+ 60 not the whole exon
        #upexc=$(echo "$us" | awk  'BEGIN {OFS="\t"} {print $1,$3,$3+60,$4,$5,$6}')
        #dsexc=$(echo "$ds" | awk 'BEGIN {OFS="\t"} {print $1,$9-60,$9,$4,$5,$6}')
				#CE scan region
				junc1=$(echo "$ds" | awk '{print $2}')
				junc2=$(echo "$ds" | awk '{print $3}')
				upex1=$(echo "$us" | awk '{print $8}')
				dsex2=$(echo "$ds" | awk '{print $9}')

				#this works better
				if [[ $usd_abs -lt $dsd_abs ]]
				then
					ce_scan1=$(echo "$ds" | awk -v c1=$junc1 -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,c1,$9,$10,$5,$6,gen}')
				else
					ce_scan1=$(echo "$us" | awk -v c1=$junc2 -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,c1,$8,$10,$5,$6,gen}')
				fi



      fi
    else
      #if [[ $dsodiff == 1 ]]
      #then
			if [[ $upexonl -lt 60 ]]
			then
        upex=$(echo "$us" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8,$9,$10,$5,$6,gen}')
			else
				upex=$(echo "$us" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$9-60,$9,$10,$5,$6,gen}')
			fi
        #ce=$(echo $ds | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$2,$8,$4,$5,$6,gen}')
        use=$(echo "$us" | awk '{print $9}')
        ce1=$(echo "$dsn" | awk -v g="$use" 'BEGIN {OFS="\t"} {print $1,g,$8-1,$10,$5,$6}')
        ce=$(echo "$ce1" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$2,$3,$10,$5,$6,gen}')


        #ce=$(echo $ds | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$2+1,$8-1,$4,$5,$6,gen}')
        #echo cc $ce
			if [[ $dsnexonl -lt 60 ]]
			then
        dsex=$(echo "$dsn" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8-1,$9,$10,$5,$6,gen}')
			else
				dsex=$(echo "$dsn" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8-1,$8+60-1,$10,$5,$6,gen}')
			fi

			#CE scan region
			junc1=$(echo "$ds" | awk '{print $2}')
			junc2=$(echo "$ds" | awk '{print $3}')
			#gave problem with single bp differences
			upex2=$(echo "$us" | awk '{print $9}')
			dsex1=$(echo "$ds" | awk '{print $8}')
			#this works better
			if [[ $usd_abs -lt $dsd_abs ]]
			then
				ce_scan1=$(echo "$ds" | awk -v c1=$junc2 -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,c1,$8,$10,$5,$6,gen}')
			else
				ce_scan1=$(echo "$us" | awk -v c1=$junc1 -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,c1,$9,$10,$5,$6,gen}')
			fi
  fi


	cat $sample | awk -v gene=$gene_name -v TX=$TxID 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene,TX}' >> res_ce_all/IGV_ce_inclusion.csv
  #echo for sample
  #cat "$sample"
  echo finally upex is "$upex"
  echo finally ce is "$ce"
  echo finally dsex is "$dsex"

	#NOW REMOVE ANY EVENT FOR WHICH START>END
	upex_st=$(echo "$upex" | awk '{print $2}')
	upex_en=$(echo "$upex" | awk '{print $3}')

	ce_st=$(echo "$ce" | awk '{print $2}')
	ce_en=$(echo "$ce" | awk '{print $3}')

	dsex_st=$(echo "$dsex" | awk '{print $2}')
	dsex_en=$(echo "$dsex" | awk '{print $3}')

	if [ $upex_st -lt $upex_en ] && [ $ce_st -lt $ce_en ] && [ $dsex_st -lt $dsex_en ]
	then

  #finally paste three segments to a bed file
  #first remove if any such file exists
  rm -f "$gene_name"_nt.bed
  echo "$upex">>"$gene_name"_nt.bed #all_junctions.bed
  echo "$ce">>"$gene_name"_nt.bed #all_junctions.bed
  echo "$dsex">>"$gene_name"_nt.bed #all_junctions.bed

  #cat "$gene_name"_nt.bed >> res_intronic_range/intronic_range.bed

	rm "$gene_name"_nt.bed

#else
	#echo "$upex">>res_intronic_range/intronic_range_star_gt_end.bed
	#echo "$ce">>res_intronic_range/intronic_range_star_gt_end.bed
	#echo "$dsex">>res_intronic_range/intronic_range_star_gt_end.bed

fi

	#also coe ce boundary range
	echo "$upex">>res_ce_all/ce_all_scan_range_junctions.bed
	echo "$ce_scan1">>res_ce_all/ce_all_scan_range_junctions.bed
	echo "$dsex">>res_ce_all/ce_all_scan_range_junctions.bed

	echo "$ce_scan1">>res_ce_all/ce_all_scan_range.bed #all_junctions.bed


fi

elif [ "$diff_exon_abs" -eq 1 ] && [ $dsovlp -le 5 ] && [ $usovlp -le 5 ] #deals with exon_joining events few bp of
then
	cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_ce_all/IGV_skiptics.csv
	#cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_intronic_range/non_ce_events.csv

	echo $line_csv | awk 'BEGIN {FS=","} {print $1,$2,$3,$4,$5,$6,$7}' | awk 'BEGIN {OFS=","} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_all/non_ce_events.csv

elif [ "$diff_exon_abs" -eq 1 ]
then
  #THESE ARE CE INTRONIC events
  #CE event with difference of one exon
  echo $gene_name is a CE event #>> res_intronic_range/intronic_range_events.txt #>>ce_intronic.txt

  #get start and end of dso
  usd=$(echo $us | awk '{print $13}')
  usd_abs=${usd#-}
  dsd=$(echo $ds | awk '{print $13}')
  dsd_abs=${dsd#-}

  #echo usd_abs $usd_abs
  #echo dsd_abs $dsd_abs
  #echo strand "$strnd"


  #record 1 (downstream exon). chr# end_ce+d       end_ce+d_nuc_sz      strand  gene_name
	#ADDED TO REMOVE ALL EXONS THAT LIE INSIDE AN INTRON
	if [ $usd_abs -le 5 ] || [ $dsd_abs -le 5 ]
	then
  if [[ "$strnd" == "+" ]]
  then
    #echo here
    #if [[ "$usd_abs" -le "$dsd_abs" ]]
    #then
    #  echo and here
		if [[ $upexonl -lt 60 ]]
		then
      upex=$(echo "$us" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8,$9,$10,$5,$6,gen}')
		else
			upex=$(echo "$us" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$9-60,$9,$10,$5,$6,gen}')
		fi

      #echo upex $upex
      #ce=$(echo "$us" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$2+1,$8+1,$4,$5,$6,gen}')
      #ce=$(echo "$ds" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$2+1,$8-1,$4,$5,$6,gen}')
      use=$(echo "$us" | awk '{print $9}')
      ce1=$(echo "$ds" | awk -v g="$use" 'BEGIN {OFS="\t"} {print $1,g,$8-1,$10,$5,$6}')
      ce=$(echo "$ce1" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$2,$3,$10,$5,$6,gen}')
      #echo ce $ce
			if [[ $dnexonl -lt 60 ]]
			then
      	dsex=$(echo "$ds" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8-1,$9,$10,$5,$6,gen}')
			else
				dsex=$(echo "$ds" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8-1,$8+60-1,$10,$5,$6,gen}')
			fi

			#CE scan region
			junc1=$(echo "$ds" | awk '{print $2}')
			junc2=$(echo "$ds" | awk '{print $3}')
			#gave problem with single bp differences
			upex2=$(echo "$us" | awk '{print $9}')
			dsex1=$(echo "$ds" | awk '{print $8}')
			#this works better
			if [[ $usd_abs -lt $dsd_abs ]]
			then
				ce_scan1=$(echo "$ds" | awk -v c1=$junc2 -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,c1,$8,$10,$5,$6,gen}')
			else
				ce_scan1=$(echo "$us" | awk -v c1=$junc1 -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,c1,$9,$10,$5,$6,gen}')
			fi

fi
if [[ "$strnd" == "-" ]]
then

    #if [[ "$usd" -lt 0 ]]
    #then
    #  echo and here
		if [[ $upexonl -lt 60 ]]
		then
      upex=$(echo "$us" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8-1,$9,$10,$5,$6,gen}')
		else
			upex=$(echo "$us" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8-1,$8+60-1,$10,$5,$6,gen}')
		fi

      #echo upex $upex
      #ce=$(echo "$us" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$2+1,$8+1,$4,$5,$6,gen}')
      #ce=$(echo "$ds" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$2+1,$8-1,$4,$5,$6,gen}')
      use=$(echo "$ds" | awk '{print $9}')
      ce1=$(echo "$us" | awk -v g="$use" 'BEGIN {OFS="\t"} {print $1,g,$8-1,$10,$5,$6}')
      ce=$(echo "$ce1" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$2,$3,$10,$5,$6,gen}')
      #echo ce $ce
			if [[ $dnexonl -lt 60 ]]
			then
      	dsex=$(echo "$ds" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8-1,$9,$10,$5,$6,gen}')
			else
				dsex=$(echo "$ds" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$9-60,$9,$10,$5,$6,gen}')
			fi
			#CE scan region
			junc1=$(echo "$ds" | awk '{print $2}')
			junc2=$(echo "$ds" | awk '{print $3}')
			upex1=$(echo "$us" | awk '{print $8}')
			dsex2=$(echo "$ds" | awk '{print $9}')

			#this works better
			if [[ $usd_abs -lt $dsd_abs ]]
			then
				ce_scan1=$(echo "$ds" | awk -v c1=$junc1 -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,c1,$9,$10,$5,$6,gen}')
			else
				ce_scan1=$(echo "$us" | awk -v c1=$junc2 -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,c1,$8,$10,$5,$6,gen}')
			fi

  fi
cat $sample | awk -v gene=$gene_name -v TX=$TxID 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene,TX}' >> res_ce_all/IGV_ce_inclusion.csv
  echo for sample
  cat $sample
  echo finally upex is "$upex"
  echo finally ce is "$ce"
  echo finally dsex is "$dsex"

	#NOW REMOVE ANY EVENT FOR WHICH START>END
	upex_st=$(echo "$upex" | awk '{print $2}')
	upex_en=$(echo "$upex" | awk '{print $3}')

	ce_st=$(echo "$ce" | awk '{print $2}')
	ce_en=$(echo "$ce" | awk '{print $3}')

	dsex_st=$(echo "$dsex" | awk '{print $2}')
	dsex_en=$(echo "$dsex" | awk '{print $3}')

	if [ $upex_st -lt $upex_en ] && [ $ce_st -lt $ce_en ] && [ $dsex_st -lt $dsex_en ]
	then

  #finally paste three segments to a bed file
  #first remove if any such file exists
  rm -f "$gene_name"_nt.bed
  echo "$upex">>"$gene_name"_nt.bed #all_junctions.bed
  echo "$ce">>"$gene_name"_nt.bed #all_junctions.bed
  echo "$dsex">>"$gene_name"_nt.bed #all_junctions.bed
  #cat "$gene_name"_nt.bed >> res_intronic_range/intronic_range.bed

	rm "$gene_name"_nt.bed
#else
	#echo "$upex">>res_intronic_range/intronic_range_star_gt_end.bed
	#echo "$ce">>res_intronic_range/intronic_range_star_gt_end.bed
	#echo "$dsex">>res_intronic_range/intronic_range_star_gt_end.bed

fi

	#also coe ce boundary range
	echo "$upex">>res_ce_all/ce_all_scan_range_junctions.bed
	echo "$ce_scan1">>res_ce_all/ce_all_scan_range_junctions.bed
	echo "$dsex">>res_ce_all/ce_all_scan_range_junctions.bed

	echo "$ce_scan1">>res_ce_all/ce_all_scan_range.bed #all_junctions.bed
else
	#BOTH ENDS OF THESE EVENTS LIE INSIDE INTRON
	echo $gene_name has both ends lying inside intron #>> res_intronic_range/intronic_range_events.txt #>>ce_intronic.txt
	#echo $gene_name is unknown event type for gene "$gene_name" >> res_intronic_range/non_ce_events.txt
	echo $gene_name is unknown event type for gene "$gene_name" >> res_ce_all/non_ce_events.txt
	#cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_intronic_range/non_ce_events.csv
	echo $line_csv | awk 'BEGIN {FS=","} {print $1,$2,$3,$4,$5,$6,$7}' | awk 'BEGIN {OFS=","} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_all/non_ce_events.csv

fi
#fi
else
  #echo $gene_name is unknown event type for gene "$gene_name" >> res_intronic_range/non_ce_events.txt
	echo $gene_name is unknown event type for gene "$gene_name" >> res_ce_all/non_ce_events.txt
	#cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_intronic_range/non_ce_events.csv
	echo $line_csv | awk 'BEGIN {FS=","} {print $1,$2,$3,$4,$5,$6,$7}' | awk 'BEGIN {OFS=","} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_all/non_ce_events.csv
fi

#rm $sample1
#rm $sample_nt.bed
#rm $sample1-3UTR.bed
#rm $sample1-5UTR.bed
#rm $allexons
#rm t$allexons
done
########## THIS SECTION WRITES NEW FILE res_ce_all/ce_all_scan_intron.bed to scan whole intron between two exons for probable ce events
nrecrds=$(cat res_ce_all/ce_all_scan_range_junctions.bed | wc -l)
readarray -t all_data < res_ce_all/ce_all_scan_range_junctions.bed
i=0
while [ $i -lt $nrecrds ]
  do
          line11=${all_data[$i]}
          line12=${all_data[$(($i+1))]}
          line13=${all_data[$(($i+2))]}
     			((i=i+3))
					strand=$(echo "$line11" | awk 'BEGIN {FS="\t"} {print $6}')
					if [ "$strand" == "+" ]
					then
						#st=$(echo $line11 | awk 'BEGIN {FS="\t"} {print $2}')
						#end=$(echo $line11 | awk 'BEGIN {FS="\t"} {print $3}')
						#off=$(($end-$st))
						#st1=$(echo $line13 | awk 'BEGIN {FS="\t"} {print $2}')
						#end1=$(echo $line13 | awk 'BEGIN {FS="\t"} {print $3}')
						#off1=$(($end1-$st1))

	          echo "$line11" "$line13" | awk -v o=$off 'BEGIN {OFS="\t"} {print $1,$3,$9,$11,$12,$13,$14}' >> res_ce_all/ce_all_scan_intron.bed
					elif [ "$strand" == "-" ]
					then
						#st=$(echo $line13 | awk 'BEGIN {FS="\t"} {print $2}')
						#end=$(echo $line13 | awk 'BEGIN {FS="\t"} {print $3}')
						#off=$(($end-$st))
						#st1=$(echo $line11 | awk 'BEGIN {FS="\t"} {print $2}')
						#end1=$(echo $line11 | awk 'BEGIN {FS="\t"} {print $3}')
						#off1=$(($end1-$st1))

						echo "$line13" "$line11" | awk 'BEGIN {OFS="\t"} {print $1,$3,$9,$11,$12,$13,$14}' >> res_ce_all/ce_all_scan_intron.bed
					fi
 done

################## END SECTION WRITES NEW FILE res_ce_all/ce_all_scan_intron.bed
#also copy for coverage calculations and avoid intronic_range calculations
cat res_ce_all/ce_all_scan_range.bed  > res_ce_all/ce_all_scan_unique_range.bed
fi

#CRYPTICS SECTION ENDS HERE

#EVENT COVERAGE CALCULATIONS STARTS HERE
#Step 3. Now calculate coverages for each ce-range from all TDP samples and aggregate them for each junction
#flag for step 3
#FILES NEEDED FOR THIS STEP: res_ce_all/ce_all_scan_unique_range.bed, bam_files/*.bam, bam_beds/*.bed for all samples in bam_files folder
coverages_bed_files_create_flg=0
#if [ $coverages_bed_files_create_flg -eq 1 ]
if [ $1 -eq 6 ] || [ $1 -eq 8 ]
then
echo NOW STARTED COVERAGE CALCULATIONS, TAKES LONG TIME OFTEN HOURS/DAYS FOR ">" 200 EVENTS
mkdir -p coverages
########NOW IDENTIFY CE coordinates
#BAMSAMPLES_ctrl="bamsamples_ctrl.txt"
#BAMSAMPLES_tdp="bamsamples_tdp.txt"
#bamsamples_tdp=$(cut -f1 "$BAMSAMPLES_tdp")
#bamsamples_tdp=$(ls bam_files/*.bam | cut -d'/' -f2 | cut -d'.' -f1)
bamsamples_tdp=$(ls *.bam | cut -d'.' -f1)
#bamsamples_ctrl=$(cut -f1 "$BAMSAMPLES_ctrl")
#bedtools coverage -a "$gene_name"_ce.bed -b /Volumes/SYEDSHAH/MichaelLab/BAM_FILES/TDP_100bp/Short-read_bam_files/TDP43-E_S5.pass2Aligned.sortedByCoord.out.bam.bed -d -s > $sample.bed.coverage
#bedtools coverage -a "$gene_name"_ce.bed -b /Volumes/SYEDSHAH/MichaelLab/BAM_FILES/TDP_100bp/Short-read_bam_files/bed/TDP43-E_S5.pass2Aligned.sortedByCoord.out.bam.bed -d -s > $sample.bed.coverage
#first do coverage calculations for each ce
i=1 #counter to use for calculations
last_genename=""
while IFS='\t' read -r line #col2 col3 col4 col5 col6 col7
do
	#now do calculations fot this ce across all replicates
	#invert the range if it is not already
	genename=$(echo $line | awk '{print $7}')
	if [[ $last_genename == $genename ]]
	then
		i=$((i+1))
		last_genename=$genename
	else
		i=1
		last_genename=$genename
	fi


	chrm_start=$(echo $line | awk '{print $2}')
	chrm_end=$(echo $line | awk '{print $3}')
	#fnm_ctrl="$genename"_"$i"_ctrl.cov.bed
	fnm_tdp="$genename"_"$i"_TDP.cov.bed
	#echo fnm_ctrl $fnm_ctrl
	echo fnm_tdp $fnm_tdp
	if [[ "$chrm_start" -lt "$chrm_end" ]]
	then
		echo $line | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' > temp_coord.bed
	else
		echo $line | awk 'BEGIN {OFS="\t"} {print $1,$3,$2,$4,$5,$6,$7}' > temp_coord.bed
	fi

	#and now for TDP43
	echo coverage calculations for tdp sampels for:
	cat temp_coord.bed
	for bamsamp in $bamsamples_tdp
	do
		bedtools coverage -sorted -g bam_beds/"$bamsamp".chromosomes.txt -a temp_coord.bed -b bam_beds/"$bamsamp"-sorted.bam.bed -d > temp_coord_"$bamsamp"_tdp.bed.cov
		#bedtools coverage -a temp_coord.bed -b bam_beds/"$bamsamp".bam.bed -d > temp_coord_"$bamsamp"_tdp.bed.cov
		#bedtools coverage -a "$gene_name"_ce.bed -b TDP43-F_S6.bam.bed -d -s > "$gene_name"_ce_TDP43-F_S6.bed.cov
		#bedtools coverage -a "$gene_name"_ce.bed -b TDP43-G_S7.bam.bed -d -s > "$gene_name"_ce_TDP43-G_S7.bed.cov
	done
	#finally sum all coverages across replicates for this ce
	samples1=""
	for bamsamp1 in $bamsamples_tdp
	do
	samples1+=" "temp_coord_"$bamsamp1"_tdp.bed.cov
	done
	echo got all coverage files for tdp $samples1
	read -a samples2 <<< "$samples1"

	paste "${samples2[@]}" | awk -v numFiles=${#samples2[@]} -v OFS='\t' '
	  {
	    row = sep = ""
	    for(i=1; i < NF/numFiles; ++i) { row = row sep $i; sep = OFS }
	    sum = $(NF/numFiles) # last header col. / (1st) data col. to sum
	    if (NR > 0) { for(i=2; i<=numFiles; ++i) sum += $(NF/numFiles * i) } # add other cols.
	    printf "%s%s%s\n", row, OFS, sum
	  }
	' > coverages/$fnm_tdp
	#also remove temp_coord_ for each sample done
	for bamsamp in $bamsamples_tdp
	do
		rm temp_coord_"$bamsamp"_tdp.bed.cov
	done

#done < res_ce_all/ce_all_scan_unique_range.bed #NEW APPROACH SCANNING WHOLE INTRON
done < res_ce_all/ce_all_scan_intron.bed #res_ce_all/ce_all_scan_unique_range.bed
fi
#fi
#EVENT COVERAGE CALCULATIONS ENDS HERE
#Step 3 ends

#CE IDENTIFICATION STARTS HERE

#Step 4. using majiq_coverages_automateV1.R script identify final ce_extension_coord.bed and ce_extension_coord.bed files
fact=3
total_splicing_events=$(cat $inpfile | wc -l)
unknow_events=0
[ -e res_ce_all/non_ce_events.csv ] && unknow_events=$(cat res_ce_all/non_ce_events.csv | wc -l)
not_found=0
[ -e res_ce_all/EnsDB_tx_not_found.csv ] && not_found=$(cat res_ce_all/EnsDB_tx_not_found.csv | wc -l)
problematic=0
[ -e res_ce_all/IGV_problematic_junctions.csv ] && problematic=$(cat res_ce_all/IGV_problematic_junctions.csv  | wc -l)
#total_events=$(($all_ce+$unknow_events+$not_found+$problematic))
unseccesful_r1=0
[ -e res_ce_all/skipped_ce.csv ] && unseccesful_r1=$(cat res_ce_all/skipped_ce.csv|wc -l)
unseccesful_r=$(($unseccesful_r1/3))
echo Total splicing events read are: $total_splicing_events >> res_ce_all/CE_all_FINAL_STATS.txt
echo Out of these $total_splicing_events total events >> res_ce_all/CE_all_FINAL_STATS.txt
echo Events not found in EnsDB are: $not_found, please see res_ce_all/EnsDB_tx_not_found.csv file >> res_ce_all/CE_all_FINAL_STATS.txt
echo Events that are not CE: $unknow_events , please see res_ce_all/non_ce_events.txt and res_ce_all/non_ce_events.csv file >> res_ce_all/CE_all_FINAL_STATS.txt
echo Events that are somewhat problematic: $problematic , please see res_ce_all/problematic_junctions.txt and res_ce_all/IGV_problematic_junctions.csv files >> res_ce_all/CE_all_FINAL_STATS.txt
ce_boundary_events=$(($total_splicing_events-$not_found-$unknow_events-$problematic))
echo '#######################################' >> res_ce_all/CE_all_FINAL_STATS.txt

if [ $1 -eq 6 ] || [ $1 -eq 7 ] || [ $1 -eq 8 ]
#if [ $ce_boundary_events -gt 0 ]
then
  echo Now starting ce_boundary calculations for a remaining total of: $ce_boundary_events events, BY INVOKING R SCRIPT  >> res_ce_all/CE_all_FINAL_STATS.txt
  #echo Identifying ce_boundary coordinates
  echo IDENTIFYING CE BOUNDARIES BY CALLING Auto_CoverV4_layered_intronV2.R script

  Rscript Auto_CoverV4_layered_intronV2.R #Auto_CoverV2.R
  echo BACK FROM CE_BOUNDARY CALCULATIONS
  echo '################################################' >> res_ce_all/CE_all_FINAL_STATS.txt
#Rscript majiq_coverages_automateV1.R #majiq_coverages_automateV1_227.R #majiq_coverages_automateV1.R
#Step 4 ENDS


#Step 5. Finally identify nt and AA sequences

#get nt and aa for ce boundaries

#First remove duplicates from ce' boundary coordinates
#NEW CODE STARTS HERE
nrecrds=0
nrecrdst=0
[ -e res_ce_all/ce_inclusion_coord.bed ] && nrecrds=$(cat res_ce_all/ce_inclusion_coord.bed | wc -l) && nrecrdst=$(($nrecrds/3))


if [ $nrecrds -ne 0 ]
then
	echo Back From R Session, Now checking "FOR DUPLICATES" for a total of ce_inclusion events: $nrecrdst  >> res_ce_all/CE_all_FINAL_STATS.txt
#echo Checking any Repeated ce_boundary events >> res_ce_all/CE_all_FINAL_STATS.txt
readarray -t all_data < res_ce_all/ce_inclusion_coord.bed
#also read ce_inclusion_coord_sashimi.bed to get TxID
readarray -t sashimi_data < res_ce_all/ce_inclusion_coord_sashimi.bed
#also read csv file
readarray -t csv_data < res_ce_all/IGV_R_returned_ce_inclusion.csv
csvi=0
#echo total ce_boundary events - excluding IR_events - read are $nrecrdst
i=0
while [ $i -lt $nrecrds ]
  do
          line11=${all_data[$i]}
          line12=${all_data[$(($i+1))]}
          line13=${all_data[$(($i+2))]}
					#TxID
					#for sashimi
					lines1=${sashimi_data[$i]}
          lines2=${sashimi_data[$(($i+1))]}
          lines3=${sashimi_data[$(($i+2))]}

					#txid_recrd=${sashimi_data[$i]}
					#echo got txid_recrd $txid_recrd
					#txid=$(echo $txid_recrd | awk 'BEGIN {OFS="\t"} {print $8}')
					#echo got txid $txid
#         echo i is $i, $line11, $line12, $line13
					csv_ln=${csv_data[$csvi]}
					((csvi=csvi+1))

     #echo "Welcome $i times"
     ((i=i+3))
     j=$i
     #now go through rest of the data
     flg=0
        while [ $j -lt $nrecrds ]
        do
                line21=${all_data[$j]}
          line22=${all_data[$(($j+1))]}
          line23=${all_data[$(($j+2))]}
          if [ "${line11[*]}" == "${line21[*]}" ] && [ "${line12[*]}" == "${line22[*]}" ] && [ "${line13[*]}" == "${line23[*]}" ]
          then
            #echo i is $i, $line11, $line12, $line13
            #echo j is $j, $line21, $line22, $line23
            flg=1
          fi
          ((j=j+3))
        done
        if [ $flg -eq 0 ]
				then
                 echo $line11 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_all/ce_inclusion_coord_uniq.bed
                 echo $line12 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_all/ce_inclusion_coord_uniq.bed
                 echo $line13 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_all/ce_inclusion_coord_uniq.bed

								 #also save bed file for sashimi plots
								 #echo $line11 | awk -v TX=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,TX}' >> res_ce_all/ce_inclusion_coord_uniq_sashimi.bed
                 #echo $line12 | awk -v TX=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,TX}' >> res_ce_all/ce_inclusion_coord_uniq_sashimi.bed
                 #echo $line13 | awk -v TX=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,TX}' >> res_ce_all/ce_inclusion_coord_uniq_sashimi.bed
								 echo $lines1 | awk  'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' >> res_ce_all/ce_inclusion_coord_uniq_sashimi.bed
                 echo $lines2 | awk  'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' >> res_ce_all/ce_inclusion_coord_uniq_sashimi.bed
								 #Also SAVING CE coordinates

								 echo $lines2 | awk  'BEGIN {OFS="\t"} {print $1,$2,$3,$3-$2,$5,$6,$7,$8}' >> res_ce_all/ce_inclusion_coord_only.bed
                 echo $lines3 | awk  'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' >> res_ce_all/ce_inclusion_coord_uniq_sashimi.bed

								 #also save csv file
								 echo $csv_ln >> res_ce_all/IGV_unique_ce_inclusion.csv
							 else
			 					echo $line11 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_all/ce_inclusion_coord_repeated.bed
			 					echo $line12 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_all/ce_inclusion_coord_repeated.bed
			 					echo $line13 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_all/ce_inclusion_coord_repeated.bed

			          fi


 done

 nrecrds_uniqt=$(cat res_ce_all/ce_inclusion_coord_uniq.bed | wc -l)
 nrecrds_uniq=$(($nrecrds_uniqt/3))
 repeated_ce_boundary_events=$(($nrecrdst-$nrecrds_uniq))
 echo Total unique ce_inclusion events are: $nrecrds_uniq, please see res_ce_all/ce_inclusion_coord_uniq.bed and res_ce_all/IGV_unique_ce_inclusion.csv >> res_ce_all/CE_all_FINAL_STATS.txt
 #echo total unique ce_boundary events read are $nrecrds_uniq
 echo Total repeated ce_inclusion events are: $repeated_ce_boundary_events, please see res_ce_all/ce_inclusion_coord_repeated.bed >> res_ce_all/CE_all_FINAL_STATS.txt
 echo total repeated ce_inclusion events were $repeated_ce_boundary_events


 #now copy back data to ce_coord
 #cat res_ce_all/ce_inclusion_coord_uniq.bed	> res_ce_all/ce_extension_coord.bed
 #uniqrecrds1=$(cat res_ce_all/ce_extension_coord.bed | wc -l)
 #uniqrecrds=$(($uniqrecrds1/3))
 #echo total records after duplicate remval are $uniqrecrds
#NEW CODE ENDS HERE
#cat res_ce_all/ce_extension_coord.bed > ce_coord_temp.bed
#cat ce_coord_temp.bed | awk '!a[$0]++'  > res_ce_all/ce_extension_coord.bed
#rm ce_coord_temp.bed

#now call bedtools getfasta function to get nt sequence from reference_genome.fa file

fasta_ce_flag=0
if [ $fasta_ce_flag -eq 1 ]
then
echo started nt seq for res_ce_all/ce_inclusion_coord_uniq.bed
bedtools getfasta -fi /Volumes/SYEDSHAH/MichaelLab/ref_genome/GRCh38.p13.genome.fa -bed res_ce_all/ce_inclusion_coord_uniq.bed -s > ce_inclusion_nt.fasta

#now remove >chr lines from the resulting file
awk '!/^>chr/' ce_inclusion_nt.fasta > ce_all_nt1.fasta
#combine the two files to get desired csv file - chrXX,start,end,genename
paste -d"\t" res_ce_all/ce_inclusion_coord_uniq.bed ce_all_nt1.fasta > res_ce_all/ce_inclusion_nt.bed
#remove temp file
rm ce_all_nt1.fasta
rm ce_inclusion_nt.fasta
#and finally transeq compatible
#awk -F "\t" '{print ">sp|"$7"_"$1"_"$2"_"$3"\n"$8}' res_ce_all/ce_inclusion_nt.bed > res_ce_all/ce_inclusion_transeq_in.fasta
#Add strand infor
awk -F "\t" '{if($6=="+") print ">sp|"$7"_"$1"_"$2"_"$3"_plus""\n"$8; else print ">sp|"$7"_"$1"_"$2"_"$3"_minus""\n"$8}' res_ce_all/ce_inclusion_nt.bed > res_ce_all/ce_inclusion_transeq_in.fasta

#awk -F "\t" '{print ">sp|"$7"-"$1":"$2"-"$3"\n"$8}' res_ce_all/ce_inclusion_nt.bed > ce_inclusion_transeq_in.fasta


echo now doing the 3-frame translation

TEMPFILE=$(echo res_ce_all/ce_inclusion_transeq_in.fasta  | perl -pe 's/\.fasta$/.temp/')
OUTPUTFILE=$(echo res_ce_all/ce_inclusion_transeq_in.fasta  | perl -pe 's/\.fasta$/.trans/')
#echo doing 3-frame translation
# Translate sequence
transeq  -sequence res_ce_all/ce_inclusion_transeq_in.fasta  -outseq "$TEMPFILE" -frame F

# Rename sequence
perl -sape 's/>/>sp|/' "$TEMPFILE" > "$OUTPUTFILE"



#also remove all newlines from the
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' "$OUTPUTFILE" > res_ce_all/CE_INCLUSION_AA.fasta

#also remove "$OUTPUTFILE" file
rm "$OUTPUTFILE"
# Remove temp file
rm "$TEMPFILE"


my_sym='>'
while read -r line
do
    a="${line[@]:0:1}"
    if [[ "$my_sym" == "$a" ]]; then
	line1="${line[@]}"
	#echo $line
    let i=0
    else
       b=$(echo "$line" | tr "*" '\n')
       while IFS= read li; #for ab in "${b[@]}"
	do
	size=${#li}
	if (( $size > 8 )); then
		let "i+=1"
		#let "j+=1"
	echo "$line1"_"$i"
	echo $li
	fi
	#echo $size
	done <<<"$b"


    fi
done < res_ce_all/CE_INCLUSION_AA.fasta > res_ce_all/FINAL_CE_INCLUSION_AA.fasta

#Also concatenate
echo NOW DOING NT AND AA TRANSLATION FOR FUSED COORDINATES FOR CE_INCLUSION
[ -e res_ce_all/ce_inclusion_fused.transeq_in.fasta ] &&  rm res_ce_all/ce_inclusion_fused.transeq_in.fasta


#rm res_ce_all/ce_inclusion_fused.transeq_in.fasta
#ALSO ADD even_number flag for sashimi plots against PEAKS results
event_i=1
cat res_ce_all/ce_inclusion_nt.bed | while read r1; read r2; read r3
do
  lstrnd=$(echo "$r1" | awk '{print $6}')
  if [[ "$lstrnd" == "+" ]]
  then
    #echo "$r1" | awk '{print $8}' #"$r2" "$r3"
    #title=$(echo "$r1" "$r2" "$r3" | awk '{print ">sp|"$23"_"$1"_"$2"_"$3"_"$10"_"$11"_"$18"_"$19"_plus"}')
		#ADD EVENT Number flag
		title=$(echo "$r1" "$r2" "$r3" | awk -v evid=$event_i '{print ">sp|"$23"_"$1"_"$2"_"$3"_"$10"_"$11"_"$18"_"$19"_"evid"_plus"}')
    #echo $title
    echo "$r1" "$r2" "$r3" | awk -v t=$title '{OFS="\t"} {print t"\n"$8$16$24}' >> res_ce_all/ce_inclusion_fused.transeq_in.fasta
  else
    #echo "$r1" "$r2" "$r3"
    #title=$(echo "$r3" "$r2" "$r1" | awk '{print ">sp|"$23"_"$1"_"$18"_"$19"_"$10"_"$11"_"$2"_"$3"_minus"}')
		#Add EVENT Number flag
		title=$(echo "$r3" "$r2" "$r1" | awk -v evid=$event_i '{print ">sp|"$23"_"$1"_"$18"_"$19"_"$10"_"$11"_"$2"_"$3"_"evid"_minus"}')
    #echo $title
    echo "$r1" "$r2" "$r3" | awk -v t=$title '{OFS="\t"} {print t"\n"$8$16$24}' >> res_ce_all/ce_inclusion_fused.transeq_in.fasta
  fi
	((event_i=event_i+1))
done

#now do the 3-frame translation

TEMPFILE=$(echo res_ce_all/ce_inclusion_fused.transeq_in.fasta | perl -pe 's/\.fasta$/.temp/')
OUTPUTFILE=$(echo res_ce_all/ce_inclusion_fused.transeq_in.fasta | perl -pe 's/\.fasta$/.trans/')

# Translate sequence
transeq  -sequence res_ce_all/ce_inclusion_fused.transeq_in.fasta  -outseq "$TEMPFILE" -frame F

# Rename sequence
perl -sape 's/>/>sp|/' "$TEMPFILE" > "$OUTPUTFILE"



#also remove all newlines from the
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' "$OUTPUTFILE" > res_ce_all/CE_INCLUSION_FUSED_AA.fasta

#also remove "$OUTPUTFILE" file
rm "$OUTPUTFILE"
# Remove temp file
rm "$TEMPFILE"


my_sym='>'
while read -r line
do
    a="${line[@]:0:1}"
    if [[ "$my_sym" == "$a" ]]; then
	line1="${line[@]}"
	#echo $line
    let i=0
    else
       b=$(echo "$line" | tr "*" '\n')
       while IFS= read li; #for ab in "${b[@]}"
	do
	size=${#li}
	if (( $size > 8 )); then
		let "i+=1"
		#let "j+=1"
	echo "$line1"_"$i"
	echo $li
	fi
	#echo $size
	done <<<"$b"


    fi
done < res_ce_all/CE_INCLUSION_FUSED_AA.fasta > res_ce_all/FINAL_CE_INCLUSION_FUSED_AA.fasta

else
	echo Back From R Session, total ce_inclusion events: $nrecrdst  >> res_ce_all/CE_all_FINAL_STATS.txt
fi

fi #END fasta_ce_flg
#and now for half IR

#First remove duplicates from IR events
#NEW CODE STARTS HERE
nrecrds_ir=0
nrecrds_irt=0
[ -e res_ce_all/ce_extension_coord.bed ] && nrecrds_ir=$(cat res_ce_all/ce_extension_coord.bed | wc -l) && nrecrds_irt=$(($nrecrds_ir/3))

if [ $nrecrds_ir -ne 0 ]
then
echo Total ce_extension events are: $nrecrds_irt, Now checking for repeated ce_extension events >> res_ce_all/CE_all_FINAL_STATS.txt
#echo checking for repeated IR events >> res_ce_all/CE_all_FINAL_STATS.txt
readarray -t all_data < res_ce_all/ce_extension_coord.bed
#also read ce_inclusion_coord_sashimi.bed to get TxID
readarray -t sashimi_data_ext < res_ce_all/ce_extension_coord_sashimi.bed

#also read csv file
readarray -t csv_data_ir < res_ce_all/IGV_R_returned_ce_extension.csv
csvi=0

#echo total IR events read are $nrecrds_irt
i=0
while [ $i -lt $nrecrds_ir ]
  do
          line11=${all_data[$i]}
          line12=${all_data[$(($i+1))]}
          line13=${all_data[$(($i+2))]}

					#for sashimi
					lines1=${sashimi_data_ext[$i]}
          lines2=${sashimi_data_ext[$(($i+1))]}
          lines3=${sashimi_data_ext[$(($i+2))]}

					#TxID
					txid_recrd=${sashimi_data_ext[$i]}
					txid=$(echo $txid_recrd | awk 'BEGIN {OFS="\t"} {print $8}')

#         echo i is $i, $line11, $line12, $line13
					csv_ln=${csv_data_ir[$csvi]}
					((csvi=csvi+1))

     #echo "Welcome $i times"
     ((i=i+3))
     j=$i
     #now go through rest of the data
     flg=0
        while [ $j -lt $nrecrds_ir ]
        do
                line21=${all_data[$j]}
          line22=${all_data[$(($j+1))]}
          line23=${all_data[$(($j+2))]}
          if [ "${line11[*]}" == "${line21[*]}" ] && [ "${line12[*]}" == "${line22[*]}" ] && [ "${line13[*]}" == "${line23[*]}" ]
          then
            #echo i is $i, $line11, $line12, $line13
            #echo j is $j, $line21, $line22, $line23
            flg=1
          fi
          ((j=j+3))
        done
        if [ $flg -eq 0 ]
				then
                 echo $line11 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_all/ce_extension_coord_uniq.bed
                 echo $line12 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_all/ce_extension_coord_uniq.bed
                 echo $line13 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_all/ce_extension_coord_uniq.bed
								 #also save bed file for with TxID for sashimi plots
								 #echo $line11 | awk -v TX=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,TX}' >> res_ce_all/ce_extension_coord_uniq_sashimi.bed
                 #echo $line12 | awk -v TX=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,TX}' >> res_ce_all/ce_extension_coord_uniq_sashimi.bed
                 #echo $line13 | awk -v TX=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,TX}' >> res_ce_all/ce_extension_coord_uniq_sashimi.bed
								 echo $lines1 | awk  'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' >> res_ce_all/ce_extension_coord_uniq_sashimi.bed
                 echo $lines2 | awk  'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' >> res_ce_all/ce_extension_coord_uniq_sashimi.bed
								 #Also save ce coordinates only
								 echo $lines2 | awk  'BEGIN {OFS="\t"} {print $1,$2,$3,$3-$2,$5,$6,$7,$8}' >> res_ce_all/ce_extension_coord_only.bed
                 echo $lines3 | awk  'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' >> res_ce_all/ce_extension_coord_uniq_sashimi.bed

								 #also save csv file
								 echo $csv_ln >> res_ce_all/IGV_unique_ce_extension.csv

				else
					echo $line11 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_all/ce_extension_coord_repeated.bed
					echo $line12 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_all/ce_extension_coord_repeated.bed
					echo $line13 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_all/ce_extension_coord_repeated.bed

         fi
 done


  #get total unique coordinates
 nrecrds_uniq_irt=$(cat res_ce_all/ce_extension_coord_uniq.bed | wc -l)
 nrecrds_uniq_ir=$(($nrecrds_uniq_irt/3))
 repeated_IR_events=$(($nrecrds_irt-$nrecrds_uniq_ir))

 echo Total unique ce_extension events are: $nrecrds_uniq_ir, please see res_ce_all/ce_extension_coord_uniq.bed and res_ce_all/IGV_unique_ce_extension.csv files >> res_ce_all/CE_all_FINAL_STATS.txt
 #echo total unique IR events read are $nrecrds_uniq_ir
 echo Total repeated ce_extension events are: $repeated_IR_events, please see res_ce_all/ce_extension_coord_repeated.bed >> res_ce_all/CE_all_FINAL_STATS.txt
 #echo THOSE WERE ALL PERTINENT STAT - Please let us know if something is missing or some more stats can be useful!!! >> res_ce_all/CE_all_FINAL_STATS.txt
 echo Total repeated ce_extension events were $repeated_IR_events
 #now copy back data to ce_coord
 #cat res_ce_all/ce_extension_coord_uniq.bed	> res_ce_all/ce_extension_coord.bed
#NEW CODE ENDS HERE - IR

fasta_ce_extension_flg=0
if [ $fasta_ce_extension_flg -eq 1 ]
then
echo now finding nt seq for res_ce_all/ce_extension_coord_uniq.bed
#now call bedtools getfasta function to get nt sequence from reference_genome.fa file

bedtools getfasta -fi /Volumes/SYEDSHAH/MichaelLab/ref_genome/GRCh38.p13.genome.fa -bed res_ce_all/ce_extension_coord_uniq.bed -s > ce_extension_nt.fasta

#now remove >chr lines from the resulting file
awk '!/^>chr/' ce_extension_nt.fasta > IR_all_nt1.fasta
#combine the two files to get desired csv file - chrXX,start,end,genename
paste -d"\t" res_ce_all/ce_extension_coord_uniq.bed IR_all_nt1.fasta > res_ce_all/ce_extension_nt.bed
#remove temp file
rm IR_all_nt1.fasta
rm ce_extension_nt.fasta
#and finally transeq compatible
#awk -F "\t" '{print ">sp|"$7"_"$1"_"$2"_"$3"\n"$8}' res_ce_all/ce_extension_nt.bed > res_ce_all/ce_extension_transeq_in.fasta
#Add strand infor
awk -F "\t" '{if($6=="+") print ">sp|"$7"_"$1"_"$2"_"$3"_plus""\n"$8; else print ">sp|"$7"_"$1"_"$2"_"$3"_minus""\n"$8}' res_ce_all/ce_extension_nt.bed > res_ce_all/ce_extension_transeq_in.fasta

#awk -F "\t" '{print ">sp|"$7"-"$1":"$2"-"$3"\n"$8}' res_ce_all/ce_inclusion_nt.bed > ce_inclusion_transeq_in.fasta

echo now doing the 3-frame translation

TEMPFILE=$(echo res_ce_all/ce_extension_transeq_in.fasta  | perl -pe 's/\.fasta$/.temp/')
OUTPUTFILE=$(echo res_ce_all/ce_extension_transeq_in.fasta  | perl -pe 's/\.fasta$/.trans/')

# Translate sequence
transeq  -sequence res_ce_all/ce_extension_transeq_in.fasta  -outseq "$TEMPFILE" -frame F

# Rename sequence
perl -sape 's/>/>sp|/' "$TEMPFILE" > "$OUTPUTFILE"



#also remove all newlines from the
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' "$OUTPUTFILE" > res_ce_all/CE_EXTENSION_AA.fasta

#also remove "$OUTPUTFILE" file
rm "$OUTPUTFILE"
# Remove temp file
rm "$TEMPFILE"


my_sym='>'
while read -r line
do
    a="${line[@]:0:1}"
    if [[ "$my_sym" == "$a" ]]; then
	line1="${line[@]}"
	#echo $line
    let i=0
    else
       b=$(echo "$line" | tr "*" '\n')
       while IFS= read li; #for ab in "${b[@]}"
	do
	size=${#li}
	if (( $size > 8 )); then
		let "i+=1"
		#let "j+=1"
	echo "$line1"_"$i"
	echo $li
	fi
	#echo $size
	done <<<"$b"


    fi
done < res_ce_all/CE_EXTENSION_AA.fasta > res_ce_all/FINAL_CE_EXTENSION_AA.fasta




#and for IR
#rm res_ce_all/ce_extension_fused.transeq_in.fasta
echo FINALLY DOING NT AND AA TRANSLATION FOR FUSED COORDINATES FOR CE_EXTENSION EVENTS
[ -e res_ce_all/ce_extension_fused.transeq_in.fasta ] && rm res_ce_all/ce_extension_fused.transeq_in.fasta
event_i=1
cat res_ce_all/ce_extension_nt.bed | while read r1; read r2; read r3
do
  lstrnd=$(echo "$r1" | awk '{print $6}')
  if [[ "$lstrnd" == "+" ]]
  then
    #echo "$r1" "$r2" "$r3"
    #title=$(echo "$r1" "$r2" "$r3" | awk '{print ">sp|"$23"_"$1"_"$2"_"$3"_"$10"_"$11"_"$18"_"$19"_plus"}')
		title=$(echo "$r1" "$r2" "$r3" | awk -v evid=$event_i '{print ">sp|"$23"_"$1"_"$2"_"$3"_"$10"_"$11"_"$18"_"$19"_"evid"_plus"}')
    #echo $title
    echo "$r1" "$r2" "$r3" | awk -v t=$title '{OFS="\t"} {print t"\n"$8$16$24}' >> res_ce_all/ce_extension_fused.transeq_in.fasta
  else
    #echo "$r1" "$r2" "$r3"
    #title=$(echo "$r1" "$r2" "$r3" | awk '{print ">sp|"$23"_"$1"_"$18"_"$19"_"$10"_"$11"_"$2"_"$3"_minus"}')
		title=$(echo "$r1" "$r2" "$r3" | awk -v evid=$event_i '{print ">sp|"$23"_"$1"_"$18"_"$19"_"$10"_"$11"_"$2"_"$3"_"evid"_minus"}')
    #echo $title
    echo "$r1" "$r2" "$r3" | awk -v t=$title '{OFS="\t"} {print t"\n"$8$16$24}' >> res_ce_all/ce_extension_fused.transeq_in.fasta
  fi
	((event_i=event_i+1))
done


#now do the 3-frame translation

TEMPFILE=$(echo res_ce_all/ce_extension_fused.transeq_in.fasta | perl -pe 's/\.fasta$/.temp/')
OUTPUTFILE=$(echo res_ce_all/ce_extension_fused.transeq_in.fasta | perl -pe 's/\.fasta$/.trans/')

# Translate sequence
transeq  -sequence res_ce_all/ce_extension_fused.transeq_in.fasta  -outseq "$TEMPFILE" -frame F

# Rename sequence
perl -sape 's/>/>sp|/' "$TEMPFILE" > "$OUTPUTFILE"



#also remove all newlines from the
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' "$OUTPUTFILE" > res_ce_all/CE_EXTENSION_FUSED_AA.fasta

#also remove "$OUTPUTFILE" file
rm "$OUTPUTFILE"
# Remove temp file
rm "$TEMPFILE"


my_sym='>'
while read -r line
do
    a="${line[@]:0:1}"
    if [[ "$my_sym" == "$a" ]]; then
	line1="${line[@]}"
	#echo $line
    let i=0
    else
       b=$(echo "$line" | tr "*" '\n')
       while IFS= read li; #for ab in "${b[@]}"
	do
	size=${#li}
	if (( $size > 8 )); then
		let "i+=1"
		#let "j+=1"
	echo "$line1"_"$i"
	echo $li
	fi
	#echo $size
	done <<<"$b"
    fi
done < res_ce_all/CE_EXTENSION_FUSED_AA.fasta > res_ce_all/FINAL_CE_EXTENSION_FUSED_AA.fasta

else
	echo Total ce_extension events are: $nrecrds_irt >> res_ce_all/CE_all_FINAL_STATS.txt
fi

fi #ENDS fasta_ce_extension_flg
#Step 5 ENDS
#NOW FOR IR EVENTS
#First remove duplicates from IR events
#NEW CODE STARTS HERE
nrecrds_ir=0
nrecrds_irt=0
[ -e res_ce_all/IR_coord.bed ] && nrecrds_ir=$(cat res_ce_all/IR_coord.bed | wc -l) #&& nrecrds_irt=$(($nrecrds_ir/3))

if [ $nrecrds_ir -ne 0 ]
then
echo TOTAL IR EVENTS ARE: $nrecrds_ir, NOW CHECKING FOR REPEATED IR EVENTS >> res_ce_all/CE_all_FINAL_STATS.txt
echo TOTAL IR EVENTS ARE: $nrecrds_ir
#echo checking for repeated IR events >> res_ce_all/CE_all_FINAL_STATS.txt
readarray -t all_data < res_ce_all/IR_coord.bed
#also read ce_inclusion_coord_sashimi.bed to get TxID
readarray -t sashimi_data_ext < res_ce_all/IR_coord_sashimi.bed

#also read csv file
readarray -t csv_data_ir < res_ce_all/IGV_R_returned_IR.csv


csvi=0

#echo total IR events read are $nrecrds_irt
i=0
while [ $i -lt $nrecrds_ir ]
  do
          line11=${all_data[$i]}
					#for sashimi
					lines1=${sashimi_data_ext[$i]}
					#TxID
					txid_recrd=${sashimi_data_ext[$i]}
					txid=$(echo $txid_recrd | awk 'BEGIN {OFS="\t"} {print $8}')

#         echo i is $i, $line11, $line12, $line13
					csv_ln=${csv_data_ir[$csvi]}
					((csvi=csvi+1))

     #echo "Welcome $i times"
     ((i=i+1))
     j=$i
     #now go through rest of the data
     flg=0
        while [ $j -lt $nrecrds_ir ]
        do
                line21=${all_data[$j]}
          if [ "${line11[*]}" == "${line21[*]}" ] #&& [ "${line12[*]}" == "${line22[*]}" ] && [ "${line13[*]}" == "${line23[*]}" ]
          then
            flg=1
          fi
          ((j=j+1))
        done
        if [ $flg -eq 0 ]
				then
                 echo $line11 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_all/IR_coord_uniq.bed
								 echo $lines1 | awk  'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' >> res_ce_all/IR_coord_uniq_sashimi.bed
                 #echo $lines2 | awk  'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' >> res_ce_all/ce_extension_coord_uniq_sashimi.bed
								 #Also save ce coordinates only
								 echo $lines1 | awk  'BEGIN {OFS="\t"} {print $1,$2,$3,$3-$2,$5,$6,$7,$8}' >> res_ce_all/IR_coord_only.bed
                 #echo $lines3 | awk  'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' >> res_ce_all/ce_extension_coord_uniq_sashimi.bed

								 #also save csv file
								 #echo $csv_ln >> res_ce_all/IGV_unique_IR.csv
								 echo $csv_ln | awk 'BEGIN {FS=","} {print $2,$3,$4,$7,$8,$9}'| awk 'BEGIN {OFS=","} {print $1,$2,$3,$4,$5,$6}' >> res_ce_all/IGV_unique_IR.csv

				else
					echo $line11 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_all/IR_coord_repeated.bed
         fi
 done
  #get total unique coordinates
 nrecrds_uniq_irt=$(cat res_ce_all/IR_coord_uniq.bed | wc -l)
 #nrecrds_uniq_ir=$(($nrecrds_uniq_irt/3))
 repeated_IR_events=$(($nrecrds_ir-$nrecrds_uniq_irt))

 echo TOTAL UNIQUE IR EVENTS ARE: $nrecrds_uniq_irt, please see res_ce_all/IR_coord_uniq.bed and res_ce_all/IGV_unique_IR.csv files >> res_ce_all/CE_all_FINAL_STATS.txt
 #echo total unique IR events read are $nrecrds_uniq_ir
 echo Total repeated IR events are: $repeated_IR_events, please see res_ce_all/IR_coord_repeated.bed >> res_ce_all/CE_all_FINAL_STATS.txt
 echo THOSE WERE ALL PERTINENT STAT - Please let us know if something is missing or some more stats can be useful!!! >> res_ce_all/CE_all_FINAL_STATS.txt
 echo Total repeated IR events were $repeated_IR_events
 #now copy back data to ce_coord
 #cat res_ce_all/ce_extension_coord_uniq.bed	> res_ce_all/ce_extension_coord.bed
#NEW CODE ENDS HERE - IR


fasta_ir_flg=0
if [ $fasta_ir_flg -eq 1 ]
then
echo now finding nt seq for res_ce_all/IR_coord_uniq.bed


#now call bedtools getfasta function to get nt sequence from reference_genome.fa file
bedtools getfasta -fi /Volumes/SYEDSHAH/MichaelLab/ref_genome/GRCh38.p13.genome.fa -bed res_ce_all/IR_coord_uniq.bed -s > res_ce_all/IR_coord_uniq_nt.fa

#now remove >chr lines from the resulting file
awk '!/^>chr/' res_ce_all/IR_coord_uniq_nt.fa > res_ce_all/IR_coord_uniq_nt1.fa


paste -d"\t" res_ce_all/IR_coord_uniq.bed res_ce_all/IR_coord_uniq_nt1.fa > res_ce_all/temp.IR_coord_uniq_nt.bed

rm res_ce_all/IR_coord_uniq_nt1.fa
#rm res_ce_all/IR_coord_uniq_nt.fa
#awk -F "\t" '{if($6=="+") print ">sp|"$7"_"$1"_"$2"-"$3"_plus""\n"$8;else print ">sp|"$7"_"$1"_"$2"-"$3"_minus""\n"$8}' res_ce_all/temp.IR_coord_uniq_nt.bed > res_ce_all/temp.IR_coord_uniq_nt.transeq_in.fasta
i=0
awk -F "\t" '{if($6=="+") print ">sp|"$7"_"$1"_"$2"_"$3"_"((i=i+1))"_plus""\n"$8;else print ">sp|"$7"_"$1"_"$2"_"$3"_"((i=i+1))"_minus""\n"$8}' res_ce_all/temp.IR_coord_uniq_nt.bed > res_ce_all/temp.IR_coord_uniq_nt.transeq_in.fasta
rm res_ce_all/temp.IR_coord_uniq_nt.bed
#also save nt sequence

TEMPFILE=$(echo res_ce_all/temp.IR_coord_uniq_nt.transeq_in.fasta | perl -pe 's/\.fasta$/.temp/')
OUTPUTFILE=$(echo res_ce_all/temp.IR_coord_uniq_nt.transeq_in.fasta | perl -pe 's/\.fasta$/.trans/')

# Translate sequence

transeq  -sequence res_ce_all/temp.IR_coord_uniq_nt.transeq_in.fasta -outseq "$TEMPFILE" -frame F


# Rename sequence
perl -sape 's/>/>sp|/' "$TEMPFILE" > "$OUTPUTFILE"


#no concatenate all files
cat "$OUTPUTFILE" >> res_ce_all/temp.IR_coord_uniq_AA.fasta
# Remove temp file
rm "$TEMPFILE"
rm "$OUTPUTFILE"

rm res_ce_all/temp.IR_coord_uniq_nt.transeq_in.fasta

#also remove all newlines from the
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' res_ce_all/temp.IR_coord_uniq_AA.fasta > res_ce_all/IR_coord_uniq_AA.fasta

#remove
rm res_ce_all/temp.IR_coord_uniq_AA.fasta

#finally break each AA sequence at stop codons (*) and remove AA sequence of length <9
my_sym='>'
while read -r line
do
    a="${line[@]:0:1}"
    if [[ "$my_sym" == "$a" ]]; then
	line1="${line[@]}"
	#echo $line
    let i=0
    else
       b=$(echo "$line" | tr "*" '\n')
       while IFS= read li; #for ab in "${b[@]}"
	do
	size=${#li}
	if (( $size > 8 )); then
		let "i+=1"
		#let "j+=1"
	echo "$line1"_"$i" #for + starnd
  #echo "$line1"_"$i"
	echo $li
	fi
	#echo $size
	done <<<"$b"


    fi
done < res_ce_all/IR_coord_uniq_AA.fasta > res_ce_all/FINAL_IR_AA.fasta


else
	echo TOTAL IR EVENTS ARE: $nrecrds_irt >> res_ce_all/CE_all_FINAL_STATS.txt
fi

fi #ENDS fasta_ir_flg
#SASHIMI PLOTS SECTION
if [ $1 -eq 7 ] || [ $1 -eq 8 ] #generate sashimi plots for each category
then

  if [ -f res_ce_all/ce_inclusion_coord_uniq_sashimi.bed ]
  then
    num_ce_incl=$(cat res_ce_all/ce_inclusion_coord_uniq_sashimi.bed | wc -l | awk '{print $1/3}')
    echo STARTING SASHIMI PLOTS FOR $num_ce_incl Events, For events > 100 It Might take several hours
    source run_sashimi.sh 3
  else
    echo File res_ce_all/ce_inclusion_coord_uniq_sashimi.bed does not exists so, NO CE_INCLUSION SASHIMI PLOTS ARE GENERATED
  fi

  if [ -f res_ce_all/ce_extension_coord_uniq_sashimi.bed ]
  then
    num_ce_ext=$(cat res_ce_all/ce_extension_coord_uniq_sashimi.bed | wc -l | awk '{print $1/3}')
    echo STARTING SASHIMI PLOTS FOR $num_ce_ext Events, For events > 100 It Might take several hours Might take several hours
    source run_sashimi.sh 4
  else
    echo File res_ce_all/ce_extension_coord_uniq_sashimi.bed does not exists so, NO CE_EXTENSION SASHIMI PLOTS ARE GENERATED

  fi

  if [ -f res_ce_all/IR_coord_uniq_sashimi.bed ]
  then
    num_IR=$(cat res_ce_all/IR_coord_uniq_sashimi.bed | wc -l)
    echo STARTING SASHIMI PLOTS FOR $num_IR Events, For events > 100 It Might take several hours Might take several hours
    source run_sashimi.sh 5
  else
    echo File res_ce_all/IR_coord_uniq_sashimi.bed does not exists so, NO CE_EXTENSION SASHIMI PLOTS ARE GENERATED

  fi
fi

else #ENDS Auto_CoverV4_layered_intronV2 and CE calculations
	echo CE_BOUNDARY_EVENTS ARE $ce_boundary_events SO ABONDONED CE_BOUNDARY CALCULATIONS
fi

######END IR EVENTS
#FINAL STATISTICS
#fi #END for 5 || 6
fi #END for $1 = 6 || 7
echo For final statistics on ce events, please see res_ce_all/CE_all_FINAL_STATS.txt file



echo ALL DONE - hopefully - Successfully

#remove all bed files now

#rm *.bed
