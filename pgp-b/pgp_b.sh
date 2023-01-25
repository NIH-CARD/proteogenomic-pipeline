#!/bin/bash

#HOW TO RUN: bash pgp_b_ce_inclusion.sh ce_extension_pgp1.csv COV_CUTOFF NUM_EVENTS_PROCESSED
#EXAMPLE: bash pgp_b_ce_inclusion.sh ce_extension_pgp1.csv .6 1



#PLEASE MAKE SURE TO CONCATENATE
#1. CE_INCLUSION_FUSED_AA.fasta FOR EACH CUTOFF (WITH 40% FOLLOWED BY NEXT CUTOFF and so on untill all hand_curated events are accounted for)
#2. ce_inclusion_fused.transeq_in.fasta
#3. FINAL_CE_INCLUSION_FUSED_AA.fasta
#4. IGV_unique_ce_inclusion.csv
###########################

#CHECK IF 3 ARGUMENTS ARE PROVIDED
if [ $# -eq 1 ]
then
	[ -e res_skiptics/Summary.Skiptics.txt ] && rm res_skiptics/Summary.Skiptics.txt
	mkdir -p res_skiptics
	if [ -f "$1" ]
	then
		splicing_events_file=$(echo "$1"|cut -d. -f1)
		sort -t$',' -k5 $1 > sorted_"$splicing_events_file".csv
	else
		#echo line 402 NO SPLICING_EVENTS csv file is provided, ABORTING!!!!
		echo PLEASE PROVIDE SPLICING_EVENTS csv file  "'"bash pgp_b.sh splicing_events.csv"'" and RERUN
		echo PLEASE PROVIDE SPLICING_EVENTS csv file  "'"bash pgp_b.sh splicing_events.csv"'" and RERUN >> res_skiptics/Summary.Skiptics.txt
		exit 1
	fi
	echo "################################ STARTED SKIPTICS IDENTIFICATION "
	echo "################################ STARTED SKIPTICS IDENTIFICATION " >> res_skiptics/Summary.Skiptics.txt
	[ -e EnsDB_tx_not_found.csv ] && rm EnsDB_tx_not_found.csv
	skiptics_flg=1
	#skiptics_sashimi_flg=0
	if [ $skiptics_flg -eq 1 ]
	then

			echo '#############################'
	    echo INVOKING SKIPTICS Script esV5_layered_CDSV3.sh FOR FILE sorted_"$splicing_events_file".csv from pgp_b.sh
			echo INVOKING SKIPTICS Script esV5_layered_CDSV3.sh FOR FILE sorted_"$splicing_events_file".csv from pgp_b.sh >> res_skiptics/Summary.Skiptics.txt
			echo #############################
			#source esV5_layered.sh sorted_$1 principal_tx2_merged_samples.csv #$2 #clean_"$splicing_events_file".csv #$1 #$2
	    source esV5_layered_CDSV3.sh sorted_"$splicing_events_file".csv principal_txs.csv #principal_tx2_merged_samples.csv #$2 #clean_majiq_events.csv #$1 #$2

			echo '#############################'
			echo BACK FROM SKIPTICS Script
			echo '#############################'
			echo "################################ DONE SKIPTICS IDENTIFICATION ", PLEASE SEE res_skiptics/Summary.Skiptics.txt for skiptics details >> res_skiptics/Summary.Skiptics.txt
			#CALL ggsashimi for SASHIMI PLOTS
			#if [ $1 -eq 4 ] || [ $1 -eq 8 ] || [ $1 -eq 3 ] #8 is for start to finish
			#then
				echo "################################ STARTED CREATING SASHIMI PLOTS FOR SKIPTICS EVENTS - MAY TAKE MANY HOURS " >> res_skiptics/Summary.Skiptics.txt
				echo "################################ STARTED CREATING SASHIMI PLOTS FOR SKIPTICS EVENTS - MAY TAKE MANY HOURS "
				#source run_sashimi.sh 1
	      source run_sashimiV1.sh res_skiptics/IGV_unique_skiptics_translated.csv res_skiptics/skiptics_uniq_sashimi.bed 1
				echo "################################ DONE - RESULTANT pdf is sashimi_plots/all_skiptics_sashimi.pdf " >> res_skiptics/Summary.Skiptics.txt
			#fi

	  	echo DONE WITH SKIPTICS CALCULATIONS, PLEASE CHECK res_skiptics/Summary.Skiptics.TXT FILE FOR ALL STATISTICS AND:
		echo IGV_skiptics_uniq.csv: contains all unique exon skip events for which nt and AA seq re calculated
		echo cds_skiptics_fused_transeq_in.fasta: contains fastaID and nt seq for each exon skip event in IGV_skiptics_uniq.csv "("except those with non-csd")"
		echo AllFrames_PEAKS_SKIPTICS_FUSED_AA.fasta: contains fastaID and 3-frame translated AA seq for each exon skip event in IGV_skiptics_uniq.csv "("except those with non-csd")"
		echo all_skiptics_sashimi.pdf: contains sashimi plot for all exon skip events in IGV_skiptics_uniq.csv
		echo There are also several other files worth looking
		echo DONE WITH SKIPTICS CALCULATIONS, PLEASE CHECK res_skiptics/Summary.Skiptics.txt FILE FOR ALL STATISTICS AND: >> res_skiptics/Summary.Skiptics.txt
		echo IGV_skiptics_uniq.csv: contains all unique exon skip events for which nt and AA seq re calculated >> res_skiptics/Summary.Skiptics.txt
		echo cds_skiptics_fused_transeq_in.fasta: contains fastaID and nt seq for each exon skip event in IGV_skiptics_uniq.csv "("except those with non-csd")" >> res_skiptics/Summary.Skiptics.txt
		echo AllFrames_PEAKS_SKIPTICS_FUSED_AA.fasta: contains fastaID and 3-frame translated AA seq for each exon skip event in IGV_skiptics_uniq.csv "("except those with non-csd")" >> res_skiptics/Summary.Skiptics.txt
		echo all_skiptics_sashimi.pdf: contains sashimi plot for all exon skip events in IGV_skiptics_uniq.csv >> res_skiptics/Summary.Skiptics.txt
		echo THERE ARE SEVERAL other files worth looking "("PLEAE UNCOMMENT BELOW TO GET TEHSE FILES GENERATED")" >> res_skiptics/Summary.Skiptics.txt

		#also delete all unnecessary files
		del_flg=1
		if [ $del_flg -eq 1 ]
		then
		mv res_skiptics/IGV_unique_skiptics_translated.csv res_skiptics/Exon_skip.csv
		rm res_skiptics/IGV_unique_skiptics_translated.csv
		rm res_skiptics/aa.fasta
		rm res_skiptics/AllFrames_PEAKS_SKIPTICS_FUSED_AA.fasta
		rm res_skiptics/AllFrames_SKIPTICS_FUSED_AA1.fasta
		rm res_skiptics/cds_IGV_unique_skiptics_translated.csv
		rm res_skiptics/cds_SKIPTICS_FUSED_AA.fasta
		rm res_skiptics/cds_SKIPTICS_FUSED_AA1.fasta

		rm res_skiptics/cds_skiptics_uniq_sashimi.bed
		rm res_skiptics/cds_skiptics_unique.bed
		rm res_skiptics/cds_successful_frames_list.csv
		rm res_skiptics/IGV_all_others.csv
		rm res_skiptics/IGV_skiptics_repeated.csv
		rm res_skiptics/IGV_skiptics_uniq.csv
		rm res_skiptics/IGV_skiptics.csv
		rm res_skiptics/IGV_unique_skiptics_translated.csv
		rm res_skiptics/non_protein_coding.bed
		rm res_skiptics/protein_coding_a.bed
		rm res_skiptics/protein_coding.bed
		rm res_skiptics/skiptics_all_coordinates.bed
		rm res_skiptics/skiptics_repeated.bed 
		rm res_skiptics/skiptics_sashimi.bed
		rm res_skiptics/skiptics_uniq_sashimi.bed
		rm res_skiptics/skiptics_unique.bed
		rm res_skiptics/skiptics_us_exons.txt
		rm res_skiptics/all_frames_final_aa.fasta
		rm res_skiptics/AllFrames_SKIPTICS_FUSED_AA2.fasta
		rm res_skiptics/cds_SKIPTICS_FUSED_AA2.fasta
		rm res_skiptics/final_aa.fasta

		fi

		exit 1
	else
		echo skiptics_flg is set to 0, Please set it to 1 and then rerun
		echo skiptics_flg is set to 0, Please set it to 1 and then rerun >> res_skiptics/Summary.Skiptics.txt
		exit 1
	fi
elif [ $# -ne 3 ]
then
	echo THIS SCRIPT REQUIRES EVENTS csv FILE, COVERAGE CUTOFF BETWEEN "> " 0 and "< " 1 and "#" OF CE_INCLUSION EVENTS PROCESSED SO FAR, 1 IF NONE SO FAR
	echo PLEAE INPUT REQUIRED PARAMETERS, EXITING NOW !!!!!!!!!!
	exit 1
else

	cryptics_flg=1
	if [ $cryptics_flg -eq 1 ]
	then

		mkdir -p res_ce_"$2"
		mkdir -p event_bedfiles

		[ -e events_to_tx_mapping_valid.csv ] && rm events_to_tx_mapping_valid.csv
		[ -e events_tx_mapping_invalid.csv ] && rm events_tx_mapping_invalid.csv

		[ "$(ls -A res_ce_"$2"/)" ] && rm res_ce_"$2"/*.*


		echo NOW STARTING CE CALCULATIONS
		echo NOW STARTING CE CALCULATIONS >> res_ce_"$2"/Summary_stats.txt
		COUNT_EVENT=$3 #THESE ARE UNIQUE CE_INCLUSION EVENTS DETECTED + 1 FOR 40% CUTOFF.
		#NEXT SORT FILE ON COLUMN 5 - THIS IS IMPORTANT AS PIPELINE DEPENDS ON IT
		echo got input file $1, Now sorting it
		echo got input file $1, Now sorting it >> res_ce_"$2"/Summary_stats.txt

		sort -t$',' -k5 $1 > sorted_"$1"


		inpfile=sorted_"$1" #all_non_skiptics.csv

		events_bed_create_flg=1
		if [ $events_bed_create_flg -eq 1 ]
		then
			[ "$(ls -A event_bedfiles/)" ] && rm event_bedfiles/*.*
			[ -e res_ce_"$2"/EnsDB_tx_not_found.csv ] && rm res_ce_"$2"/EnsDB_tx_not_found.csv
			echo Now Calling TxEnsDB103_layeredV6.R >> res_ce_"$2"/Summary_stats.txt
			Rscript TxEnsDB103_layeredV6.R "$inpfile" principal_txs.csv res_ce_"$2"/Summary_stats.txt #$tx_lst
			echo Back from TxEnsDB103_layeredV6.R >> res_ce_"$2"/Summary_stats.txt
			#Rscript TxEnsDB103_layeredV6.R sorted_ce_inclusion_pgp1.csv principal_txs.csv res_ce_"$2"/Summary_stats.txt #$tx_lst
			
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
			allexons=$(echo "$sample" | cut -d'/' -f2 |cut -d'_' -f2)
			gene_name=$(echo "$allexons" |rev| cut -d'.' -f2-|cut -d'-' -f2-|rev)

			sortBed -i event_bedfiles/$allexons > event_bedfiles/t$allexons
			#Also read Tx Files to retrieve selected Tx - should find better ways
			TxID=$(head -1 event_bedfiles/TxID$allexons | awk '{print $7}')

			#echo got TxID $TxID

			strnd=$(cat "$sample" | awk '{print $6}')
			ds=$(bedtools closest -a $sample -b event_bedfiles/t$allexons -s -D a -iu -d -t first )
			us=$(bedtools closest -a $sample -b event_bedfiles/t$allexons -s  -D a -id -d -t first)

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
			then
				echo $gene_name is a CE event >> res_ce_"$2"/ex1_ex2_ce.txt  #>> res_intronic_range/intronic_range_events.txt
				echo ds $ds >> res_ce_"$2"/ex1_ex2_ce.txt
				echo us $us >> res_ce_"$2"/ex1_ex2_ce.txt

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
				echo now dsn $dsn >> res_ce_"$2"/ex1_ex2_ce.txt
				echo now usn $usn >> res_ce_"$2"/ex1_ex2_ce.txt
				echo dsodiff $dsodiff >> res_ce_"$2"/ex1_ex2_ce.txt
				echo up exon length $usnexonl and dsnexonl $dsnexonl

				if [ $usnexonl == "." ] || [ $dsnexonl == "." ]
				then
					#echo $gene_name is ce with up exon length $usnexonl and dsnexonl $dsnexonl , LEAVING it for now PLEASE MAKE SURE IT IS PROPERLY IDENTIFIED >> res_intronic_range/problematic_junctions.txt
					echo $gene_name is ce with up exon length $usnexonl and dsnexonl $dsnexonl , LEAVING it for now PLEASE MAKE SURE IT IS PROPERLY IDENTIFIED >> res_ce_"$2"/problematic_junctions.txt
					#cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_intronic_range/IGV_problematic_junctions.csv
					cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_ce_"$2"/IGV_problematic_junctions.csv

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

					gene_id=$(echo "$line_csv"| awk 'BEGIN {FS=","} {print $6}')
					o_gene_id=$(echo "$line_csv"| awk 'BEGIN {FS=","} {print $5}')
					cat $sample | awk -v gene=$gene_name -v TX=$TxID -v og=$o_gene_id -v gid=$gene_id 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene,TX,og,gid}' >> res_ce_"$2"/IGV_ce_inclusion.csv
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

					fi

					#also coe ce boundary range
					echo "$upex">>res_ce_"$2"/ce_all_scan_range_junctions.bed
					echo "$ce_scan1">>res_ce_"$2"/ce_all_scan_range_junctions.bed
					echo "$dsex">>res_ce_"$2"/ce_all_scan_range_junctions.bed

					echo "$ce_scan1">>res_ce_"$2"/ce_all_scan_range.bed #all_junctions.bed


				fi

			elif [ "$diff_exon_abs" -eq 1 ] && [ $dsovlp -le 5 ] && [ $usovlp -le 5 ] #deals with exon_joining events few bp of
			then
				cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_ce_"$2"/IGV_skiptics.csv
				#cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_intronic_range/non_ce_events.csv

				echo $line_csv | awk 'BEGIN {FS=","} {print $1,$2,$3,$4,$5,$6,$7}' | awk 'BEGIN {OFS=","} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_"$2"/non_ce_events.csv

			elif [ "$diff_exon_abs" -eq 1 ]
			then
				#THESE ARE CE INTRONIC events
				#CE event with difference of one exon

				#echo $gene_name is a CE event #>> res_intronic_range/intronic_range_events.txt #>>ce_intronic.txt

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

				gene_id=$(echo "$line_csv"| awk 'BEGIN {FS=","} {print $6}')
				o_gene_id=$(echo "$line_csv"| awk 'BEGIN {FS=","} {print $5}')

				cat $sample | awk -v gene=$gene_name -v TX=$TxID -v og=$o_gene_id -v gid=$gene_id 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene,TX,og,gid}' >> res_ce_"$2"/IGV_ce_inclusion.csv


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

				fi

				#also coe ce boundary range
				echo "$upex">>res_ce_"$2"/ce_all_scan_range_junctions.bed
				echo "$ce_scan1">>res_ce_"$2"/ce_all_scan_range_junctions.bed
				echo "$dsex">>res_ce_"$2"/ce_all_scan_range_junctions.bed

				echo "$ce_scan1">>res_ce_"$2"/ce_all_scan_range.bed #all_junctions.bed
			else
				#BOTH ENDS OF THESE EVENTS LIE INSIDE INTRON
				echo $gene_name has both ends lying inside intron #>> res_intronic_range/intronic_range_events.txt #>>ce_intronic.txt
				#echo $gene_name is unknown event type for gene "$gene_name" >> res_intronic_range/non_ce_events.txt
				echo $gene_name is unknown event type for gene "$gene_name" >> res_ce_"$2"/non_ce_events.txt
				#cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_intronic_range/non_ce_events.csv
				echo $line_csv | awk 'BEGIN {FS=","} {print $1,$2,$3,$4,$5,$6,$7}' | awk 'BEGIN {OFS=","} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_"$2"/non_ce_events.csv

			fi

			else
				#echo $gene_name is unknown event type for gene "$gene_name" >> res_intronic_range/non_ce_events.txt
				echo $gene_name is unknown event type for gene "$gene_name" >> res_ce_"$2"/non_ce_events.txt
				#cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_intronic_range/non_ce_events.csv
				echo $line_csv | awk 'BEGIN {FS=","} {print $1,$2,$3,$4,$5,$6,$7}' | awk 'BEGIN {OFS=","} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_"$2"/non_ce_events.csv
			fi

			#rm $sample1
			#rm $sample_nt.bed
			#rm $sample1-3UTR.bed
			#rm $sample1-5UTR.bed
			#rm $allexons
			#rm t$allexons
		done
		########## THIS SECTION WRITES NEW FILE res_ce_"$2"/ce_all_scan_intron.bed to scan whole intron between two exons for probable ce events
		nrecrds=$(cat res_ce_"$2"/ce_all_scan_range_junctions.bed | wc -l)
		readarray -t all_data < res_ce_"$2"/ce_all_scan_range_junctions.bed
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

					echo "$line11" "$line13" | awk -v o=$off 'BEGIN {OFS="\t"} {print $1,$3,$9,$11,$12,$13,$14}' >> res_ce_"$2"/ce_all_scan_intron.bed
				elif [ "$strand" == "-" ]
				then
					#st=$(echo $line13 | awk 'BEGIN {FS="\t"} {print $2}')
					#end=$(echo $line13 | awk 'BEGIN {FS="\t"} {print $3}')
					#off=$(($end-$st))
					#st1=$(echo $line11 | awk 'BEGIN {FS="\t"} {print $2}')
					#end1=$(echo $line11 | awk 'BEGIN {FS="\t"} {print $3}')
					#off1=$(($end1-$st1))

					echo "$line13" "$line11" | awk 'BEGIN {OFS="\t"} {print $1,$3,$9,$11,$12,$13,$14}' >> res_ce_"$2"/ce_all_scan_intron.bed
				fi
		done

		################## END SECTION WRITES NEW FILE res_ce_"$2"/ce_all_scan_intron.bed
		#also copy for coverage calculations and avoid intronic_range calculations
		cat res_ce_"$2"/ce_all_scan_range.bed  > res_ce_"$2"/ce_all_scan_unique_range.bed
	fi

	######### Decide ce_extension or ce_inclusion

	if [ "$2" = "ext" ];
	then
		mv res_ce_"$2"/Summary_stats.txt res_ce_"$2"/Summary_stats.txt
		[ -e res_ce_"$2"/ce_extension_coord.bed ] && rm res_ce_"$2"/ce_extension_coord.bed
		nrecrds=0
		#nrecrdst=0
		[ -e res_ce_"$2"/ce_all_scan_range_junctions.bed ] && nrecrds=$(cat res_ce_"$2"/ce_all_scan_range_junctions.bed | wc -l) && nrecrdst=$(($nrecrds/3))
		echo GOT $nrecrdst CE_EXTENSION events
		echo GOT $nrecrdst CE_EXTENSION events >> res_ce_"$2"/Summary_stats.txt
		if [ $nrecrds -ne 0 ]
		then
			#GET up/dn exon coordinates
			readarray -t all_data < res_ce_"$2"/ce_all_scan_range_junctions.bed
			i=0
			while [ $i -lt $nrecrds ]
			do
				line11=${all_data[$i]}
				line12=${all_data[$(($i+1))]}
				line13=${all_data[$(($i+2))]}
				#echo "Welcome $i times"
				((i=i+3))
				#strand=$(echo $line11 | awk '{print $6}')
				echo $line11 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_"$2"/ce_extension_coord.bed

				scn1=$(echo $line12 | awk '{print $2}')
				scn2=$(echo $line12 | awk '{print $3}')
				strnd=$(echo $line12 | awk '{print $6}')
				if [ $scn1 -lt $scn2 ]
				then
					#if [ "$strnd" == "+" ] #As this is start of range, so we subtract 1 (as it is the start of the feature) for proper nt selection
					#then
					#AS THIS IS START OF RANGE, SO WE SUBTRACT 1 (AS IT IS THE START OF THE FEATURE) FOR PROPER NT SELECTION
					echo $line12 | awk 'BEGIN {OFS="\t"} {print $1,$2-1,$3-1,$4,$5,$6,$7}' >> res_ce_"$2"/ce_extension_coord.bed

				else
					#echo $line12 | awk 'BEGIN {OFS="\t"} {print $1,$3,$2-1,$4,$5,$6,$7}' >> res_ce_"$2"/ce_extension_coord.bed
					echo $line12 | awk 'BEGIN {OFS="\t"} {print $1,$3,$2,$4,$5,$6,$7}' >> res_ce_"$2"/ce_extension_coord.bed
				fi
				echo $line13 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_"$2"/ce_extension_coord.bed
			done
			nrecrds_ir=0
			nrecrds_irt=0
			[ -e res_ce_"$2"/ce_extension_coord.bed ] && nrecrds_ir=$(cat res_ce_"$2"/ce_extension_coord.bed | wc -l) && nrecrds_irt=$(($nrecrds_ir/3))
			if [ $nrecrds_ir -ne 0 ]
			then

				echo Now checking for repeated ce_extension events >> res_ce_"$2"/Summary_stats.txt
				#echo checking for repeated IR events >> res_ce_"$2"/Summary_stats.txt
				readarray -t all_data < res_ce_"$2"/ce_extension_coord.bed
				#also read IGV_ce_inclusion.csv to get TxID - It is safe as we consider all events as ce_extension here
				readarray -t get_txid < res_ce_"$2"/IGV_ce_inclusion.csv

				#also read csv file
				readarray -t csv_data_ir < sorted_$1
				csvi=0

				#echo total IR events read are $nrecrds_irt
				i=0
				while [ $i -lt $nrecrds_ir ]
				do
					line11=${all_data[$i]}
					line12=${all_data[$(($i+1))]}
					line13=${all_data[$(($i+2))]}

					csv_ln=${csv_data_ir[$csvi]}
					txid_rec=${get_txid[$csvi]}
					txid=$(echo $txid_rec | awk 'BEGIN {FS=","} {print $9}')

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
						#MAKE COLUMN 4 as 1 to compare when removing non-cds lines
						echo $line11 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,1,$5,$6,$7}' >> res_ce_"$2"/ce_extension_coord_uniq.bed
						echo $line12 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,1,$5,$6,$7}' >> res_ce_"$2"/ce_extension_coord_uniq.bed
						echo $line13 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,1,$5,$6,$7}' >> res_ce_"$2"/ce_extension_coord_uniq.bed
						#Also save ce coordinates only
						echo $line12 | awk  -v tx=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$3-$2,$5,$6,$7,tx}' >> res_ce_"$2"/ce_extension_coord_only.txt
						###### Also save ce_extension_coord_uniq_sashimi.bed for sashimi plots
						# Also save for sashimi
						echo $line11 | awk -v tx=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,tx}' >> res_ce_"$2"/ce_extension_coord_uniq_sashimi.bed
						echo $line12 | awk -v tx=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,1,$5,$6,$7,tx}' >> res_ce_"$2"/ce_extension_coord_uniq_sashimi.bed
						echo $line13 | awk -v tx=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,tx}' >> res_ce_"$2"/ce_extension_coord_uniq_sashimi.bed
						#####

						######## End saving ce_extension_coord_uniq_sashimi.bed

						#also save csv file
						echo $csv_ln | awk 'BEGIN {FS=","} {print $1,$2,$3,$4,$5,$6}' | awk -v tx=$txid '{OFS=","} {print $1":"$2"-"$3,$1,$2,$3,1,0,$4,$5,tx,$5,$6}' >> res_ce_"$2"/IGV_unique_ce_extension.csv


					else
						echo $csv_ln awk 'BEGIN {FS=","} {print $1,$2,$3,$4,$5,$6}' | awk '{OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6}' >> res_ce_"$2"/ce_extension_coord_repeated.csv
					fi
				done
			fi
			#get total unique coordinates
			nrecrds_uniq_irt=$(cat res_ce_"$2"/ce_extension_coord_uniq.bed | wc -l)
			nrecrds_uniq_ir=$(($nrecrds_uniq_irt/3))
			repeated_IR_events=$(($nrecrds_irt-$nrecrds_uniq_ir))

			echo Total unique ce_extension events are: $nrecrds_uniq_ir, please see res_ce_"$2"/ce_extension_coord_uniq.bed and res_ce_"$2"/IGV_unique_ce_extension.csv files >> res_ce_"$2"/Summary_stats.txt
			#echo total unique IR events read are $nrecrds_uniq_ir
			if [ $repeated_IR_events -gt 0 ]
			then
				echo Total repeated ce_extension events are: $repeated_IR_events, please see res_ce_"$2"/ce_extension_coord_repeated.csv >> res_ce_"$2"/Summary_stats.txt
			fi
			#echo THOSE WERE ALL PERTINENT STAT - Please let us know if something is missing or some more stats can be useful!!! >> res_ce_"$2"/Summary_stats.txt
			echo Total repeated ce_extension events were $repeated_IR_events
			#now copy back data to ce_coord
			#cat res_ce_"$2"/ce_extension_coord_uniq.bed	> res_ce_"$2"/ce_extension_coord.bed

			##########CHECK FOR orf
			orf_flg=1
			if [ $orf_flg -eq 1 ]
			then

				echo NOW CALLING get_orf_cds.R TO IDENTIFY PROTEIN CODING GENES FROM THE LIST OF IDENTIFIED CRYPTICS AND CDS PHASE
				echo NOW CALLING get_orf_cds.R TO IDENTIFY PROTEIN CODING GENES FROM THE LIST OF IDENTIFIED CRYPTICS AND CDS PHASE >>res_ce_"$2"/Summary_stats.txt
				Rscript get_orf_cds.R res_ce_"$2"/ce_extension_coord_only.txt res_ce_"$2"/Summary_stats.txt
				echo DONE WITH CDS LIST and PHASE FOR CE_EXTENSION list, Please see res_ce_"$2"/protein_coding.bed
				echo DONE WITH CDS LIST and PHASE FOR CE_EXTENSION list, Please see res_ce_"$2"/protein_coding.bed >>res_ce_"$2"/Summary_stats.txt
			else
				echo SKIPPING ORF CALCULATIONS, PLEASE SET orf_flg flag to 1 to repeat these calculations
				echo SKIPPING ORF CALCULATIONS, PLEASE SET orf_flg flag to 1 to repeat these calculations >> res_ce_"$2"/Summary_stats.txt
			fi
			#get difference between two files
			grep -v -f res_ce_"$2"/protein_coding.bed res_ce_"$2"/ce_extension_coord_only.txt > res_ce_"$2"/non_protein_coding.bed

			#NOW UPDATE LISTS BY REMOVING NON CODING TXS
			echo NOW UPDATING CE_EXTENSION LIST WITH CDS, Please see non-protein coding list for non cds events res_ce_"$2"/non_protein_coding.bed
			echo NOW UPDATING CE_EXTENSION LIST WITH CDS, Please see non-protein coding list for non cds events res_ce_"$2"/non_protein_coding.bed >>res_ce_"$2"/Summary_stats.txt


			#Also modify skiptics_unique.bed list by removing non-coding genes
			############
			[ -e res_ce_"$2"/cds_ce_extension_coord_uniq.bed ] && rm res_ce_"$2"/cds_ce_extension_coord_uniq.bed

			nrecrds=$(cat res_ce_"$2"/protein_coding.bed | wc -l)
			nrecrds_u=$(cat res_ce_"$2"/ce_extension_coord_uniq.bed | wc -l)

			readarray -t cds_data < res_ce_"$2"/protein_coding_a.bed

			readarray -t uniq_data < res_ce_"$2"/ce_extension_coord_uniq.bed
			#Also MODIFY

			[ -e res_ce_"$2"/cds_IGV_unique_ce_extension.csv ] && rm res_ce_"$2"/cds_IGV_unique_ce_extension.csv
			[ -e res_ce_"$2"/cds_ce_extension_coord_only.txt ] && rm res_ce_"$2"/cds_ce_extension_coord_only.txt
			readarray -t all_sashimi_csv < res_ce_"$2"/IGV_unique_ce_extension.csv

			i=0
			eventn=0
			while [ $i -lt $nrecrds_u ]
			do
				line21=${uniq_data[$i]}
				line22=${uniq_data[$(($i+1))]}
				line23=${uniq_data[$(($i+2))]}

				#echo got cds line $cds
				j=0
				while [ $j -lt $nrecrds ]
				do
					cds=${cds_data[$j]}

					#line23=${all_data[$(($j+2))]}
					if [ "${cds[*]}" == "${line22[*]}" ] #|| [ "${cds[*]}" == "${line22[*]}" ]
					then
						echo $line21 | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_"$2"/cds_ce_extension_coord_uniq.bed
						echo $line22 | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_"$2"/cds_ce_extension_coord_uniq.bed
						echo $line23 | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_"$2"/cds_ce_extension_coord_uniq.bed
						#also copy sashimi bed and csv files
						#bed21=${all_sashimi_bed[$i]}
						#bed22=${all_sashimi_bed[$(($i+1))]}
						#echo $bed21 | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' >> res_skiptics/cds_skiptics_uniq_sashimi.bed
						#echo $bed22 | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' >> res_skiptics/cds_skiptics_uniq_sashimi.bed
						csv21=${all_sashimi_csv[$j]}
						txid=$(echo $csv21|awk 'BEGIN {FS=","} {print $9}')

						echo $line22 | awk -v tx=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$3-$2,$5,$6,$7,tx}' >> res_ce_"$2"/cds_ce_extension_coord_only.txt

						echo $csv21 >> res_ce_"$2"/cds_IGV_unique_ce_extension.csv
						break
					fi
					((j=j+1))
				done
				((i=i+3))
			done
			echo DONE WITH GENERATION OF VALID CDS bed and csv files for FINAL NT AND AA TRANSLATION, Please see res_ce_"$2"/cds_ce_extension_coord_uniq.bed and res_ce_"$2"/cds_IGV_unique_ce_extension.csv
			echo DONE WITH GENERATION OF VALID CDS bed and csv files for FINAL NT AND AA TRANSLATION, Please see res_ce_"$2"/cds_ce_extension_coord_uniq.bed and res_ce_"$2"/cds_IGV_unique_ce_extension.csv >> res_ce_"$2"/Summary_stats.txt

			##########END CHECK FOR ORF
			echo NOW FINDING NT SEQ FOR res_ce_"$2"/cds_ce_extension_coord_uniq.bed
			#now call bedtools getfasta function to get nt sequence from reference_genome.fa file

			bedtools getfasta -fi /Volumes/SYEDSHAH/MichaelLab/ref_genome/GRCh38.p13.genome.fa -bed res_ce_"$2"/cds_ce_extension_coord_uniq.bed -s > ce_extension_nt.fasta

			#now remove >chr lines from the resulting file
			awk '!/^>chr/' ce_extension_nt.fasta > IR_all_nt1.fasta
			#combine the two files to get desired csv file - chrXX,start,end,genename
			paste -d"\t" res_ce_"$2"/cds_ce_extension_coord_uniq.bed IR_all_nt1.fasta > res_ce_"$2"/cds_ce_extension_nt.bed
			#remove temp file
			rm IR_all_nt1.fasta
			rm ce_extension_nt.fasta
			#and finally transeq compatible
			#awk -F "\t" '{print ">sp|"$7"_"$1"_"$2"_"$3"\n"$8}' res_ce_"$2"/cds_ce_extension_nt.bed > res_ce_"$2"/ce_extension_transeq_in.fasta
			#Add strand infor
			fasta_each_flg=0
			if [ $fasta_each_flg -eq 1 ] #start fasta_each_flg
			then

			awk -F "\t" '{if($6=="+") print ">sp|"$7"_"$1"_"$2"_"$3"_plus""\n"$8; else print ">sp|"$7"_"$1"_"$2"_"$3"_minus""\n"$8}' res_ce_"$2"/cds_ce_extension_nt.bed > res_ce_"$2"/ce_extension_transeq_in.fasta

			#awk -F "\t" '{print ">sp|"$7"-"$1":"$2"-"$3"\n"$8}' res_ce_"$2"/ce_inclusion_nt.bed > ce_inclusion_transeq_in.fasta


			echo now doing the 3-frame translation

			TEMPFILE=$(echo res_ce_"$2"/ce_extension_transeq_in.fasta  | perl -pe 's/\.fasta$/.temp/')
			OUTPUTFILE=$(echo res_ce_"$2"/ce_extension_transeq_in.fasta  | perl -pe 's/\.fasta$/.trans/')

			# Translate sequence
			transeq  -sequence res_ce_"$2"/ce_extension_transeq_in.fasta  -outseq "$TEMPFILE" -frame F

			# Rename sequence
			perl -sape 's/>/>sp|/' "$TEMPFILE" > "$OUTPUTFILE"



			#also remove all newlines from the
			awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' "$OUTPUTFILE" > res_ce_"$2"/CE_EXTENSION_AA.fasta

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
			done < res_ce_"$2"/CE_EXTENSION_AA.fasta > res_ce_"$2"/FINAL_CE_EXTENSION_AA.fasta
		fi #END fasta_each_flg

		echo FINALLY DOING NT AND AA TRANSLATION FOR FUSED COORDINATES FOR CE_EXTENSION EVENTS
		echo FINALLY DOING NT AND AA TRANSLATION FOR FUSED COORDINATES FOR CE_EXTENSION EVENTS >> res_ce_"$2"/Summary_stats.txt
		[ -e res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA.fasta ] && rm res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA.fasta
		[ -e res_ce_"$2"/cds_ce_extension_fused.transeq_in.fasta ] && rm res_ce_"$2"/cds_ce_extension_fused.transeq_in.fasta
		event_i=$COUNT_EVENT
		#Also read unsuccessful_frames_list.csv file to decide the frame of translation
		i=0
		readarray -t cds_data < res_ce_"$2"/cds_successful_frames_list.csv

		cat res_ce_"$2"/cds_ce_extension_nt.bed | while read r1; read r2; read r3
		do
			#echo TRANSLATING EVENT $event_i
			lstrnd=$(echo "$r1" | awk '{print $6}')
			if [[ "$lstrnd" == "+" ]]
			then
			#echo "$r1" "$r2" "$r3"
			#title=$(echo "$r1" "$r2" "$r3" | awk '{print ">sp|"$23"_"$1"_"$2"_"$3"_"$10"_"$11"_"$18"_"$19"_plus"}')
				title=$(echo "$r1" "$r2" "$r3" | awk -v evid=$event_i '{print ">sp|"$23"_"$1"_"$2"_"$3"_"$10"_"$11"_"$18"_"$19"_"evid"_plus"}')
			#echo $title
			echo "$r1" "$r2" "$r3" | awk -v t=$title '{OFS="\t"} {print t"\n"$8$16$24}' > res_ce_"$2"/cds_ce_extension_fused.transeq_in1.fasta
			else
			#echo "$r1" "$r2" "$r3"
			#title=$(echo "$r1" "$r2" "$r3" | awk '{print ">sp|"$23"_"$1"_"$18"_"$19"_"$10"_"$11"_"$2"_"$3"_minus"}')
				title=$(echo "$r1" "$r2" "$r3" | awk -v evid=$event_i '{print ">sp|"$23"_"$1"_"$18"_"$19"_"$10"_"$11"_"$2"_"$3"_"evid"_minus"}')
			#echo $title
			echo "$r1" "$r2" "$r3" | awk -v t=$title '{OFS="\t"} {print t"\n"$8$16$24}' > res_ce_"$2"/cds_ce_extension_fused.transeq_in1.fasta
			fi
			cat res_ce_"$2"/cds_ce_extension_fused.transeq_in1.fasta >> res_ce_"$2"/cds_ce_extension_fused.transeq_in.fasta
			#also get translation fram
			frame_line=${cds_data[$i]}
			framen=$(echo $frame_line | awk 'BEGIN {FS=","} {print $9}')
		
			#now do the 3-frame translation

			TEMPFILE=$(echo res_ce_"$2"/cds_ce_extension_fused.transeq_in1.fasta | perl -pe 's/\.fasta$/.temp/')
			OUTPUTFILE=$(echo res_ce_"$2"/cds_ce_extension_fused.transeq_in1.fasta | perl -pe 's/\.fasta$/.trans/')

			# Translate sequence
			#transeq  -sequence res_ce_"$2"/cds_ce_extension_fused.transeq_in1.fasta  -outseq "$TEMPFILE" -frame=$framen
			transeq  -sequence res_ce_"$2"/cds_ce_extension_fused.transeq_in1.fasta  -outseq "$TEMPFILE" -frame F

			# Rename sequence
			perl -sape 's/>/>sp|/' "$TEMPFILE" > "$OUTPUTFILE"



			#also remove all newlines from the
			awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' "$OUTPUTFILE" >> res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA.fasta

			#also remove "$OUTPUTFILE" file
			rm "$OUTPUTFILE"
			# Remove temp file
			rm "$TEMPFILE"
			((event_i=event_i+1))
			((i=i+1))

		done
		rm res_ce_"$2"/cds_ce_extension_fused.transeq_in1.fasta
		echo GOT FUSED NT and AA FASTA FILE FOR CE_EXTENSION EVENTS, PLEASE SEE res_ce_"$2"/cds_ce_extension_fused.transeq_in.fasta and res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA.fasta file
		echo GOT FUSED NT and AA FASTA FILE FOR CE_EXTENSION EVENTS, PLEASE SEE res_ce_"$2"/cds_ce_extension_fused.transeq_in.fasta and res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA.fasta file >> res_ce_"$2"/Summary_stats.txt

		#NOW CALL R SCRIPT TO MATCH cannonical frame to get final_aa.fasta
		echo NOW CALLING check_aaV4_allFrames TO MATCH cannonical frame to get final_aa.fasta  >> res_ce_"$2"/Summary_stats.txt 
		Rscript check_aaV4_allFrames.R res_ce_"$2"/aa.fasta res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA.fasta res_ce_"$2"/Summary_stats.txt res_ce_"$2"/cds_ce_extension_coord_only.txt res_ce_"$2"/cds_IGV_unique_ce_extension.csv
		#cp final_aa.fasta res_ce_"$2"/final_aa.fasta
		stop_codon_flg=1
		if [ $stop_codon_flg -eq 1 ] && [ -f res_ce_"$2"/final_aa.fasta ]
		then

			#now truncate lines at *
			#sed "s/*.*//" res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA.fasta > res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA1.fasta
			sed "s/*.*//" res_ce_"$2"/final_aa.fasta > res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA1.fasta
			echo DONE WITH TRUNCATION OF FASTA LINES @ STOP CODONS FOR CDS FASTA, PLEASe SEE res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA1.fasta file
			echo DONE WITH TRUNCATION OF FASTA LINES @ STOP CODONS FOR CDS FASTA, PLEASe SEE res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA1.fasta file >> res_ce_"$2"/Summary_stats.txt
			#Also truncate lines at * for allframes files
			sed "s/*.*//" res_ce_"$2"/all_frames_final_aa.fasta > res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA1.fasta
			echo DONE WITH TRUNCATION OF FASTA LINES @ STOP CODONS FOR ALL FRAMES FASTA, PLEASe SEE res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA.fasta file
			echo DONE WITH TRUNCATION OF FASTA LINES @ STOP CODONS FOR ALL FRAMES FASTA, PLEASe SEE res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA.fasta file >> res_ce_"$2"/Summary_stats.txt
			else
				cat res_ce_"$2"/final_aa.fasta > res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA1.fasta
				cat res_ce_"$2"/all_frames_final_aa.fasta > res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA1.fasta
			fi

			#delete empty line if any
			#echo NOW REMOVING EMPTY LINES LEFT BY EVENTS WITH STOP CODONS IN THE BEGINNING ""ARHGAP22_chr10"" is one EXAMPLE
			echo NOW REMOVING EMPTY LINES LEFT BY EVENTS WITH STOP CODONS IN THE BEGINNING ""ARHGAP22_chr10"" is one EXAMPLE >> res_ce_"$2"/Summary_stats.txt
			[ -e res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA2.fasta ] && rm res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA2.fasta
			[ -e res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA_Empty.fasta ] && rm res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA_Empty.fasta

			i=0
			cat res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA1.fasta | while read r1; read r2
			do

				if [ "$r2" != "" ]; then
						echo "$r1" >> res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA2.fasta
						echo "$r2" >> res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA2.fasta
					else
						echo "$r1" >> res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA_Empty.fasta
						((i=i+1))
				fi

			done
			if [ $(echo res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA_Empty.fasta|wc -l) -gt 0 ]
			then
				echo DONE WITH REMOVING EMPTY LINES GOT $(echo res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA_Empty.fasta|wc -l) CE_EXTENSION HAVING STOP CODON AT THE BEGINNING, PLEASE SEE res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA_Empty.fasta file
				echo DONE WITH REMOVING EMPTY LINES GOT $(echo res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA_Empty.fasta|wc -l) CE_EXTENSION HAVING STOP CODON AT THE BEGINNING, PLEASE SEE res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA_Empty.fasta file >> res_ce_"$2"/Summary_stats.txt
			else
				cat res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA1.fasta > res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA2.fasta
			fi

			##########ALLFRAMES FASTA FILE
			echo NOW REMOVING EMPTY LINES LEFT BY EVENTS WITH STOP CODONS IN THE BEGINNING ""ARHGAP22_chr10"" is one EXAMPLE >> res_ce_"$2"/Summary_stats.txt
			[ -e res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA2.fasta ] && rm res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA2.fasta
			[ -e res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA_Empty.fasta ] && rm res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA_Empty.fasta

			i=0
			cat res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA1.fasta | while read r1; read r2
			do

				if [ "$r2" != "" ]; then
						echo "$r1" >> res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA2.fasta
						echo "$r2" >> res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA2.fasta
					else
						echo "$r1" >> res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA_Empty.fasta
						((i=i+1))
				fi

			done
			if [ $(echo res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA_Empty.fasta|wc -l) -gt 0 ]
			then
				echo DONE WITH REMOVING EMPTY LINES GOT $(echo res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA_Empty.fasta|wc -l) CE_EXTENSION HAVING STOP CODON AT THE BEGINNING, PLEASE SEE res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA_Empty.fasta file
				echo DONE WITH REMOVING EMPTY LINES GOT $(echo res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA_Empty.fasta|wc -l) CE_EXTENSION HAVING STOP CODON AT THE BEGINNING, PLEASE SEE res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA_Empty.fasta file >> res_ce_"$2"/Summary_stats.txt
			else
				cat res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA1.fasta > res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA2.fasta
			fi


		##########



			#empty_lines=$(grep -n '^$' res_skiptics/new_fasta.txt | cut -d: -f1)
			echo STARTED REMOVING AA "<" 8
			echo STARTED REMOVING AA "<" 8  >> res_ce_"$2"/Summary_stats.txt


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
			done < res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA2.fasta > res_ce_"$2"/cds_PEAKS_CE_EXTENSION_FUSED_AA1.fasta

			#also remove trailing X in fasta files
			awk 'NR%2==0{sub(/X$/,"")}1' res_ce_"$2"/cds_PEAKS_CE_EXTENSION_FUSED_AA1.fasta > res_ce_"$2"/cds_PEAKS_CE_EXTENSION_FUSED_AA.fasta
			[ -e res_ce_"$2"/cds_PEAKS_CE_EXTENSION_FUSED_AA1.fasta ] && rm res_ce_"$2"/cds_PEAKS_CE_EXTENSION_FUSED_AA1.fasta

			#########ALLFRAMES
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
			done < res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA2.fasta > res_ce_"$2"/AllFrames_PEAKS_CE_EXTENSION_FUSED_AA1.fasta

			#also remove trailing X in fasta files
			awk 'NR%2==0{sub(/X$/,"")}1' res_ce_"$2"/AllFrames_PEAKS_CE_EXTENSION_FUSED_AA1.fasta > res_ce_"$2"/AllFrames_PEAKS_CE_EXTENSION_FUSED_AA.fasta
			[ -e res_ce_"$2"/AllFrames_PEAKS_CE_EXTENSION_FUSED_AA1.fasta ] && rm res_ce_"$2"/AllFrames_PEAKS_CE_EXTENSION_FUSED_AA1.fasta

		###################

			echo DONE REMOVING AA "<" 8, PLEASE SEE res_ce_"$2"/cds_PEAKS_CE_EXTENSION_FUSED_AA.fasta and res_ce_"$2"/AllFrames_PEAKS_CE_EXTENSION_FUSED_AA.fasta
			echo DONE REMOVING AA "<" 8, PLEASE SEE res_ce_"$2"/cds_PEAKS_CE_EXTENSION_FUSED_AA.fasta and res_ce_"$2"/AllFrames_PEAKS_CE_EXTENSION_FUSED_AA.fasta >> res_ce_"$2"/Summary_stats.txt

			echo FOR FINAL STATISTICS ON CE EXTENSION EVENTS, PLEASE SEE res_ce_"$2"/Summary_stats.txt file
			echo ALL DONE - HOPEFULLY - SUCCESSFULLY >> res_ce_"$2"/Summary_stats.txt

		else
				#echo Total ce_extension events are: $nrecrds_ir >> res_ce_"$2"/Summary_stats.txt
				echo CAUTION!  No final_aa.fasta FILE WAS GENERATED, PLEASE TAKE A CLOSER LOOK   >> res_ce_"$2"/Summary_stats.txt
				#echo FILE res_ce_"$2"/ce_extension_coord.bed HAS ZERO RECORDS SO EXITING >> res_ce_"$2"/Summary_stats.txt
		fi

		echo FOR FINAL STATISTICS ON CE EVENTS, PLEASE SEE res_ce_"$2"/Summary_stats.txt file

		echo now GENERATING SASHIMI PLOTS >> res_ce_"$2"/Summary_stats.txt
		echo now GENERATING SASHIMI PLOTS
		if [ -f res_ce_"$2"/ce_extension_coord_uniq_sashimi.bed ]
		then
			num_ce_ext=$(cat res_ce_"$2"/ce_extension_coord_uniq_sashimi.bed | wc -l | awk '{print $1/3}')
			echo STARTING SASHIMI PLOTS FOR $num_ce_ext Events, For events ">" 100 It Might take several hours >> res_ce_"$2"/Summary_stats.txt
			#source run_sashimi.sh 4
			source run_sashimiV1.sh res_ce_"$2"/IGV_unique_ce_extension.csv res_ce_"$2"/ce_extension_coord_uniq_sashimi.bed 4
		else
			echo File res_ce_"$2"/ce_extension_coord_uniq_sashimi.bed does not exists so, NO CE_EXTENSION SASHIMI PLOTS ARE GENERATED
			echo File res_ce_"$2"/ce_extension_coord_uniq_sashimi.bed does not exists so, NO CE_EXTENSION SASHIMI PLOTS ARE GENERATED >> res_ce_"$2"/Summary_stats.txt

		fi

		echo CLEANING UP >> res_ce_"$2"/Summary_stats.txt

		rm res_ce_"$2"/IGV_ce_inclusion.csv
		rm res_ce_"$2"/IGV_unique_ce_extension.csv
		mv res_ce_"$2"/cds_IGV_unique_ce_extension.csv res_ce_"$2"/CE_extension.csv
		rm res_ce_"$2"/cds_ce_extension_coord_uniq.bed
		rm res_ce_"$2"/cds_ce_extension_nt.bed
		rm res_ce_"$2"/cds_successful_frames_list.csv
		rm res_ce_"$2"/cds_unsuccessful_frames_list.csv
		rm res_ce_"$2"/ce_all_scan_intron.bed
		rm res_ce_"$2"/ce_all_scan_range.bed
		rm res_ce_"$2"/ce_all_scan_range_junctions.bed
		rm res_ce_"$2"/ce_all_scan_unique_range.bed
		rm res_ce_"$2"/ce_extension_coord.bed
		rm res_ce_"$2"/ce_extension_coord_uniq.bed
		rm res_ce_"$2"/ce_extension_coord_uniq_sashimi.bed
		rm res_ce_"$2"/non_protein_coding.bed
		rm res_ce_"$2"/protein_coding.bed
		rm res_ce_"$2"/protein_coding_a.bed

		rm res_ce_"$2"/AllFrames_CE_EXTENSION_FUSED_AA1.fasta
		rm res_ce_"$2"/AllFrames_PEAKS_CE_EXTENSION_FUSED_AA.fasta
		rm res_ce_"$2"/aa.fasta
		rm res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA.fasta
		rm res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA1.fasta
		#cds_PEAKS_CE_EXTENSION_FUSED_AA.fasta
		#cds_ce_extension_fused.transeq_in.fasta


		echo ALL DONE - HOPEFULLY - SUCCESSFULLY >> res_ce_"$2"/Summary_stats.txt
		echo Done CE_EXTENSION, Now Exitting >> res_ce_"$2"/Summary_stats.txt
		exit 1

	fi 
	######### End Decide ce_extension


	#Step 4. using majiq_coverages_automateV1.R script identify final ce_extension_coord.bed and ce_extension_coord.bed files
	fact=3
	total_majiq_events=$(cat $inpfile | wc -l)
	unknow_events=0
	[ -e res_ce_"$2"/non_ce_events.csv ] && unknow_events=$(cat res_ce_"$2"/non_ce_events.csv | wc -l)
	not_found=0
	[ -e res_ce_"$2"/EnsDB_tx_not_found.csv ] && not_found=$(cat res_ce_"$2"/EnsDB_tx_not_found.csv | wc -l)
	problematic=0
	[ -e res_ce_"$2"/IGV_problematic_junctions.csv ] && problematic=$(cat res_ce_"$2"/IGV_problematic_junctions.csv  | wc -l)
	#total_events=$(($all_ce+$unknow_events+$not_found+$problematic))
	unseccesful_r1=0
	[ -e res_ce_"$2"/skipped_ce.csv ] && unseccesful_r1=$(cat res_ce_"$2"/skipped_ce.csv|wc -l)
	unseccesful_r=$(($unseccesful_r1/3))
	echo Total events read are: $total_majiq_events >> res_ce_"$2"/Summary_stats.txt
	echo Out of these $total_majiq_events total events >> res_ce_"$2"/Summary_stats.txt
	echo Events not found in EnsDB are: $not_found, please see res_ce_"$2"/EnsDB_tx_not_found.csv file >> res_ce_"$2"/Summary_stats.txt
	echo Events that are not CE: $unknow_events , please see res_ce_"$2"/non_ce_events.txt and res_ce_"$2"/non_ce_events.csv file >> res_ce_"$2"/Summary_stats.txt
	echo Events that are somewhat problematic: $problematic , please see res_ce_"$2"/problematic_junctions.txt and res_ce_"$2"/IGV_problematic_junctions.csv files >> res_ce_"$2"/Summary_stats.txt
	ce_boundary_events=$(($total_majiq_events-$not_found-$unknow_events-$problematic))
	echo '#######################################' >> res_ce_"$2"/Summary_stats.txt
	echo Now starting ce_boundary calculations for a remaining total of: $ce_boundary_events events, BY INVOKING R SCRIPT  >> res_ce_"$2"/Summary_stats.txt
	#echo Identifying ce_boundary coordinates
	echo identifying ce boundaries by calling Auto_CoverV4_layered_intronV3.R script
	if [ $ce_boundary_events -gt 0 ]
	then
		Rscript Auto_CoverV4_layered_intronV3.R res_ce_"$2"/IGV_ce_inclusion.csv res_ce_"$2"/ce_all_scan_range_junctions.bed $2 res_ce_"$2"/ce_all_scan_intron.bed #Auto_CoverV4_layered_intron.R #Auto_CoverV2.R
	else
		echo CE_BOUNDARY_EVENTS ARE $ce_boundary_events SO ABONDONED CE_BOUNDARY CALCULATIONS
	fi
	echo BACK FROM CE_BOUNDARY CALCULATIONS
	echo '################################################' >> res_ce_"$2"/Summary_stats.txt
	#Rscript majiq_coverages_automateV1.R #majiq_coverages_automateV1_227.R #majiq_coverages_automateV1.R
	#Step 4 ENDS


	#Step 5. Finally identify nt and AA sequences

	#get nt and aa for ce boundaries

	#First remove duplicates from ce' boundary coordinates
	#NEW CODE STARTS HERE
	nrecrds=0
	nrecrdst=0
	[ -e res_ce_"$2"/ce_inclusion_coord.bed ] && nrecrds=$(cat res_ce_"$2"/ce_inclusion_coord.bed | wc -l) && nrecrdst=$(($nrecrds/3))


	if [ $nrecrds -ne 0 ]
	then
		echo Back From R Session, Now checking "FOR DUPLICATES" for a total of ce_inclusion events: $nrecrdst  >> res_ce_"$2"/Summary_stats.txt
		#echo Checking any Repeated ce_boundary events >> res_ce_"$2"/Summary_stats.txt
		readarray -t all_data < res_ce_"$2"/ce_inclusion_coord.bed
		#also read ce_inclusion_coord_sashimi.bed to get TxID
		readarray -t sashimi_data < res_ce_"$2"/ce_inclusion_coord_sashimi.bed
		#also read csv file
		readarray -t csv_data < res_ce_"$2"/IGV_R_returned_ce_inclusion.csv
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
				echo $line11 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_"$2"/ce_inclusion_coord_uniq.bed
				echo $line12 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_"$2"/ce_inclusion_coord_uniq.bed
				echo $line13 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_"$2"/ce_inclusion_coord_uniq.bed

				#also save bed file for sashimi plots
				#echo $line11 | awk -v TX=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,TX}' >> res_ce_"$2"/ce_inclusion_coord_uniq_sashimi.bed
				#echo $line12 | awk -v TX=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,TX}' >> res_ce_"$2"/ce_inclusion_coord_uniq_sashimi.bed
				#echo $line13 | awk -v TX=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,TX}' >> res_ce_"$2"/ce_inclusion_coord_uniq_sashimi.bed
				echo $lines1 | awk  'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' >> res_ce_"$2"/ce_inclusion_coord_uniq_sashimi.bed
				echo $lines2 | awk  'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' >> res_ce_"$2"/ce_inclusion_coord_uniq_sashimi.bed
				#Also SAVING CE coordinates

				echo $lines2 | awk  'BEGIN {OFS="\t"} {print $1,$2,$3,$3-$2,$5,$6,$7,$8}' >> res_ce_"$2"/ce_inclusion_coord_only.txt
				echo $lines3 | awk  'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' >> res_ce_"$2"/ce_inclusion_coord_uniq_sashimi.bed

				#also save csv file
				echo $csv_ln >> res_ce_"$2"/IGV_unique_ce_inclusion.csv
				echo $csv_ln | awk 'BEGIN {FS=","} {print $2,$3,$4,$7,$10,$11}' | awk '{OFS=","} {print $1,$2,$3,$4,$5,$6}' >> res_ce_"$2"/CE_inclusion.csv
			else
				echo $line11 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_"$2"/ce_inclusion_coord_repeated.bed
				echo $line12 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_"$2"/ce_inclusion_coord_repeated.bed
				echo $line13 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_"$2"/ce_inclusion_coord_repeated.bed

			fi


		done

		nrecrds_uniqt=$(cat res_ce_"$2"/ce_inclusion_coord_uniq.bed | wc -l)
		nrecrds_uniq=$(($nrecrds_uniqt/3))
		repeated_ce_boundary_events=$(($nrecrdst-$nrecrds_uniq))
		echo Total unique ce_inclusion events are: $nrecrds_uniq, please see res_ce_"$2"/ce_inclusion_coord_uniq.bed and res_ce_"$2"/IGV_unique_ce_inclusion.csv >> res_ce_"$2"/Summary_stats.txt
		#echo total unique ce_boundary events read are $nrecrds_uniq
		echo Total repeated ce_inclusion events are: $repeated_ce_boundary_events, please see res_ce_"$2"/ce_inclusion_coord_repeated.bed >> res_ce_"$2"/Summary_stats.txt
		echo total repeated ce_inclusion events were $repeated_ce_boundary_events

		######## NEW SECTION FOR SELECTING CANONICAL FRAME ONLY
		##########CHECK FOR orf
		orf_flg=1
		if [ $orf_flg -eq 1 ]
		then

			echo NOW CALLING get_orf_cds.R TO IDENTIFY PROTEIN CODING GENES FROM THE LIST OF IDENTIFIED CRYPTICS AND CDS PHASE
			echo NOW CALLING get_orf_cds.R TO IDENTIFY PROTEIN CODING GENES FROM THE LIST OF IDENTIFIED CRYPTICS AND CDS PHASE >>res_ce_"$2"/Summary_stats.txt

			#Rscript get_orf_cds.R res_ce_"$2"/ce_inclusion_coord_only.txt res_ce_"$2"/Summary_stats.txt

			Rscript get_orf_cds.R res_ce_"$2"/ce_inclusion_coord_only.txt res_ce_"$2"/Summary_stats.txt cryptics

			#[ -e protein_coding.bed ] && mv protein_coding.bed res_ce_"$2"/protein_coding.bed
			#[ -e cds_unsuccessful_frames_list.csv ] && mv cds_unsuccessful_frames_list.csv res_ce_"$2"/cds_unsuccessful_frames_list.csv
			#[ -e cds_successful_frames_list.csv ] && mv cds_successful_frames_list.csv res_ce_"$2"/cds_successful_frames_list.csv
			#[ -e protein_coding_a.bed ] && mv protein_coding_a.bed res_ce_"$2"/protein_coding_a.bed #this is to compare with skiptics_unique.bed

			echo DONE WITH CDS LIST and PHASE FOR CE_INCLUSION list, Please see res_ce_"$2"/protein_coding.bed
			echo DONE WITH CDS LIST and PHASE FOR CE_INCLUSION list, Please see res_ce_"$2"/protein_coding.bed >>res_ce_"$2"/Summary_stats.txt
		else
			echo SKIPPING ORF CALCULATIONS, PLEASE SET orf_flg flag to 1 to repeat these calculations
			echo SKIPPING ORF CALCULATIONS, PLEASE SET orf_flg flag to 1 to repeat these calculations >> res_ce_"$2"/Summary_stats.txt
		fi
		#get difference between two files
		grep -v -f res_ce_"$2"/protein_coding.bed res_ce_"$2"/ce_inclusion_coord_only.txt > res_ce_"$2"/non_protein_coding.bed


		#NOW UPDATE LISTS BY REMOVING NON CODING TXS
		echo FOUND ??????? $(cat res_ce_"$2"/non_protein_coding.bed | wc -l) ????????? EVENTS WHICH ARE EITHER non-coding OR HAVE 0 LENGTH CDS
		echo FOUND ??????? $(cat res_ce_"$2"/non_protein_coding.bed | wc -l) ????????? EVENTS WHICH ARE EITHER non-coding OR HAVE 0 LENGTH CDS >> res_ce_"$2"/Summary_stats.txt
		echo NOW UPDATING CE_INCLUSION LIST WITH CDS, Please see non-protein coding list for non cds events res_ce_"$2"/non_protein_coding.bed
		echo NOW UPDATING CE_INCLUSION LIST WITH CDS, Please see non-protein coding list for non cds events res_ce_"$2"/non_protein_coding.bed >>res_ce_"$2"/Summary_stats.txt


		#Also modify skiptics_unique.bed list by removing non-coding genes
		############
		[ -e res_ce_"$2"/cds_ce_inclusion_coord_uniq.bed ] && rm res_ce_"$2"/cds_ce_inclusion_coord_uniq.bed
		[ -e res_ce_"$2"/cds_ce_inclusion_coord_only.txt ] && rm res_ce_"$2"/cds_ce_inclusion_coord_only.txt
		nrecrds=$(cat res_ce_"$2"/protein_coding.bed | wc -l)
		nrecrds_u=$(cat res_ce_"$2"/ce_inclusion_coord_uniq.bed | wc -l)

		readarray -t cds_data < res_ce_"$2"/protein_coding_a.bed

		readarray -t uniq_data < res_ce_"$2"/ce_inclusion_coord_uniq.bed
		#Also MODIFY

		[ -e res_ce_"$2"/cds_IGV_unique_ce_inclusion.csv ] && rm res_ce_"$2"/cds_IGV_unique_ce_inclusion.csv
		readarray -t all_sashimi_csv < res_ce_"$2"/IGV_unique_ce_inclusion.csv

		i=0
		eventn=0
		while [ $i -lt $nrecrds_u ]
		do
			line21=${uniq_data[$i]}
			line22=${uniq_data[$(($i+1))]}
			line23=${uniq_data[$(($i+2))]}

				#echo got cds line $cds
				j=0
			while [ $j -lt $nrecrds ]
			do
				cds=${cds_data[$j]}

				#line23=${all_data[$(($j+2))]}
				if [ "${cds[*]}" == "${line22[*]}" ] #|| [ "${cds[*]}" == "${line22[*]}" ]
				then

					csv21=${all_sashimi_csv[$eventn]}
					#get TxID
					txid=$(echo $csv21|awk 'BEGIN {FS=","} {print $9}')

					echo $line21 | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_"$2"/cds_ce_inclusion_coord_uniq.bed
					echo $line22 | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_"$2"/cds_ce_inclusion_coord_uniq.bed
					echo $line23 | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_ce_"$2"/cds_ce_inclusion_coord_uniq.bed

					echo $line22 | awk -v tx=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$3-$2,$5,$6,$7,tx}' >> res_ce_"$2"/cds_ce_inclusion_coord_only.txt
					#also copy sashimi bed and csv files
					#bed21=${all_sashimi_bed[$i]}
					#bed22=${all_sashimi_bed[$(($i+1))]}
					#echo $bed21 | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' >> res_skiptics/cds_skiptics_uniq_sashimi.bed
					#echo $bed22 | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' >> res_skiptics/cds_skiptics_uniq_sashimi.bed
					#csv21=${all_sashimi_csv[$j]}

					echo $csv21 >> res_ce_"$2"/cds_IGV_unique_ce_inclusion.csv
					break
				fi
				((j=j+1))
			done
			((i=i+3))
			((eventn=eventn+1))
		done
		echo AFTER REMOVING non-cds/zero-cds events, WE HAVE "*****" $(cat res_ce_"$2"/cds_IGV_unique_ce_inclusion.csv|wc -l) EVENTS "******" FOR FINAL PROCESSING, PLEASE SEE res_ce_"$2"/cds_IGV_unique_ce_inclusion.csv and res_ce_"$2"/cds_ce_inclusion_coord_uniq.bed FILES FOR DETAILS
		echo AFTER REMOVING non-cds/zero-cds events, WE HAVE "*****" $(cat res_ce_"$2"/cds_IGV_unique_ce_inclusion.csv|wc -l) EVENTS "******" FOR FINAL PROCESSING, PLEASE SEE res_ce_"$2"/cds_IGV_unique_ce_inclusion.csv and res_ce_"$2"/cds_ce_inclusion_coord_uniq.bed FILES FOR DETAILS >> res_ce_"$2"/Summary_stats.txt
		echo DONE WITH GENERATION OF VALID CDS bed and csv files for FINAL NT AND AA TRANSLATION, Please see res_ce_"$2"/cds_ce_inclusion_coord_uniq.bed and res_ce_"$2"/cds_IGV_unique_ce_inclusion.csv
		echo DONE WITH GENERATION OF VALID CDS bed and csv files for FINAL NT AND AA TRANSLATION, Please see res_ce_"$2"/cds_ce_inclusion_coord_uniq.bed and res_ce_"$2"/cds_IGV_unique_ce_inclusion.csv >> res_ce_"$2"/Summary_stats.txt

		echo                                           >> res_ce_"$2"/Summary_stats.txt
		echo "---------------------------------------" >> res_ce_"$2"/Summary_stats.txt


		########END NEW SECTION FOR SELECTING CANONICAL FRAME ONLY

		#now copy back data to ce_coord
		#cat res_ce_"$2"/ce_inclusion_coord_uniq.bed	> res_ce_"$2"/ce_extension_coord.bed
		#uniqrecrds1=$(cat res_ce_"$2"/ce_extension_coord.bed | wc -l)
		#uniqrecrds=$(($uniqrecrds1/3))
		#echo total records after duplicate remval are $uniqrecrds
		#NEW CODE ENDS HERE
		#cat res_ce_"$2"/ce_extension_coord.bed > ce_coord_temp.bed
		#cat ce_coord_temp.bed | awk '!a[$0]++'  > res_ce_"$2"/ce_extension_coord.bed
		#rm ce_coord_temp.bed

		#now call bedtools getfasta function to get nt sequence from reference_genome.fa file


		echo STARTED NT SEQ FOR res_ce_"$2"/cds_ce_inclusion_coord_uniq.bed
		bedtools getfasta -fi /Volumes/SYEDSHAH/MichaelLab/ref_genome/GRCh38.p13.genome.fa -bed res_ce_"$2"/cds_ce_inclusion_coord_uniq.bed -s > ce_inclusion_nt.fasta

		#now remove >chr lines from the resulting file
		awk '!/^>chr/' ce_inclusion_nt.fasta > ce_all_nt1.fasta
		#combine the two files to get desired csv file - chrXX,start,end,genename
		paste -d"\t" res_ce_"$2"/cds_ce_inclusion_coord_uniq.bed ce_all_nt1.fasta > res_ce_"$2"/ce_inclusion_nt.bed
		#remove temp file
		rm ce_all_nt1.fasta
		rm ce_inclusion_nt.fasta
		#and finally transeq compatible
		#awk -F "\t" '{print ">sp|"$7"_"$1"_"$2"_"$3"\n"$8}' res_ce_"$2"/ce_inclusion_nt.bed > res_ce_"$2"/ce_inclusion_transeq_in.fasta
		#Add strand infor
		awk -F "\t" '{if($6=="+") print ">sp|"$7"_"$1"_"$2"_"$3"_plus""\n"$8; else print ">sp|"$7"_"$1"_"$2"_"$3"_minus""\n"$8}' res_ce_"$2"/ce_inclusion_nt.bed > res_ce_"$2"/ce_inclusion_transeq_in.fasta

		#awk -F "\t" '{print ">sp|"$7"-"$1":"$2"-"$3"\n"$8}' res_ce_"$2"/ce_inclusion_nt.bed > ce_inclusion_transeq_in.fasta

		fasta_each_flg=0
		if [ $fasta_each_flg -eq 1 ] #start fasta_each_flg
		then

			echo now doing the 3-frame translation

			TEMPFILE=$(echo res_ce_"$2"/ce_inclusion_transeq_in.fasta  | perl -pe 's/\.fasta$/.temp/')
			OUTPUTFILE=$(echo res_ce_"$2"/ce_inclusion_transeq_in.fasta  | perl -pe 's/\.fasta$/.trans/')
			#echo doing 3-frame translation
			# Translate sequence
			transeq  -sequence res_ce_"$2"/ce_inclusion_transeq_in.fasta  -outseq "$TEMPFILE" -frame F

			# Rename sequence
			perl -sape 's/>/>sp|/' "$TEMPFILE" > "$OUTPUTFILE"



			#also remove all newlines from the
			awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' "$OUTPUTFILE" > res_ce_"$2"/CE_INCLUSION_AA.fasta

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
			done < res_ce_"$2"/CE_INCLUSION_AA.fasta > res_ce_"$2"/FINAL_CE_INCLUSION_AA.fasta
		fi ##end fasta_each_flg
		#Also concatenate
		echo NOW DOING NT AND AA TRANSLATION FOR FUSED COORDINATES FOR CE_INCLUSION
		echo NOW DOING NT AND AA TRANSLATION FOR FUSED COORDINATES FOR CE_INCLUSION >> res_ce_"$2"/Summary_stats.txt
		#[ -e res_ce_"$2"/ce_inclusion_fused.transeq_in.fasta ] &&  rm res_ce_"$2"/ce_inclusion_fused.transeq_in.fasta

		[ -e res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA.fasta ] && rm res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA.fasta
		[ -e res_ce_"$2"/cds_ce_inclusion_fused.transeq_in.fasta ] && rm res_ce_"$2"/cds_ce_inclusion_fused.transeq_in.fasta


		#rm res_ce_"$2"/ce_inclusion_fused.transeq_in.fasta
		#ALSO ADD even_number flag for sashimi plots against PEAKS results
		event_i=$COUNT_EVENT
		cat res_ce_"$2"/ce_inclusion_nt.bed | while read r1; read r2; read r3
		do
			lstrnd=$(echo "$r1" | awk '{print $6}')
			if [[ "$lstrnd" == "+" ]]
			then
				#echo "$r1" | awk '{print $8}' #"$r2" "$r3"
				#title=$(echo "$r1" "$r2" "$r3" | awk '{print ">sp|"$23"_"$1"_"$2"_"$3"_"$10"_"$11"_"$18"_"$19"_plus"}')
					#ADD EVENT Number flag
					title=$(echo "$r1" "$r2" "$r3" | awk -v evid=$event_i '{print ">sp|"$23"_"$1"_"$2"_"$3"_"$10"_"$11"_"$18"_"$19"_"evid"_plus"}')
				#echo $title
				echo "$r1" "$r2" "$r3" | awk -v t=$title '{OFS="\t"} {print t"\n"$8$16$24}' > res_ce_"$2"/cds_ce_inclusion_fused.transeq_in1.fasta
			else
				#echo "$r1" "$r2" "$r3"
				#title=$(echo "$r3" "$r2" "$r1" | awk '{print ">sp|"$23"_"$1"_"$18"_"$19"_"$10"_"$11"_"$2"_"$3"_minus"}')
					#Add EVENT Number flag
					title=$(echo "$r3" "$r2" "$r1" | awk -v evid=$event_i '{print ">sp|"$23"_"$1"_"$18"_"$19"_"$10"_"$11"_"$2"_"$3"_"evid"_minus"}')
				#echo $title
				echo "$r1" "$r2" "$r3" | awk -v t=$title '{OFS="\t"} {print t"\n"$8$16$24}' > res_ce_"$2"/cds_ce_inclusion_fused.transeq_in1.fasta
			fi
			cat res_ce_"$2"/cds_ce_inclusion_fused.transeq_in1.fasta >> res_ce_"$2"/cds_ce_inclusion_fused.transeq_in.fasta

			TEMPFILE=$(echo res_ce_"$2"/cds_ce_inclusion_fused.transeq_in1.fasta | perl -pe 's/\.fasta$/.temp/')
			OUTPUTFILE=$(echo res_ce_"$2"/cds_ce_inclusion_fused.transeq_in1.fasta | perl -pe 's/\.fasta$/.trans/')

			# Translate sequence
			transeq  -sequence res_ce_"$2"/cds_ce_inclusion_fused.transeq_in1.fasta  -outseq "$TEMPFILE" -frame F

			# Rename sequence
			perl -sape 's/>/>sp|/' "$TEMPFILE" > "$OUTPUTFILE"



			#also remove all newlines from the
			awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' "$OUTPUTFILE" >> res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA.fasta

			#also remove "$OUTPUTFILE" file
			rm "$OUTPUTFILE"
			# Remove temp file
			rm "$TEMPFILE"

			((event_i=event_i+1))
		done
		rm res_ce_"$2"/cds_ce_inclusion_fused.transeq_in1.fasta

		echo GOT FUSED NT and AA FASTA FILE FOR CE_INCLUSION EVENTS, PLEASE SEE res_ce_"$2"/cds_ce_inclusion_fused.transeq_in.fasta and res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA.fasta file
		echo GOT FUSED NT and AA FASTA FILE FOR CE_INCLUSION EVENTS, PLEASE SEE res_ce_"$2"/cds_ce_inclusion_fused.transeq_in.fasta and res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA.fasta file >> res_ce_"$2"/Summary_stats.txt

		#NOW CALL R SCRIPT TO MATCH cannonical frame to get final_aa.fasta
		echo NOW CALL check_aaV4_allFrames.R TO MATCH cannonical frame to get final_aa.fasta >>res_ce_"$2"/Summary_stats.txt
		Rscript check_aaV4_allFrames.R res_ce_"$2"/aa.fasta res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA.fasta res_ce_"$2"/Summary_stats.txt res_ce_"$2"/cds_ce_inclusion_coord_only.txt res_ce_"$2"/cds_IGV_unique_ce_inclusion.csv
		echo Back from check_aaV4_allFrames.R TO MATCH cannonical frame to get final_aa.fasta >>res_ce_"$2"/Summary_stats.txt

			if [ -f res_ce_"$2"/final_aa.fasta ] #final_aa.fasta exists
			then

				echo AFTER REMOVING NON-CDS AA FROM res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA.fasta FILES, WE FINALLY HAVE "***" $(cat res_ce_"$2"/final_aa.fasta|wc -l|awk '{print $1/2}') EVENTS "***" THAT ARE TRANSLATED, PLEASE SEE res_ce_"$2"/final_aa.fasta FILE
				echo AFTER REMOVING NON-CDS AA FROM res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA.fasta FILES, WE FINALLY HAVE "***" $(cat res_ce_"$2"/final_aa.fasta|wc -l|awk '{print $1/2}') EVENTS "***" THAT ARE TRANSLATED, PLEASE SEE res_ce_"$2"/final_aa.fasta FILE >> res_ce_"$2"/Summary_stats.txt
				stop_codon_flg=1
				if [ $stop_codon_flg -eq 1 ]
				then

					#now truncate lines at *
					#sed "s/*.*//" res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA.fasta > res_ce_"$2"/cds_CE_EXTENSION_FUSED_AA1.fasta
					sed "s/*.*//" res_ce_"$2"/final_aa.fasta > res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA1.fasta
					echo DONE WITH TRUNCATION OF FASTA LINES @ STOP CODONS, PLEASe SEE res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA1.fasta file
					echo DONE WITH TRUNCATION OF FASTA LINES @ STOP CODONS, PLEASe SEE res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA1.fasta file >> res_ce_"$2"/Summary_stats.txt

					#Also truncate lines at * for allframes files
					sed "s/*.*//" res_ce_"$2"/all_frames_final_aa.fasta > res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA1.fasta
					echo DONE WITH TRUNCATION OF FASTA LINES @ STOP CODONS FOR ALL FRAMES FASTA, PLEASe SEE res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA1.fasta file
					echo DONE WITH TRUNCATION OF FASTA LINES @ STOP CODONS FOR ALL FRAMES FASTA, PLEASe SEE res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA1.fasta file >> res_ce_"$2"/Summary_stats.txt
				else
					cat res_ce_"$2"/final_aa.fasta > res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA1.fasta
					cat res_ce_"$2"/all_frames_final_aa.fasta > res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA1.fasta
				fi

				#delete empty line if any
				#echo NOW REMOVING EMPTY LINES LEFT BY EVENTS WITH STOP CODONS IN THE BEGINNING ""ARHGAP22_chr10"" is one EXAMPLE
				echo NOW REMOVING EMPTY LINES LEFT BY EVENTS WITH STOP CODONS IN THE BEGINNING ""ARHGAP22_chr10"" is one EXAMPLE >> res_ce_"$2"/Summary_stats.txt
				[ -e res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA2.fasta ] && rm res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA2.fasta
				[ -e res_ce_"$2"/cds_INCLUSION_FUSED_AA_Empty.fasta ] && rm res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA_Empty.fasta

				i=0
				cat res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA1.fasta | while read r1; read r2
				do
					if [ "$r2" != "" ]; then
						echo "$r1" >> res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA2.fasta
						echo "$r2" >> res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA2.fasta
					else
						echo GOT EMPTY LINE FOR FASTA_ID: "$r1"
						echo "$r1" >> res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA_Empty.fasta
						((i=i+1))
					fi

				done
				if [ $(echo res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA_Empty.fasta|wc -l) -gt 0 ]
				then
					echo DONE WITH REMOVING EMPTY LINES GOT $(echo res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA_Empty.fasta|wc -l) CE_INCLUSION HAVING STOP CODON AT THE BEGINNING, PLEASE SEE res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA_Empty.fasta file
					echo DONE WITH REMOVING EMPTY LINES GOT $(echo res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA_Empty.fasta|wc -l) CE_INCLUSION HAVING STOP CODON AT THE BEGINNING, PLEASE SEE res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA_Empty.fasta file >> res_ce_"$2"/Summary_stats.txt
				else
					cat res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA1.fasta > res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA2.fasta
				fi

				##########ALLFRAMES FASTA FILE
				echo NOW REMOVING EMPTY LINES LEFT BY EVENTS WITH STOP CODONS IN THE BEGINNING ""ARHGAP22_chr10"" is one EXAMPLE >> res_ce_"$2"/Summary_stats.txt
				[ -e res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA2.fasta ] && rm res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA2.fasta
				[ -e res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA_Empty.fasta ] && rm res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA_Empty.fasta

				i=0
				cat res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA1.fasta | while read r1; read r2
				do
					if [ "$r2" != "" ]; then
						echo "$r1" >> res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA2.fasta
						echo "$r2" >> res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA2.fasta
					else
						echo "$r1" >> res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA_Empty.fasta
						((i=i+1))
					fi

				done
				if [ $(echo res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA_Empty.fasta|wc -l) -gt 0 ]
				then
					echo DONE WITH REMOVING EMPTY LINES GOT $(echo res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA_Empty.fasta|wc -l) CE_INCLUSION HAVING STOP CODON AT THE BEGINNING, PLEASE SEE res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA_Empty.fasta file
					echo DONE WITH REMOVING EMPTY LINES GOT $(echo res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA_Empty.fasta|wc -l) CE_INCLUSION HAVING STOP CODON AT THE BEGINNING, PLEASE SEE res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA_Empty.fasta file >> res_ce_"$2"/Summary_stats.txt
				else
					cat res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA1.fasta > res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA2.fasta
				fi


				##########

				truncate_flg=1
				if [ $truncate_flg -eq 1 ]
				then
					echo STARTED REMOVING AA "<" 8
					echo STARTED REMOVING AA "<" 8  >> res_ce_"$2"/Summary_stats.txt

					#now do the 3-frame translation



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
					done < res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA2.fasta > res_ce_"$2"/cds_PEAKS_CE_INCLUSION_FUSED_AA1.fasta

					#also remove trailing X in fasta files
					awk 'NR%2==0{sub(/X$/,"")}1' res_ce_"$2"/cds_PEAKS_CE_INCLUSION_FUSED_AA1.fasta > res_ce_"$2"/cds_PEAKS_CE_INCLUSION_FUSED_AA.fasta
					[ -e res_ce_"$2"/cds_PEAKS_CE_INCLUSION_FUSED_AA1.fasta ] && rm res_ce_"$2"/cds_PEAKS_CE_INCLUSION_FUSED_AA1.fasta
				else
					echo skipping removing AA "<=8" lines from FASTA file
					echo skipping removing AA "<=8" lines from FASTA file >> res_ce_"$2"/Summary_stats.txt
					awk 'NR%2==0{sub(/X$/,"")}1' res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA2.fasta > res_ce_"$2"/cds_PEAKS_CE_INCLUSION_FUSED_AA.fasta
					[ -e res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA2.fasta ] && rm res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA2.fasta

				fi
				#statements
				#########ALLFRAMES
				if [ $truncate_flg -eq 1 ]
				then
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
					done < res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA2.fasta > res_ce_"$2"/AllFrames_PEAKS_CE_INCLUSION_FUSED_AA1.fasta

					#also remove trailing X in fasta files
					awk 'NR%2==0{sub(/X$/,"")}1' res_ce_"$2"/AllFrames_PEAKS_CE_INCLUSION_FUSED_AA1.fasta > res_ce_"$2"/AllFrames_PEAKS_CE_INCLUSION_FUSED_AA.fasta
					[ -e res_ce_"$2"/AllFrames_PEAKS_CE_INCLUSION_FUSED_AA1.fasta ] && rm res_ce_"$2"/AllFrames_PEAKS_CE_INCLUSION_FUSED_AA1.fasta
				else
					echo skipping removing AA "<=8" lines from FASTA file res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA2.fasta
					echo skipping removing AA "<=8" lines from FASTA file res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA2.fasta >> res_ce_"$2"/Summary_stats.txt

					#also remove trailing X in fasta files
					awk 'NR%2==0{sub(/X$/,"")}1' res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA2.fasta > res_ce_"$2"/AllFrames_PEAKS_CE_INCLUSION_FUSED_AA.fasta
					[ -e res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA2.fasta ] && rm res_ce_"$2"/AllFrames_CE_INCLUSION_FUSED_AA2.fasta

				fi
			###################
			else
				echo CAUTION!  No final_aa.fasta FILE WAS GENERATED, PLEASE TAKE A CLOSER LOOK   >> res_ce_"$2"/Summary_stats.txt
				#echo BACK FROM R SESSION, TOTAL CE_INCLUSION EVENTS: $nrecrdst  >> res_ce_"$2"/Summary_stats.txt
			fi




		#Also concatenating ce_extension and IR events for next iteration
		echo ALSO CONCATENATING CE_EXTENSION AND IR EVENTS - IF ANY - into REMAINING_EVENTS.csv file FOR NEXT ITERATION
		echo ALSO CONCATENATING CE_EXTENSION AND IR EVENTS - IF ANY - into REMAINING_EVENTS.csv file FOR NEXT ITERATION >> res_ce_"$2"/Summary_stats.txt

		[ -e res_ce_"$2"/REMAINING_EVENTS.csv ] && rm res_ce_"$2"/REMAINING_EVENTS.csv

		nrecrds1=$(cat res_ce_"$2"/IGV_R_returned_ce_extension.csv | wc -l )

		if [ $nrecrds1 -gt 0 ]
		then
			echo COPYING REMAINING $nrecrds1 ce_extension EVENTS INTO res_ce_"$2"/REMAINING_EVENTS.csv FILE >> res_ce_"$2"/Summary_stats.txt
			readarray -t all_csv < res_ce_"$2"/IGV_R_returned_ce_extension.csv
			i=0
			while [ $i -lt $nrecrds1 ]
			do
				line1=${all_csv[$i]}
				echo $line1 | awk 'BEGIN {FS=","} {print $2,$3,$4,$7,$8,$11}' | awk 'BEGIN {OFS=","} {print $1,$2,$3,$4,$5,$6}' >> res_ce_"$2"/REMAINING_EVENTS.csv
				((i=i+1))
			done
		fi

		nrecrds2=$(cat res_ce_"$2"/IGV_R_returned_IR.csv | wc -l )

		if [ $nrecrds2 -gt 0 ]
		then
			echo COPYING REMAINING $nrecrds2 ir EVENTS INTO res_ce_"$2"/REMAINING_EVENTS.csv FILE >> res_ce_"$2"/Summary_stats.txt
			readarray -t all_csv < res_ce_"$2"/IGV_R_returned_IR.csv
			i=0
			while [ $i -lt $nrecrds2 ]
			do
				line1=${all_csv[$i]}
				echo $line1 | awk 'BEGIN {FS=","} {print $2,$3,$4,$7,$8,$11}' |awk 'BEGIN {OFS=","} {print $1,$2,$3,$4,$5,$6}' >> res_ce_"$2"/REMAINING_EVENTS.csv
				((i=i+1))
			done
		fi

		#sort -t$',' -k5 res_ce_"$2"/REMAINING_EVENTS.csv > REMAINING_EVENTS.csv

		remaining_events=$(($nrecrds1+$nrecrds2))
		processed=$(cat res_ce_"$2"/cds_IGV_unique_ce_inclusion.csv|wc -l)
		next_count=$(($processed+$COUNT_EVENT))
		if [  $remaining_events -gt 0 ]
		then
			echo "{{{{{{{{{{{{{  IMPORTANT - WE STILL HAVE " $remaining_events events that are either ce_extension/IR "}}}}}}}}}}}}}"
			echo "{{{{{{{{{{{{{  IMPORTANT - WE STILL HAVE " $remaining_events events that are either ce_extension/IR "}}}}}}}}}}}}}" >> res_ce_"$2"/Summary_stats.txt
			echo "Current total events processed are: " $processed, PLEASE START NEW RUN WITH THIS NUMBER FOR PROPER FASTAID IN AA FASTA FILE FOR "***** PORPER BACKMAPPING PEAKS SEARCH RESULTS"
			echo "Current total events processed are: " $processed, PLEASE START NEW RUN WITH THIS NUMBER FOR PROPER FASTAID IN AA FASTA FILE FOR "***** PORPER BACKMAPPING PEAKS SEARCH RESULTS" >> res_ce_"$2"/Summary_stats.txt
			echo TO CLASSIFY THESE EVENTS, PLEASE COPY FILE res_ce_"$2"/REMAINING_EVENTS.csv in current folder and RE-RUN THIS SCRIPT WITH FOLLOWING COMMONAD LINE "("by adjusting coverage_cutoff value")":
			echo TO CLASSIFY THESE EVENTS, PLEASE COPY FILE res_ce_"$2"/REMAINING_EVENTS.csv in current folder and RE-RUN THIS SCRIPT WITH FOLLOWING COMMONAD LINE "("by adjusting coverage_cutoff value")": >> res_ce_"$2"/Summary_stats.txt
			echo bash pgp_b.sh REMAINING_EVENTS.csv coverage_cutoff $processed+1
			echo bash pgp_b.sh REMAINING_EVENTS.csv coverage_cutoff $processed+1 >> res_ce_"$2"/Summary_stats.txt


		else
			echo DONE PROCESSING ALL CE_INCLUSION EVENTS
			echo DONE PROCESSING ALL CE_INCLUSION EVENTS >> res_ce_"$2"/Summary_stats.txt

			echo "!!!!!!!!!!" PLEASE NOTE THAT THIS WAS FINAL ITERATION OF CE_INCLUSION, PLEASE RUN CE_EXTENSION script with starting NUMBER OF EVENTS SET TO :$next_count FOR PROPER BACKMAPPING "!!!!!!!!!!!"
			echo "!!!!!!!!!!" PLEASE NOTE THAT THIS WAS FINAL ITERATION OF CE_INCLUSION, PLEASE RUN CE_EXTENSION script with starting NUMBER OF EVENTS SET TO :$next_count FOR PROPER BACKMAPPING "!!!!!!!!!!!" >> res_ce_"$2"/Summary_stats.txt


		fi
		echo FOR FINAL STATISTICS ON CE EVENTS, PLEASE SEE res_ce_"$2"/Summary_stats.txt file

	else




	[ -e res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA1.fasta ] && rm res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA1.fasta
	[ -e res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA2.fasta ] && rm res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA2.fasta
	[ -e res_ce_"$2"/cds_PEAKS_CE_INCLUSION_FUSED_AA1.fasta ] && rm res_ce_"$2"/cds_PEAKS_CE_INCLUSION_FUSED_AA1.fasta
	[ -e res_ce_"$2"/cds_PEAKS_CE_INCLUSION_FUSED_AA.fasta ] && rm res_ce_"$2"/cds_PEAKS_CE_INCLUSION_FUSED_AA.fasta
	[ -e res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA_Empty.fasta ] && rm res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA_Empty.fasta
	echo check_aaV4_allFrames.R DID NOT RETURN ANY final_aa.fasta file, so exiting
	echo check_aaV4_allFrames.R DID NOT RETURN ANY final_aa.fasta file, so exiting >> res_ce_"$2"/Summary_stats.txt
	echo FOR FINAL STATISTICS ON CE EVENTS, PLEASE SEE res_ce_"$2"/Summary_stats.txt file

	fi

	echo Finally run sashimi plots
	if [ -f res_ce_"$2"/ce_inclusion_coord_uniq_sashimi.bed ]
	then
		num_ce_incl=$(cat res_ce_"$2"/ce_inclusion_coord_uniq_sashimi.bed | wc -l | awk '{print $1/3}')
		echo STARTING SASHIMI PLOTS FOR $num_ce_incl Events, For events ">" 100 It Might take several hours
		#source run_sashimi.sh 3
		source run_sashimiV1.sh res_ce_"$2"/IGV_unique_ce_inclusion.csv res_ce_"$2"/ce_inclusion_coord_uniq_sashimi.bed 3
	else
		echo File res_ce_"$2"/ce_inclusion_coord_uniq_sashimi.bed does not exists so, NO CE_INCLUSION SASHIMI PLOTS ARE GENERATED
	fi

	echo ALL DONE - hopefully - Successfully

	echo FINAL CLEANUP - Removing intermediate/unnecessary files
	rm res_skiptics/*sashimi*
	rm res_ce_all/*sashimi*
	rm *.bed
	rm all_events_bed_sashimi.tab
	rm dump.fasta
	rm res_ce_"$2"/ce_all_scan_intron.bed
	rm res_ce_"$2"/ce_all_scan_range_junctions.bed
	rm res_ce_"$2"/ce_all_scan_range.bed
	rm res_ce_"$2"/ce_all_scan_unique_range.bed
	rm res_ce_"$2"/ce_extension_coord_only.bed
	rm res_ce_"$2"/ce_extension_coord_uniq.bed
	rm res_ce_"$2"/ce_extension_coord.bed
	rm res_ce_"$2"/ce_inclusion_coord_only.bed
	rm res_ce_"$2"/ce_inclusion_coord_repeated.bed
	rm res_ce_"$2"/ce_inclusion_coord_uniq.bed
	rm res_ce_"$2"/ce_inclusion_coord.bed
	rm res_ce_"$2"/IR_coord_only.bed
	rm res_ce_"$2"/IR_coord_repeated.bed
	rm res_ce_"$2"/IR_coord_uniq.bed
	rm res_ce_"$2"/IR_coord.bed
	rm res_ce_"$2"/IGV_R_returned_ce_extension.csv
	rm res_ce_"$2"/IGV_R_returned_ce_inclusion.csv
	rm res_ce_"$2"/IGV_R_returned_IR.csv
	rm res_ce_"$2"/IGV_ce_inclusion.csv
	#rm res_ce_"$2"/IGV_coverage_file_notfound_ce.csv
	rm res_ce_"$2"/IGV_unique_ce_inclusion.csv
	rm res_ce_"$2"/cds_IGV_unique_ce_inclusion.csv
	rm res_ce_"$2"/cds_successful_frames_list.csv
	rm res_ce_"$2"/cds_unsuccessful_frames_list.csv
	rm res_ce_"$2"/aa.fasta
	rm res_ce_"$2"/cds_ce_inclusion_coord_uniq.bed
	rm res_ce_"$2"/cds_CE_INCLUSION_FUSED_AA.fasta
	rm res_ce_"$2"/ce_inclusion_coord_only.txt
	rm res_ce_"$2"/ce_inclusion_coord_sashimi.bed
	rm res_ce_"$2"/ce_inclusion_coord_uniq_sashimi.bed
	rm res_ce_"$2"/ce_inclusion_transeq_in.fasta
	rm res_ce_"$2"/IR_coord_sashimi.bed
	rm res_ce_"$2"/non_protein_coding.bed
	rm res_ce_"$2"/protein_coding_a.bed
	rm res_ce_"$2"/protein_coding.bed


	#remove all bed files now

	#rm *.bed
	#fi #ce_extension block
fi # opening if
