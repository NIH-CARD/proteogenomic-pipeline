#!/bin/bash
#HOW TO RUN: pgp_b_ce_ir.sh IR_pgp1.csv
#NEXT SORT FILE ON COLUMN 5 - THIS IS IMPORTANT AS PIPELINE DEPENDS ON IT


sort -t$',' -k5 $1 > sorted_$1

mkdir -p res_IR


#Step 2 - STARTING ce calculations
echo "############### NOW STARTING CE IDENTIFICATION #########" > res_IR/Summary_stats.txt

echo NOW STARTING CE IDENTIFICATION
	#inpfile=sorted_$1
	inpfile=sorted_$1 #all_non_skiptics.csv
#tx_lst=$2
#flag to run this section
cryptics_flg=1
if [ $cryptics_flg -eq 1 ]
then

	#mkdir -p res_intronic_range
	mkdir -p event_bedfiles


	[ -e events_to_tx_mapping_valid.csv ] && rm events_to_tx_mapping_valid.csv
	[ -e events_tx_mapping_invalid.csv ] && rm events_tx_mapping_invalid.csv

	########IR SECTION

	[ -e res_IR/IR_coord.bed ] && rm res_IR/IR_coord.bed
	[ -e res_IR/IR_coord_sashimi.bed ] && rm res_IR/IR_coord_sashimi.bed
	[ -e res_IR/IGV_R_returned_IR.csv ] && rm res_IR/IGV_R_returned_IR.csv
	[ -e res_IR/IR_coord_uniq.bed ] && rm res_IR/IR_coord_uniq.bed
	[ -e res_IR/IR_coord_uniq_sashimi.bed ] && rm res_IR/IR_coord_uniq_sashimi.bed
	[ -e res_IR/IR_coord_only.bed ] && rm res_IR/IR_coord_only.bed
	[ -e res_IR/IGV_unique_IR.csv ] && rm res_IR/IGV_unique_IR.csv

	[ -e res_IR/IGV_repeated_IR.csv ] && rm res_IR/IGV_repeated_IR.csv



	######END IR SECTION

	#[ "$(ls -A event_bedfiles/)" ] && rm event_bedfiles/*.*

	[ -e res_IR/unique_IR.csv ] && rm res_IR/unique_IR.csv

	[ -e res_IR/unique_ce_extension.csv ] && rm res_IR/unique_ce_extension.csv

	[ -e res_IR/unique_ce_inclusion.csv ] && rm res_IR/unique_ce_inclusion.csv

	[ -e res_IR/ce_all_scan_intron.bed ] && rm res_IR/ce_all_scan_intron.bed

	[ -e res_IR/ce_inclusion_coord_only.bed ] && rm res_IR/ce_inclusion_coord_only.bed
	[ -e res_IR/ce_extension_coord_only.bed ] && rm res_IR/ce_extension_coord_only.bed
	[ -e res_IR/IGV_R_returned_ce_extension.csv ] && rm res_IR/IGV_R_returned_ce_extension.csv
	[ -e res_IR/IGV_R_returned_ce_inclusion.csv ] && rm res_IR/IGV_R_returned_ce_inclusion.csv
	[ -e res_IR/IGV_skiptics.csv ] && rm res_IR/IGV_skiptics.csv
	[ -e res_IR/IGV_unique_ce_extension.csv ] && rm res_IR/IGV_unique_ce_extension.csv
	[ -e res_IR/IGV_unique_ce_inclusion.csv ] && rm res_IR/IGV_unique_ce_inclusion.csv
	[ -e res_IR/non_ce_events.csv ] && rm res_IR/non_ce_events.csv

	[ -e res_IR/ce_extension_coord_sashimi.bed ] && rm res_IR/ce_extension_coord_sashimi.bed
	[ -e res_IR/ce_extension_coord_uniq_sashimi.bed ] && rm res_IR/ce_extension_coord_uniq_sashimi.bed
	[ -e res_IR/ce_extension_coord.bed ] && rm res_IR/ce_extension_coord.bed
	[ -e res_IR/ce_inclusion_coord.bed ] && rm res_IR/ce_inclusion_coord.bed

	[ -e res_IR/ce_all_scan_unique_range.bed ] && rm res_IR/ce_all_scan_unique_range.bed
	[ -e res_IR/ce_inclusion_coord_uniq_sashimi.bed ] && rm res_IR/ce_inclusion_coord_uniq_sashimi.bed
	[ -e res_IR/ce_extension_coord_uniq_sashimi.bed ] && rm res_IR/ce_extension_coord_uniq_sashimi.bed
	[ -e res_IR/IGV_skiptics.csv ] && rm res_IR/IGV_skiptics.csv
	[ -e res_IR/IGV_unique_ce_inclusion.csv ] && rm res_IR/IGV_unique_ce_inclusion.csv
	[ -e res_IR/IGV_unique_ce_extension.csv ] && rm res_IR/IGV_unique_ce_extension.csv
	[ -e res_IR/ex1_ex2_ce.txt ] && rm res_IR/ex1_ex2_ce.txt

	[ -e res_IR/ce_all_scan_range_junctions.bed ] && rm res_IR/ce_all_scan_range_junctions.bed
	[ -e res_IR/ce_all_scan_range.bed ] && rm res_IR/ce_all_scan_range.bed


	[ -e res_IR/ce_extension_coord_uniq.bed ] && rm res_IR/ce_extension_coord_uniq.bed
	[ -e res_IR/ce_extension_coord_repeated.bed ] && rm res_IR/ce_extension_coord_repeated.bed


	[ -e res_IR/Summary_stats.txt ] && rm res_IR/Summary_stats.txt

	[ -e res_IR/IGV_ce_inclusion.csv ] && rm res_IR/IGV_ce_inclusion.csv


	[ -e res_IR/ce_inclusion_coord_uniq.bed ] && rm res_IR/ce_inclusion_coord_uniq.bed
	[ -e res_IR/ce_inclusion_coord_repeated.bed ] && rm res_IR/ce_inclusion_coord_repeated.bed

	[ -e res_IR/non_ce_events.txt ] && rm res_IR/non_ce_events.txt
	[ -e res_IR/non_ce_events.csv ] && rm res_IR/non_ce_events.csv

	[ -e res_IR/IGV_problematic_junctions.csv ] && rm res_IR/IGV_problematic_junctions.csv
	[ -e res_IR/problematic_junctions.txt ] && rm res_IR/problematic_junctions.txt


	events_bed_create_flg=1
	if [ $events_bed_create_flg -eq 1 ]
	then
		[ "$(ls -A event_bedfiles/)" ] && rm event_bedfiles/*.*
		[ -e res_IR/EnsDB_tx_not_found.csv ] && rm res_IR/EnsDB_tx_not_found.csv
		Rscript TxEnsDB103_layeredV6.R $inpfile principal_txs.csv res_IR/Summary_stats.txt
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
		#echo processing "$sample"
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
	#cat "$sample"

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
			echo $gene_name is a CE event >> res_IR/ex1_ex2_ce.txt  #>> res_intronic_range/intronic_range_events.txt
			echo ds $ds >> res_IR/ex1_ex2_ce.txt
			echo us $us >> res_IR/ex1_ex2_ce.txt

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
		echo now dsn $dsn >> res_IR/ex1_ex2_ce.txt
			echo now usn $usn >> res_IR/ex1_ex2_ce.txt
		echo dsodiff $dsodiff >> res_IR/ex1_ex2_ce.txt
			echo up exon length $usnexonl and dsnexonl $dsnexonl

			if [ $usnexonl == "." ] || [ $dsnexonl == "." ]
			then
				#echo $gene_name is ce with up exon length $usnexonl and dsnexonl $dsnexonl , LEAVING it for now PLEASE MAKE SURE IT IS PROPERLY IDENTIFIED >> res_intronic_range/problematic_junctions.txt
				echo $gene_name is ce with up exon length $usnexonl and dsnexonl $dsnexonl , LEAVING it for now PLEASE MAKE SURE IT IS PROPERLY IDENTIFIED >> res_IR/problematic_junctions.txt
				#cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_intronic_range/IGV_problematic_junctions.csv
				cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_IR/IGV_problematic_junctions.csv

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

		gene_id=$(echo "$line_csv"| awk 'BEGIN {FS=","} {print $6}')
		o_gene_id=$(echo "$line_csv"| awk 'BEGIN {FS=","} {print $5}')
		cat $sample | awk -v gene=$gene_name -v TX=$TxID -v og=$o_gene_id -v gid=$gene_id 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene,TX,og,gid}' >> res_IR/IGV_ce_inclusion.csv
	#echo for sample
	#cat "$sample"
	#echo finally upex is "$upex"
	#echo finally ce is "$ce"
	#echo finally dsex is "$dsex"

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
		echo "$upex">>res_IR/ce_all_scan_range_junctions.bed
		echo "$ce_scan1">>res_IR/ce_all_scan_range_junctions.bed
		echo "$dsex">>res_IR/ce_all_scan_range_junctions.bed

		echo "$ce_scan1">>res_IR/ce_all_scan_range.bed #all_junctions.bed


	fi

	elif [ "$diff_exon_abs" -eq 1 ] && [ $dsovlp -le 5 ] && [ $usovlp -le 5 ] #deals with exon_joining events few bp of
	then
		cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_IR/IGV_skiptics.csv
		#cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_intronic_range/non_ce_events.csv

		echo $line_csv | awk 'BEGIN {FS=","} {print $1,$2,$3,$4,$5,$6,$7}' | awk 'BEGIN {OFS=","} {print $1,$2,$3,$4,$5,$6,$7}' >> res_IR/non_ce_events.csv

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

		gene_id=$(echo "$line_csv"| awk 'BEGIN {FS=","} {print $6}')
		o_gene_id=$(echo "$line_csv"| awk 'BEGIN {FS=","} {print $5}')

		cat $sample | awk -v gene=$gene_name -v TX=$TxID -v og=$o_gene_id -v gid=$gene_id 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene,TX,og,gid}' >> res_IR/IGV_ce_inclusion.csv
	#echo for sample
	#cat $sample
	#echo finally upex is "$upex"
	#echo finally ce is "$ce"
	#echo finally dsex is "$dsex"

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
		echo "$upex">>res_IR/ce_all_scan_range_junctions.bed
		echo "$ce_scan1">>res_IR/ce_all_scan_range_junctions.bed
		echo "$dsex">>res_IR/ce_all_scan_range_junctions.bed

		echo "$ce_scan1">>res_IR/ce_all_scan_range.bed #all_junctions.bed
	else
		#BOTH ENDS OF THESE EVENTS LIE INSIDE INTRON
		echo $gene_name has both ends lying inside intron #>> res_intronic_range/intronic_range_events.txt #>>ce_intronic.txt
		#echo $gene_name is unknown event type for gene "$gene_name" >> res_intronic_range/non_ce_events.txt
		echo $gene_name is unknown event type for gene "$gene_name" >> res_IR/non_ce_events.txt
		#cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_intronic_range/non_ce_events.csv
		echo $line_csv | awk 'BEGIN {FS=","} {print $1,$2,$3,$4,$5,$6,$7}' | awk 'BEGIN {OFS=","} {print $1,$2,$3,$4,$5,$6,$7}' >> res_IR/non_ce_events.csv

	fi
	#fi
	else
	#echo $gene_name is unknown event type for gene "$gene_name" >> res_intronic_range/non_ce_events.txt
		echo $gene_name is unknown event type for gene "$gene_name" >> res_IR/non_ce_events.txt
		#cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_intronic_range/non_ce_events.csv
		echo $line_csv | awk 'BEGIN {FS=","} {print $1,$2,$3,$4,$5,$6,$7}' | awk 'BEGIN {OFS=","} {print $1,$2,$3,$4,$5,$6,$7}' >> res_IR/non_ce_events.csv
	fi

	#rm $sample1
	#rm $sample_nt.bed
	#rm $sample1-3UTR.bed
	#rm $sample1-5UTR.bed
	#rm $allexons
	#rm t$allexons
	done
	########## THIS SECTION WRITES NEW FILE res_IR/ce_all_scan_intron.bed to scan whole intron between two exons for probable ce events
	nrecrds=$(cat res_IR/ce_all_scan_range_junctions.bed | wc -l)
	readarray -t all_data < res_IR/ce_all_scan_range_junctions.bed
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

				echo "$line11" "$line13" | awk -v o=$off 'BEGIN {OFS="\t"} {print $1,$3,$9,$11,$12,$13,$14}' >> res_IR/ce_all_scan_intron.bed
						elif [ "$strand" == "-" ]
						then
							#st=$(echo $line13 | awk 'BEGIN {FS="\t"} {print $2}')
							#end=$(echo $line13 | awk 'BEGIN {FS="\t"} {print $3}')
							#off=$(($end-$st))
							#st1=$(echo $line11 | awk 'BEGIN {FS="\t"} {print $2}')
							#end1=$(echo $line11 | awk 'BEGIN {FS="\t"} {print $3}')
							#off1=$(($end1-$st1))

							echo "$line13" "$line11" | awk 'BEGIN {OFS="\t"} {print $1,$3,$9,$11,$12,$13,$14}' >> res_IR/ce_all_scan_intron.bed
						fi
	done

	################## END SECTION WRITES NEW FILE res_IR/ce_all_scan_intron.bed
	#also copy for coverage calculations and avoid intronic_range calculations
	cat res_IR/ce_all_scan_range.bed  > res_IR/ce_all_scan_unique_range.bed
fi


#CRYPTICS SECTION ENDS HERE

#INPUT: IR ONLY
#GET UP/DN stream exons and Translate intervening intron
IR_flg=1
if [ "$IR_flg" -eq 1 ]
then
	[ -e res_IR/IR_coord.bed ] && rm res_IR/IR_coord.bed
	nrecrds=0
	#nrecrdst=0
	[ -e res_IR/ce_all_scan_range_junctions.bed ] && nrecrds=$(cat res_IR/ce_all_scan_range_junctions.bed | wc -l) && nrecrdst=$(($nrecrds/3))
	echo GOT $nrecrdst IR events
	if [ $nrecrds -ne 0 ]
	then

#GET up/dn exon coordinates
	readarray -t all_data < res_IR/ce_all_scan_range_junctions.bed
	i=0

	while [ $i -lt $nrecrds ]
	  do
	          line11=${all_data[$i]}
	          line12=${all_data[$(($i+1))]}
	          line13=${all_data[$(($i+2))]}



	     #echo "Welcome $i times"
	     ((i=i+3))
			 strand=$(echo $line11 | awk '{print $6}')
			 if [ "$strand" == "+" ]
			 then
				 echo $line11 $line13| awk 'BEGIN {OFS="\t"} {print $1,$3,$9,$11,$12,$13,$14}' >> res_IR/IR_coord.bed
			 else
				 echo $line13 $line11| awk 'BEGIN {OFS="\t"} {print $1,$3,$9,$11,$12,$13,$14}' >> res_IR/IR_coord.bed
			 fi
	 done
	 nrecrds_ir=0

	 [ -e res_IR/IR_coord.bed ] && nrecrds_ir=$(cat res_IR/IR_coord.bed | wc -l) #&& nrecrds_irt=$(($nrecrds_ir/3))

	 if [ $nrecrds_ir -ne 0 ]
	 then
	 echo TOTAL IR EVENTS ARE: $nrecrds_ir, NOW CHECKING FOR REPEATED IR EVENTS >> res_IR/Summary_stats.txt
	 echo TOTAL IR EVENTS ARE: $nrecrds_ir
	 #echo checking for repeated IR events >> res_IR/Summary_stats.txt
	 readarray -t all_data < res_IR/IR_coord.bed
	 #also read ce_inclusion_coord_sashimi.bed to get TxID
	 #readarray -t sashimi_data_ext < res_IR/IR_coord_sashimi.bed

	 #also read IGV_ce_inclusion.csv to get TxID - It is safe as we consider all events as ce_extension here
	 readarray -t get_txid < res_IR/IGV_ce_inclusion.csv

	 #also read csv file
	 readarray -t csv_data_ir < sorted_$1


	 #csvi=0

	 #echo total IR events read are $nrecrds_irt
	 i=0
	csvi=0
	 while [ $i -lt $nrecrds_ir ]
	   do
	          line11=${all_data[$i]}
						csv_ln=${csv_data_ir[$i]}

						txid_rec=${get_txid[$csvi]}
						txid=$(echo $txid_rec | awk 'BEGIN {FS=","} {print $9}')
						((csvi=csvi+1))


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
	                 echo $line11 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_IR/IR_coord_uniq.bed
									 #Also for sashimi plots
									 echo $line11 | awk -v tx=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,1,$5,$6,$7,tx}' >> res_IR/IR_coord_uniq_sashimi.bed

	 								 #Also save ce coordinates only
	 								 echo $line11 | awk  -v tx=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$3-$2,$5,$6,$7,tx}' >> res_IR/IR_coord_only.bed
									 #echo $csv_ln | awk 'BEGIN {FS=","} {print $1,$2,$3,$4,$5,$6}'| awk 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6}' >> res_IR/IGV_unique_IR.csv
									 echo $csv_ln | awk 'BEGIN {FS=","} {print $1,$2,$3,$4,$5,$6}'| awk 'BEGIN {OFS=","} {print $1,$2,$3,$4,$5,$6}' >> res_IR/IGV_unique_IR.csv
	 				else
	 					echo $line11 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_IR/IR_coord_repeated.bed
						echo $csv_ln | awk 'BEGIN {FS=","} {print $1,$2,$3,$4,$5,$6}'| awk 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6}' >> res_IR/IGV_repeated_IR.csv
	          fi
	  done
	   #get total unique coordinates
	  nrecrds_uniq_irt=$(cat res_IR/IR_coord_uniq.bed | wc -l)
	  #nrecrds_uniq_ir=$(($nrecrds_uniq_irt/3))
	  repeated_IR_events=$(($nrecrds_ir-$nrecrds_uniq_irt))

	  echo TOTAL UNIQUE IR EVENTS ARE: $nrecrds_uniq_irt, please see res_IR/IR_coord_uniq.bed and res_IR/IGV_unique_IR.csv files >> res_IR/Summary_stats.txt
	  #echo total unique IR events read are $nrecrds_uniq_ir
	  echo Total repeated IR events are: $repeated_IR_events, please see res_IR/IR_coord_repeated.bed >> res_IR/Summary_stats.txt
	  echo THOSE WERE ALL PERTINENT STAT - Please let us know if something is missing or some more stats can be useful!!! >> res_IR/Summary_stats.txt
	  echo Total repeated IR events were $repeated_IR_events
	  #now copy back data to ce_coord
	  #cat res_IR/ce_extension_coord_uniq.bed	> res_IR/ce_extension_coord.bed
	 #NEW CODE ENDS HERE - IR


	 echo now finding nt seq for res_IR/IR_coord_uniq.bed
	 #now call bedtools getfasta function to get nt sequence from reference_genome.fa file
	 #bedtools getfasta -fi /Volumes/SYEDSHAH/MichaelLab/ref_genome/GRCh38.p13.genome.fa -bed res_IR/IR_coord_uniq.bed -s > res_IR/IR_coord_uniq_nt.fa
	 bedtools getfasta -fi GRCh38.p13.genome.fa -bed res_IR/IR_coord_uniq.bed -s > res_IR/IR_coord_uniq_nt.fa

	 #now remove >chr lines from the resulting file
	 awk '!/^>chr/' res_IR/IR_coord_uniq_nt.fa > res_IR/IR_coord_uniq_nt1.fa


	 paste -d"\t" res_IR/IR_coord_uniq.bed res_IR/IR_coord_uniq_nt1.fa > res_IR/temp.IR_coord_uniq_nt.bed

	 rm res_IR/IR_coord_uniq_nt1.fa
	 #rm res_IR/IR_coord_uniq_nt.fa
	 #awk -F "\t" '{if($6=="+") print ">sp|"$7"_"$1"_"$2"-"$3"_plus""\n"$8;else print ">sp|"$7"_"$1"_"$2"-"$3"_minus""\n"$8}' res_IR/temp.IR_coord_uniq_nt.bed > res_IR/temp.IR_coord_uniq_nt.transeq_in.fasta
	 i=0
	 awk -F "\t" '{if($6=="+") print ">sp|"$7"_"$1"_"$2"_"$3"_"((i=i+1))"_plus""\n"$8;else print ">sp|"$7"_"$1"_"$2"_"$3"_"((i=i+1))"_minus""\n"$8}' res_IR/temp.IR_coord_uniq_nt.bed > res_IR/IR_coord_uniq_nt.transeq_in.fasta
	 rm res_IR/temp.IR_coord_uniq_nt.bed
	 #also save nt sequence

	 TEMPFILE=$(echo res_IR/IR_coord_uniq_nt.transeq_in.fasta | perl -pe 's/\.fasta$/.temp/')
	 OUTPUTFILE=$(echo res_IR/IR_coord_uniq_nt.transeq_in.fasta | perl -pe 's/\.fasta$/.trans/')

	 # Translate sequence

	 transeq  -sequence res_IR/IR_coord_uniq_nt.transeq_in.fasta -outseq "$TEMPFILE" -frame F


	 # Rename sequence
	 perl -sape 's/>/>sp|/' "$TEMPFILE" > "$OUTPUTFILE"


	 #no concatenate all files
	 cat "$OUTPUTFILE" >> res_IR/temp.IR_AA.fasta
	 # Remove temp file
	 rm "$TEMPFILE"
	 rm "$OUTPUTFILE"

	 #rm res_IR/temp.IR_coord_uniq_nt.transeq_in.fasta

	 #also remove all newlines from the
	 awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' res_IR/temp.IR_AA.fasta > res_IR/IR_AA.fasta

	 #remove
	 rm res_IR/temp.IR_AA.fasta

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
	 done < res_IR/IR_AA.fasta > res_IR/FINAL_IR_AA1.fasta

	 #also remove trailing X from res_IR/FINAL_IR_AA1.fasta file
	 awk 'NR%2==0{sub(/X$/,"")}1' res_IR/FINAL_IR_AA1.fasta > res_IR/FINAL_IR_AA.fasta
	 [ -e res_IR/FINAL_IR_AA1.fasta ] && rm res_IR/FINAL_IR_AA1.fasta

	 else
	 	echo TOTAL IR EVENTS ARE: $nrecrds_irt >> res_IR/Summary_stats.txt
	 fi

	 ######END IR EVENTS
	 #FINAL STATISTICS

	 echo For final statistics on ce events, please see res_IR/Summary_stats.txt file

	 echo NOW sashimi plots
	 if [ -f res_IR/IR_coord_uniq_sashimi.bed ]
   then
     num_IR=$(cat res_IR/IR_coord_uniq_sashimi.bed | wc -l)
     echo STARTING SASHIMI PLOTS FOR $num_IR Events, For events ">" 100 It Might take several hours Might take several hours
     source run_sashimiV1.sh res_IR/IGV_unique_IR.csv res_IR/IR_coord_uniq_sashimi.bed 5
   else
     echo File res_IR/IR_coord_uniq_sashimi.bed does not exists so, NO CE_EXTENSION SASHIMI PLOTS ARE GENERATED

   fi
   echo CLEANING UP >> res_IR/Summary_stats.txt
   
	rm res_IR/IGV_ce_inclusion.csv
	mv res_IR/IGV_unique_IR.csv res_IR/IR.csv
	rm res_IR/IR_AA.fasta


	#rm res_IR/IR_coord_only.bed
	rm res_IR/IR_coord_uniq.bed
	rm res_IR/IR_coord_uniq_nt.fa
	#IR_coord_uniq_nt.transeq_in.fasta
	rm res_IR/IR_coord_uniq_sashimi.bed
	
	rm res_IR/ce_all_scan_intron.bed
	rm res_IR/ce_all_scan_range.bed
	rm res_IR/ce_all_scan_range_junctions.bed
	rm res_IR/ce_all_scan_unique_range.bed
	rm res_IR/IR_coord.bed

	 echo ALL DONE - hopefully - Successfully
else
	echo file res_IR/ce_all_scan_range_junctions.bed has zero records so EXITING
fi


fi

#DONE IR ONLY


#remove all bed files now

#rm *.bed
