#THIS VERSION NOW ADDS EVENT_ID FLAG TO FINAL AA FILE FOR FUSED EVENTS TO MAP BACK SASHIMI PLOTS FOR THE PEAKS EVENTS

#STARTED USING TxEnsDB103_layeredV2.R which takes care of the Txs seleceted not encapsulating the event - 01-13-2022
#Adding up/dn stream exon length for sashimi plots
#above change INTERFERE with removal of repeated entries, so CREATING NEW bed file (skiptics_uniq_sashimi.bed) that conatin exon lengths
#STARTING NEW CODE TO TIGHTLY ONLY EXON_SKIP EVENTS - Dec2, 2021
#ONLY EXON_SKIP EVENTS - USE

mkdir -p res_skiptics
mkdir -p event_bedfiles

#[ -e event_bedfiles/*.bed ] &&
#rm event_bedfiles/*.bed
[ "$(ls -A event_bedfiles/)" ] && rm event_bedfiles/*.*
[ "$(ls -A es_event_bedfiles/)" ] && rm es_event_bedfiles/*.*


[ -e events_to_tx_mapping_valid.csv ] && rm events_to_tx_mapping_valid.csv
[ -e events_tx_mapping_invalid.csv ] && rm events_tx_mapping_invalid.csv

echo CLEARING UP res_skiptics folder from previous run if any. SASHIMI FOLDERS WILL ONLY BE DELETED BEFORE CREATING NEW SASHIMI PLOTS.

[ -e res_skiptics/EnsDB_tx_not_found.csv ] && rm res_skiptics/EnsDB_tx_not_found.csv
[ -e res_skiptics/skiptics_sashimi.bed ] && rm res_skiptics/skiptics_sashimi.bed
[ -e res_skiptics/skiptics_uniq_sashimi.bed ] && rm res_skiptics/skiptics_uniq_sashimi.bed

[ -e res_skiptics/all_processed.txt ] && rm res_skiptics/all_processed.txt
[ -e res_skiptics/skiptics_verbose.txt ] && rm res_skiptics/skiptics_verbose.txt
[ -e res_skiptics/IGV_skiptics.csv ] && rm res_skiptics/IGV_skiptics.csv

[ -e res_skiptics/IGV_skiptics.csv ] && rm res_skiptics/IGV_skiptics.csv
[ -e res_skiptics/IGV_all_others.csv ] && rm res_skiptics/IGV_all_others.csv

[ -e all_non_skiptics.csv ] && rm all_non_skiptics.csv

[ -e res_skiptics/IGV_unique_skiptics_translated.csv ] && rm res_skiptics/IGV_unique_skiptics_translated.csv
[ -e res_skiptics/skiptics_all_coordinates.bed ] && rm res_skiptics/skiptics_all_coordinates.bed
[ -e res_skiptics/skiptics_unique.bed ] && rm res_skiptics/skiptics_unique.bed
[ -e res_skiptics/skiptics_repeated.bed ] && rm res_skiptics/skiptics_repeated.bed
[ -e res_skiptics/skiptics_fused_transeq_in.fasta ] && rm res_skiptics/skiptics_fused_transeq_in.fasta
[ -e res_skiptics/IGV_multi_skiptics.csv ] && rm res_skiptics/IGV_multi_skiptics.csv
[ -e res_skiptics/IGV_skiptics_repeated.csv ] && rm res_skiptics/IGV_skiptics_repeated.csv

[ -e res_skiptics/SKIPTICS_FINAL_STATS.txt ] && rm res_skiptics/SKIPTICS_FINAL_STATS.txt
#[ -e still_all_not_found.csv ] && mv still_all_not_foun.csv res_skiptics/.
#[ -e res/output.txt ] && rm res/output.txt
inpfile=$1 #RBP_all.csv
txfile=$2
#tx_lst=$2
#cat all_majiq_events.csv | awk 'BEGIN {FS="\,"}($7 !~ "novel_exon_skip")' > all_except_exon_skip.csv
echo "######################################" >> res_skiptics/SKIPTICS_FINAL_STATS.txt
echo CALLING TxEnsDB103_layeredV5.R FROM SKIPTICS SCRITP TO GENERATE BED FILES FOR EACH EVENT >> res_skiptics/SKIPTICS_FINAL_STATS.txt
Rscript TxEnsDB103_layeredV5.R $inpfile $txfile #$tx_lst
[ -e EnsDB_tx_not_found.csv ] && cp EnsDB_tx_not_found.csv res_skiptics/.
echo BACK FROM TxEnsDB103_layeredV5.R script >> res_skiptics/SKIPTICS_FINAL_STATS.txt
echo "#######################################" >> res_skiptics/SKIPTICS_FINAL_STATS.txt
#Rscript es_Rgeneid_transcriptV1_EnsDB103.R $inpfile #all_majiq_events.csv #exon_skip.csv
#also read csv file to filter out non_es events for ce script
readarray -t csv_data < all_tx_events.csv
csvi=0

flg_bed=1
if [ $flg_bed -eq 1 ]
then
samples=$(ls event_bedfiles/temp_*.bed)

for sample in $samples
do
	#read csv entry
	csv_ln=${csv_data[$csvi]}
	((csvi=csvi+1))

	#echo processing "$sample"
	echo processing "$sample" >> res_skiptics/all_processed.txt
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
	us=$(bedtools closest -a $sample -b event_bedfiles/t$allexons -s  -D a -id -d -t last)

	#echo ds $ds
	#echo us $us

	#get up/dn exon lengths
	upexonl=$(echo "$us" | awk '{print $10}')
	dnexonl=$(echo "$ds" | awk '{print $10}')

  #get up and down stream exon numbers
	upexon=$(echo "$us" | awk '{print $11}')
	dnexon=$(echo "$ds" | awk '{print $11}')

	#get overlap with up/dn exon - 0 means complete overlap which is assumed as exon skip event
	dsovlp=$(echo "$ds" | awk '{print $13}')
	usovlp=$(echo "$us" | awk '{print $13}')


  diff_exon=$(($upexon-$dnexon))
	#take absolute value
	diff_exon_abs=${diff_exon#-}

  if [ "$upexon" == "$dnexon" ]  && [ $upexonl -eq $dnexonl ] #|| [ $diff_exon_abs -gt 1 ]
	then
		#THIS SECTION IS FOR EXON_SKIP events
		dsi=$(bedtools intersect -a $sample -b event_bedfiles/t$allexons -s -wo)
		usi=$(bedtools intersect -a $sample -b event_bedfiles/t$allexons -s -wo)

		dsovlp=$(echo "$dsi" | awk '{print $13}')
		usovlp=$(echo "$usi" | awk '{print $13}')
		#get up/dn exon lengths
		upexonl=$(echo "$dsi" | awk '{print $10}')
		dnexonl=$(echo "$dsi" | awk '{print $10}')


		diff1=$(($dnexonl-$dsovlp))
		diff2=$(($upexonl-$usovlp)) #if single exon is intersected

		if [ $dsovlp -eq $usovlp ] && [ $upexonl -eq $dnexonl ] && [ $diff1 -eq 1 ] && [ $diff2 -eq 1 ]
		then

			#As algorithm returns skipped exon for exon_skip events, so determine coordinates of actual exons involved in the event
			#echo "$us" | awk 'BEGIN {OFS="\t"} {print $1,$2-10,$3,$4,$5,$6}' > dump_us.bed
			echo "$us" | awk 'BEGIN {OFS="\t"} {print $1,$2-10,$2,$4,$5,$6}' > dump_us.bed
			#echo "$ds" | awk 'BEGIN {OFS="\t"} {print $1,$2,$3+10,$4,$5,$6}' > dump_ds.bed
			echo "$ds" | awk 'BEGIN {OFS="\t"} {print $1,$3,$3+10,$4,$5,$6}' > dump_ds.bed

      dsn=$(bedtools closest -a dump_ds.bed -b event_bedfiles/t$allexons -s -D a -iu -d -t last )
      usn=$(bedtools closest -a dump_us.bed -b event_bedfiles/t$allexons -s  -D a -id -d -t first)

      #get up and down stream exon numbers
      upexon=$(echo "$usn" | awk '{print $11}')
      dnexon=$(echo "$dsn" | awk '{print $11}')

      dsovlp=$(echo "$dsn" | awk '{print $13}')
  		usovlp=$(echo "$usn" | awk '{print $13}')


      if [ $upexon -eq -1 ] || [ $dnexon -eq -1 ] || [ $dsovlp -ne 0 ] || [ $usovlp -ne 0 ]
      then

	        cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_skiptics/IGV_all_others.csv
					#also save csv file
					echo $csv_ln >> all_non_skiptics.csv
					echo PROBABLY THIS NEVER REACHES

      else
	      echo $gene_name is an exon_skip event >> res_skiptics/skiptics_verbose.txt
				echo dsi is $dsi >> res_skiptics/skiptics_verbose.txt
				echo usi is $usi >> res_skiptics/skiptics_verbose.txt

				echo us exon is "$us" >> res_skiptics/skiptics_verbose.txt
				echo ds exon is "$ds" >> res_skiptics/skiptics_verbose.txt
				echo dump_us exon is >> res_skiptics/skiptics_verbose.txt
				cat  dump_us.bed >> res_skiptics/skiptics_verbose.txt
				echo dump_ds exon is  >> res_skiptics/skiptics_verbose.txt
				cat dump_ds.bed >> res_skiptics/skiptics_verbose.txt

				echo dsn $dsn >> res_skiptics/skiptics_verbose.txt
				echo usn $usn >> res_skiptics/skiptics_verbose.txt
				#usn=$(bedtools closest -a dump_us.bed -b t$allexons -s  -D a -id -d -t last)
				upex=$(echo "$usn" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8-1,$9,$4,$5,$6,gen}') #-1 is to comply with bedtools getfasta
				dsex=$(echo "$dsn" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8-1,$9,$4,$5,$6,gen}') #-1 is to comply with bedtools getfasta
				#saving exon size instead for sashimi plots
				upex1=$(echo "$usn" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8-1,$9,$10,$5,$6,gen}') #-1 is to comply with bedtools getfasta
				dsex1=$(echo "$dsn" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8-1,$9,$10,$5,$6,gen}') #-1 is to comply with bedtools getfasta


				#fi
				echo finally upex is "$upex" >> res_skiptics/skiptics_verbose.txt
			  #echo finally ce is "$ce"
			  echo finally dsex is "$dsex" >> res_skiptics/skiptics_verbose.txt
			  #finally paste three segments to a bed file
			  #first remove if any such file exists
			  rm -f "$gene_name"_nt.bed
			  echo "$upex">>"$gene_name"_nt.bed #all_junctions.bed
			  #echo "$ce">>"$gene_name"_nt.bed #all_junctions.bed
			  echo "$dsex">>"$gene_name"_nt.bed #all_junctions.bed

			  cat "$gene_name"_nt.bed >> res_skiptics/skiptics_all_coordinates.bed
				#and also save for sashimi plote
				echo "$upex1">>"$gene_name"_nt1.bed #all_junctions.bed
			  echo "$dsex1">>"$gene_name"_nt1.bed #all_junctions.bed

				cat "$gene_name"_nt1.bed >> res_skiptics/skiptics_sashimi.bed

	      rm "$gene_name"_nt.bed
				rm "$gene_name"_nt1.bed
				cat $sample | awk -v gene=$gene_name -v TX=$TxID 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene,TX}' >> res_skiptics/IGV_skiptics.csv
    fi
	else
		cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_skiptics/IGV_all_others.csv
		echo $csv_ln >> all_non_skiptics.csv
		fi
	elif [ $diff_exon_abs -ge 1 ]  && [ $dsovlp -eq 0 ] && [ $usovlp -eq 0 ] #THIS IS FOR MULTI-SKIP
	then
		#As algorithm detects middle exon for exon_skip events, so determine coordinates of actual exons involved in the event
		echo "$us" | awk 'BEGIN {OFS="\t"} {print $1,$2-10,$3,$4,$5,$6}' > dump_us.bed
		echo "$ds" | awk 'BEGIN {OFS="\t"} {print $1,$2,$3+10,$4,$5,$6}' > dump_ds.bed

			dsn=$(bedtools closest -a dump_ds.bed -b event_bedfiles/t$allexons -s -D a -iu -d -t last )
			usn=$(bedtools closest -a dump_us.bed -b event_bedfiles/t$allexons -s  -D a -id -d -t first)
			upex=$(echo "$usn" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8-1,$9,$4,$5,$6,gen}')
			dsex=$(echo "$dsn" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8-1,$9,$4,$5,$6,gen}')
			#saving exon size instead for sashimi plots
			upex1=$(echo "$usn" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8-1,$9,$10,$5,$6,gen}')
			dsex1=$(echo "$dsn" | awk -v gen="$gene_name" 'BEGIN {OFS="\t"} {print $1,$8-1,$9,$10,$5,$6,gen}')


		#echo finally upex is "$upex"
		#echo finally ce is "$ce"
		#echo finally dsex is "$dsex"
		#finally paste three segments to a bed file
		#first remove if any such file exists
		rm -f "$gene_name"_nt.bed
		echo "$upex">>"$gene_name"_nt.bed #all_junctions.bed
		#echo "$ce">>"$gene_name"_nt.bed #all_junctions.bed
		echo "$dsex">>"$gene_name"_nt.bed #all_junctions.bed

		cat "$gene_name"_nt.bed >> res_skiptics/skiptics_all_coordinates.bed

		#and also save for sashimi plote
		echo "$upex1">>"$gene_name"_nt1.bed #all_junctions.bed
		echo "$dsex1">>"$gene_name"_nt1.bed #all_junctions.bed

		cat "$gene_name"_nt1.bed >> res_skiptics/skiptics_sashimi.bed

		rm "$gene_name"_nt.bed
		rm "$gene_name"_nt1.bed

		cat $sample | awk -v gene=$gene_name -v TX=$TxID 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene,TX}' >> res_skiptics/IGV_skiptics.csv

		cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_skiptics/IGV_multi_skiptics.csv
else
  cat $sample | awk -v gene=$gene_name 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,gene}' >> res_skiptics/IGV_all_others.csv
	echo $csv_ln >> all_non_skiptics.csv
	fi

done
fi


#######NEW CODE
nrecrds=$(cat res_skiptics/skiptics_all_coordinates.bed | wc -l)
nrecrdst=$(($nrecrds/2))
#echo total intronic_range events - at the start - are $nrecrdst


#PRINT STATS SO FAR
total_events=$(cat $inpfile | wc -l)
echo All events read are: $total_events >> res_skiptics/SKIPTICS_FINAL_STATS.txt
ensdb_notfound=0
[ -e res_skiptics/EnsDB_tx_not_found.csv ] && ensdb_notfound=$(cat res_skiptics/EnsDB_tx_not_found.csv| wc -l)
echo Events not found in EnsDB are: $ensdb_notfound , PLEASE SEE res_skiptics/EnsDB_tx_not_found.csv file >> res_skiptics/SKIPTICS_FINAL_STATS.txt
all_others=0
[ -e res_skiptics/IGV_all_others.csv ] && all_others=$(cat res_skiptics/IGV_all_others.csv | wc -l)
remaining_events_processed=$(($total_events-$ensdb_notfound))
echo Remaining Events Processed are: $remaining_events_processed >> res_skiptics/SKIPTICS_FINAL_STATS.txt
echo Out of these $remaining_events_processed events, non_skiptics are: $all_others, PLEASE SEE res_skiptics/IGV_all_others.csv >> res_skiptics/SKIPTICS_FINAL_STATS.txt
echo Total skiptic events are: $nrecrdst >> res_skiptics/SKIPTICS_FINAL_STATS.txt
if [ $nrecrdst -gt 0 ]
then
multi_es=0
[ -e res_skiptics/IGV_multi_skiptics.csv ] && multi_es=$(cat res_skiptics/IGV_multi_skiptics.csv | wc -l)
echo Of these $nrecrdst skiptic events, multi_skiptic events are: $multi_es, PLEASE SEE res_skiptics/IGV_multi_skiptics.csv >> res_skiptics/SKIPTICS_FINAL_STATS.txt



readarray -t all_data < res_skiptics/skiptics_all_coordinates.bed
#also read csv for TxID
readarray -t all_csv_data < res_skiptics/IGV_skiptics.csv
#also read sashimi bed file to get unique entries
readarray -t all_sashimi_data < res_skiptics/skiptics_sashimi.bed

i=0
eventn=0
while [ $i -lt $nrecrds ]
  do
          line11=${all_data[$i]}
          line12=${all_data[$(($i+1))]}
					#sashimi records
					lines1=${all_sashimi_data[$i]}
          lines2=${all_sashimi_data[$(($i+1))]}

          #line13=${all_data[$(($i+2))]}
#         echo i is $i, $line11, $line12, $line13
				txid_recrd=${all_csv_data[$eventn]}
				txid=$(echo $txid_recrd | awk 'BEGIN {FS=","} {print $9}')
        ((eventn=eventn+1))
     #echo "Welcome $i times"
     ((i=i+2))
     j=$i
     #now go through rest of the data
     flg=0
        while [ $j -lt $nrecrds ]
        do
                line21=${all_data[$j]}
          line22=${all_data[$(($j+1))]}
          #line23=${all_data[$(($j+2))]}
          if [ "${line11[*]}" == "${line21[*]}" ] && [ "${line12[*]}" == "${line22[*]}" ]
          then
            #echo i is $i, $line11, $line12, $line13
            #echo j is $j, $line21, $line22, $line23
            flg=1
          fi
          ((j=j+2))
        done
        if [ $flg -eq 0 ]
				then
                 echo $line11 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_skiptics/skiptics_unique.bed
                 echo $line12 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_skiptics/skiptics_unique.bed
								 #This is used of sashimi plots, should do better

								 echo $lines1 | awk -v tx=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,tx}' >> res_skiptics/skiptics_uniq_sashimi.bed
								 echo $lines2 | awk -v tx=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,tx}' >> res_skiptics/skiptics_uniq_sashimi.bed

                 #echo $line13 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_skiptics_intronic_range/intronic_range_all_unique.bed
				else
					echo $line11 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_skiptics/skiptics_repeated.bed
					echo $line12 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_skiptics/skiptics_repeated.bed
					#echo $line13 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' >> res_skiptics_intronic_range/intronic_range_repeated.bed
					echo $line11 | awk 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,$7}' >> res_skiptics/IGV_skiptics_repeated.csv
					echo $line12 | awk 'BEGIN {OFS=","} {print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6,$7}' >> res_skiptics/IGV_skiptics_repeated.csv

         fi
 done
 uniq_eventst=$(cat res_skiptics/skiptics_unique.bed | wc -l)
 uniq_events=$(($uniq_eventst/2))
 rep_eventst=$(cat res_skiptics/skiptics_repeated.bed | wc -l)
 rep_events=$(($rep_eventst/2))

echo Out of total $nrecrdst skiptic events, unique skiptic events are: $uniq_events, PLEASE SEE res_skiptics/skiptics_unique.bed >> res_skiptics/SKIPTICS_FINAL_STATS.txt
echo And repeated skiptic events are: $rep_events, PLEASE SEE res_skiptics/skiptics_repeated.bed and res_skiptics/IGV_skiptics_repeated.csv >> res_skiptics/SKIPTICS_FINAL_STATS.txt
echo THOSE ARE ALL PERTINENT STATISTICS, PLEASE LET US KNOW IF ANY OTHER INFORMATION MIGHT BE USEFUL!!! >> res_skiptics/SKIPTICS_FINAL_STATS.txt

#also read sashimi bed file to get TxID
readarray -t all_sashimi_data < res_skiptics/skiptics_uniq_sashimi.bed
i=0
cat res_skiptics/skiptics_unique.bed | while read r1; read r2
#write for IGV
do
	lines1=${all_sashimi_data[$i]}
	txid=$(echo "$lines1" | awk 'BEGIN {FS="\t"} {print $8}')
	#echo got txid $txid at line $i
	((i=i+2))
		echo "$r1" "$r2" | awk -v tx=$txid 'BEGIN {OFS=","} {print $1":"$3"-"$9,$1,$3,$9,$4,$5,$6,$7,tx}' >> res_skiptics/IGV_unique_skiptics_translated.csv
done
##########NEW CODE ENDS
#fasta_flg=0
#if [ $fasta_flg -eq 1 ]
#then
echo STARTING NT and AA TRANSLATION
#now call bedtools getfasta function to get nt sequence from reference_genome.fa file
bedtools getfasta -fi GRCh38.p13.genome.fa -bed res_skiptics/skiptics_unique.bed -s > res_skiptics/skiptics_unique_nt.fasta

#now remove >chr lines from the res_skipticsulting file
awk '!/^>chr/' res_skiptics/skiptics_unique_nt.fasta > skiptics_unique_nt1.fasta
#combine the two files to get desired csv file - chrXX,start,end,genename
paste -d"\t" res_skiptics/skiptics_unique.bed skiptics_unique_nt1.fasta > res_skiptics/skiptics_unique_nt.bed
#remove temp file
rm skiptics_unique_nt1.fasta
#and finally transeq compatible
awk -F "\t" '{if($6=="+") print ">sp|"$7"_"$1"_"$2"_"$3"_plus""\n"$8; else print ">sp|"$7"_"$1"_"$2"_"$3"_minus""\n"$8}' res_skiptics/skiptics_unique_nt.bed > res_skiptics/skiptics_unique_transeq_in.fasta
rm res_skiptics/skiptics_unique_nt.bed
#awk -F "\t" '{print ">sp|"$7"-"$1":"$2"-"$3"\n"$8}' ce_all_unique_coord_nt.bed > ce_all_unique_coord_nt.transeq_in.fasta

#now do the 3-frame translation

TEMPFILE=$(echo res_skiptics/skiptics_unique_transeq_in.fasta  | perl -pe 's/\.fasta$/.temp/')
OUTPUTFILE=$(echo res_skiptics/skiptics_unique_transeq_in.fasta  | perl -pe 's/\.fasta$/.trans/')

# Translate sequence
transeq -sequence res_skiptics/skiptics_unique_transeq_in.fasta  -outseq "$TEMPFILE" -frame F


# Rename sequence
perl -sape 's/>/>sp|/' "$TEMPFILE" > "$OUTPUTFILE"



#also remove all newlines from the
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' "$OUTPUTFILE" > res_skiptics/SKIPTICS_AA.fasta
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
done < res_skiptics/SKIPTICS_AA.fasta > res_skiptics/FINAL_SKIPTICS_AA.fasta

echo AND NOW STARTING NT and AA TRANSLATION FOR CONCATENATED EXON RANGES
#NOW concatenate nt sequence for exon_skip
#now form concatenated coordinates
event_i=1
cat res_skiptics/skiptics_unique.bed | while read r1; read r2
do
		strnd1=$(echo "$r1" | awk '{print $6}')
		if [[ $strnd1 == "+" ]]
		then

			echo "$r1" "$r2" | awk '{OFS="\t"} {print}' | tr ' ' '\n' > dump.bed
		else
			echo "$r2" "$r1" | awk '{OFS="\t"} {print}' | tr ' ' '\n' > dump.bed
		fi

		bedtools getfasta -fi GRCh38.p13.genome.fa -bed dump.bed -s > dump.fasta
		#now remove >chr lines from the resulting file
		awk '!/^>chr/' dump.fasta > dump1.fasta

		#combine the two files to get desired bed file - chrXX,start,end,genename
		#paste -d"\t" dump.bed dump1.fasta > dump.genome_seq.bed
		if [[ $strnd1 == "+" ]]
		then
			#title=$(cat dump.bed | tr -s '\n' ' ' | awk '{print ">sp|"$1"_"$2"_"$3"_"$9"_"$10"_"$14"_plus"}')
			title=$(cat dump.bed | tr -s '\n' ' ' | awk -v evid=$event_i '{print ">sp|"$14"_"$1"_"$2"_"$3"_"$9"_"$10"_"evid"_plus"}')
			awk '!/^>chr/' dump.fasta | tr -d '\n' | awk -v p=$title '{print p,$1}' | tr ' ' '\n' >> res_skiptics/skiptics_fused_transeq_in.fasta
		else
			#title=$(cat dump.bed | tr -s '\n' ' ' | awk '{print ">sp|"$1"_"$9"_"$10"_"$2"_"$3"_"$14"_minus"}')
			#NOW CHANGIN
			title=$(cat dump.bed | tr -s '\n' ' ' | awk -v evid=$event_i '{print ">sp|"$14"_"$1"_"$2"_"$3"_"$9"_"$10"_"evid"_minus"}')
			#title=$(cat dump.bed | tr -s '\n' ' ' | awk -v evid=$event_i '{print ">sp|"$14"_"$1"_"$9"_"$10"_"$2"_"$3"_"evid"_minus"}')
			awk '!/^>chr/' dump.fasta | tr -d '\n' | awk -v p=$title '{print p,$1}' | tr ' ' '\n' >> res_skiptics/skiptics_fused_transeq_in.fasta

		fi
		((event_i=event_i+1))
done

#cleanup
#rm dump.bed
#rm dump.fasta
rm dump1.fasta
#rm dump.genome_seq.bed

#now do the 3-frame translation

TEMPFILE=$(echo res_skiptics/skiptics_fused_transeq_in.fasta  | perl -pe 's/\.fasta$/.temp/')
OUTPUTFILE=$(echo res_skiptics/skiptics_fused_transeq_in.fasta  | perl -pe 's/\.fasta$/.trans/')

# Translate sequence
transeq -sequence res_skiptics/skiptics_fused_transeq_in.fasta  -outseq "$TEMPFILE" -frame F

# Rename sequence
perl -sape 's/>/>sp|/' "$TEMPFILE" > "$OUTPUTFILE"

#also remove all newlines from the
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' "$OUTPUTFILE" > res_skiptics/SKIPTICS_FUSED_AA.fasta
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
done < res_skiptics/SKIPTICS_FUSED_AA.fasta > res_skiptics/FINAL_SKIPTICS_FUSED_AA.fasta
#copy all event bed files to es_bedfiles folder for debug purpose
#mkdir -p es_event_bedfiles
[ "$(ls -A es_event_bedfiles/)" ] && rm es_event_bedfiles/*.* && mkdir -p es_event_bedfiles && mv event_bedfiles/*.bed es_event_bedfiles/.
#[ "$(ls -A event_bedfiles/)" ] && mkdir -p es_event_bedfiles && mv event_bedfiles/*.bed es_event_bedfiles/.
#final clean up
#rm *.bed

else

[ "$(ls -A es_event_bedfiles/)" ] && rm es_event_bedfiles/*.*
[ "$(ls -A event_bedfiles/)" ] && mkdir -p es_event_bedfiles && mv event_bedfiles/*.bed es_event_bedfiles/.
#final clean up
#rm *.bed
echo AS SKIPTIC EVENTS ARE: $nrecrdst , SO EXITING >> res_skiptics/SKIPTICS_FINAL_STATS.txt
fi
#fi #FASTA FILE CREATION FLAG
