#!/bin/bash
#declare encodings
#-*- coding: utf-8 -*-

#THIS SI FINAL SCRIPT FOR SASHIMI PLOTS FOR ALL PARTS OF PGP

#THIS SCRIPT CONTAINS SASHIMI PLOT CODE FOR
#0. PLEASE NOTE THAT FOLLOWING 2 FILES (OR SOFT LINKS) SHOULD BE IN CURRENT FOLDER
	#01: ggsashimi_txV4.py
	#02: Homo_sapiens.GRCh38.103.chr.sorted_new.gtf
#1. Skiptic Events
#2. ALL MAJIQ EVENTS
#3. CE (INCLUDING INCLUSION, EXTENSION AND IR) events


#First CHECK IF Called from pgp-a/b or pgp-c
#CHECK IF 3 ARGUMENTS ARE PROVIDED
if [ $# -eq 3 ]
then


	######### NEW - Now read input csv and bed files and flag

	inp_csv=$1
	inp_bed=$2

	#get folder
	inp_prefix=$(echo $1|cut -d'/' -f1)

	if [ $3 -eq 1 ]
	then
		[ "$(ls -A "$inp_prefix"/sashimi_plots/)" ] && rm "$inp_prefix"/sashimi_plots/*.*
		mkdir -p "$inp_prefix"/sashimi_plots
		bed="$inp_bed" #res_skiptics/skiptics_uniq_sashimi.bed  #skiptics_sashimi.bed #"sam_junc1.txt" #"junctions.txt"
		readarray -t all_bed_data < $bed
		nrecrds=$(cat $bed | wc -l)
		nrecrdst=$((nrecrds/2))
		echo read $nrecrdst records
		csv="$inp_csv" #res_skiptics/IGV_unique_skiptics_translated.csv #"sam_junc1.txt" #"junctions.txt"
		readarray -t all_csv_data < $csv
		#nrecrdst=$(($nrecrds/3))
		#while IFS= read -r line;
		i=0
		eventn=0

		while [ $i -lt $nrecrds ]
		do
			#construct string for ggsashimi
			line1=${all_bed_data[$i]}
		line2=${all_bed_data[$(($i+1))]}
			#line3=${all_bed_data[$(($i+2))]}
			#echo line1 $line1 line 2 $line2
		((i=i+2))
			#((eventn=eventn+1))
			#read strand
			strnd=$(echo $line1 | awk '{print $6}')
			#also read TxID
			TxID=$(echo $line1 | awk '{print $8}')

			#us exon length
			exon1=0 #$(echo $line1 | awk '{print $4}')
			#ds exon length
			exon2=0 #$(echo $line2 | awk '{print $4}')

			#if [ $strnd == "+" ]
			#then

			line=$(echo $line1 $line2 | awk '{print $1":"$2-50"-"$11+50}')



			#else
			#	line=$(echo $line1 $line2 | awk '{print $1":"$2"-"$11}')
			#fi
			#echo line $line
			event=${all_csv_data[$eventn]}
			gene_name=$(echo $event | awk 'BEGIN {FS=","} {print $8}')
			#fn=$(echo $line1 $line2 $line3 | awk '{print $7"-"$1"_"$2"-"$9}')
			fn=$(echo $event | awk 'BEGIN {FS=","} {print $8"-"$2"_"$3"-"$4}')
			#string for majiq event
			chr_name=$(echo $event | awk 'BEGIN {FS=","} {print $2}')
			start=$(echo $event | awk 'BEGIN {FS=","} {print $3}')
			end=$(echo $event | awk 'BEGIN {FS=","} {print $4}')
			#gene_name=$(echo $event | awk 'BEGIN {FS=","} {print $8}')
			#intron=$(($end-$start))
			#gene_name-start-intron-end #for now using strand from ggsashimi
			majiq_event=$(echo $chr_name $start $end | awk '{print $1"-"$2"-"$3}')
			#also get actual event identified
			event_identified=$(echo $line1 $line2 | awk '{print "PGPEvent-"$3"-"$10}')


			((eventn=eventn+1))
			echo processing event num $eventn and event $fn #and majiq event is $majiq_event and event identified is $event_identified

			#./ggsashimi_txV5.py -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -GeneName "$gene_name" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o res_skiptics/sashimi_plots/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
			#Here removed -PGPTx flag
			#./ggsashimi_txV5.py -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -GeneName "$gene_name" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o res_skiptics/sashimi_plots/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
			./ggsashimi_txV3.py -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -GeneName "$gene_name" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o "$inp_prefix"/sashimi_plots/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=2.5 --height=2.5 --width=18 -P palette.txt
					##rsvg-convert -f pdf -o "$inp_prefix"/sashimi_plots/"$fn".pdf "$inp_prefix"/sashimi_plots/"$fn".svg
				#[ -e res_skiptics/sashimi_plots/"$fn".svg ] && rm res_skiptics/sashimi_plots/"$fn".svg

		######
		#./ggsashimi_tx.py -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o "$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
		########


		done
		#now merge all pdf's

		#pdfunite "$inp_prefix"/sashimi_plots/*.pdf "$inp_prefix"/all_skiptics_sashimi.pdf
		python merge_sashimis.py "$inp_prefix"/sashimi_plots/

	fi

	######THIS IS FOR ALL MAJIQ EVENTS

	if [ $3 -eq 2 ]
	then
		[ "$(ls -A all_events_sashimi/)" ] && rm all_events_sashimi/*.*
		mkdir -p all_events_sashimi
		bed=$2_all_sashimi.bed  #skiptics_sashimi.bed #"sam_junc1.txt" #"junctions.txt"
		readarray -t all_bed_data < $bed
		nrecrds=$(cat $bed | wc -l)
		echo read $nrecrds records
		csv=$2_all_sashimi.csv #"sam_junc1.txt" #"junctions.txt"
		readarray -t all_csv_data < $csv

		i=0
		eventn=0
		while [ $i -lt $nrecrds ]
		do
			#construct string for ggsashimi
			line1=${all_bed_data[$i]}
		((i=i+1))
			#read strand
			#strnd=$(echo $line1 | awk '{print $6}')
			strnd='+'
			#also read TxID
			TxID=$(echo $line1 | awk '{print $8}')
			line=$(echo $line1 | awk '{print $1":"$2-50"-"$3+50}') #this is the actual event
			#get majiq event
			event=${all_csv_data[$eventn]}
			((eventn=eventn+1))
			#fn=$(echo $line1 $line2 $line3 | awk '{print $7"-"$1"_"$2"-"$9}')
			#THIS IS TO KEEP DIFFERENT FILE NAMES FOR EVENTS WITH IDENTICAL COORDINATES AND GENE_NAMES
				gene_name=$(echo $event | awk 'BEGIN {FS=","} {print $5}')
						if [ $i -eq 1 ]
						then
							temp_gene=$(echo $event | awk 'BEGIN {FS=","} {print $5}')
							trackj=1
						elif [ "$temp_gene" == "$gene_name" ]
						then
							((trackj=trackj+1))
						else
							temp_gene=$(echo $event | awk 'BEGIN {FS=","} {print $5}')
							trackj=1
						fi
			#SHOULD NOT HAVE THIS BUT NEED TO KEEP AS PER REQUEST
			fn=$(echo $event | awk -v tr=$trackj 'BEGIN {FS=","} {print $5"-"$1"-"$2"-"$3"-"tr}') #also gene_id to avoid same file names for repeated events
			#string for majiq event
			#echo event $event
			chr_name=$(echo $event | awk 'BEGIN {FS=","} {print $1}')
			start=$(echo $event | awk 'BEGIN {FS=","} {print $2}')
			end=$(echo $event | awk 'BEGIN {FS=","} {print $3}')

			#majiq_event=$(echo $chr_name $start $end | awk '{print $1"-"$2"-"$3}')
			majiq_event=$(echo $chr_name $start $end | awk '{print $1"-"$2"-"$3}')
			exon1=0
			exon2=0
			#also get actual event identified
			event_identified=$(echo $line1 | awk '{print "None ""-"$2"-"$3}')


			echo processing event num $eventn and event $fn #and majiq event is $majiq_event

		#	./ggsashimi_txV4.py -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -GeneName "$gene_name" -MajiqStrnd "$strnd" -ORIG 0 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o all_events_sashimi/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
		#             ##rsvg-convert --unlimited  -f pdf -o all_events_sashimi/"$fn".pdf all_events_sashimi/"$fn".svg
			#./ggsashimi_txV5.py -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -GeneName "$gene_name" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o all_events_sashimi/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
			#Here removed -PGPTx flag
			#./ggsashimi_txV5.py -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -GeneName "$gene_name" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o all_events_sashimi/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
			./ggsashimi_txV3.py -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -GeneName "$gene_name" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o all_events_sashimi/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=2.5 --height=2.5 --width=18 -P palette.txt
			#./ggsashimi_txV32.py -S "minus" -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -GeneName "$gene_name" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o all_events_sashimi/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=2.5 --height=2.5 --width=18 -P palette.txt

										##rsvg-convert --unlimited  -f pdf -o all_events_sashimi/"$fn".pdf all_events_sashimi/"$fn".svg

				#[ -e all_events_sashimi/"$fn".svg ] && rm all_events_sashimi/"$fn".svg

		done
		#now merge all pdf's

		#pdfunite all_events_sashimi/*.pdf all_events_sashimi/all_events_sashimi.pdf
		python merge_sashimis.py all_events_sashimi/
	fi


	#THIS SECTION IS FOR CE_INCLUSION EVENTS
	#WILL MERGE INCLUSION AND EXTENSION EVENTS
	if [ $3 -eq 3 ]
	then
		#[ "$(ls -A res_ce_all/ce_incl_sashimi_plots/)" ] && rm res_ce_all/ce_incl_sashimi_plots/*.*
		[ "$(ls -A "$inp_prefix"/ce_incl_sashimi_plots/)" ] && rm "$inp_prefix"/ce_incl_sashimi_plots/*.*
		mkdir -p "$inp_prefix"/ce_incl_sashimi_plots
		#cat ce_inclusion_coord_uniq_sashimi.bed ce_extension_coord_uniq_sashimi.bed > all_ce.bed
		#cat IGV_unique_ce_inclusion.csv IGV_unique_ce_extension.csv > all_ce.csv
		bed=$inp_bed #res_ce_all/ce_inclusion_coord_uniq_sashimi.bed #"sam_junc1.txt" #"junctions.txt"
		#readarray -t all_bed_data < $bed

		readarray -t all_bed_data < $bed
		nrecrds=$(cat $bed | wc -l)
		nrecrdst=$((nrecrds/3))
		echo read $nrecrdst records
		#csv=res_ce_all/IGV_unique_ce_inclusion.csv #"sam_junc1.txt" #"junctions.txt"
		csv=$inp_csv #"sam_junc1.txt" #"junctions.txt"
		readarray -t all_csv_data < $csv

		# bed=ce_extension_coord_uniq_sashimi.bed #"sam_junc1.txt" #"junctions.txt"
		# readarray -t all_bed_data < $bed
		# nrecrds=$(cat $bed | wc -l)
		# nrecrdst=$((nrecrds/3))
		# echo read $nrecrdst records
		# csv=IGV_unique_ce_extension.csv #"sam_junc1.txt" #"junctions.txt"
		# readarray -t all_csv_data < $csv


		#nrecrdst=$(($nrecrds/3))
		#while IFS= read -r line;
		i=0
		eventn=0
		while [ $i -lt $nrecrds ]
		do
			#construct string for ggsashimi
			line1=${all_bed_data[$i]}
		line2=${all_bed_data[$(($i+1))]}
			line3=${all_bed_data[$(($i+2))]}

			#echo line1 $line1 line 2 $line2
		((i=i+3))
			#((eventn=eventn+1))
			#read strand
			strnd=$(echo $line1 | awk '{print $6}')
			#also read TxID
			TxID=$(echo $line1 | awk '{print $8}')

			if [ $strnd == "+" ]
			then
				strndflg=plus
				#us exon length
				exon1=$(echo $line1 | awk '{print $4}')
				#ds exon length
				exon2=$(echo $line3 | awk '{print $4}')
				#now modify to reflect whole up/dn exons
				up=$(echo $line1 | awk '{print $2}')
				upn=$(($up-$exon1))
				#ds exon
				dn=$(echo $line3 | awk '{print $3}')
				dnn=$(($dn+$exon2))

				fulltitle=$(echo $line1 $line2 $line3 | awk '{print $1"-"$2":"$3"-"$10":"$11"-"$18":"$19}')
				#line=$(echo $line1 $line2 $line3 | awk -v ex1="$exon1" ex2="$exon2" '{print $1":"($2-50-$ex1)"-"($19+50+$ex2)}')
				#line=$(echo $line1 $upn $dnn | awk '{print $1":"$9-50"-"$10+50}')
				line=$(echo $line1 $upn $dnn | awk '{print $1":"$9"-"$10}')
			else
				strndflg=minus
				#us exon length
				exon1=$(echo $line3 | awk '{print $4}')
				#ds exon length
				exon2=$(echo $line1 | awk '{print $4}')
				#now modify to reflect whole up/dn exons
				dn=$(echo $line3 | awk '{print $2}')
				dnn=$(($dn-$exon1))
				#ds exon
				up=$(echo $line1 | awk '{print $3}')
				upn=$(($up+$exon2))

				fulltitle=$(echo $line1 $line2 $line3 | awk '{print $1"-"$18":"$19"-"$10":"$11"-"$2":"$3}')
				#line=$(echo $line1 $line2 $line3 | awk -v ex1="$exon1" ex2="$exon2" '{print $1":"($18-50-$ex1)"-"($2+50+$ex2)}')
				#line=$(echo $line1 $dnn $upn | awk '{print $1":"$9-50"-"$10+50}')
				line=$(echo $line1 $dnn $upn | awk '{print $1":"$9"-"$10}') #60 nt are already off
			fi
			#echo line $line
			event=${all_csv_data[$eventn]}
			#fn=$(echo $line1 $line2 $line3 | awk '{print $7"-"$1"_"$2"-"$9}')
			fn=$(echo $event | awk 'BEGIN {FS=","} {print $8"-"$2"_"$3"-"$4}')

			#string for majiq event
			chr_name=$(echo $event | awk 'BEGIN {FS=","} {print $2}')
			start=$(echo $event | awk 'BEGIN {FS=","} {print $3}')
			end=$(echo $event | awk 'BEGIN {FS=","} {print $4}')
			gene_name=$(echo $event | awk 'BEGIN {FS=","} {print $8}')
			intron=$(($end-$start))
			#gene_name-start-intron-end #for now using strand from ggsashimi
			#majiq_event=$(echo $gene_name $start $intron $end | awk '{print $1"-"$2"-"$3"-"$4}')
			majiq_event=$(echo $chr_name $start $end | awk '{print $1"-"$2"-"$3}')
			#also get actual event identified
			#event_identified=$(echo $line2 | awk '{print "PGPEvent""-"$2"-"$3}')
			event_identified=$(echo $line2 | awk '{print $1"-"$2"-"$3}')

			((eventn=eventn+1))
			echo processing event num $eventn and event $fn #and line is $line


			#./ggsashimi_txV5.py  -PGPTx "$event_identified" -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -FullTitle "$fulltitle" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o res_ce_all/ce_incl_sashimi_plots/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
		#Here removed -PGPTx
			./ggsashimi_txV3.py -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -GeneName "$gene_name" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -FullTitle "$fulltitle" -MajiqTx "$majiq_event" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o "$inp_prefix"/ce_incl_sashimi_plots/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=2.5 --height=2.5 --width=18 -P palette.txt
					##rsvg-convert -f pdf -o "$inp_prefix"/ce_incl_sashimi_plots/"$fn".pdf "$inp_prefix"/ce_incl_sashimi_plots/"$fn".svg
				#[ -e res_ce_all/ce_incl_sashimi_plots/"$fn".svg ] && rm res_ce_all/ce_incl_sashimi_plots/"$fn".svg

		done
		#now merge all pdf's

		#pdfunite "$inp_prefix"/ce_incl_sashimi_plots/*.pdf "$inp_prefix"/ce_incl_all_sashimi_plots.pdf

		python merge_sashimis.py "$inp_prefix"/ce_incl_sashimi_plots/

	fi
	#CE_EXTENSION

	if [ $3 -eq 4 ]
	then
		[ "$(ls -A "$inp_prefix"/ce_ext_sashimi_plots/)" ] && rm "$inp_prefix"/ce_ext_sashimi_plots/*.*
		mkdir -p "$inp_prefix"/ce_ext_sashimi_plots
		bed=$inp_bed #res_ce_all/ce_extension_coord_uniq_sashimi.bed #"sam_junc1.txt" #"junctions.txt"
		readarray -t all_bed_data < $bed
		nrecrds=$(cat $bed | wc -l)
		nrecrdst=$((nrecrds/3))
		echo read $nrecrdst records
		csv=$inp_csv #res_ce_all/IGV_unique_ce_extension.csv #"sam_junc1.txt" #"junctions.txt"
		readarray -t all_csv_data < $csv

		#nrecrdst=$(($nrecrds/3))
		#while IFS= read -r line;
		i=0
		eventn=0
		while [ $i -lt $nrecrds ]
		do
			#construct string for ggsashimi
			line1=${all_bed_data[$i]}
		line2=${all_bed_data[$(($i+1))]}
			line3=${all_bed_data[$(($i+2))]}

			#echo line1 $line1 line 2 $line2
		((i=i+3))
			#((eventn=eventn+1))
			#read strand
			strnd=$(echo $line1 | awk '{print $6}')
			#also read TxID
			TxID=$(echo $line1 | awk '{print $8}')

			if [ $strnd == "+" ]
			then
				#strndflg=plus
				#us exon length
				exon1=$(echo $line1 | awk '{print $4}')
				#ds exon length
				exon2=$(echo $line3 | awk '{print $4}')
				#now modify to reflect whole up/dn exons
				up=$(echo $line1 | awk '{print $2}')
				upn=$(($up-$exon1))
				#ds exon
				dn=$(echo $line3 | awk '{print $3}')
				dnn=$(($dn+$exon2))

				fulltitle=$(echo $line1 $line2 $line3 | awk '{print $1"-"$2":"$3"-"$10":"$11"-"$18":"$19}')
				#line=$(echo $line1 $line2 $line3 | awk -v ex1="$exon1" ex2="$exon2" '{print $1":"($2-50-$ex1)"-"($19+50+$ex2)}')
				line=$(echo $line1 $upn $dnn | awk '{print $1":"$9"-"$10}')
			else
				#strndflg=minus
				#us exon length
				exon1=$(echo $line3 | awk '{print $4}')
				#ds exon length
				exon2=$(echo $line1 | awk '{print $4}')
				#now modify to reflect whole up/dn exons
				dn=$(echo $line3 | awk '{print $2}')
				dnn=$(($dn-$exon1))
				#ds exon
				up=$(echo $line1 | awk '{print $3}')
				upn=$(($up+$exon2))

				fulltitle=$(echo $line1 $line2 $line3 | awk '{print $1"-"$18":"$19"-"$10":"$11"-"$2":"$3}')
				#line=$(echo $line1 $line2 $line3 | awk -v ex1="$exon1" ex2="$exon2" '{print $1":"($18-50-$ex1)"-"($2+50+$ex2)}')
				line=$(echo $line1 $dnn $upn | awk '{print $1":"$9"-"$10}')
			fi
			#echo line $line
			event=${all_csv_data[$eventn]}
			#fn=$(echo $line1 $line2 $line3 | awk '{print $7"-"$1"_"$2"-"$9}')
			fn=$(echo $event | awk 'BEGIN {FS=","} {print $8"-"$2"_"$3"-"$4}')

			#string for majiq event
			chr_name=$(echo $event | awk 'BEGIN {FS=","} {print $2}')
			start=$(echo $event | awk 'BEGIN {FS=","} {print $3}')
			end=$(echo $event | awk 'BEGIN {FS=","} {print $4}')
			gene_name=$(echo $event | awk 'BEGIN {FS=","} {print $8}')
			intron=$(($end-$start))
			#gene_name-start-intron-end #for now using strand from ggsashimi
			#majiq_event=$(echo $gene_name $start $intron $end | awk '{print $1"-"$2"-"$3"-"$4}')
			majiq_event=$(echo $chr_name $start $end | awk '{print $1"-"$2"-"$3}')
			#also get actual event identified
			#event_identified=$(echo $line2 | awk '{print "PGPEvent""-"$2"-"$3}')
			event_identified=$(echo $line2 | awk '{print $1"-"$2"-"$3}')

			((eventn=eventn+1))
			echo processing event num $eventn and event $fn #, pgp_identified event $event_identified and line is $line


			#./ggsashimi_txV5.py -PGPTx "$event_identified" -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -FullTitle "$fulltitle" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o res_ce_all/ce_ext_sashimi_plots/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
			#Here removed -PGPTx
			./ggsashimi_txV3.py -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -GeneName "$gene_name" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -FullTitle "$fulltitle" -MajiqTx "$majiq_event" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o "$inp_prefix"/ce_ext_sashimi_plots/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=2.5 --height=2.5 --width=18 -P palette.txt
				##rsvg-convert -f pdf -o "$inp_prefix"/ce_ext_sashimi_plots/"$fn".pdf "$inp_prefix"/ce_ext_sashimi_plots/"$fn".svg
			#[ -e res_ce_all/ce_ext_sashimi_plots/"$fn".svg ] && rm res_ce_all/ce_ext_sashimi_plots/"$fn".svg

		done

		#now merge all pdf's

		#pdfunite "$inp_prefix"/ce_ext_sashimi_plots/*.pdf "$inp_prefix"/ce_exten_all_sashimi_plots.pdf
		python merge_sashimis.py "$inp_prefix"/ce_ext_sashimi_plots/
	fi

	######THIS IS FOR ALL IR EVENTS

	if [ $3 -eq 5 ]
	then
		[ "$(ls -A "$inp_prefix"/ir_sashimi_plots/)" ] && rm "$inp_prefix"/ir_sashimi_plots/*.*
		mkdir -p "$inp_prefix"/ir_sashimi_plots
		bed=$inp_bed #res_ce_all/IR_coord_uniq_sashimi.bed  #skiptics_sashimi.bed #"sam_junc1.txt" #"junctions.txt"
		readarray -t all_bed_data < $bed
		nrecrds=$(cat $bed | wc -l)
		echo read $nrecrds records
		csv=$inp_csv #res_ce_all/IGV_unique_IR.csv #"sam_junc1.txt" #"junctions.txt"
		readarray -t all_csv_data < $csv

		i=0
		eventn=0

		while [ $i -lt $nrecrds ]
		do
				#construct string for ggsashimi
				line1=${all_bed_data[$i]}
		((i=i+1))
				#read strand
				strnd=$(echo $line1 | awk '{print $6}')
				#also read TxID
				TxID=$(echo $line1 | awk '{print $8}')
				line=$(echo $line1 | awk '{print $1":"$2-50"-"$3+50}') #this is the actual event
				#get majiq event
				event=${all_csv_data[$eventn]}
				((eventn=eventn+1))
				#fn=$(echo $line1 $line2 $line3 | awk '{print $7"-"$1"_"$2"-"$9}')
				#THIS IS TO KEEP DIFFERENT FILE NAMES FOR EVENTS WITH IDENTICAL COORDINATES AND GENE_NAMES
					gene_name=$(echo $event | awk 'BEGIN {FS=","} {print $5}')
										if [ $i -eq 1 ]
										then
												temp_gene=$(echo $event | awk 'BEGIN {FS=","} {print $5}')
												trackj=1
										elif [ "$temp_gene" == "$gene_name" ]
										then
												((trackj=trackj+1))
										else
												temp_gene=$(echo $event | awk 'BEGIN {FS=","} {print $5}')
												trackj=1
										fi
				#SHOULD NOT HAVE THIS BUT NEED TO KEEP AS PER REQUEST
				fn=$(echo $event | awk -v tr=$trackj 'BEGIN {FS=","} {print $5"-"$1"-"$2"-"$3"-"tr}') #also gene_id to avoid same file names for repeated events
				#string for majiq event
				#echo event $event
				chr_name=$(echo $event | awk 'BEGIN {FS=","} {print $1}')
				start=$(echo $event | awk 'BEGIN {FS=","} {print $2}')
				end=$(echo $event | awk 'BEGIN {FS=","} {print $3}')

				#majiq_event=$(echo $chr_name $start $end | awk '{print $1"-"$2"-"$3}')
				majiq_event=$(echo $chr_name $start $end | awk '{print $1"-"$2"-"$3}')
				exon1=0
				exon2=0
				#also get actual event identified
				event_identified=$(echo $line1 | awk '{print $1"-"$2"-"$3}')


				echo processing event num $eventn and event $fn #and majiq event is $majiq_event
				#./ggsashimi_txV5.py -PGPTx "$event_identified"  -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -GeneName "$gene_name" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o res_ce_all/ir_sashimi_plots/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
						./ggsashimi_txV3.py -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -GeneName "$gene_name" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o "$inp_prefix"/ir_sashimi_plots/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=2.5 --height=2.5 --width=18 -P palette.txt
					##rsvg-convert --unlimited  -f pdf -o "$inp_prefix"/ir_sashimi_plots/"$fn".pdf "$inp_prefix"/ir_sashimi_plots/"$fn".svg
				#[ -e res_ce_all/ir_sashimi_plots/"$fn".svg ] && rm res_ce_all/ir_sashimi_plots/"$fn".svg

		done
		#now merge all pdf's

		#pdfunite "$inp_prefix"/ir_sashimi_plots/*.pdf "$inp_prefix"/IR_all_sashimi_plots.pdf
		python merge_sashimis.py "$inp_prefix"/ir_sashimi_plots/
	fi

else

	echo came in for peaksbackmapping
	######THIS PART IS FOR BACKMAPPING
	flg_sashimi_files_only=1
	if [ $flg_sashimi_files_only -eq 1 ]
	then
		step01_flg=1
		if [ $step01_flg -eq 1 ]
		then
			echo "################################ GENERATING BED FILE FOR SASHIMI PLOTS FOR ALL EVENTS"
			#echo "################################ GENERATING BED FILE FOR SASHIMI PLOTS FOR ALL EVENTS" >> ALL_STATS.txt
			#IMPORTANT - THIS CODE RELIES ON THE OUTPUT OF BEDTOOLS CLOSEST FUNCTION
			#GIVEN A MAJIQ EVENT (as a bed file) and BED FILE FOR THE TRANSCRIPT IT BELIEVES TO BE PART OF
			#BEDTOOLS CLOSEST FUNCTION USES INPUT RANGE (AMJIQ JUNTION) and finds closest exons (and their ranges) from the TRANSCRIPT BED FILE
			#BEDTOOLS CLOSEST RETURNS A BED FILE WITH FOLLOWING OUTPUT
			#INPUT: chr,start,end,1,0,+ (majiq event with strand)
			#OUTPUT: chr,start,end,1,0,+,chr,start,end,size,exon_num,+,distance
			#here first 6 entries are the original majiq input and next 7 entries are the resulting closest exon coordinates, its size (in bp),exon_num,strand and distance from reference
			#SHOULD ADD CHECKS ON EXON_NUM READING FROM FILE


			#dir=$(echo $1|cut -d'/' -f1)
			inp=$(echo $1|cut -d'.' -f1)

			echo run_sashimiV1.sh, got inp as $inp

			#also get folder
			folder=$(echo $inp|cut -d'/' -f1)
			echo run_sashimiV1.sh, got folder as $folder
			#also get event TYPE
			eventyp=$(echo $2|cut -d'/' -f2|cut -d'.' -f1)
			echo run_sashimiV1.sh, got eventype as $eventyp


			[ -e "$inp"_majiq.bed ] && rm "$inp"_majiq.bed #for abs(ex1-ex2)>1
			[ -e "$inp"_majiq.csv ] && rm "$inp"_majiq.csv
			[ -e all_tx_events.csv ] && rm all_tx_events.csv
			[ -e majiq_events_sashimi2.bed ] && rm majiq_events_sashimi2.bed #for abs(ex1-ex2)>1
			[ -e majiq_events_sashimi2.csv ] && rm majiq_events_sashimi2.csv
			[ -e majiq_events_sashimi01.bed ] && rm majiq_events_sashimi01.bed #+ strand for abs(ex1-ex2)=0
			[ -e majiq_events_sashimi01.csv ] && rm majiq_events_sashimi01.csv
			[ -e majiq_events_sashimi02.bed ] && rm majiq_events_sashimi02.bed #- strand for abs(ex1-ex2)=0
			[ -e majiq_events_sashimi02.csv ] && rm majiq_events_sashimi02.csv
			[ -e majiq_events_progress2.txt ] && rm majiq_events_progress2.txt
			[ -e majiq_events_progress01.txt ] && rm majiq_events_progress01.txt
			[ -e majiq_events_progress02.txt ] && rm majiq_events_progress02.txt
			[ -e majiq_events_progress_all.txt ] && rm majiq_events_progressall.txt
			[ -e majiq_events_progress1.txt ] && rm majiq_events_progress1.txt

			flg=1
			if [ $flg -eq 1 ]
			then
					mkdir -p event_bedfiles
					[ "$(ls -A event_bedfiles/)" ] && rm event_bedfiles/*.*

					#Rscript pgp-c_TxEnsDB103_layeredV4_PEAKS.R $1 #$tx_lst
					Rscript TxEnsDB103_layeredV6.R $1 principal_txs.csv

			fi

			readarray -t csv_data < all_tx_events.csv
			csvi=0

			samples=$(ls event_bedfiles/temp_*.bed)

			for sample in $samples
			do
				#read csv entry
				csv_ln=${csv_data[$csvi]}
				((csvi=csvi+1))

				echo processing "$sample"
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
							#echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> majiq_events_sashimi2.bed
							echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> "$inp"_majiq.bed
							#echo $csv_ln >> majiq_events_sashimi2.csv
							echo $csv_ln >> "$inp"_majiq.csv
						else
							#echo ds $ds >> majiq_events_progress2.txt
							#echo us $us >> majiq_events_progress2.txt
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
									#echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> majiq_events_sashimi2.bed
									echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> "$inp"_majiq.bed
									#echo $csv_ln >> majiq_events_sashimi2.csv
									echo $csv_ln >> "$inp"_majiq.csv
								else
									echo ds 1 $ds >> majiq_events_progress2.txt
									echo us 1 $us >> majiq_events_progress2.txt

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
									#echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> majiq_events_sashimi2.bed
									echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> "$inp"_majiq.bed
									#echo $csv_ln >> majiq_events_sashimi2.csv
									echo $csv_ln >> "$inp"_majiq.csv
								else
									echo ds 2 $ds >> majiq_events_progress2.txt
									echo us 2 $us >> majiq_events_progress2.txt

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
								#echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> majiq_events_sashimi01.bed
								echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> "$inp"_majiq.bed

								#echo $csv_ln >> majiq_events_sashimi01.csv
								echo $csv_ln >> "$inp"_majiq.csv
							else
								echo ds $ds >> majiq_events_progress01.txt
								echo us $us >> majiq_events_progress01.txt

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
								#echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> majiq_events_sashimi02.bed
								echo $us $start $end $strnd $gene_name $TxID | awk 'BEGIN {OFS="\t"} {print $1,$14,$15,1,0,$16,$17,$18}' >> "$inp"_majiq.bed
								#echo $csv_ln >> majiq_events_sashimi02.csv
								echo $csv_ln >> "$inp"_majiq.csv
							else
								echo ds $ds >> majiq_events_progress02.txt
								echo us $us >> majiq_events_progress02.txt

								echo diff_exon_abs is $diff_exon_abs selected event $sample has event_st $event_st selected start $start event end $event_end selected end $end - please check
							fi
						fi
						#first check if star > event_start, then select upstream exon
					else
						echo ds $ds >>  majiq_events_progress_all.txt
						echo us $us >>  majiq_events_progress_all.txt
				fi


				#echo $lines1 | awk -v tx=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,tx}' >> res_skiptics/skiptics_uniq_sashimi.bed
				#echo $lines2 | awk -v tx=$txid 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,tx}' >> res_skiptics/skiptics_uniq_sashimi.bed

			done
		fi
	
		flag=$4

		if [ $flag -eq 1 ]
		then

			mkdir -p "$folder"/sashimi_plots/"$eventyp"
			rm -f "$folder"/sashimi_plots/"$eventyp"/*.pdf
			rm -f "$folder"/sashimi_plots/"$eventyp"/*.svg

			bed="$inp"_majiq.bed #majiq_events_sashimi_all.bed  #skiptics_sashimi.bed #"sam_junc1.txt" #"junctions.txt"
			readarray -t all_bed_data < $bed
			nrecrds=$(cat $bed | wc -l)
			echo read $nrecrds records
			csv="$inp"_majiq.csv #majiq_events_sashimi_all.csv #"sam_junc1.txt" #"junctions.txt"
			readarray -t all_csv_data < $csv
			#read peaks bed file
			peak_bed=$2
			readarray -t all_peak_data < $peak_bed

			#also read peaks csv file (to get AA seq)
			readarray -t all_peak_csv < $3

			i=0
			eventn=0
			rm -f temp_pdfs/*.pdf
			while [ $i -lt $nrecrds ]
			do
				#construct string for ggsashimi
				line1=${all_bed_data[$i]}
				AA=${all_peak_csv[$i]}
			((i=i+1))
				#read strand
				strnd=$(echo $line1 | awk '{print $6}')
				#also read TxID
				TxID=$(echo $line1 | awk '{print $8}')
				line=$(echo $line1 | awk '{print $1":"$2-50"-"$3+50}') #this is the actual event
				#get majiq event
				event=${all_csv_data[$eventn]}
				peak_event=${all_peak_data[$eventn]}
				((eventn=eventn+1))
				#fn=$(echo $line1 $line2 $line3 | awk '{print $7"-"$1"_"$2"-"$9}')
				#THIS IS TO KEEP DIFFERENT FILE NAMES FOR EVENTS WITH IDENTICAL COORDINATES AND GENE_NAMES
					gene_name=$(echo $event | awk 'BEGIN {FS=","} {print $5}')
							if [ $i -eq 1 ]
							then
								temp_gene=$(echo $event | awk 'BEGIN {FS=","} {print $5}')
								trackj=1
							elif [ "$temp_gene" == "$gene_name" ]
							then
								((trackj=trackj+1))
							else
								temp_gene=$(echo $event | awk 'BEGIN {FS=","} {print $5}')
								trackj=1
							fi
				#SHOULD NOT HAVE THIS BUT NEED TO KEEP AS PER REQUEST
				fn=$(echo $event | awk -v tr=$trackj 'BEGIN {FS=","} {print $5"-"$1"-"$2"-"$3"-"tr}') #also gene_id to avoid same file names for repeated events
				#string for majiq event
				#echo event $event
				chr_name=$(echo $event | awk 'BEGIN {FS=","} {print $1}')
				start=$(echo $event | awk 'BEGIN {FS=","} {print $2}')
				end=$(echo $event | awk 'BEGIN {FS=","} {print $3}')

				#majiq_event=$(echo $chr_name $start $end | awk '{print $1"-"$2"-"$3}')
				majiq_event=$(echo $chr_name $start $end | awk '{print $1"-"$2"-"$3}')

				#ALSO ADD PEAKS EVENT
				pchr_name=$(echo $peak_event | awk 'BEGIN {FS="\t"} {print $1}')
				pstart=$(echo $peak_event | awk 'BEGIN {FS="\t"} {print $2}')
				pend=$(echo $peak_event | awk 'BEGIN {FS="\t"} {print $3}')


				pstart1=$(echo "$peak_event" | awk 'BEGIN {FS="\t"} {print $2}')
				pend1=$(echo "$peak_event" | awk 'BEGIN {FS="\t"} {print $3}')

				#Also get length of event range
				sz1=$(($pend1-$pstart1))
				#echo sz1 is "$sz1"


				peaks_event=$(echo $pchr_name $pstart $pend | awk '{print $1"-"$2"-"$3}')

				#Also get AA seq
				AAseq1=$(echo "$AA" | awk 'BEGIN {FS=","} {print $1}')
				#now add new line if string is larger than 100 characters
				length_AA_flg=0
				if [ "${#AAseq1}" -gt 100 ]
				then
					AAseq2=$(echo "$AAseq1" | cut -c 0-100 )
					AAseq3=$(echo "$AAseq1" | cut -c 101-"${#AAseq1}" )
					#now concatenate
					AAseq="$AAseq2\n$AAseq3"
					length_AA_flg=1
				else
					AAseq="$AAseq1"
				fi
				#also get nt seq
				NTseq1=$(echo "$AA" | awk 'BEGIN {FS=","} {print $7}')
				length_nt_flg=0
				if [ "$((${#NTseq1}-2))" -gt 100 ]
				then
					NTseq2=$(echo "$NTseq1" | cut -c 1-100 )
					NTseq3=$(echo "$NTseq1" | cut -c 101-"$((${#NTseq1}))")
					#now concatenate
					NTseq="$NTseq2""\n""$NTseq3"
					length_nt_flg=1
				else
					NTseq="$NTseq1"
				fi
				#NOw concatenate gene_name, NTseq and AA seq - here -4 is due to _num and two characters for new line
				if [ $length_nt_flg -eq 1 ] && [ $length_AA_flg -eq 0 ]
				then
					#Added additional 2 empty lines to squeeze x-axis
					#gene_name1=$(echo "$gene_name" "$AAseq" "${#AAseq}" "$NTseq"  "$((${#NTseq}-4))" | awk '{print "GeneName:"$1"\nAA:"$2"("$3")\nNT:"$4"("$5")"}')
					#gene_name1=$(echo "$gene_name" "$AAseq" "${#AAseq}" "$NTseq"  "$((${#NTseq}-4))" | awk '{print "Gene: "$1"\nAA:"$2"("$3")\nNT:"$4"("$5")\n\n"}')
					#gene_name1=$(echo "$gene_name" "$AAseq" "${#AAseq}" "$NTseq"  "$((${#NTseq}-4))" | awk '{print "Gene: "$1" - AA:"$2"("$3") - NT:"$4"("$5")"}')
					gene_name1=$(echo "$gene_name" "$AAseq" "${#AAseq}" "$NTseq"  "$((${#NTseq}-4))" | awk '{print "Gene: "$1}') # - AA:"$2"("$3") - NT:"$4"("$5")"}')
				elif [ $length_nt_flg -eq 1 ] && [ $length_AA_flg -eq 1 ]
				then
					#gene_name1=$(echo "$gene_name" "$AAseq" "$((${#AAseq}-2))" "$NTseq"  "$((${#NTseq}-4))" | awk '{print "GeneName:"$1"\nAA:"$2"("$3")\nNT:"$4"("$5")"}')
					#gene_name1=$(echo "$gene_name" "$AAseq" "$((${#AAseq}-2))" "$NTseq"  "$((${#NTseq}-4))" | awk '{print "Gene: "$1"\nAA:"$2"("$3")\nNT:"$4"("$5")\n\n"}')
					#gene_name1=$(echo "$gene_name" "$AAseq" "$((${#AAseq}-2))" "$NTseq"  "$((${#NTseq}-4))" | awk '{print "Gene: "$1" - AA:"$2"("$3") - NT:"$4"("$5")"}')
					gene_name1=$(echo "$gene_name" "$AAseq" "$((${#AAseq}-2))" "$NTseq"  "$((${#NTseq}-4))" | awk '{print "Gene: "$1}') #'" - AA:"$2"("$3") - NT:"$4"("$5")"}')

				elif [ $length_nt_flg -eq 0 ] && [ $length_AA_flg -eq 1 ]
				then
					#gene_name1=$(echo "$gene_name" "$AAseq" "$((${#AAseq}-2))" "$NTseq"  "$((${#NTseq}-2))" | awk '{print "GeneName:"$1"\nAA:"$2"("$3")\nNT:"$4"("$5")"}')
					#gene_name1=$(echo "$gene_name" "$AAseq" "$((${#AAseq}-2))" "$NTseq"  "$((${#NTseq}-2))" | awk '{print "Gene: "$1"\nAA:"$2"("$3")\nNT:"$4"("$5")\n\n"}')
					#gene_name1=$(echo "$gene_name" "$AAseq" "$((${#AAseq}-2))" "$NTseq"  "$((${#NTseq}-2))" | awk '{print "Gene: "$1" - AA:"$2"("$3") - NT:"$4"("$5")"}')
					gene_name1=$(echo "$gene_name" "$AAseq" "$((${#AAseq}-2))" "$NTseq"  "$((${#NTseq}-2))" | awk '{print "Gene: "$1}') #'" - AA:"$2"("$3") - NT:"$4"("$5")"}')
				else
					#gene_name1=$(echo "$gene_name" "$AAseq" "${#AAseq}" "$NTseq"  "$((${#NTseq}-2))" | awk '{print "GeneName:"$1"\nAA:"$2"("$3")\nNT:"$4"("$5")"}')
					#gene_name1=$(echo "$gene_name" "$AAseq" "${#AAseq}" "$NTseq"  "$((${#NTseq}-2))" | awk '{print "Gene: "$1"\nAA:"$2"("$3")\nNT:"$4"("$5")\n\n"}')
					#gene_name1=$(echo "$gene_name" "$AAseq" "${#AAseq}" "$NTseq"  "$((${#NTseq}-2))" | awk '{print "Gene: "$1" - AA:"$2"("$3") - NT:"$4"("$5")"}')
					gene_name1=$(echo "$gene_name" "$AAseq" "${#AAseq}" "$NTseq"  "$((${#NTseq}-2))" | awk '{print "Gene: "$1}') #'" - AA:"$2"("$3") - NT:"$4"("$5")"}')
				fi

				########END PEAKS EVENT

				exon1=0
				exon2=0
				#also get actual event identified
				event_identified=$(echo $line1 | awk '{print "None ""-"$2"-"$3}')
				echo processing event num $eventn and event $fn #and majiq event is $majiq_event
				#./ggsashimi_txV2.py -b all_bams.tsv -c "$line" -g /Volumes/SYEDSHAH/MichaelLab/ref_genome/Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -PeaksAA "$AAseq" -PeaksFlg 1 -PeaksTx "$peaks_event" -GeneName "$gene_name1" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o temp_pdfs/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
			#	./ggsashimi_txV3.py -A "median_j" -b all_bams.tsv -c "$line" -g /Volumes/SYEDSHAH/MichaelLab/ref_genome/Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -PeaksAA "$sz1" -PeaksFlg 1 -PeaksTx "$peaks_event" -GeneName "$gene_name1" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o temp_pdfs/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
			#	##rsvg-convert -f pdf -o temp_pdfs/"$fn".pdf temp_pdfs/"$fn".svg
				./ggsashimi_txV3.py -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -PeaksAA "$sz1" -PeaksFlg 1 -PeaksTx "$peaks_event" -GeneName "$gene_name1" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o "$folder"/sashimi_plots/"$eventyp"/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
				##rsvg-convert -f pdf -o "$folder"/sashimi_plots/"$eventyp"/"$fn".pdf "$folder"/sashimi_plots/"$eventyp"/"$fn".svg

				[ -e "$folder"/sashimi_plots/"$eventyp"/"$fn".svg ] && rm "$folder"/sashimi_plots/"$eventyp"/"$fn".svg

			done
			#now merge all pdf's
			#fname=$(echo $2 |cut -d'.' -f1)
			#pdfunite temp_pdfs/*.pdf sashimi_plots/$fname.pdf
			#pdfunite *.pdf all_events_sashimi.pdf
			#rm -f temp_pdfs/*.pdf
			#pdfunite "$folder"/sashimi_plots/"$eventyp"/*.pdf "$folder"/sashimi_plots/"$eventyp"/all_"$eventyp".pdf
			python merge_sashimis.py "$folder"/sashimi_plots/"$eventyp"/
		fi


		flag=$4
		#THIS ADDS TWO EXONS IN PEAKS EVENTS SHARING EXON_CE OR CE_EXON
		#USES PeaksFlg=2 for ggsashimi_txV2.py
		#THIS SECTION DEALS WITH UEX_CE, CE_DEX for CE events (flag=0 or 1) and UEX_DEX for SKIPTICS (flag=2)
		if [ $flag -eq 2 ]
		then

			mkdir -p "$folder"/sashimi_plots/"$eventyp"
			rm -f "$folder"/sashimi_plots/"$eventyp"/*.pdf
			rm -f "$folder"/sashimi_plots/"$eventyp"/*.svg

			echo current wkd is "$pwd"
			echo came in "$flag -eq 2" area with bed file $inp.bed
			bed="$inp"_majiq.bed #majiq_events_sashimi_all.bed  #skiptics_sashimi.bed #"sam_junc1.txt" #"junctions.txt"
			readarray -t all_bed_data < $bed
			nrecrds=$(cat $bed | wc -l)
			echo read $nrecrds records
			csv="$inp"_majiq.csv #majiq_events_sashimi_all.csv #"sam_junc1.txt" #"junctions.txt"
			readarray -t all_csv_data < $csv
			#ALSO READ PEAKS EVENTS BED FILE, HERE IT HAS 2 EXONS PER EVENT
			peak_bed=$2
			readarray -t all_peak_data < $peak_bed

			#also read peaks csv file (to get AA seq)
			readarray -t all_peak_csv < $3

			i=0
			eventn=0
			rm -f temp_pdfs/*.pdf
			while [ $i -lt $nrecrds ]
			do
				#construct string for ggsashimi
				line1=${all_bed_data[$i]}
				#get majiq event
				event=${all_csv_data[$i]}
				AA=${all_peak_csv[$i]}
			((i=i+1))
				#read strand
				strnd=$(echo $line1 | awk '{print $6}')
				#also read TxID
				TxID=$(echo $line1 | awk '{print $8}')
				line=$(echo $line1 | awk '{print $1":"$2-50"-"$3+50}') #this is the actual event
				#echo eventn is $eventn
				peak_event=${all_peak_data[$eventn]}
				((eventn=eventn+1))
				#echo now eventn is $eventn
				peak_event1=${all_peak_data[$eventn]}
				((eventn=eventn+1))
				#fn=$(echo $line1 $line2 $line3 | awk '{print $7"-"$1"_"$2"-"$9}')
				#THIS IS TO KEEP DIFFERENT FILE NAMES FOR EVENTS WITH IDENTICAL COORDINATES AND GENE_NAMES
					gene_name=$(echo $event | awk 'BEGIN {FS=","} {print $5}')
							if [ $i -eq 1 ]
							then
								temp_gene=$(echo $event | awk 'BEGIN {FS=","} {print $5}')
								trackj=1
							elif [ "$temp_gene" == "$gene_name" ]
							then
								((trackj=trackj+1))
							else
								temp_gene=$(echo $event | awk 'BEGIN {FS=","} {print $5}')
								trackj=1
							fi
				#SHOULD NOT HAVE THIS BUT NEED TO KEEP AS PER REQUEST
				fn=$(echo $event | awk -v tr=$trackj 'BEGIN {FS=","} {print $5"-"$1"-"$2"-"$3"-"tr}') #also gene_id to avoid same file names for repeated events
				#string for majiq event
				#echo event $event
				chr_name=$(echo $event | awk 'BEGIN {FS=","} {print $1}')
				start=$(echo $event | awk 'BEGIN {FS=","} {print $2}')
				end=$(echo $event | awk 'BEGIN {FS=","} {print $3}')

				#majiq_event=$(echo $chr_name $start $end | awk '{print $1"-"$2"-"$3}')
				majiq_event=$(echo $chr_name $start $end | awk '{print $1"-"$2"-"$3}')

				#ALSO ADD 2 PEAKS EVENT - one for us/ds and other for ce
				pchr_name=$(echo "$peak_event" | awk 'BEGIN {FS="\t"} {print $1}')
				pstart=$(echo "$peak_event" | awk 'BEGIN {FS="\t"} {print $2}')
				pend=$(echo "$peak_event" | awk 'BEGIN {FS="\t"} {print $3}')

				#pchr_name=$(echo $peak_event | awk 'BEGIN {FS="\t"} {print $1}')
				#echo peak_event is $peak_event
				#echo peak_event1 is $peak_event1
				pchr_name1=$(echo "$peak_event1" | awk 'BEGIN {FS="\t"} {print $1}')
				pstart1=$(echo "$peak_event1" | awk 'BEGIN {FS="\t"} {print $2}')
				pend1=$(echo "$peak_event1" | awk 'BEGIN {FS="\t"} {print $3}')

				#echo pchr_name1 is $pchr_name1 pstart1 is $pstart1 and pend1 is $pend1


				#peaks_event=$(echo $pchr_name $pstart $pend | awk '{print $1"-"$2"-"$3}')
				peaks_event=$(echo "$pchr_name" "$pstart" "$pend" "$pstart1" "$pend1" | awk '{print $1"-"$2"-"$3"-"$4"-"$5}')

				#Also get AA seq
				#AAseq=$(echo "$AA" | awk 'BEGIN {FS=","} {print $1}')
				#also get nt seq
				NTseq1=$(echo "$AA" | awk 'BEGIN {FS=","} {print $7}')
				length_nt_flg=0
				if [ "$((${#NTseq1}-2))" -gt 100 ]
				then
					NTseq2=$(echo "$NTseq1" | cut -c 1-100 )
					NTseq3=$(echo "$NTseq1" | cut -c 101-"$((${#NTseq1}))")
					#now concatenate
					NTseq="$NTseq2""\n""$NTseq3"
					length_nt_flg=1
				else
					NTseq="$NTseq1"
				fi

				AAseq1=$(echo "$AA" | awk 'BEGIN {FS=","} {print $1}')
				length_AA_flg=0
				if [ "${#AAseq1}" -gt 100 ]
				then
					AAseq2=$(echo "$AAseq1" | cut -c 0-100 )
					AAseq3=$(echo "$AAseq1" | cut -c 101-"${#AAseq1}" )
					#now concatenate
					AAseq="$AAseq2\n$AAseq3"
					length_AA_flg=1
				else
					AAseq="$AAseq1"
				fi

				#NOw concatenate gene_name and NTseq
				###########
				#gene_name1=$(echo "$gene_name" "$NTseq" | awk '{print $1"\n"$2}')
				if [ $length_nt_flg -eq 1 ] && [ $length_AA_flg -eq 0 ]
				then
					#gene_name1=$(echo "$gene_name" "$AAseq" "${#AAseq}" "$NTseq"  "$((${#NTseq}-4))" | awk '{print "Gene: "$1" - AA:"$2"("$3") - NT:"$4"("$5")"}')
					gene_name1=$(echo "$gene_name" "$AAseq" "${#AAseq}" "$NTseq"  "$((${#NTseq}-4))" | awk '{print "Gene: "$1}') #'" - AA:"$2"("$3") - NT:"$4"("$5")"}')
				elif [ $length_nt_flg -eq 1 ] && [ $length_AA_flg -eq 1 ]
				then
					#gene_name1=$(echo "$gene_name" "$AAseq" "$((${#AAseq}-2))" "$NTseq"  "$((${#NTseq}-4))" | awk '{print "Gene: "$1" - AA:"$2"("$3") - NT:"$4"("$5")"}')
					gene_name1=$(echo "$gene_name" "$AAseq" "$((${#AAseq}-2))" "$NTseq"  "$((${#NTseq}-4))" | awk '{print "Gene: "$1}') #'" - AA:"$2"("$3") - NT:"$4"("$5")"}')

				elif [ $length_nt_flg -eq 0 ] && [ $length_AA_flg -eq 1 ]
				then
					#gene_name1=$(echo "$gene_name" "$AAseq" "$((${#AAseq}-2))" "$NTseq"  "$((${#NTseq}-2))" | awk '{print "Gene: "$1" - AA:"$2"("$3") - NT:"$4"("$5")"}')
					gene_name1=$(echo "$gene_name" "$AAseq" "$((${#AAseq}-2))" "$NTseq"  "$((${#NTseq}-2))" | awk '{print "Gene: "$1}') #'" - AA:"$2"("$3") - NT:"$4"("$5")"}')
				else
					#gene_name1=$(echo "$gene_name" "$AAseq" "${#AAseq}" "$NTseq"  "$((${#NTseq}-2))" | awk '{print "Gene: "$1" - AA:"$2"("$3") - NT:"$4"("$5")"}')
					gene_name1=$(echo "$gene_name" "$AAseq" "${#AAseq}" "$NTseq"  "$((${#NTseq}-2))" | awk '{print "Gene: "$1}') #'" - AA:"$2"("$3") - NT:"$4"("$5")"}')
				fi

				############
				#now addind upexon and ce part or ce and dnexon part separately

				if [ "$5" -eq 0 ]
				then
					NTseque1=$(echo "$AA" | awk 'BEGIN {FS=","} {print $8}')
					if [ "$((${#NTseque1}-2))" -gt 100 ]
					then
						NTseque2=$(echo "$NTseque1" | cut -c 1-100 )
						NTseque3=$(echo "$NTseque1" | cut -c 101-"$((${#NTseque1}))")
						#now concatenate
						NTseque="$NTseque2""\n""$NTseque3"
					else
						NTseque="$NTseque1"
					fi

					NTseqce1=$(echo "$AA" | awk 'BEGIN {FS=","} {print $9}')
					if [ "$((${#NTseqce1}-2))" -gt 100 ]
					then
						NTseqce2=$(echo "$NTseqce1" | cut -c 1-100 )
						NTseqce3=$(echo "$NTseqce1" | cut -c 101-"$((${#NTseqce1}))")
						#now concatenate
						NTseqce="$NTseqce2""\n""$NTseqce3"
					else
						NTseqce="$NTseqce1"
					fi

					AAseque1=$(echo "$AA" | awk 'BEGIN {FS=","} {print $10}')
					if [ "${#AAseque1}" -gt 100 ]
					then
						AAseque2=$(echo "$AAseque1" | cut -c 0-100 )
						AAseque3=$(echo "$AAseque1" | cut -c 101-"${#AAseque1}" )
						#now concatenate
						AAseque="$AAseque2\n$AAseque3"
					else
						AAseque="$AAseque1"
					fi

					AAseqce1=$(echo "$AA" | awk 'BEGIN {FS=","} {print $11}')
					if [ "${#AAseqce1}" -gt 100 ]
					then
						AAseqce2=$(echo "$AAseqce1" | cut -c 0-100 )
						AAseqce3=$(echo "$AAseqce1" | cut -c 101-"${#AAseqce1}" )
						#now concatenate
						AAseqce="$AAseqce2\n$AAseqce3"
					else
						AAseqce="$AAseqce1"
					fi

					#gene_name2=$(echo "$gene_name" "$AAseque" "${#AAseque}" "$NTseque" "$((${#NTseque}-2))" "$AAseqce" "${#AAseqce}" "$NTseqce" "$((${#NTseqce}-2))" "$AAseq" "${#AAseq}" "$NTseq" "$((${#NTseq}-2))" | awk '{print "Gene: "$1"\nUEX:"$2"("$3")-"$4"("$5")\nCE:"$6"("$7")-"$8"("$9")\nAA_all:"$10"("$11")\nNt_all:"$12"("$13")"}')
					#removing AA_All
					#gene_name2=$(echo "$gene_name" "$AAseque" "${#AAseque}" "$NTseque" "$((${#NTseque}-2))" "$AAseqce" "${#AAseqce}" "$NTseqce" "$((${#NTseqce}-2))" "$AAseq" "${#AAseq}" "$NTseq" "$((${#NTseq}-2))" | awk '{print "Gene: "$1"\nUEX:"$2"("$3")-"$4"("$5")\nCE:"$6"("$7")-"$8"("$9")"}') #'"\nAA_all:"$10"("$11")\nNt_all:"$12"("$13")"}')
					#gene_name2=$(echo "$gene_name" "$AAseque" "${#AAseque}" "$NTseque" "$((${#NTseque}-2))" "$AAseqce" "${#AAseqce}" "$NTseqce" "$((${#NTseqce}-2))" "$AAseq" "${#AAseq}" "$NTseq" "$((${#NTseq}-2))" | awk '{print "Gene: "$1" - UEX:"$2"("$3")-"$4"("$5") - CE:"$6"("$7")-"$8"("$9")"}') #'"\nAA_all:"$10"("$11")\nNt_all:"$12"("$13")"}')
					gene_name2=$(echo "$gene_name" "$AAseque" "${#AAseque}" "$NTseque" "$((${#NTseque}-2))" "$AAseqce" "${#AAseqce}" "$NTseqce" "$((${#NTseqce}-2))" "$AAseq" "${#AAseq}" "$NTseq" "$((${#NTseq}-2))" | awk '{print "Gene: "$1}') #'" - UEX:"$2"("$3")-"$4"("$5") - CE:"$6"("$7")-"$8"("$9")"}') #'"\nAA_all:"$10"("$11")\nNt_all:"$12"("$13")"}')
					#Also get length of each range
					len1=$(echo "$AA" | awk 'BEGIN {FS=","} {print $13}')
					len2=$(echo "$AA" | awk 'BEGIN {FS=","} {print $14}')
					sz1=$(($len2-$len1))
					len1=$(echo "$AA" | awk 'BEGIN {FS=","} {print $15}')
					len2=$(echo "$AA" | awk 'BEGIN {FS=","} {print $16}')
					sz2=$(($len2-$len1))

					#Also add UEX_CE coordinates
					#JuncSpanning_Coord=$(echo "$AA" "$sz1" "$sz2" | awk 'BEGIN {FS=","} {print $12":"$13"-"$14":"$17"\n"$12":"$15"-"$16":"$18}')
					JuncSpanning_Coord=$(echo "$AA", $sz1, $sz2 | awk 'BEGIN {FS=","} {print $12":"$13"-"$14":"$17"\n"$12":"$15"-"$16":"$18}')
				elif [ "$5" -eq 1 ]
				then
					NTseque=$(echo "$AA" | awk 'BEGIN {FS=","} {print $8}')
					NTseqce=$(echo "$AA" | awk 'BEGIN {FS=","} {print $9}')
					AAseque=$(echo "$AA" | awk 'BEGIN {FS=","} {print $10}')
					AAseqce=$(echo "$AA" | awk 'BEGIN {FS=","} {print $11}')
					#gene_name2=$(echo "$gene_name" "$NTseque" "$NTseqce" "$NTseq" | awk '{print "GeneName:"$1"\nNt_ce:"$2"\nNt_dex:"$3"\nNt_all:"$4}')

					#gene_name2=$(echo "$gene_name" "$AAseque" "$NTseque" "$AAseqce" "$NTseqce" "$NTseq" "$AAseq" "${#AAseq}" | awk '{print "GeneName:"$1"\nCE:"$2"-"$3"\nDEX:"$4"-"$5"\nNt_all:"$6"\nAA_all:"$7"("$8")"}')
					#gene_name2=$(echo "$gene_name" "$AAseque" "${#AAseque}" "$NTseque" "$((${#NTseque}-2))" "$AAseqce" "${#AAseqce}" "$NTseqce" "$((${#NTseqce}-2))" "$AAseq" "${#AAseq}" "$NTseq" "$((${#NTseq}-2))" | awk '{print "Gene: "$1"\nCE:"$2"("$3")-"$4"("$5")\nDEX:"$6"("$7")-"$8"("$9")\nAA_all:"$10"("$11")\nNt_all:"$12"("$13")"}')
					#removing AA_All
					#gene_name2=$(echo "$gene_name" "$AAseque" "${#AAseque}" "$NTseque" "$((${#NTseque}-2))" "$AAseqce" "${#AAseqce}" "$NTseqce" "$((${#NTseqce}-2))" "$AAseq" "${#AAseq}" "$NTseq" "$((${#NTseq}-2))" | awk '{print "Gene: "$1"\nCE:"$2"("$3")-"$4"("$5")\nDEX:"$6"("$7")-"$8"("$9")"}') #'"\nAA_all:"$10"("$11")\nNt_all:"$12"("$13")"}')
					#gene_name2=$(echo "$gene_name" "$AAseque" "${#AAseque}" "$NTseque" "$((${#NTseque}-2))" "$AAseqce" "${#AAseqce}" "$NTseqce" "$((${#NTseqce}-2))" "$AAseq" "${#AAseq}" "$NTseq" "$((${#NTseq}-2))" | awk '{print "Gene: "$1" - CE:"$2"("$3")-"$4"("$5") - DEX:"$6"("$7")-"$8"("$9")"}') #'"\nAA_all:"$10"("$11")\nNt_all:"$12"("$13")"}')
					gene_name2=$(echo "$gene_name" "$AAseque" "${#AAseque}" "$NTseque" "$((${#NTseque}-2))" "$AAseqce" "${#AAseqce}" "$NTseqce" "$((${#NTseqce}-2))" "$AAseq" "${#AAseq}" "$NTseq" "$((${#NTseq}-2))" | awk '{print "Gene: "$1}') #'" - CE:"$2"("$3")-"$4"("$5") - DEX:"$6"("$7")-"$8"("$9")"}') #'"\nAA_all:"$10"("$11")\nNt_all:"$12"("$13")"}')
					#Also add UEX_CE coordinates
					#Also get length of each range
					len1=$(echo "$AA" | awk 'BEGIN {FS=","} {print $13}')
					len2=$(echo "$AA" | awk 'BEGIN {FS=","} {print $14}')
					sz1=$(($len2-$len1))
					len1=$(echo "$AA" | awk 'BEGIN {FS=","} {print $15}')
					len2=$(echo "$AA" | awk 'BEGIN {FS=","} {print $16}')
					sz2=$(($len2-$len1))
					#JuncSpanning_Coord=$(echo "$AA", $sz1, $sz2 | awk 'BEGIN {FS=","} {print $12":"$13"-"$14"\n"$12":"$15"-"$16"\n"$17"-"$18}')
					JuncSpanning_Coord=$(echo "$AA", $sz1, $sz2 | awk  'BEGIN {FS=","} {print $12":"$13"-"$14":"$17"\n"$12":"$15"-"$16":"$18}')
					#JuncSpanning_Coord=$(echo "$AA" | awk 'BEGIN {FS=","} {print $12":"$13"-"$14"\n"$12":"$15"-"$16}')
				else
					NTseque=$(echo "$AA" | awk 'BEGIN {FS=","} {print $8}')
					NTseqce=$(echo "$AA" | awk 'BEGIN {FS=","} {print $9}')
					AAseque=$(echo "$AA" | awk 'BEGIN {FS=","} {print $10}')
					AAseqce=$(echo "$AA" | awk 'BEGIN {FS=","} {print $11}')

					#gene_name2=$(echo "$gene_name" "$AAseque" "${#AAseque}" "$NTseque" "$((${#NTseque}-2))" "$AAseqce" "${#AAseqce}" "$NTseqce" "$((${#NTseqce}-2))" "$AAseq" "${#AAseq}" "$NTseq" "$((${#NTseq}-2))" | awk '{print "Gene: "$1"\nUEX:"$2"("$3")-"$4"("$5")\nDEX:"$6"("$7")-"$8"("$9")\nAA_all:"$10"("$11")\nNt_all:"$12"("$13")"}')
					#removing AA_All
					#gene_name2=$(echo "$gene_name" "$AAseque" "${#AAseque}" "$NTseque" "$((${#NTseque}-2))" "$AAseqce" "${#AAseqce}" "$NTseqce" "$((${#NTseqce}-2))" "$AAseq" "${#AAseq}" "$NTseq" "$((${#NTseq}-2))" | awk '{print "Gene: "$1" - UEX:"$2"("$3")-"$4"("$5") - DEX:"$6"("$7")-"$8"("$9")"}') #'"\nAA_all:"$10"("$11")\nNt_all:"$12"("$13")"}')
					gene_name2=$(echo "$gene_name" "$AAseque" "${#AAseque}" "$NTseque" "$((${#NTseque}-2))" "$AAseqce" "${#AAseqce}" "$NTseqce" "$((${#NTseqce}-2))" "$AAseq" "${#AAseq}" "$NTseq" "$((${#NTseq}-2))" | awk '{print "Gene: "$1}') #'" - UEX:"$2"("$3")-"$4"("$5") - DEX:"$6"("$7")-"$8"("$9")"}') #'"\nAA_all:"$10"("$11")\nNt_all:"$12"("$13")"}')
					#Also get length of each range
					len1=$(echo "$AA" | awk 'BEGIN {FS=","} {print $13}')
					len2=$(echo "$AA" | awk 'BEGIN {FS=","} {print $14}')
					sz1=$(($len2-$len1))
					len1=$(echo "$AA" | awk 'BEGIN {FS=","} {print $15}')
					len2=$(echo "$AA" | awk 'BEGIN {FS=","} {print $16}')
					sz2=$(($len2-$len1))

					#Also add UEX_CE coordinates
					#JuncSpanning_Coord=$(echo "$AA" "$sz1" "$sz2" | awk 'BEGIN {FS=","} {print $12":"$13"-"$14":"$17"\n"$12":"$15"-"$16":"$18}')
					JuncSpanning_Coord=$(echo "$AA", $sz1, $sz2 | awk 'BEGIN {FS=","} {print $12":"$13"-"$14":"$17"\n"$12":"$15"-"$16":"$18}')

				fi


				########END PEAKS EVENT

				exon1=0
				exon2=0
				#also get actual event identified
				event_identified=$(echo $line1 | awk '{print "None ""-"$2"-"$3}')
				echo processing event num $i and event $fn #and majiq event is $majiq_event and peaks event is $peaks_event

				###should remove this
				#./ggsashimi_txV3.py -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -PeaksAA "$sz1" -PeaksFlg 1 -PeaksTx "$peaks_event" -GeneName "$gene_name1" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o "$folder"/sashimi_plots/"$eventyp"/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
				####

				./ggsashimi_txV3.py -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -PeaksJuncSpanningCoord "$JuncSpanning_Coord" -PeaksAA "$AAseq" -PeaksFlg 2 -PeaksTx "$peaks_event" -GeneName "$gene_name2" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o "$folder"/sashimi_plots/"$eventyp"/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
				##rsvg-convert -f pdf -o "$folder"/sashimi_plots/"$eventyp"/"$fn".pdf "$folder"/sashimi_plots/"$eventyp"/"$fn".svg

				[ -e "$folder"/sashimi_plots/"$eventyp"/"$fn".svg ] && rm "$folder"/sashimi_plots/"$eventyp"/"$fn".svg

			done
			#now merge all pdf's
			#fname=$(echo $2 |cut -d'.' -f1)
			#pdfunite temp_pdfs/*.pdf sashimi_plots/$fname.pdf
			#pdfunite *.pdf all_events_sashimi.pdf
			#rm -f temp_pdfs/*.pdf
			#pdfunite "$folder"/sashimi_plots/"$eventyp"/*.pdf "$folder"/sashimi_plots/"$eventyp"/all_$fnam"$eventyp".pdf
			python merge_sashimis.py "$folder"/sashimi_plots/"$eventyp"/
		fi
	fi

fi