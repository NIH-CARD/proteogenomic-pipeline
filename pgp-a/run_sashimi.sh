#!/bin/bash
#declare encodings
#-*- coding: utf-8 -*-

#THIS SCRIPT CONTAINS SASHIMI PLOT CODE FOR
#0. PLEASE NOTE THAT FOLLOWING 2 FILES (OR SOFT LINKS) SHOULD BE IN CURRENT FOLDER
	#01: ggsashimi_txV4.py
	#02: Homo_sapiens.GRCh38.103.chr.sorted_new.gtf
#1. Skiptic Events
#2. ALL MAJIQ EVENTS
#3. CE (INCLUDING INCLUSION, EXTENSION AND IR) events
if [ $1 -eq 1 ]
then
	[ "$(ls -A res_skiptics/sashimi_plots/)" ] && rm res_skiptics/sashimi_plots/*.*
	mkdir -p res_skiptics/sashimi_plots
bed=res_skiptics/skiptics_uniq_sashimi.bed  #skiptics_sashimi.bed #"sam_junc1.txt" #"junctions.txt"
readarray -t all_bed_data < $bed
nrecrds=$(cat $bed | wc -l)
nrecrdst=$((nrecrds/2))
echo read $nrecrdst records
csv=res_skiptics/IGV_unique_skiptics_translated.csv #"sam_junc1.txt" #"junctions.txt"
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
	echo processing event num $eventn event $fn and majiq event is $majiq_event and event identified is $event_identified

	./ggsashimi_txV4.py -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -GeneName "$gene_name" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o res_skiptics/sashimi_plots/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
             rsvg-convert -f pdf -o res_skiptics/sashimi_plots/"$fn".pdf res_skiptics/sashimi_plots/"$fn".svg
        #[ -e res_skiptics/sashimi_plots/"$fn".svg ] && rm res_skiptics/sashimi_plots/"$fn".svg

######
#./ggsashimi_tx.py -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o "$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
########


done
#now merge all pdf's

pdfunite res_skiptics/sashimi_plots/*.pdf res_skiptics/all_skiptics_sashimi.pdf

fi

######THIS IS FOR ALL MAJIQ EVENTS

if [ $1 -eq 2 ]
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
	event_identified=$(echo $line1 | awk '{print "None ""-"$2"-"$3}')


	echo processing event num $eventn event $fn and majiq event is $majiq_event

#	./ggsashimi_txV4.py -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -GeneName "$gene_name" -MajiqStrnd "$strnd" -ORIG 0 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o all_events_sashimi/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
#             rsvg-convert --unlimited  -f pdf -o all_events_sashimi/"$fn".pdf all_events_sashimi/"$fn".svg
	./ggsashimi_txV4.py -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -GeneName "$gene_name" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o all_events_sashimi/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
					              rsvg-convert --unlimited  -f pdf -o all_events_sashimi/"$fn".pdf all_events_sashimi/"$fn".svg

        #[ -e all_events_sashimi/"$fn".svg ] && rm all_events_sashimi/"$fn".svg

done
#now merge all pdf's

pdfunite all_events_sashimi/*.pdf all_events_sashimi/all_events_sashimi.pdf
fi


#THIS SECTION IS FOR CE_INCLUSION EVENTS
#WILL MERGE INCLUSION AND EXTENSION EVENTS
if [ $1 -eq 3 ]
then
	[ "$(ls -A res_ce_all/ce_incl_sashimi_plots/)" ] && rm res_ce_all/ce_incl_sashimi_plots/*.*
	mkdir -p res_ce_all/ce_incl_sashimi_plots
#cat ce_inclusion_coord_uniq_sashimi.bed ce_extension_coord_uniq_sashimi.bed > all_ce.bed
#cat IGV_unique_ce_inclusion.csv IGV_unique_ce_extension.csv > all_ce.csv
bed=res_ce_all/ce_inclusion_coord_uniq_sashimi.bed #"sam_junc1.txt" #"junctions.txt"
readarray -t all_bed_data < $bed
nrecrds=$(cat $bed | wc -l)
nrecrdst=$((nrecrds/3))
echo read $nrecrdst records
csv=res_ce_all/IGV_unique_ce_inclusion.csv #"sam_junc1.txt" #"junctions.txt"
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
	echo processing event num $eventn event $fn and line is $line


	./ggsashimi_txV4.py  -PGPTx "$event_identified" -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -FullTitle "$fulltitle" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o res_ce_all/ce_incl_sashimi_plots/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
             rsvg-convert -f pdf -o res_ce_all/ce_incl_sashimi_plots/"$fn".pdf res_ce_all/ce_incl_sashimi_plots/"$fn".svg
        #[ -e res_ce_all/ce_incl_sashimi_plots/"$fn".svg ] && rm res_ce_all/ce_incl_sashimi_plots/"$fn".svg

done
#now merge all pdf's

pdfunite res_ce_all/ce_incl_sashimi_plots/*.pdf res_ce_all/ce_incl_all_sashimi_plots.pdf

fi
#CE_EXTENSION

if [ $1 -eq 4 ]
then
	[ "$(ls -A res_ce_all/ce_ext_sashimi_plots/)" ] && rm res_ce_all/ce_ext_sashimi_plots/*.*
	mkdir -p res_ce_all/ce_ext_sashimi_plots
bed=res_ce_all/ce_extension_coord_uniq_sashimi.bed #"sam_junc1.txt" #"junctions.txt"
readarray -t all_bed_data < $bed
nrecrds=$(cat $bed | wc -l)
nrecrdst=$((nrecrds/3))
echo read $nrecrdst records
csv=res_ce_all/IGV_unique_ce_extension.csv #"sam_junc1.txt" #"junctions.txt"
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
	echo processing event num $eventn, majiq event $fn, pgp_identified event $event_identified and line is $line


	./ggsashimi_txV4.py -PGPTx "$event_identified" -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -FullTitle "$fulltitle" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o res_ce_all/ce_ext_sashimi_plots/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
		rsvg-convert -f pdf -o res_ce_all/ce_ext_sashimi_plots/"$fn".pdf res_ce_all/ce_ext_sashimi_plots/"$fn".svg
	#[ -e res_ce_all/ce_ext_sashimi_plots/"$fn".svg ] && rm res_ce_all/ce_ext_sashimi_plots/"$fn".svg

done

#now merge all pdf's

pdfunite res_ce_all/ce_ext_sashimi_plots/*.pdf res_ce_all/ce_exten_all_sashimi_plots.pdf
fi

######THIS IS FOR ALL IR EVENTS

if [ $1 -eq 5 ]
then
	[ "$(ls -A res_ce_all/ir_sashimi_plots/)" ] && rm res_ce_all/ir_sashimi_plots/*.*
	mkdir -p res_ce_all/ir_sashimi_plots
bed=res_ce_all/IR_coord_uniq_sashimi.bed  #skiptics_sashimi.bed #"sam_junc1.txt" #"junctions.txt"
readarray -t all_bed_data < $bed
nrecrds=$(cat $bed | wc -l)
echo read $nrecrds records
csv=res_ce_all/IGV_unique_IR.csv #"sam_junc1.txt" #"junctions.txt"
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


        echo processing event num $eventn event $fn and majiq event is $majiq_event
        ./ggsashimi_txV4.py -PGPTx "$event_identified"  -A "median_j" -b all_bams.tsv -c "$line" -g Homo_sapiens.GRCh38.103.chr.sorted_new.gtf -GeneName "$gene_name" -MajiqStrnd "$strnd" -ORIG 1 -UEX "$exon1" -DEX "$exon2" -MajiqTx "$majiq_event" -PGPTx "$event_identified" -Majiq "$fn" -Tx "$TxID" -M 1 -C 3 -o res_ce_all/ir_sashimi_plots/"$fn" -O 3 --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
             rsvg-convert --unlimited  -f pdf -o res_ce_all/ir_sashimi_plots/"$fn".pdf res_ce_all/ir_sashimi_plots/"$fn".svg
        #[ -e res_ce_all/ir_sashimi_plots/"$fn".svg ] && rm res_ce_all/ir_sashimi_plots/"$fn".svg

done
#now merge all pdf's

pdfunite res_ce_all/ir_sashimi_plots/*.pdf res_ce_all/IR_all_sashimi_plots.pdf
fi
