#!/bash/bin

#PLEASE MAKE SURE THAT CONDITIONS.TXT IS PRESENT IN CURRENT FOLDER
#Step 0. Run stringtie to generate gtf files from each sample bam file
#Read sample files (using bam sample names here) and create gtf files for each sample

#samples=$(ls *.bam | cut -d'.' -f 1 | awk '{ print length, $0 }' | sort -n | cut -d" " -f2-)
samples=$(ls *.bam)
flg=0
if [ $flg -eq 1 ]
then
mkdir -p assembly
[ -e samples_nt_tdp.txt ] && rm samples_nt_tdp.txt
#samples=$(ls *.bam | cut -d'.' -f 1 | awk '{ print length, $0 }' | sort -n | cut -d" " -f2-)
for sample in $samples
do
    echo STARTED PROCESSING SAMPLE: $sample FOR GTF CREATIONS
    sam=$(echo $sample | cut -d. -f1)
    echo $sam
    stringtie $sample -l "$sam" -p 8 -A -C -B -G /home/ext_syed_datatecnica_com/mountfolder/ref_genome/gencode.v38.annotation.gtf -o assembly/"$sam".gtf
    echo assembly/"$sam".gtf >> samples_nt_tdp.txt
done


#cat samples_nt_tdp_exc.txt > samples_nt_tdp.txt
#FILE CONTAINING FULL PATH TO SAMPLE NAMES
cat samples_nt_tdp.txt

#FILENAME="samples_nt_tdp.txt"
#step1 - merge all gtf files - please update locations of genome gtf file and stringtie generated gtf files (in samples_nt_tdp.txt).
echo STARTED MERGING .gtf FILES
stringtie --merge -G /home/ext_syed_datatecnica_com/mountfolder/ref_genome/gencode.v38.annotation.gtf -p 8 -o stringtie_merged_nttdp.gtf samples_nt_tdp.txt
#step2 -get sample names and remove trailing gtf from samples file used in previous step
#COMMENTED THIS - ADDED
#samples=$(cut -d/ -f6 "$FILENAME" | cut -d. -f1)
#ADDED THIS (in /data/seddighis/multi_cell/JCM6188-14_S14.gtf , we gtf file is after 4th "/", so cut -d/ f4 will pick the JCM6188-14_S14.gtf file)- ADDED
samples=$(cut -d/ -f2 samples_nt_tdp.txt | cut -d. -f1)
#PLEASE MAKE SURE THAT BAM files are stored in folder /data/seddighis/multi_cell/ otherwise it will not work - ADDED
[ -e file1.txt ] && rm file1.txt
[ -e file2.txt ] && rm file2.txt
#RESTIMATION
for sample in $samples
do
  echo PROCESSING RE-ESTIMATION FOR SAMPLE: "$sample"
  #please update paths for bam files for each sample.
  stringtie "$sample".pass2Aligned.sortedByCoord.out.bam -G stringtie_merged_nttdp.gtf -e -p 8 -o assembly/"$sample".merged.gtf
  echo $sample.merged >> file1.txt
  echo assembly/$sample.merged.gtf >>file2.txt
done

[ -e input_prepDE.txt ] && rm input_prepDE.txt
#create input file for prepDE.py script
paste -d'\t' file1.txt file2.txt > input_prepDE.txt
#create input file for run_deseq2.R
#ADD string "sample" in first row
sed '1 s/^/samples\n/' file1.txt > file11.txt
fi
#Also create conditions.txt file
yes "ctrl" | head -n 16 > conditions.txt
yes "tdp43kd" | head -n 9 >> conditions.txt
#ADD string "condition" in the first row
sed '1 s/^/condition\n/' conditions.txt > conditions1.txt
#PLEASE NOTE THAT COLUMN BATCH SHOULD CONTAIN B (FOR n=4, kd=3) or A (FOR n=12, kd=6). PLEASE ADJUST IF YOUR CONDITIONS ARE DIFFERENT
yes "B" | head -n 4 > file3.txt
yes "A" | head -n 12 >> file3.txt
yes "B" | head -n 3 >> file3.txt
yes "A" | head -n 6 >> file3.txt
#ADD string "batch" in the first row
sed '1 s/^/batch\n/' file3.txt > file33.txt

#CREATE DESIGN MATRIX - HERE WE ARE ADDING BATCH COLUMN TO COUNTER BATCH EFFECTS
#paste -d'\t' file11.txt conditions1.txt > new_samples.csv
paste -d'\t' file11.txt conditions1.txt file33.txt > new_samples.csv
#step3 - Finally read input_prepDE.txt file and generate counts table for input to DESeq2
#INPUTFNAME="input_prepDE.txt"
#echo "$INPUTFNAME"
cat input_prepDE.txt
prepDE.py -i input_prepDE.txt
#step4 - use generated gene_count_matrix.csv (containing count values for all genes in all conditions, it is number of genes*number of samples (here 18 samples, 12 control and 6 TDP43)) file alongwith new_samples.csv () file to run DESeq2.
#Now run deseq2
#Rscript run_deseq2.R
