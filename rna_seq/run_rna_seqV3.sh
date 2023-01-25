#!/bin/bash -l
#tools required: cat, zcat, cutadapt, fastqc, STAR, samtools
#Step 1: Create STAR INDEX - Please make sure that GRCh38.p13.genome.fa and gencode.v38.annotation.gtf (or files for desired geneome must be present in current directory)

#PLEASE COMMENT THESE BELOW TWO LINES, PROBABLY WILL NOT WORK ON BIOWULF
#Request  24 processors and at least 98 GB of RAM
##$ -pe omp 24
##$ -l mem_per_core=4G


#NOTE: Please uncomment following two line if creating STAR INDEX from scratch. Otherwise make sure that star_index folder exists in the current directory
#otherwise please refere to correct star_index directory

#set this path if using already generated STAR Index
#star_index=../../RawData/star_index

#PLEASE MODIFY BELOW LINE TO EXTRACT UNIQUE SAMPLE NAMES FROM FASTQ FILES
SAMPLES=$(ls -1 *R1*.gz | awk -F '_' '{print $1"_"$2}' | sort | uniq )

#PLEASE MAKE flg=1 IF YOU WANT TO CREATE STAR INDEX FILES AS WELL AS CONCATENATE fastq files
star_flg=0
if [ $star_flg -eq 1 ]
then

#PLEASE NOTE THAT REFerence creation take lot of memory, please increase RAM allocation using --limitBAMsortRAM flag in below line if code crashes at this stage - It is proven to be working with RAM>500GB for GRCh38.p13 genome
#PLEASE UNCOMMENT THIS LINE FOR STAR REFERENCE FILES CREATION
mkdir -p star_index
STAR --runMode genomeGenerate   --runThreadN 32 --limitBAMsortRAM=644245094400  --genomeDir star_index --genomeFastaFiles GRCh38.p13.genome.fa --sjdbGTFfile gencode.v38.annotation.gtf --sjdbOverhang 149
fi
#Step 2-4. (concatenaton of samples from Lane{1,2}, Adapter Trimming for Sample_R{1,2} to get trimmed samples, qc for each sample and star pass1
concat_flg=0
if [ $concat_flg -eq 1 ]
mkdir -p concatenated
echo starting sample concatenation for L{1,2} for all samples

for i in $SAMPLES
do
  echo processing sample $i to concatenate L001 and L002 files for R1 and R2
  #concatenate each sample from two lanes for R1
  cat $i\_L001_R1_001.fastq.gz $i\_L002_R1_001.fastq.gz > concatenated/$i\.R1.fastq.gz
  #concatenate each sample from two lanes for R2
  cat $i\_L001_R2_001.fastq.gz $i\_L002_R2_001.fastq.gz > concatenated/$i\.R2.fastq.gz
done
echo done sample concatenation for L{1,2} for all samples and starting adapter trimming and fastqc step
fi

#PLEAE CHANGE BELOLW LINE TO flag_qc=1 IF YOU WANT TO CUT ADAPTERS FROM raw fastq files and GENERATE QC for each sample
#COMMENT fastqc lines if only trimming is required
flag_qc=0
if [ $flag_qc -eq 1 ]
then
mkdir -p trimmed
mkdir -p qc_raw
for i in $SAMPLES
do
  #trimm adapters from Illumina Truseq
  echo processing sample $i for cutadapt
  cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --overlap 6 -q 20 -j 16 --minimum-length 25 concatenated/$i\.R1.fastq.gz  concatenated/$i\.R2.fastq.gz -o trimmed/$i\.R1.trim.fastq.gz -p trimmed/$i\.R2.trim.fastq.gz #&> trimmed/$i\.trim.log
  echo processing sample $i for fastqc
  fastqc -t 16 trimmed/$i\.R1.trim.fastq.gz trimmed/$i\.R2.trim.fastq.gz -o qc_raw/
done
fi
#####THIS SECTION DOES 2 pass START ALIGNMENT and indexing
mkdir -p starpass1
echo completed concatenation and starting star pass1
for i in $SAMPLES
do
  #star pass1
  echo processing sample $i for star pass1
  STAR  --outTmpDir starTmp --runThreadN 16 --limitBAMsortRAM=322122547200 --genomeDir star_index/ --runMode alignReads --readFilesIn <(zcat trimmed/$i\.R1.trim.fastq.gz trimmed/$i\.R2.trim.fastq.gz) --outFileNamePrefix  starpass1/$i\.pass1 --outSAMattrRGline ID:$i PL:Illumina SM:$i --outSAMstrandField intronMotif -–outSAMattributes XS --outSAMtype BAM SortedByCoordinate
  rm -fr starTmp
done
echo completed star pass1 for all samples, moving on to samtools for indexing
for i in $SAMPLES
do
  #samtools index generation
  echo generating index file for sample $i afte star pass1
  samtools index -@ 16 starpass1/$i\.pass1Aligned.sortedByCoord.out.bam
done
echo completed indexing for star pass1 for all samples, moving on to star pass2
#star pass 2
mkdir -p starpass2
for i in $SAMPLES
do
  echo started processing sample $i for star pass2
  cat starpass1/$i\.pass1SJ.out.tab | awk -v OFS='\t' '{print $1,$2,$3,$4}' > starpass1/temp_$i\.sjdb.txt
  STAR  --outTmpDir starTmp --runThreadN 16 --limitBAMsortRAM=322122547200 --genomeDir star_index/ --runMode alignReads --sjdbFileChrStartEnd starpass1/temp_$i\.sjdb.txt --readFilesIn <(zcat trimmed/$i\.R1.trim.fastq.gz trimmed/$i\.R2.trim.fastq.gz) --outFileNamePrefix  starpass2/$i\.pass2 --outSAMattrRGline ID:$i PL:Illumina SM:$i --outSAMstrandField intronMotif -–outSAMattributes XS --outSAMtype BAM SortedByCoordinate
  rm starpass1/temp_$i\.sjdb.txt
  rm -fr starTmp
done
echo completed star pass2 for all samples, moving on to samtools for indexing
for i in $SAMPLES
do
  echo processing sample $i for samtools
  samtools index -@ 16 starpass2/$i\.pass2Aligned.sortedByCoord.out.bam
done

echo ALL DONE -hopefully- successfully
