# proteogenomic-pipeline
This repository comprises of three pipelines for biomarker discovery project. Each pipeline make use of a combination of bash (and R) scripts to accomplish specific task.

1. RNA-seq analysis pipeline	
	This pipeline performs typical rna-seq analysis steps including 1) merging fastq files (from different lanes and forward/reverse reads), 2) trimming of fastq files, 3) quality control (fastq), 4) STAR index creation and 5) STAR alignment.
2. Differential expression analysis
	This pipeline performs differential gene analysis by 1) generating gene and transcript abundance file for each sample using StringTie2, 2) generates necessary input files for DESeq2 analysis (please see DESeq2 README file).  

3. proteogenomic-pipeline
	This pipeline identifies potential cryptics (cryptic-exon(ce) and ce-extension events), skiptics (single and multiple exon skipping events) and intron (IR) retention events, creates sashimi plots and generates files containing aminoacid (AA) and nucleotide (nt) sequences from a list (generated from MAJIQ software) of miss-spliced events. AA file of the potential cryptic, skiptics and IR events is used as input to PEAKS search tools for identification of potential biomarkers. It comprises of various bash and R scripts and is divided into 3 parts.
	A.1) Identification of potential splicing events (cryptics, skiptics and IR) along with sashimi plot for each event for easy visualization of all events to visually prune miss-categorized events.
	A.2) Read coverage files for each of the cryptics are calculated.
	B. Final list (after manual modification from part A above) of each event type is used to geenrate AA and nt fasta files alongwith sashimi plots for each category. 
	C) BACKMAPPING: After final AA fasta file (generated from all event types in part B) is search using PEAKS search algorithm for potential biomarker events, output of PEAKS search (alongside files from part B) is used to map back selected peptides (to the event they are selected from) as sashimi plots to validate each event type (cryptics, skiptics and IR). 





