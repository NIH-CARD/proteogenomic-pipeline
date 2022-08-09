inp=$1
inp1=$(echo $1 | cut -d'.' -f1)


#now call bedtools getfasta function to get nt sequence from reference_genome.fa file
bedtools getfasta -fi /Volumes/SYEDSHAH/MichaelLab/ref_genome/GRCh38.p13.genome.fa -bed $inp -s > "$inp1"_nt.fa

#now remove >chr lines from the resulting file
awk '!/^>chr/' "$inp1"_nt.fa > "$inp1"_nt1.fa


paste -d"\t" $inp "$inp1"_nt1.fa > temp."$inp1"_nt.bed

rm "$inp1"_nt1.fa
#rm "$inp1"_nt.fa
awk -F "\t" '{if($6=="+") print ">sp|"$7"_"$1"_"$2"-"$3"_plus""\n"$8;else print ">sp|"$7"_"$1"_"$2"-"$3"_minus""\n"$8}' temp."$inp1"_nt.bed > temp."$inp1"_nt.transeq_in.fasta
rm temp."$inp1"_nt.bed
#also save nt sequence

TEMPFILE=$(echo temp."$inp1"_nt.transeq_in.fasta | perl -pe 's/\.fasta$/.temp/')
OUTPUTFILE=$(echo temp."$inp1"_nt.transeq_in.fasta | perl -pe 's/\.fasta$/.trans/')

# Translate sequence

transeq -sequence temp."$inp1"_nt.transeq_in.fasta -outseq "$TEMPFILE" -frame F


# Rename sequence
perl -sape 's/>/>sp|/' "$TEMPFILE" > "$OUTPUTFILE"


#no concatenate all files
cat "$OUTPUTFILE" >> temp."$inp1"_AA.fasta
# Remove temp file
rm "$TEMPFILE"
rm "$OUTPUTFILE"

rm temp."$inp1"_nt.transeq_in.fasta

#also remove all newlines from the
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' temp."$inp1"_AA.fasta > "$inp1"_AA.fasta

#remove
rm temp."$inp1"_AA.fasta
