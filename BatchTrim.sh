#!/bin/bash
set -e

# script by ellisrichardj
# This will perform quality timming on all samples within a directory
# Requires trimmomatic-0.30.jar
# Assumes 'adapter.fasta' is in ~/ReferenceSequences/

# Version 0.1.1 23/10/14
# Version 0.1.2 06/11/14 - Run gzip (recompress trimmed fastq) in the background to increase speed

# check for mandatory positional parameters
if [ $# -lt 1 ]; then
  echo "
Usage: $0 {path to Data Folder} "
exit 1
fi

DataFolder="$1"
Count=0
Start=$(date +%s)

for file in "$DataFolder"/*_*_*_R1_*.gz
do
	((Count=Count+1))
	fname=$(basename "$file")
	samplename=${fname%%_*}
	echo "Quality trimming sample "$Count": "$samplename""
trimmomatic-0.30.jar PE -threads 4 -phred33 "$DataFolder"/"$samplename"_*_*_R1_*.gz "$DataFolder"/"$samplename"_*_*_R2_*.gz "$samplename"_R1_trim_paired.fastq "$samplename"_R1_unpaired.fastq "$samplename"_R2_trim_paired.fastq "$samplename"_R2_unpaired.fastq ILLUMINACLIP:/home/sequence/ReferenceSequences/adapter.fasta:2:30:10 SLIDINGWINDOW:10:15 MINLEN:80
	rm "$samplename"_R1_unpaired.fastq
	rm "$samplename"_R2_unpaired.fastq
	gzip "$samplename"_R1_trim_paired.fastq &
	gzip "$samplename"_R2_trim_paired.fastq	&
done

End=$(date +%s)
TimeTaken=$((End-Start))

echo  | awk -v D=$TimeTaken '{printf "Trimmed '$Count' samples in: %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'
