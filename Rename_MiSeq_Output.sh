#!/bin/bash
set -e

# script by ellisrichardj
# This will rename contig files that are output from MiSeqReporter with the sample name, 
# rather than a numerical identifier.  Assumes that config files are in an 'Alignment' subdirectory
# of the directory that contains the raw fastq

# Version 0.1 10/02/15

# check for mandatory positional parameters
if [ $# != 1 ]; then
  echo "Usage: $0 {path to Illumina MiSeq Data Folder} "
exit 1
fi

DataFolder="$1"

for file in "$DataFolder"/*_*_*_R1_*.gz
do
	fname=$(basename "$file")
	samplename=${fname%%_*}
	fileend=${fname#*_S}
	samplenumber=${fileend%%_*}
	echo "$samplename"=Sample"$samplenumber"
	mv "$DataFolder"/Alignment/Contigs."$samplenumber".fa "$DataFolder"/Alignment/"$samplename"_contigs.fa
done


