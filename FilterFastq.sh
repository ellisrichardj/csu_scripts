#!/bin/bash
set -e

# script by ellisrichardj
# This will use kraken to pre-filter fastq files to contain just the organism of interest

# Requirements:
# Trimmomatic
# kraken
# seqtk

# Version 0.0.1 27/07/17 Initial version

# set defaults for the options

taxon=*

# parse the options
while getopts 't:' opt ; do
  case $opt in
    t) taxon=$OPTARG ;;
  esac
done

# skip over the processed options
shift $((OPTIND-1)) 


# check for mandatory positional parameters
if [ $# -lt 2 ]; then
  echo "
Usage: $0 [-t taxon of interest] <fastq_R1>  <fastq_R2>

Taxon of interest can be any name that can be filtered by regex
e.g. Mycobact will select all Mycobacterium sequences
"
exit 1
fi

# Set variables
R1=$(readlink -f "$1")
R2=$(readlink -f "$2")
now=$(date '+%x %R')
Start=$(date +%s)

	sfile1=$(basename "$R1")
	sfile2=$(basename "$R2")
	samplename=${sfile1%%_*}

threads=$(grep -c ^processor /proc/cpuinfo)

#Quality trim
trimmomatic-0.30.jar PE -threads $threads -phred33 "$R1" "$R2" "$samplename"_R1_trim_paired.fastq "$samplename"_R1_unpaired.fastq "$samplename"_R2_trim_paired.fastq "$samplename"_R2_unpaired.fastq ILLUMINACLIP:/home/richard/ReferenceSequences/adapter.fasta:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60
rm "$samplename"_R1_unpaired.fastq
rm "$samplename"_R2_unpaired.fastq

#Filter for specified taxon
echo "Now running kraken on the trimmed reads"
kraken --threads 4 --quick --preload --paired --fastq-input "$samplename"_R1_trim_paired.fastq "$samplename"_R2_trim_paired.fastq | kraken-translate > "$samplename".trans
echo "Identifying "$taxon" reads"
grep $taxon "$samplename".trans > "$samplename"_"$taxon".trans
echo "Generating filtered fastq files"
~/seqtk/seqtk subseq "$samplename"_R1_trim_paired.fastq "$samplename"_"$taxon".trans | gzip > "$samplename"_"$taxon"_R1.fastq.gz &
~/seqtk/seqtk subseq "$samplename"_R1_trim_paired.fastq "$samplename"_"$taxon".trans | gzip > "$samplename"_"$taxon"_R2.fastq.gz &
wait

End=$(date +%s)
TimeTaken=$((End-Start))
echo "Filtered fastqs are "$samplename"_"$taxon"_R1.fastq.gz and "$samplename"_"$taxon"_R1.fastq.gz"
echo  | awk -v D=$TimeTaken '{printf "Filtered dataset in: %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'


