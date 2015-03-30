#!/bin/bash

# Script by ellisrichardj for denovo assembly
# Starting with raw illumina output (fastq.gz) this script will perform adapter removal
# and quality trimming of the reads before using Velvet for the assembly. 

# Version 0.1.1 10/09/14
# Version 0.1.2 14/01/15 Change commands to assume that velveth and velvetg are in PATH
# Version 0.1.3 16/01/15 Added k-mer value to file names
# Version 0.1.4 30/03/15 Remove options for read tracking and generation of AMOS file to improve speed
#			and minimize disk space requirements

# set our defaults for the options
KVALUE=101
CUTOFF=auto

# parse the options
while getopts 'c:k:' opt ; do
  case $opt in
    c) CUTOFF=$OPTARG ;;
    k) KVALUE=$OPTARG ;;
  esac
done

# skip over the processed options
shift $((OPTIND-1)) 

# check for mandatory positional parameters
if [ $# -lt 3 ]; then
  echo "Usage: $0 [options] outputFolder Read1 Read2"
  echo "Options: -k kmervalue (def: $KVALUE) | -c cov_cutoff (def: $CUTOFF)"
  exit 1
fi
DIR="$1"
LEFT="$2"
RIGHT="$3"
# do the quality trimming and assembly
trimmomatic-0.30.jar PE -threads 6 -phred33 "$LEFT" "$RIGHT" "$DIR"_R1_trim_paired.fastq "$DIR"_R1_unpaired.fastq "$DIR"_R2_trim_paired.fastq "$DIR"_R2_unpaired.fastq ILLUMINACLIP:/home/richard/ReferenceSequences/adapter.fasta:2:30:10 MINLEN:"$KVALUE"
rm "$DIR"_R1_unpaired.fastq
rm "$DIR"_R2_unpaired.fastq
echo "Assembling with K=$KVALUE and cutoff=$CUTOFF"
velveth "$DIR"_"$KVALUE" "$KVALUE" -shortPaired -fmtAuto -separate "$DIR"_R1_trim_paired.fastq "$DIR"_R2_trim_paired.fastq
velvetg "$DIR"_"$KVALUE" -exp_cov auto -cov_cutoff "$CUTOFF" -clean yes #-read_trkg yes -amos_file yes

mv "$DIR"_"$KVALUE"/{Log,"$DIR"_Log} 
mv "$DIR"_"$KVALUE"/{contigs.fa,"$DIR"_"$KVALUE"_contigs.fa}
# mv "$DIR"_"$KVALUE"/{velvet_asm.afg,"$DIR"_"$KVALUE"_velvet_asm.afg}
rm "$DIR"_R1_trim_paired.fastq
rm "$DIR"_R2_trim_paired.fastq

echo "Results are in: "$DIR"_"$KVALUE""
