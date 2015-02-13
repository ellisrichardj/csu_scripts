#!/bin/bash
set -e

# script by ellisrichardj
# This will merge all files with matching names in each directory
# Useful for merging data from 2 (or possible more) Illumina runs of the same samples
# Merged files will be written to current working directory

# Version 0.0.1 06/11/14
# Version 0.1.1 22/01/15 - Altered to match sample names only rather than full file name
# Version 0.1.2 10/02/15 - Allows two 'cat' processes to run concurrently

# check for mandatory positional parameters
if [ $# -lt 2 ]; then
  echo "
Usage: $0 {path to Folder 1}  {path to Folder 2}"
exit 1
fi

Dir1="$1"
Dir2="$2"

for file in "$Dir1"/*_*R1*.gz
do
   	fname=$(basename "$file")
	samplename=${fname%%_*}
	D1R1="$Dir1"/"$samplename"_*R1*.gz
	D2R1="$Dir2"/"$samplename"_*R1*.gz
	D1R2="$Dir1"/"$samplename"_*R2*.gz
	D2R2="$Dir2"/"$samplename"_*R2*.gz
   if [ -f "$Dir2"/"$samplename"_*R1*.gz ]; then
	echo "Merging "$samplename""
	cat $D1R1 $D2R1 > $PWD/"$samplename"_R1.fastq.gz &
	cat $D1R2 $D2R2 > $PWD/"$samplename"_R2.fastq.gz &
	wait
   else
	echo "Matching file "$samplename" does not exist in both directories: skipping"
   fi
done

