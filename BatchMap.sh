#!/bin/bash
set -e

# script by ellisrichardj
# This will perform mapping of all samples within a directory
# Requires BWA and samtools

# Version 0.1.1 21/10/14

# 

# check for mandatory positional parameters
if [ $# -lt 2 ]; then
  echo "
Usage: $0 <path to reference> <path to Data Folder> "
exit 1
fi

Start=$(date +%s)

Ref=$1
DataFolder="$2"
Count=0
threads=$(grep -c ^processor /proc/cpuinfo)

	ref=$(basename "$Ref")
	refname=${ref%%_*}
	reffile=${ref%%.*}

homedir=$(pwd)
cp "$Ref" "$homedir"
bwa index "$ref"

for file in "$DataFolder"/*R1*.gz
do
	((Count=Count+1))
	fname=$(basename "$file")
	samplename=${fname%%_*}
	mkdir "$samplename"_mapto_"$refname"
	cd "$samplename"_mapto_"$refname"	
	echo "Mapping sample "$Count": "$samplename" to "$ref""
	bwa mem -t "$threads" "$homedir"/"$ref" "$DataFolder"/"$samplename"*R1*.gz "$DataFolder"/"$samplename"*R2*.gz | samtools view -Su - | samtools sort - "$samplename"-"$refname"_map_sorted
	samtools index "$samplename"-"$refname"_map_sorted.bam
	samtools flagstat "$samplename"-"$refname"_map_sorted.bam > "$samplename"-"$refname"_MappingStats.txt
	cd "$homedir"	
done

End=$(date +%s)
TimeTaken=$((End-Start))

echo  | awk -v D=$TimeTaken '{printf "Mapped '$Count' samples in: %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'
