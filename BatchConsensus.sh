#!/bin/bash
set -e

# script by ellisrichardj
# This will perform consensus calling of bam files within all subdirectories of a given directory
# Requires samtools, vcf2consensus.pl

# Version 0.1.1 14/01/15

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
refbase=$(basename "$Ref")
refname=${refbase%%_*}

for file in "$DataFolder"/*/*.bam
do
	((Count=Count+1))
	dir=$(dirname "$file")
	cd $dir
	fname=$(basename "$file")
	samplename=${fname%%_*}

	# remove PCR duplicates from mapping file and re-index
	samtools rmdup "$file" "$samplename"_clean.bam
	samtools index "$samplename"_clean.bam

	# generate and correctly label consensus
	samtools mpileup -L 2000 -AEuf "$Ref" "$samplename"_clean.bam | bcftools view -cg - > "$samplename".vcf
	vcf2consensus.pl consensus -f "$Ref" "$samplename".vcf  | sed '1s/.*/>'"$samplename"_mappedto_"$refname"'/g' - > "$samplename"_mappedto_"$refname"_consensus.fas
	cd ..	
	
done

End=$(date +%s)
TimeTaken=$((End-Start))

echo  | awk -v D=$TimeTaken '{printf "Generated consensus from '$Count' samples in: %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'
