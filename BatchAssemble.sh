#!/bin/bash
set -e

# script by ellisrichardj
# This will perform denovo assembly on all samples within a directory
# Requires Velvet_assemble.sh, Velvet and trimmomatic
# Velvet must have been compliped to allow longer k-mers and multi-threading

# Version 0.1.1 06/10/14

# Changelog 0.1 to 0.1.1 - added out output for number of samples and time taken

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
	echo "Assembling sample "$Count": "$samplename""
	Velvet_assemble.sh "$samplename" "$samplename"_*_*_R1_*.gz "$samplename"_*_*_R2_*.gz	
done

End=$(date +%s)
TimeTaken=$((End-Start))

echo  | awk -v D=$TimeTaken '{printf "Assembled '$Count' samples in: %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'
