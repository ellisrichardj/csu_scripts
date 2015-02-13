#!/bin/bash
set -e

# script by ellisrichardj
# This will perform denovo assembly on all samples within a directory
# It works by running velvet with 2 different k-mer values and then merging these to produce a higher quality hybrid assembly
# Requires Velvet_assemble.sh, Velvet and trimmomatic
# Velvet must have been compliped to allow longer k-mers and multi-threading

# Version 0.1
# Version 0.1.1 06/10/14 - added output for number of samples and time taken
# Version 0.1.2 16/01/15 - added options to change velvet parameters (kmer and cutoff)
#			 - performs second assembly with lower kmer (-30)
# Version 0.1.3 21/01/15 - Collates assembly statistics (N50, largest contig, total size, etc) for all samples and k-mers
# Version 0.1.4 10/02/15 - Altered wildcard/filename pattern to allow for data merged with MergeFQ.sh 

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

LowerK=$((KVALUE-30))

# skip over the processed options
shift $((OPTIND-1)) 

# check for mandatory positional parameters
if [ $# -lt 1 ]; then
  echo "Usage: $0 {path to Data Folder} "
  echo "Options: -k kmervalue (def: $KVALUE) | -c cov_cutoff (def: $CUTOFF)"
exit 1
fi

DataFolder="$1"
Count=0
Start=$(date +%s)

echo "Sample_k-mer, Nodes, N50, MaxContig, TotalSize," >> AssemblyStats.csv
Statsfile=$(readlink -f AssemblyStats.csv)

for file in "$DataFolder"/*_R1*.gz
do
	((Count=Count+1))
	fname=$(basename "$file")
	samplename=${fname%%_*}
	mkdir "$samplename"_assembly
	ln -s "$(readlink -f "$DataFolder"/"$samplename"*R1*.gz)" "$samplename"_assembly/"$samplename"_R1.fastq.gz
	ln -s "$(readlink -f "$DataFolder"/"$samplename"*R1*.gz)" "$samplename"_assembly/"$samplename"_R2.fastq.gz
	cd "$samplename"_assembly

	echo "Assembling sample "$Count": "$samplename" with K=$KVALUE and cutoff=$CUTOFF"
	Velvet_assemble.sh -k "$KVALUE" -c "$CUTOFF" "$samplename" "$samplename"_R1.fastq.gz "$samplename"_R2.fastq.gz
	Assembly1=$(readlink -f "$samplename"_"$KVALUE"/*_contigs.fa)
	echo -e "$Assembly1""\t"$samplename"_Velvet"$KVALUE"" > "$samplename"_assemblies.txt
	awk '{ORS=""}/velvetg/{print $2","}/max/{print $4",",$9,$11,$13"\n"}' "$samplename"_"$KVALUE"/*Log >> $Statsfile
	Velvet_assemble.sh -k "$LowerK" "$samplename" "$samplename"_R1.fastq.gz "$samplename"_R2.fastq.gz

	Assembly2=$(readlink -f "$samplename"_"$LowerK"/*_contigs.fa)
	echo -e "$Assembly2""\t"$samplename"_Velvet"$LowerK"" >> "$samplename"_assemblies.txt
	awk '{ORS=""}/velvetg/{print $2","}/max/{print $4",",$9,$11,$13"\n"}' "$samplename"_"$LowerK"/*Log >> $Statsfile

	rm *.gz
	cd ..
done

End=$(date +%s)
TimeTaken=$((End-Start))

echo  | awk -v D=$TimeTaken '{printf "Assembled '$Count' samples in: %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'
