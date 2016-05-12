#!/bin/bash
set -e

# script by ellisrichardj
# This will perform denovo assembly on all samples within a directory
# It has been written with prokaryotic genomes in mind
# It works by running velvet with 3 different k-mer values and then merging these to produce a higher quality
# hybrid assembly
# Requires Velvet_assemble.sh, Velvet, CISA and trimmomatic
# Velvet must have been (re)compiled to allow longer k-mers and multi-threading (cd to velvet directory, then
# 	[sudo] make 'MAXKMERLENGTH=301' 'OPENMP=1' 'BIGASSEMBLY=1' 'LONGSEQUENCES=1'

# Version 0.1 initial version
# Version 0.1.1 06/10/14 - added output for number of samples and time taken
# Version 0.1.2 16/01/15 - added options to change velvet parameters (kmer and cutoff)
#			 - performs second assembly with lower kmer (-30)
# Version 0.1.3 21/01/15 - Collates assembly statistics (N50, largest contig, total size, etc) for all samples
#			and k-mers
# Version 0.1.4 10/02/15 - Altered wildcard/filename pattern to allow for data merged with MergeFQ.sh
# Version 0.2.0 03/03/15 - Performs assemblies with three k-mer values and now uses CISA to merge assemblies 
#			generated with multiple k-mers, added comments throughout script to help understand
#			what the script is doing
# Version 0.2.1 17/03/15 - Correct small error in on-screen messages, bug in symlinks to raw data files,
#			move final assemblies into single directory
# Version 0.2.2 16/12/15 - Extract assembly statistics for Final assembly
# Version 0.2.3 10/05/16 - Add CISA assembly stats table to Output directory, clean up intermediate files
#			after each sample


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

MidK=$((KVALUE-18))
LowerK=$((KVALUE-36))

# skip over the processed options
shift $((OPTIND-1)) 

# check for mandatory positional parameters
if [ $# -lt 1 ]; then
  echo "Usage: $0 {path to Data Folder} "
  echo "Options: -k kmervalue (default: $KVALUE; minimum value:57)"
  echo "         -c cov_cutoff (def: $CUTOFF)"
exit 1
fi

DataFolder="$1"
Count=0
Start=$(date +%s)

# Prepare stats file and final output folder
echo "Sample_k-mer, Nodes, N50, MaxContig, TotalSize," >> AssemblyStats.csv
Statsfile=$(readlink -f AssemblyStats.csv)
mkdir FinalAssemblies
OutputDir=$(pwd)/"FinalAssemblies"
echo "Sample_k-mer, Nodes, N50, MaxContig, TotalSize," >> "$OutputDir"/FinalAssemblyStats.csv
FinalStats=$(readlink -f "$OutputDir"/FinalAssemblyStats.csv)

# Run assembly for each sample
for file in "$DataFolder"/*_R1*.gz
do
	((Count=Count+1))
	fname=$(basename "$file")
	samplename=${fname%%_*}
# Create new directory for each sample
	mkdir "$samplename"_assembly
#	ln -s "$(readlink -f "$DataFolder"/"$samplename"*R1*.gz)" "$samplename"_assembly/"$samplename"_R1.fastq.gz
#	ln -s "$(readlink -f "$DataFolder"/"$samplename"*R2*.gz)" "$samplename"_assembly/"$samplename"_R2.fastq.gz
	cd "$samplename"_assembly

# Assemble sample using highest k-mer value
	echo "Assembling sample "$Count": "$samplename" with K=$KVALUE and cutoff=$CUTOFF"

	Velvet_assemble.sh -k "$KVALUE" -c "$CUTOFF" "$samplename" "$DataFolder"/"$samplename"*R1*.gz "$DataFolder"/"$samplename"*R2*.gz

	echo -e count=3 > "$samplename"_assemblies.txt
	Assembly1=$(readlink -f "$samplename"_"$KVALUE"/*_contigs.fa)
	echo -e data="$Assembly1"",title="$samplename"_Velvet"$KVALUE"" >> "$samplename"_assemblies.txt
	awk '{ORS=""}/velvetg/{print $2","}/max/{print $4",",$9,$11,$13"\n"}' "$samplename"_"$KVALUE"/*Log >> $Statsfile

# Assemble sample using mid k-mer value
echo "Assembling sample "$Count": "$samplename" with K=$MidK and cutoff=$CUTOFF"

	Velvet_assemble.sh -k "$MidK" -c "$CUTOFF" "$samplename" "$DataFolder"/"$samplename"*R1*.gz "$DataFolder"/"$samplename"*R2*.gz

	Assembly2=$(readlink -f "$samplename"_"$MidK"/*_contigs.fa)
	echo -e data="$Assembly2"",title="$samplename"_Velvet"$MidK"" >> "$samplename"_assemblies.txt
	awk '{ORS=""}/velvetg/{print $2","}/max/{print $4",",$9,$11,$13"\n"}' "$samplename"_"$MidK"/*Log >> $Statsfile

# Assemble sample using lowest k-mer value
echo "Assembling sample "$Count": "$samplename" with K=$LowerK and cutoff=$CUTOFF"

	Velvet_assemble.sh -k "$LowerK" -c "$CUTOFF" "$samplename" "$DataFolder"/"$samplename"*R1*.gz "$DataFolder"/"$samplename"*R2*.gz

	Assembly3=$(readlink -f "$samplename"_"$LowerK"/*_contigs.fa)
	echo -e data="$Assembly3"",title="$samplename"_Velvet"$LowerK"" >> "$samplename"_assemblies.txt
	echo -e Master_file="$samplename"_contigs.fa >> "$samplename"_assemblies.txt
	awk '{ORS=""}/velvetg/{print $2","}/max/{print $4",",$9,$11,$13"\n"}' "$samplename"_"$LowerK"/*Log >> $Statsfile

# run CISA to merge different k-mer assemblies
	echo "Merging Assemblies for sample "$Count": "$samplename""
	Merge.py "$samplename"_assemblies.txt
	genome=$(awk '/Whole/{print $4}' Merge_info | awk '$0>x{x=$0};END{print x}')
	echo genome="$genome" > "$samplename"_CISA.config
	echo infile="$samplename"_contigs.fa >> "$samplename"_CISA.config
	echo outfile="$samplename"_CISA.fa >> "$samplename"_CISA.config
	echo CISA=/usr/local/bioinf/CISA1.3 >> "$samplename"_CISA.config
	CISA.py "$samplename"_CISA.config

# extract CISA output assembly metrics
	FinalSize=$(cat "$samplename"_CISA.fa | awk '$0 !~ ">" {c+=length($0);} END { print c; }')
	FinalContigs=$(grep -c '>' "$samplename"_CISA.fa)
	FinalMax=$(wc -L < "$samplename"_CISA.fa)
	FinalN50=$(cat "$samplename"_CISA.fa | awk '$0 ~ ">" {print c; c=0 "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | sort -n | awk '{len[i++]=$1;sum+=$1} END {for (j=0;j<i+1;j++) {csum+=len[j]; if (csum>sum/2) {print len[j];break}}}')
	echo "$samplename"_CISA,"$FinalContigs","$FinalN50","$FinalMax","$FinalSize" >> $Statsfile
	echo "$samplename"_CISA,"$FinalContigs","$FinalN50","$FinalMax","$FinalSize" >> $FinalStats

	cp "$samplename"_CISA.fa $OutputDir
	cd ..
	rm -rf "$samplename"_assembly
done

echo "Final Assemblies for "$Count" samples are in "$OutputDir""

End=$(date +%s)
TimeTaken=$((End-Start))

echo  | awk -v D=$TimeTaken '{printf "Assembled '$Count' samples in: %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'
