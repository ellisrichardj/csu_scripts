#!/bin/bash
set -e

# Pipeline to generate the consensus sequence for an unknown influenza virus genome.  Based on Sherlock.sh and adapted to work on genome segments independently.

# Analysis steps required to find unknown infectious agent using Illumina shotgun sequence data from DNA/RNA extracted from host tissue.  Initially the host sequences are removed by mapping to the host (or closely related) genome.  

# Pipeline by ellisrichardj making use of a variety of standard bioinformatic tools

# Software required:
	# bwa - 0.7.5a or above
	# samtools
	# velvet (must be compiled for multithreading and use of k-mers of at least 101) 
	# Trimmomatic
	# BLAST+
	# IterMap.sh (available at https://github.com/ellisrichardj/csu_scripts)
	# vcf2consensus.pl (also available at https://github.com/ellisrichardj/csu_scripts)

# Ensure that each of these are in PATH or symbolic links are in a folder that is in your PATH

# Version 0.2.0	01/06/05	Initial version adapted from Sherlock version 0.2.0

# set defaults for the options
KVALUE=101
CUTOFF=auto
Blast_e_value=0.0001

Start=$(date +%s)

# parse the options
while getopts 'H:e:r:c:k:' opt ; do
  case $opt in
    H) PathToHOST=$OPTARG ;;
    e) Blast_e_value=$OPTARG ;;
    r) PathToSearchData=$OPTARG ;;
    c) CUTOFF=$OPTARG ;;
    k) KVALUE=$OPTARG ;;
  esac
done
# skip over the processed options
shift $((OPTIND-1)) 
# check for mandatory positional parameters
if [ $# -lt 2 ]; then
  echo "
Usage: $0 [options] -H path to HOST genome -r path to Reference Database Directory Read1.fastq.gz Read2.fastq.gz"
  echo "Reference database should be a directory containing a separate fasta file for each virus segment"
  echo "Options: De Novo assembly: -k Velvet kmervalue (default = 101) | -c Velvet cov_cutoff (default = auto)"
  echo "Options: Blast Output: -e Blast e value (default = 0.0001)
"
  exit 1
fi
LEFT="$1"
RIGHT="$2"

threads=$(grep -c ^processor /proc/cpuinfo)

# Get sample/reference/host names from file names
	sfile1=$(basename "$LEFT")
	sfile2=$(basename "$RIGHT")
	samplename=${sfile1%%_*}

	ref=$(basename "$PathToSearchData")

	host=$(basename "$PathToHOST")
	hostname=${host%%.*}

# Create output OutputDirectory
OutputDir="$samplename"_Sherlock_"$ref"
mkdir "$OutputDir"
echo "Getting host genome"

# Generate bwa index of host genome if it doesn't already exist
if [ -f "$PathToHOST".sa ]
then 	echo "Using existing host index"
else 	echo "Generating local index of host genome" 
	ln -s "$(readlink -f "$PathToHOST")" "$OutputDir"/HostGenome.fa
	bwa index "$OutputDir"/HostGenome.fa
	PathToHOST="$OutputDir"/HostGenome.fa
fi

# Map raw data to host genome
echo "Mapping raw reads to host reference genome"
bwa mem -t "$threads" "$PathToHOST" "$LEFT" "$RIGHT" | samtools view -Su - | samtools sort - "$OutputDir"/"$samplename"_"$hostname"_map_sorted

# Extract sequence reads that do not map to the host genome
samtools view -b -f 4 "$OutputDir"/"$samplename"_"$hostname"_map_sorted.bam > "$OutputDir"/"$samplename"_nonHost.bam

# Denovo assembly of non-host reads
velveth "$OutputDir"/"$samplename"_nonHost_"$KVALUE" "$KVALUE" -shortPaired -bam "$OutputDir"/"$samplename"_nonHost.bam
velvetg "$OutputDir"/"$samplename"_nonHost_"$KVALUE" -exp_cov auto -cov_cutoff "$CUTOFF" -ins_length 300 -clean yes -unused_reads yes

# Now run Blast of output contigs (potentially including unassembled singleton reads) against database of choice.  If you have some idea of what you are looking for you can limit the search to a particular taxonomic group (e.g. Rhabdoviridae) providing you have a database for it.


cd "$OutputDir"
echo "BLASTing contigs against selected database(s)"
for file in "$PathToSearchData"/*.fas
do
	seg=$(basename "$file")
	segname=${seg%%.*}
	ln -s "$(readlink -f "$file")" "$segname".fasta
	makeblastdb -in "$segname".fasta -parse_seqids -dbtype nucl
	blastall -b 5 -p blastn -d "$segname".fasta -i "$samplename"_nonHost_"$KVALUE"/contigs.fa -o "$samplename"_"$segname"_crunch.txt -m 8 -e "$Blast_e_value" &
done
wait

# Extract the single sequence from blast output corresponding to the highest scoring match

for crunch in *.txt
do	
	topseg=${crunch%%.*}
	part=${topseg#*_}
	segname=${part%_*}
	sort -k12,12 -rn "$crunch" | head -1 - | awk '{print $2}' - > "$topseg"_match.txt
	AccNo=$(sed -e 's/.*\([A-Z][A-Z][0-9][0-9][0-9][0-9][0-9][0-9]\).*/\1/g'  "$topseg"_match.txt)
	blastdbcmd -db "$segname".fasta -entry_batch "$topseg"_match.txt | sed '1s/.*/>'"$AccNo"-"$segname"'/g' -> "$segname"_match.fas
done

cat *_match.fas > top_matches.fas

# Map original data to selected reference sequences and generate new consensus
echo "Now mapping to generate new consensus"
IterMap.sh -i4 top_matches.fas "$LEFT" "$RIGHT"


End=$(date +%s)
TimeTaken=$((End-Start))
echo "Results are in: "$OutputDir""
echo  | awk -v D=$TimeTaken '{printf "Performed Sherlock Analysis in: %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'
