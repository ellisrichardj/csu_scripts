#!/bin/bash
set -e

# Analysis steps required to find unknown infectious agent using Illumina shotgun sequence data from DNA/RNA extracted
# from host tissue.  Initially the host sequences are removed by mapping to the host (or closely related) genome.  

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

# Version 0.1.1 May 2014 Initial version
# Version 0.1.2 25/11/14 Improved Usage guidance, added function to indicate time taken for analysis
# Version 0.1.3 04/03/15 Perform assembly directly from bam file rather than converting back fastq
# Version 0.1.4 16/04/15 symlink to host genome instead of copying to local directory
# Version 0.1.5 19/05/15 change underlying viralDetContigs to avoid blast segmentation fault (limit output hits to 5)
#			symlinks to reference database and query contigs/reads
# Version 0.1.6 20/05/15 uses exisiting bwa index if it exists
# Version 0.1.7 29/05/15 blast directly from this script without calling additional python script
# Version 0.2.0 01/06/15 Revised to incorporate IterMap; mapping raw data to best hit(s) from blast search
#			automatically generates output directory and file names based on input file names and 
#			accession numbers
# Version 0.2.1 04/08/15 Remove symlink to search data
# Version 0.2.2 06/10/15 Bugfix for location of search database; using existing blast database if available

# set defaults for the options
KVALUE=101
CUTOFF=auto
OutputDir=Sherlock_output
#BlastlngthThresh=0.5
#BlastSimThresh=0.2
Blast_e_value=0.0001

Start=$(date +%s)

# parse the options
while getopts 'H:e:r:c:k:' opt ; do
  case $opt in
#    s) BlastSimThresh=$OPTARG ;;
    H) PathToHOST=$OPTARG ;;
    e) Blast_e_value=$OPTARG ;;
    r) PathToSearchData=$OPTARG ;;
#    l) BlastlngthThresh=$OPTARG ;;
#    o) OutputDir=$OPTARG ;;
    c) CUTOFF=$OPTARG ;;
    k) KVALUE=$OPTARG ;;
  esac
done
# skip over the processed options
shift $((OPTIND-1)) 
# check for mandatory positional parameters
if [ $# -lt 2 ]; then
  echo "
Usage: $0 [options] -H path to HOST genome -r path to Reference Database Read1.fastq.gz Read2.fastq.gz"
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
	refname=${ref%%.*}

	host=$(basename "$PathToHOST")
	hostname=${host%%.*}

# Create output OutputDirectory
OutputDir="$samplename"_Sherlock_"$refname"
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
#cat "$OutputDir"/"$samplename"_nonHost_"$KVALUE"/*.fa > "$OutputDir"/"$samplename"_nonHost_"$KVALUE"/contigs_singletons.fa

# Now run Blast of output contigs (potentially including unassembled singleton reads) against database of choice.  If you have some idea of what you are looking for you can limit the search to a particular taxonomic group (e.g. Rhabdoviridae) providing you have a database for it.

cd "$OutputDir"
#ln -s "$(readlink -f "$samplename"_nonHost_"$KVALUE"/contigs_singletons.fa)" contigs_singletons.fa
#ln -s "$(readlink -f "$PathToSearchData")" Reference.fasta
if [ -f "$PathToSearchData".nhr ]
then 	echo "Using existing BLAST database"
	SearchData="$(readlink -f "$PathToSearchData")"
else 	echo "Generating local BLAST datbase" 
	ln -s "$(readlink -f "$PathToSearchData")" SearchData.fa
	makeblastdb -in "$SearchData" -parse_seqids -dbtype nucl
	SearchData="$(readlink -f "SearchData.fa")"
fi

echo "BLASTing contigs against selected database"
blastall -b 5 -p blastn -d "$SearchData" -i "$samplename"_nonHost_"$KVALUE"/contigs.fa -o "$samplename"_"$refname"_crunch.txt -m 8 -e "$Blast_e_value"
# -a "$threads"

# Extract the single sequence from blast output corresponding to the highest scoring match
sort -k12,12 -rn "$samplename"_"$refname"_crunch.txt | head -1 - | awk '{print $2}' - > top_match.txt
AccNo=$(sed -e 's/.*\([A-Z][A-Z][0-9][0-9][0-9][0-9][0-9][0-9]\).*/\1/g'  top_match.txt)
blastdbcmd -db "$SearchData" -entry_batch top_match.txt > "$AccNo".fas

# Map original data to selected reference sequences and generate new consensus
IterMap.sh -i4 "$AccNo".fas "$LEFT" "$RIGHT"


End=$(date +%s)
TimeTaken=$((End-Start))
echo "Results are in: "$OutputDir""
echo  | awk -v D=$TimeTaken '{printf "Performed Sherlock Analysis in: %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'
