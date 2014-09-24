#!/bin/bash
set -e

# Analysis steps required to find unknown infectious agent using Illumina shotgun sequence data from DNA/RNA extracted from host tissue

# Software required:
	# bwa - 0.7.5a or above
	# samtools
	# bam2fastx (part of bamUtil package)
	# velvet
	# fastx-toolkit
	# Mira (includes fastqselect tool)
	# Trimmomatic
	# viralDetContigs_150514.py
# Ensure these are in PATH or symbolic links are in a folder that is in your PATH


# set defaults for the options
KVALUE=101
CUTOFF=auto
OutputDir=Sherlock_output
BlastlngthThresh=0.5
BlastSimThresh=0.2
Blast_e_value=0.001

# parse the options
while getopts 's:h:e:r:l:o:c:k:' opt ; do
  case $opt in
    s) BlastSimThresh=$OPTARG ;;
    h) PathToHOST=$OPTARG ;;
    e) Blast_e_value=$OPTARG ;;
    r) PathToSearchData=$OPTARG ;;
    l) BlastlngthThresh=$OPTARG ;;
    o) OutputDir=$OPTARG ;;
    c) CUTOFF=$OPTARG ;;
    k) KVALUE=$OPTARG ;;
  esac
done
# skip over the processed options
shift $((OPTIND-1)) 
# check for mandatory positional parameters
if [ $# -lt 2 ]; then
  echo "
Usage: $0 [options] -h path to HOST genome -r path to Reference Database -o Output OutputDirectory name Read1 Read2"
  echo "Options: De Novo assembly: -k Velvet kmervalue | -c Velvet cov_cutoff"
  echo "Options: Blast Output: -l minimum proportion of contig length (default = 0.8) | -s  minimum similarity threshold for output (default = 0.95)
"
  exit 1
fi
LEFT="$1"
RIGHT="$2"

# Create output OutputDirectory
mkdir "$OutputDir"
echo "Getting host genome"
cp "$PathToHOST" "$OutputDir"/HostGenome.fa

# Generate bwa index of host genome
echo "Generating index of host genome" 
/usr/local/bioinf/bwa-0.7.5a/bwa index "$OutputDir"/HostGenome.fa

# Map raw data to host genome
echo "Mapping raw reads to host reference genome"
/usr/local/bioinf/bwa-0.7.5a/bwa mem -t 6 "$OutputDir"/HostGenome.fa "$LEFT" "$RIGHT" | samtools view -Su - | samtools sort - "$OutputDir"/"$OutputDir"_Host_map_sorted

# Extract sequence reads that do not map to the host genome
samtools view -b -f 4 "$OutputDir"/"$OutputDir"_Host_map_sorted.bam > "$OutputDir"/"$OutputDir"_nonHost.bam

# Convert non-mapping reads back to fastq format and extract all pairs
echo "Extracting non-host raw data"
bam2fastx -a -o "$OutputDir"/"$OutputDir"_nonHost.fasta "$OutputDir"/"$OutputDir"_nonHost.bam
fasta_formatter -t -i "$OutputDir"/"$OutputDir"_nonHost.fasta -o "$OutputDir"/"$OutputDir"_Sequences_tab.txt
awk < "$OutputDir"/"$OutputDir"_Sequences_tab.txt '{print $1}' > "$OutputDir"/"$OutputDir"_Sequences.lst
rm "$OutputDir"/"$OutputDir"_nonHost.fasta
rm "$OutputDir"/"$OutputDir"_Sequences_tab.txt
rm "$OutputDir"/HostGenome.*

gunzip -c "$LEFT" > All_R1.fastq
fastqselect -infile All_R1.fastq -name "$OutputDir"/"$OutputDir"_Sequences.lst -outfile "$OutputDir"/"$OutputDir"_R1_nonHost.fastq
gunzip -c "$RIGHT" > All_R2.fastq
fastqselect -infile All_R2.fastq -name "$OutputDir"/"$OutputDir"_Sequences.lst -outfile "$OutputDir"/"$OutputDir"_R2_nonHost.fastq
rm All_R1.fastq
rm All_R2.fastq

# Quality trimming and Denovo assembly of non-host reads
let Trimlgth="$KVALUE"+10
trimmomatic-0.30.jar PE -threads 6 -phred33 "$OutputDir"/"$OutputDir"_R1_nonHost.fastq "$OutputDir"/"$OutputDir"_R2_nonHost.fastq  "$OutputDir"/R1_nonhost_trim_paired.fastq "$OutputDir"/R1_unpaired.fastq "$OutputDir"/R2_nonhost_trim_paired.fastq "$OutputDir"/R2_unpaired.fastq ILLUMINACLIP:/home/sequence/ReferenceSequences/adapter.fasta:1:30:10:5:true MAXINFO:"$Trimlgth":0.5 MINLEN:"$Trimlgth"
rm "$OutputDir"/R1_unpaired.fastq
rm "$OutputDir"/R2_unpaired.fastq
/opt/velvet/velveth "$OutputDir"/"$OutputDir"_"$KVALUE" "$KVALUE" -shortPaired -fmtAuto "$OutputDir"/R1_nonhost_trim_paired.fastq "$OutputDir"/R2_nonhost_trim_paired.fastq
/opt/velvet/velvetg "$OutputDir"/"$OutputDir"_"$KVALUE" -exp_cov auto -cov_cutoff "$CUTOFF" -ins_length 300 -clean yes -read_trkg yes -amos_file yes -unused_reads yes
cat "$OutputDir"/"$OutputDir"_"$KVALUE"/*.fa > "$OutputDir"/"$OutputDir"_"$KVALUE"/contigs_singletons.fa

# Now run Blast of output contigs and singleton reads against database of choice.  If you have some idea of what you are looking for you can limit the search to a particular taxonomic group (e.g. Rhabdoviridae).

cd "$OutputDir"
cp "$OutputDir"_"$KVALUE"/contigs_singletons.fa .
cp "$PathToSearchData" ./Reference.fasta
viralDetContigs_150514.py "$PWD" Reference.fasta "$Blast_e_value" contigs_singletons.fa "$BlastlngthThresh" "$BlastSimThresh"
echo "Results are in: "$OutputDir""
