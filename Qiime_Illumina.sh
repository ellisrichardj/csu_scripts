#!/bin/bash

Start=$(date +%s)

# set defaults for the options
OutputDir="$PWD"/preprocessed

# parse the options
while getopts 'o:' opt ; do
  case $opt in
    o) OutputDir=$OPTARG ;;
  esac
done

# skip over the processed options
shift $((OPTIND-1)) 
# check for mandatory positional parameters
if [ $# -lt 1 ]; then
  echo "
Usage: $0 [options] path to Data Directory [required]
	Options: -o Path to Output directory (default = $PWD/preprocessed)
	"
  exit 1
fi

PathToData="$(readlink -f "$1")"
PathToOutput="$(readlink -f "$OutputDir")"

threads=$(grep -c ^processor /proc/cpuinfo)

~/microbiome_helper/run_pear.pl -p $threads -o $PathToOutput/assembled $PathToData/*.fastq.gz

rm $PathToOutput/assembled/*.discarded.fastq
rm $PathToOutput/assembled/*.unassembled.*.fastq

~/microbiome_helper/readFilter.pl -f GTGYCAGCMGCCGCGGTA -pc forward -q 25 -p 95 -l 400 -t $threads -b /usr/local/bioinf/BBMap/bbmap -o $PathToOutput/filtered_reads $PathToOutput/assembled/*.assembled.fastq

~/microbiome_helper/run_fastq_to_fasta.pl -p $threads -o $PathToOutput/fasta $PathToOutput/filtered_reads/*.fastq

rm -rf $PathToOutput/assembled/
rm -rf $PathToOutput/filtered_reads/

~/microbiome_helper/create_qiime_map.pl $PathToOutput/fasta/* > $PathToOutput/map.txt

qiime add_qiime_labels.py -i $PathToOutput/fasta -m $PathToOutput/map.txt -c FileInput -o $PathToOutput/combined_fasta

qiime pick_open_reference_otus.py -i $PathToOutput/combined_fasta/combined_seqs.fna -o $PathToOutput/OTUs_out -aO $threads -p ~/Documents/Qiime_parameters.txt

End=$(date +%s)
TimeTaken=$((End-Start))
echo "Results are in: "$PathToOutput"
echo  | awk -v D=$TimeTaken '{printf "Performed Qiime Illumina processing in: %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'
