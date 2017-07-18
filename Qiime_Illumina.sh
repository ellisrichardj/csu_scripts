#!/bin/bash

# script by ellisrichardj
# Automated processing of Illumina 16S data
# This will process a directory of Illumina data accoring to the microbiome_helper SOP
# The output will be ready for further QIIME analysis

# Requirements:
#	QIIME
#	microbiome_helper (latest version available at: https://github.com/LangilleLab/microbiome_helper)
#	Pear (latest version available at: https://github.com/frederic-mahe/PEAR)

# Version 0.1.0 18/07/17 - Comments added and minor bug fixes


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

echo "Created map file"
echo ""
echo "Creating combined fastq"
echo ""

qiime add_qiime_labels.py -i $PathToOutput/fasta -m $PathToOutput/map.txt -c FileInput -o $PathToOutput/combined_fasta

echo "Now Using Qiime to pick otus."
echo "This may take some time!"
echo ""

qiime pick_open_reference_otus.py -i $PathToOutput/combined_fasta/combined_seqs.fna -o $PathToOutput/OTUs_out -p ~/Documents/Qiime_parameters.txt

End=$(date +%s)
TimeTaken=$((End-Start))
echo "Results are in: "$PathToOutput"
echo  | awk -v D=$TimeTaken '{printf "Performed Illumina pre-processing and OTU picking in: %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'"
