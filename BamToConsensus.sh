#!/bin/bash
set -e

# check for mandatory positional parameters
if [ $# -lt 2 ]; then
  echo "
Usage: $0 {path to Reference}  {path to Bam}"

  exit 1
fi
REF="$1"
BAM="$2"

refname="${REF%.*}"
samplename="${BAM%_*}"

# remove PCR duplicates from mapping file
samtools rmdup "$BAM" "$samplename"_clean.bam

# generate and correctly label consensus
samtools mpileup -E -uf "$REF" "$samplename"_clean.bam | bcftools view -cg - | /usr/share/samtools/vcfutils.pl vcf2fq > "$samplename"_mappedto_"$refname"_consensus.fastq
fastq_to_fasta.py "$samplename"_mappedto_"$refname"_consensus.fastq
sed '1s/.*/>'"$samplename"_mappedto_"$refname"'/g' "$samplename"_mappedto_"$refname"_consensus.fasta > "$samplename"_mappedto_"$refname"_consensus.fas

# delete unwanted files
rm "$samplename"_mappedto_"$refname"_consensus.fastq
rm "$samplename"_mappedto_"$refname"_consensus.qual
rm "$samplename"_mappedto_"$refname"_consensus.fasta

# mapping statistics
samtools flagstat "$samplename"_clean.bam > "$BAM"_stats.txt

