#!/bin/bash
set -e

# script by ellisrichardj to generate a new consensus in fasta format from a bam (mapping) file
# requires adapted vcfutils script (vcf2consensus.pl) to properly call indels (original version 
# just ouputs lower quality consensus sequence in indel regions and this information was lost when
# converting from fastq to fasta)

# Currently only labels single-segment fasta file correctly.  Needs more work to use with multi
# fasta reference files

# **requires vcf2consensus.pl, samtools/bcftools and clustalw** 

# Version 0.1.5 27/11/14

# Change Log
# v0.1.1 first version 11/09/14
# v0.1.2 added line to re-index bam after duplicate removal
# v0.1.3 tidy up generation of unnecessary intermeditate files by piping
# v0.1.4 added '-A' switch to samtools mpileup command to allow inclusion of anomalous read pairs
#	 increased the read depth which inhibits indel calling to 2000 (added -L 1000 to samtools mpileup command)
# v0.1.5 Added -E switch to samtools mpileup (alternate to BAQ which appears to lead to missed SNPs)
#	In testing the -E option provided a concensus which reflected visualization of the bam file

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

# remove PCR duplicates from mapping file and re-index
samtools rmdup "$BAM" "$samplename"_clean.bam
samtools index "$samplename"_clean.bam

# generate and correctly label consensus
samtools mpileup -L 2000 -AEuf "$REF" "$samplename"_clean.bam | bcftools view -cg - > "$samplename"_"$refname".vcf
vcf2consensus.pl consensus -f "$REF" "$samplename"_"$refname".vcf  | sed '1s/.*/>'"$samplename"_mappedto_"$refname"'/g' - > "$samplename"_mappedto_"$refname"_consensus.fas

# mapping statistics
samtools flagstat "$samplename"_clean.bam > "$samplename"_MappingStats.txt

# generate pairwise alignment of reference and new concensus
cat "$REF" "$samplename"_mappedto_"$refname"_consensus.fas > unaligned.fas
clustalw -infile=unaligned.fas -outfile="$refname"_"$samplename"_aligned.fas -output=FASTA

rm unaligned.fas
rm *.dnd

