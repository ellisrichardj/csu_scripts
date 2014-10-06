#!/bin/bash
set -e

# script by ellisrichardj to iteratively map and call consensus, and then remap to
# new consensus, thereby theoretically increasing coverage and accuracy of consensus
# especially when consensus is likely to be divergent from original reference

# Version 0.1.1 06/10/14

# set defaults for the options
iter=1

# parse the options
while getopts 'i:' opt ; do
  case $opt in
    i) iter=$OPTARG ;;
  esac
done
# skip over the processed options
shift $((OPTIND-1)) 

# check for mandatory positional parameters
if [ $# -lt 3 ]; then
  echo "
Usage: $0 [-i # iterations] <path to Reference> <path to R1 fastq> <path to R1 fastq> "
exit 1
fi

Ref=$1
R1=$2
R2=$3
Start=$(date +%s)

	sfile1=$(basename "$R1")
	sfile2=$(basename "$R2")
	samplename=${sfile1%%_*}

	rfile=$(basename "$Ref")
	refname=${rfile%%_*}

mkdir "$samplename"_IterMap"$iter"
cp "$Ref" "$samplename"_IterMap"$iter"
cp "$R1" "$samplename"_IterMap"$iter"
cp "$R2" "$samplename"_IterMap"$iter"
cd "$samplename"_IterMap"$iter"

count=1
while (($count <= $iter))
do
	# mapping to original reference or most recently generated consensus
	bwa index "$rfile"
	bwa mem -t 6 "$rfile" "$R1" "$R2" | samtools view -Su - | samtools sort - "$samplename"-"$refname"-iter"$count"_map_sorted

if [ $count == $iter ]; then
	# generate and correctly label consensus using cleaned bam on final iteration
	samtools rmdup "$samplename"-"$refname"-iter"$count"_map_sorted.bam "$samplename"-"$refname"-iter"$count"_clean.bam
	samtools index "$samplename"-"$refname"-iter"$count"_clean.bam

samtools mpileup -uf "$rfile" "$samplename"-"$refname"-iter"$count"_clean.bam | bcftools view -cg - > "$samplename"-"$refname"-iter"$count".vcf
	vcf2consensus.pl consensus -f "$rfile" "$samplename"-"$refname"-iter"$count".vcf | sed '1s/.*/>'"$samplename"-"$refname"-iter"$count"'/g' - > "$samplename"-"$refname"-iter"$count"_consensus.fas

	# mapping statistics
	samtools flagstat "$samplename"-"$refname"-iter"$count"_clean.bam > "$samplename"-"$refname"-iter"$count"_MappingStats.txt
	rfile="$samplename"-"$refname"-iter"$count"_consensus.fas

else

	samtools mpileup -uf "$rfile" "$samplename"-"$refname"-iter"$count"_map_sorted.bam | bcftools view -cg - > "$samplename"-"$refname"-iter"$count".vcf
	vcf2consensus.pl consensus -f "$rfile" "$samplename"-"$refname"-iter"$count".vcf | sed '1s/.*/>'"$samplename"-"$refname"-iter"$count"'/g' - > "$samplename"-"$refname"-iter"$count"_consensus.fas

	# mapping statistics
	samtools flagstat "$samplename"-"$refname"-iter"$count"_map_sorted.bam > "$samplename"-"$refname"-iter"$count"_MappingStats.txt
	rfile="$samplename"-"$refname"-iter"$count"_consensus.fas

fi
	((count=count+1))
	echo "$rfile"
done

# generate pairwise alignment of reference and each new concensus
cat "$Ref" *.fas > unaligned.fas
clustalw -infile=unaligned.fas -outfile=Increments_aligned.fas -output=FASTA

rm unaligned.fas
rm *.dnd
rm *.gz

End=$(date +%s)
TimeTaken=$((End-Start))
echo "Results are in "$samplename"_IterMap"$iter""
echo  | awk -v D=$TimeTaken '{printf "Performed '$iter' mapping iterations in: %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'
