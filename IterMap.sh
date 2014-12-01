#!/bin/bash
set -e

# script by ellisrichardj to iteratively map and call consensus, and then remap to
# new consensus, thereby theoretically increasing coverage and accuracy of consensus
# especially when consensus is likely to be divergent from original reference

# Version 0.1.1 06/10/14 First verison
# Version 0.1.2 14/10/14 Reduced mapping stringency for first iteration (k, B and O options for bwa mem), allowed 
#	inclusion of anomalous read pairs in variant calling
# Version 0.1.4 24/11/14 increased the read depth which inhibits indel calling to 10000 from default of 250 (added -L
#	2000 to samtools mpileup command)
# Version 0.1.5 26/11/14 Added -E switch to samtools mpileup (alternate to BAQ which appears to lead to missed SNPs)
#	In testing the -E option provided a concensus which reflected visualization of the bam file

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
Usage: $0 [-i # iterations] <path to Reference> <path to R1 fastq> <path to R1 fastq> 
"
exit 1
fi

Ref=$1
R1=$2
R2=$3
Start=$(date +%s)

	sfile1=$(basename "$R1")
	sfile2=$(basename "$R2")
	samplename=${sfile1%%_*}

	ref=$(basename "$Ref")
	refname=${ref%%_*}
	reffile=${ref%%.*}

mkdir "$samplename"_IterMap"$iter"
cp "$Ref" "$samplename"_IterMap"$iter"/"$reffile".fas
cp "$R1" "$samplename"_IterMap"$iter"
cp "$R2" "$samplename"_IterMap"$iter"
cd "$samplename"_IterMap"$iter"

rfile="$reffile".fas
count=1
threads=$(grep -c ^processor /proc/cpuinfo)

while (($count <= $iter))
do
# Set reduced mapping stringency for first iteration
if [ $count == 1 ]; then
	mem=16
	mmpen=2
	gappen=4
else
	mem=19
	mmpen=4
	gappen=6
fi

	# mapping to original reference or most recently generated consensus
	bwa index "$rfile"
	bwa mem -t "$threads" -k "$mem" -B "$mmpen" -O "$gappen" "$rfile" "$R1" "$R2" | samtools view -Su - | samtools sort - "$samplename"-"$reffile"-iter"$count"_map_sorted

if [ $count == $iter ]; then
	# generate and correctly label consensus using cleaned bam on final iteration
	samtools rmdup "$samplename"-"$reffile"-iter"$count"_map_sorted.bam "$samplename"-"$reffile"-iter"$count"_clean.bam
	samtools index "$samplename"-"$reffile"-iter"$count"_clean.bam

samtools mpileup -L 10000 -AEuf "$rfile" "$samplename"-"$reffile"-iter"$count"_clean.bam | bcftools view -cg - > "$samplename"-"$reffile"-iter"$count".vcf
	vcf2consensus.pl consensus -f "$rfile" "$samplename"-"$reffile"-iter"$count".vcf | sed '1s/.*/>'"$samplename"-"$reffile"-iter"$count"'/g' - > "$samplename"-"$reffile"-iter"$count"_consensus.fas

	# mapping statistics
	samtools flagstat "$samplename"-"$reffile"-iter"$count"_clean.bam > "$samplename"-"$reffile"-iter"$count"_MappingStats.txt
	rfile="$samplename"-"$reffile"-iter"$count"_consensus.fas

else

	samtools mpileup -L 10000  -AEuf "$rfile" "$samplename"-"$reffile"-iter"$count"_map_sorted.bam | bcftools view -cg - > "$samplename"-"$reffile"-iter"$count".vcf
	vcf2consensus.pl consensus -f "$rfile" "$samplename"-"$reffile"-iter"$count".vcf | sed '1s/.*/>'"$samplename"-"$reffile"-iter"$count"'/g' - > "$samplename"-"$reffile"-iter"$count"_consensus.fas

	# mapping statistics
	samtools flagstat "$samplename"-"$reffile"-iter"$count"_map_sorted.bam > "$samplename"-"$reffile"-iter"$count"_MappingStats.txt
	rfile="$samplename"-"$reffile"-iter"$count"_consensus.fas

fi
	((count=count+1))
	echo "New Consensus: "$rfile""
done

# generate pairwise alignment of reference and each new concensus
cat *.fas > unaligned.fas
clustalw -infile=unaligned.fas -outfile=Increments_aligned.fas -output=FASTA

rm unaligned.fas
rm *.dnd
rm *.gz

End=$(date +%s)
TimeTaken=$((End-Start))
echo "Results are in "$samplename"_IterMap"$iter""
echo "New consensus after "$iter" iterations: "$rfile""
echo  | awk -v D=$TimeTaken '{printf "Performed '$iter' mapping iterations in: %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'
