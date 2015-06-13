#!/bin/bash
set -e

# script by ellisrichardj to iteratively map and call consensus, and then remap to
# new consensus, thereby theoretically increasing coverage and accuracy of consensus,
# especially when the consensus is likely to be divergent from original reference

# Version 0.1.1 06/10/14 Initial verison
# Version 0.1.2 14/10/14 Reduced mapping stringency for first iteration (k, B and O options for bwa mem), allowed 
#	inclusion of anomalous read pairs in variant calling
# Version 0.1.4 24/11/14 increased the read depth which inhibits indel calling to 10000 from default of 250 (added -L
#	10000 to samtools mpileup command)
# Version 0.1.5 26/11/14 Added -E switch to samtools mpileup (alternate to BAQ which appears to lead to missed SNPs)
#	In testing the -E option provided a concensus which reflected visualization of the bam file
# Version 0.1.6 06/02/15 Create symlinks for reference and data files rather than copyping data into new directory and 
#	use default bwa parameters (normal stringency) if performing just a single iteration
# Version 0.1.7 05/03/15 Bug fix for symlinks
# Version 0.1.8 20/03/15 Added option to specify minimum expected coverage (this will alter depth variable when calling
#	consensus from vcf); generates log file to record commands and paramaters used
# Version 0.1.9 24/03/15 added -A option to bcftools view command which allows all possible variants to be retained in
#	output; reduced bwa -T value to 10 to allow output of lower quality alignments; added -Q0 to mpileup (min map
#	quality of bases); added -p 0.2 to bcftools view
# Version 0.2.1	14/04/15 only remove unaligned reads from final cleaned BAM - this ensures statistics reflect proportion
#	of mapped reads
# Version 0.2.2 17/04/15 Add -c to bcftools view command (SNP calling)
# Version 0.2.3 22/04/15 Changed mpileup to -Q0, removed bcftools view -c and -p 0.2, added bcftools -cg back in
# Version 0.2.4 02/06/15 Appends iteration count to end of each fasta header
# Version 0.2.5 12/06/15 Cleaned up logging of commands by defining LogRun function rather than using echo

# set defaults for the options

iter=1
minexpcov=5
minQ=20

# parse the options
while getopts 'i:c:q:' opt ; do
  case $opt in
    i) iter=$OPTARG ;;
    c) minexpcov=$OPTARG ;;
    q) minQ=$OPTARG ;;
  esac
done
# skip over the processed options
shift $((OPTIND-1)) 

# check for mandatory positional parameters
if [ $# -lt 3 ]; then
  echo "
Usage: $0 [-i # iterations] [-c minimum expected coverage] [-q minimum RMS quality] <path to Reference> <path to R1 fastq> <path to R1 fastq> 
"
exit 1
fi

Ref=$1
R1=$2
R2=$3
now=$(date '+%x %R')
Start=$(date +%s)

	sfile1=$(basename "$R1")
	sfile2=$(basename "$R2")
	samplename=${sfile1%%_*}

	ref=$(basename "$Ref")
	refname=${ref%%_*}
	reffile=${ref%%.*}

mkdir "$samplename"_IterMap"$iter"
ln -s "$(readlink -f "$Ref")" "$samplename"_IterMap"$iter"/"$reffile".fas
ln -s "$(readlink -f "$R1")" "$samplename"_IterMap"$iter"/R1.fastq.gz
ln -s "$(readlink -f "$R2")" "$samplename"_IterMap"$iter"/R2.fastq.gz
cd "$samplename"_IterMap"$iter"

rfile="$reffile".fas
count=1
if [ $minexpcov -lt 5 ]; then depth=1; else depth=10; fi
threads=$(grep -c ^processor /proc/cpuinfo)

# Log commands that are run
echo -e "$now\n\tItermap v0.2.4 running with $threads cores\n\tThe following commands were run:\n"  > "$samplename"_IterMap"$iter".log
# Define function to log commands and then run them
LogRun(){
echo "$@" >> "$samplename"_IterMap"$iter".log
eval "$@"
}

	while (($count <= $iter))
	do
	# Set reduced mapping stringency for first iteration
	if [ $count == 1 ] && [ $iter != 1 ]; then
		mem=16
		mmpen=2
		gappen=4
	else
		mem=18
		mmpen=4
		gappen=6
	fi

	# mapping to original reference or most recently generated consensus
	LogRun bwa index "$rfile"
	LogRun bwa mem -T10 -t "$threads" -k "$mem" -B "$mmpen" -O "$gappen" "$rfile" R1.fastq.gz R2.fastq.gz |
	 samtools view -Su - |
	 samtools sort - "$samplename"-"$reffile"-iter"$count"_map_sorted

if [ $count == $iter ]; then
	# generate and correctly label consensus using cleaned bam on final iteration
	LogRun samtools rmdup "$samplename"-"$reffile"-iter"$count"_map_sorted.bam "$samplename"-"$reffile"-iter"$count"_clean.bam
	LogRun samtools view -bF4 -o "$samplename"-"$reffile"-iter"$count"_clean_mapOnly.bam "$samplename"-"$reffile"-iter"$count"_clean.bam
	LogRun samtools index "$samplename"-"$reffile"-iter"$count"_clean_mapOnly.bam

	LogRun samtools mpileup -L 10000 -Q0 -AEupf "$rfile" "$samplename"-"$reffile"-iter"$count"_clean.bam |
	 bcftools view -cg - > "$samplename"-"$reffile"-iter"$count".vcf
	LogRun vcf2consensus.pl consensus -d "$minexpcov" -Q "$minQ" -f "$rfile" "$samplename"-"$reffile"-iter"$count".vcf |
	 sed '/^>/ s/-iter[0-9]//;/^>/ s/$/'-iter"$count"'/' - > "$samplename"-"$reffile"-iter"$count"_consensus.fas

	# mapping statistics
	LogRun samtools flagstat "$samplename"-"$reffile"-iter"$count"_clean.bam > "$samplename"-"$reffile"-iter"$count"_MappingStats.txt
	rfile="$samplename"-"$reffile"-iter"$count"_consensus.fas

else

	LogRun samtools mpileup -L 10000 -Q0 -AEupf "$rfile" "$samplename"-"$reffile"-iter"$count"_map_sorted.bam |
	 bcftools view -cg - > "$samplename"-"$reffile"-iter"$count".vcf
	LogRun vcf2consensus.pl consensus -d "$minexpcov" -Q "$minQ" -f "$rfile" "$samplename"-"$reffile"-iter"$count".vcf |
	 sed '/^>/ s/-iter[0-9]//;/^>/ s/$/'-iter"$count"'/' - > "$samplename"-"$reffile"-iter"$count"_consensus.fas

	# mapping statistics
	LogRun samtools flagstat "$samplename"-"$reffile"-iter"$count"_map_sorted.bam > "$samplename"-"$reffile"-iter"$count"_MappingStats.txt
	rfile="$samplename"-"$reffile"-iter"$count"_consensus.fas

fi
	((count=count+1))
	echo "New Consensus: "$rfile""
done

complete=$(date '+%x %R')
echo -e "Completed processing $complete"  >> "$samplename"_IterMap"$iter".log

# generate pairwise alignment of reference and each iteration of the concensus
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
