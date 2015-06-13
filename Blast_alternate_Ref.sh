#!/bin/bash
set -e
# set our defaults for the options
BlastlngthThresh=0.5
BlastSimThresh=0.8
BlastEvalue=0.001
# parse the options
while getopts 'r:l:s:e:' opt ; do
  case $opt in
    r) PathToSearchData=$OPTARG ;;
    l) BlastlngthThresh=$OPTARG ;;
    s) BlastSimThresh=$OPTARG ;;
    e) BlastEvalue=$OPTARG ;;
  esac
done
# skip over the processed options
shift $((OPTIND-1)) 
# check for mandatory positional parameters
if [ $# -lt 1 ]; then
  echo "
Usage: $0 [options] -r path to Reference Database outputFolder"
    echo "Options: Blast Output: -l minimum proportion of contig length (default = 0.5) | -s  minimum similarity threshold for output (default = 0.8)
"
  exit 1
fi
DIR="$1"

# Now run Blast of output contigs file against database of choice - looking through the matches for each contig should reveal something interesting.  If you have some idea of what you are looking for you can limit the search to a particular taxonomic group (e.g. viruses).  If this draws a blank try Blasting the complete non-host fasta file, but this will take much longer to trawl through the results
echo "Getting Reference Data"
mkdir "$DIR"
cp contigs_singletons.fa "$DIR"
cd "$DIR"
cp "$PathToSearchData" ./Reference.fasta
viralDetContigs_new.py "$PWD" Reference.fasta "$BlastEvalue" contigs_singletons.fa  "$BlastlngthThresh" "$BlastSimThresh"



