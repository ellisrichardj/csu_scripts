#!/bin/bash
set -e

# script by ellisrichardj
# generates tab de-limited file indicating samples and paired end sequence files available in a given directory
# primarily generates input for Nullabor (https://github.com/tseemann/nullarbor), but may be useful for other things....

# skip over the processed options
shift $((OPTIND-1)) 

# check for mandatory positional parameters
if [ $# -lt 1 ]; then
  echo "Usage: $0 {path to Data Folder} "
exit 1
fi


DataFolder="$1"

for file in "$DataFolder"/*_R1*.gz
do
	fname=$(basename "$file")
	samplename=${fname%%_*}
	file2=${file/R1/R2}
	echo -e "$samplename""\t""$file""\t"$file2"" >> samples.tab
done
