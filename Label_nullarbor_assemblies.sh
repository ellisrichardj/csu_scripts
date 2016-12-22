#!/bin/bash
set -e

# script by ellisrichardj
# This will rename contig files that are output from nullabor with the sample name
# rather than a numerical identifier.  Assumes that each assembly is in a subdirectory
# with the sample name

# Version 0.1 16/09/16

# check for mandatory positional parameters
if [ $# != 1 ]; then
  echo "Usage: $0 {path to Nullarbor Results Folder} "
exit 1
fi

DataFolder="$(readlink -f "$1")"
mkdir "$DataFolder"/Assemblies

for file in "$DataFolder"/*/contigs.fa
do
	samdir=$(dirname "$file")
	samplename=${samdir##*/}
	
	echo "$samplename"
	mv "$DataFolder"/"$samplename"/contigs.fa "$DataFolder"/Assemblies/"$samplename"_contigs.fa
done


