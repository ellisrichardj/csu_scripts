#!/bin/bash
set -e

# script by ellisrichardj
# This will merge all files with matching names in each directory
# Useful for merging data from 2 (or possible more) Illumina runs of the same samples
# Merged files will be written to current working directory

# check for mandatory positional parameters
if [ $# -lt 2 ]; then
  echo "
Usage: $0 {path to Folder 1}  {path to Folder 2}"
exit 1
fi

Folder1="$1"
Folder2="$2"

for file in "$Folder1"/*.gz
do
	fname=$(basename "$file")
	echo "Merging "$fname""
	cat "$Folder1"/"$fname" "$Folder2"/"$fname" > $PWD/$fname
done

