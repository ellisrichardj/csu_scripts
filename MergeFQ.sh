#!/bin/bash
set -e

# script by ellisrichardj
# This will merge all files with matching names in each directory
# Useful for merging data from 2 (or possible more) Illumina runs of the same samples
# Merged files will be written to current working directory

# Version 0.0.1 06/11/14

# check for mandatory positional parameters
if [ $# -lt 2 ]; then
  echo "
Usage: $0 {path to Folder 1}  {path to Folder 2}"
exit 1
fi

Dir1="$1"
Dir2="$2"

for file in "$Dir1"/*.gz
do
   	fname=$(basename "$file")
   if [ -f "$Dir2"/"$fname" ]; then
	echo "Merging "$fname""
	cat "$Dir1"/"$fname" "$Dir2"/"$fname" > $PWD/$fname
   else
	echo "Matching file "$fname" does not exist in both directories: skipping"
   fi
done

