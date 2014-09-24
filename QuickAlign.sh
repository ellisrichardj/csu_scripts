#!/bin/bash
set -e

# script by ellisrichardj to generate alignment of a new concensus with the original reference

# Version 0.1 10/09/14

# check for mandatory positional parameters
if [ $# -lt 2 ]; then
  echo "
Usage: $0 {path to Reference}  {path to new consensus}"

  exit 1
fi

REF="$1"
NEW="$2"

refname="${REF%.*}"
samplename="${NEW%_*}"

cat "$REF" "$NEW" > unaligned.fas
clustalw -infile=unaligned.fas -outfile="$refname"_"$samplename"_aligned.fas

rm unaligned.fas
rm *.dnd
