#!/bin/sh
#    oncreate.sh /PATH/TO/CREATED/FILE
#
#    Handler for files found by illumina-cron.sh
#    Converts .bcl data generated by illumina NextSeq machines into .fastq
#    writing to /PATH/TO/CREATED/
#    Output from the bcl2fastq command is stored in /PATH/TO/CREATED/bcl2fastq.log
#    bcl2fastq runs when the CREATED directory name contains a sequencer id starting with N

FILEPATH=$1
DIRECTORY=$(dirname "$FILEPATH")

# Parse. Expected format: <job_id>_<sequencer>_<unknown>_<unknown>
BASEDIR=$(echo $DIRECTORY | rev | cut -d '/' -f1 | rev)
SEQUENCER=$(echo $BASEDIR | cut -d _ -f2)

# Skip non-nextseq runs
if [[ ${SEQUENCER::1} != N ]]; then
    echo non-NextSeq machine, skipping $BASEDIR
    exit 1
fi

# Print job info
echo "*** NextSeq Run Detected ***"
date
echo Run Directory: $DIRECTORY 
echo Sequencer: $SEQUENCER
echo Converting bcl to fastq...
echo

# Convert .bcl to .fastq
cd $DIRECTORY

/usr/local/bin/bcl2fastq --sample-sheet ./*.csv --no-lane-splitting --output-dir $PWD >> $PWD/bcl2fastq.log 2>> $PWD/bcl2fastq.err.log

echo exit $DIRECTORY
