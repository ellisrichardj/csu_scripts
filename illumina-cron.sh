#!/bin/sh
#     illumina-cronjob.sh
#
#     Prints a cron command for automated conversion of bcl into fastq for NextSeq machines
#     The command schedules a find on WATCHPATH every minute and calls `oncreate.sh <filepath>` on files
#     created within the last minute. 
#     Output from oncreate.sh is stored in ./illumina-cron.log
#     
#     Usage:
#       (1) call `sh illumina_cronjob.sh`
#       (2) copy output to `crontab -e`

WATCHPATH=/illumina/runs/
FILENAME=CopyComplete.txt

DIRNAME=$PWD/$(dirname "$0")
echo "* * * * * find $WATCHPATH/ -maxdepth 2 -cmin 1 -name $FILENAME -exec sh $DIRNAME/oncreate.sh {} \; >> $DIRNAME/illumina-cron.log 2>> $DIRNAME/illumina-cron.err.log"
