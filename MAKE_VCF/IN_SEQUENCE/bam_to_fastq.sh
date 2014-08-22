#!/bin/bash

# this script takes a TCGA .bam (aligned to hg18) and converts it to a pair of unaligned .fastq files

module load samtools
module load picard

# these are command line arguments; pfx is the sample name (the TCGA barcode)
pfx=$1
wrk_dir=$2
queue="normal"

# 1. pfx.bam to pfx.fastq conversion
# cd /scratch/00001/mattcowp/dakota
cd $wrk_dir
cd $pfx
echo "$pfx:" `date` > ../$pfx.log

echo "Running SamToFastQ on $pfx...."  `date` >> ../$pfx.log
(java -d64 -Xmx4g -jar /opt/apps/picard/1.98/SamToFastq.jar I=$pfx.bam F=$pfx.1.fastq F2=$pfx.2.fastq 2>&1) >> $pfx\_SamToFastQ.log
echo "SamToFastQ done: " `date` >> ../$pfx.log