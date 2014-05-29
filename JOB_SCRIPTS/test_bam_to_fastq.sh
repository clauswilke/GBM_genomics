#! /bin/bash

## here, I am submitting a job that will do two things
## first, it will test the piccard command, and 
## secondly it will create a fastq of one of the files 26 tumor samples in /scratch/00001/mattcowp/dakota/
## so that I can use the fastq to test the rest of the script in the development queue (this is the slow step)
## also, this is meant to run on stampede

#!/bin/bash

#$ -V
#$ -N rc         			# Job Name
#$ -j n                     # Combine stderr and stdout
#$ -o rc.err     			# error
#$ -o rc.out     			# Name of the output file (eg. myMPI.oJobID)
#$ -q normal                # Queue name "normal"
#$ -pe 12way 12           	# one node, 12 independent cores
#$ -A A-bio7
#$ -l h_rt=24:00:00         # Run time (hh:mm:ss) - 1.5 hours

set -x                      # Echo commands

echo `hostname`

# unset any MPI/OpenMP thread affinities
export MV2_USE_AFFINITY=0
export MV2_ENABLE_AFFINITY=0
export VIADEV_USE_AFFINITY=0
export VIADEV_ENABLE_AFFINITY=0

# Run the executable
cd /scratch/00001/mattcowp/dakota
bash $HOME/GBM_genomics/MAKE_VCF/12_batchFastExon.sh C484.TCGA-02-2486-01A-01D-1494-08.5 splitSize batchSize coresPerRun