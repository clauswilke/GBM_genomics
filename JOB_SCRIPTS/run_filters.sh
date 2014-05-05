#!/bin/bash

# SGE options
#$ -S/bin/bash
#$ -N filter0-4
#$ -m beas
#$ -M dakotaz@utexas.edu

# Create a working directory and go into it
WDIR=/state/partition1/$USER/$JOB_NAME-$JOB_ID
mkdir -p $WDIR
if [ ! -d $WDIR ]
then
  echo $WDIR not created
  exit
fi
cd $WDIR

# Copy data and config files
cp $HOME/
cp 

# run the python scripts
python 

# Copy results back to the home directory

# Cleanup
rm -rf $WDIR