#! /bin/bash

#SBATCH -J SNPcalling
#SBATCH -o SNPcalling.out
#SBATCH -e SNP calling.err
#SBATCH -n 16
#SBATCH -p normal
#SBATCH -A mattcowp
#SBATCH -t 24:00:00

set -x                      # Echo commands

echo `hostname`

## set up input variables
## pfx is the prefix (sample name) that is common between the tumor and blood files for the sample
## blood is the end of the blood file name, such that the blood file is called $pfx$blood
## tumor is the end of the tumor file name, such that the tumor file is called $pfx$tumor
## SS, MT, ST, and VS refer to each of the four possible SNP callers
##		the value for each of these variables should be "yes" or "no," depending whether you want that SNP-caller to run
pfx=""
tumor=""
blood=""
SS=""
MT=""
ST=""
VS=""

## Run the executable
bash $HOME/GBM_genomics/MAKE_VCF/IN_SEQUENCE/SNPcalling.sh $pfx $tumor $blood $SS $MT $ST $VS
