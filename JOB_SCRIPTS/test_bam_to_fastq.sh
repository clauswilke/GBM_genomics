#! /bin/bash

#SBATCH -J t1_3.1
#SBATCH -o t1_3.1.out
#SBATCH -e t1_3.1.err
#SBATCH -n 96
#SBATCH -p normal
#SBATCH -A mattcowp
#SBATCH -t 24:00:00

set -x                      # Echo commands

echo `hostname`

# Run the executable
cd /scratch/00001/mattcowp/dakota
bash $HOME/GBM_genomics/MAKE_VCF/12_batchFastExon.sh C484.TCGA-06-0190-01A-01D-1491-08.11 2000000 32 4