#! /bin/bash

#SBATCH -J t1_3.1
#SBATCH -o t1_3.1.out
#SBATCH -e t1_3.1.err
#SBATCH -n 16
#SBATCH -p normal
#SBATCH -A mattcowp
#SBATCH -t 24:00:00

set -x                      # Echo commands

echo `hostname`

# Run the executable
bash $HOME/GBM_genomics/MAKE_VCF/IN_SEQUENCE/align_fastq.sh C484.TCGA-02-2486-01A-01D-1494-08.5