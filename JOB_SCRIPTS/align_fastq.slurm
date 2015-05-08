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

data_dir="/scratch/00001/mattcowp/dakota/"

# Run the executable (tumor)
cd $data_dir/$pfx
bash $HOME/GBM_genomics/MAKE_VCF/IN_SEQUENCE/align_fastq.sh C484.TCGA-06-0125-01A-01D-1490-08.6 tumor

# run the executable (normal)
# cd $data_dir/$pfx/normal
# bash $HOME/GBM_genomics/MAKE_VCF/IN_SEQUENCE/align_fastq.sh C484.TCGA-06-0125-10A-01D-1490-08.6 blood