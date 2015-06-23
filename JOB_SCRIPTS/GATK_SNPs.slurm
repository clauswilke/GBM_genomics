#! /bin/bash

#SBATCH -J t1_3.1
#SBATCH -o t1_3.1.out
#SBATCH -e t1_3.1.err
#SBATCH -n 16
#SBATCH -p normal
#SBATCH -A mattcowp
#SBATCH -t 20:00:00

set -x                      # Echo commands

echo `hostname`

# Run the executables
cd /scratch/00001/mattcowp/dakota
bash /home1/01839/dakotaz/GBM_genomics/MAKE_VCF/SNPs_GATK.sh C484.TCGA-06-0125-02A-11D-2280-08.1 C484.TCGA-06-0125-10A-01D-1490-08.6

#bash /home1/01839/dakotaz/GBM_genomics/MAKE_VCF/IN_SEQUENCE/SNPs_SomaticSniper.sh C484.TCGA-06-0125-02A-11D-2280-08.1 C484.TCGA-06-0125-10A-01D-1490-08.6