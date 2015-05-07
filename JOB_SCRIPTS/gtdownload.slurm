#! /bin/bash

#SBATCH -J TCGAman
#SBATCH -o TCGAman.out
#SBATCH -e TCGAman.err
#SBATCH -n 16
#SBATCH -p normal
#SBATCH -A mattcowp
#SBATCH -t 5:00:00

set -x                      # Echo commands

echo `hostname`

# run the script...
cd /work/01839/dakotaz/
bash /work/01839/dakotaz/gtdownload.sh