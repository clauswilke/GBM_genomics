#!/bin/bash

for r in $@
do
    echo '#!/bin/bash' > $r.sbat
    echo '' >> $r.sbat
    echo '#SBATCH -J' $r >> $r.sbat
    echo '#SBATCH -o' $r.out >> $r.sbat
    echo '#SBATCH -e' $r.err >> $r.sbat
    echo '#SBATCH -n 16' >> $r.sbat
    echo '#SBATCH -p normal' >> $r.sbat
    echo '#SBATCH -A mattcowp' >> $r.sbat
    echo '#SBATCH -t 24:00:00' >> $r.sbat

    echo '' >> $r.sbat
    echo 'cd $SCRATCH/TCGA_NextGen_Analysis/hg19_bams' >> $r.sbat
    echo 'echo $PWD' >> $r.sbat
    echo '' >> $r.sbat

    echo 'echo 'Run:' `date`' >> $r.sbat
    echo "bash $r" >> $r.sbat
    echo 'echo 'Done: ' `date`' >> $r.sbat

done;