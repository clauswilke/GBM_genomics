#!/bin/bash

module load bwa/0.7.7
module load samtools
module load picard
module load python
module load pylauncher
module load gatk/2.5.2

# these are command line arguments; pfx is the sample name (the TCGA barcode; everything before .bam or .fastq)
pfx=$1
splitSize=$2
batchSize=$3
coresPerRun=$4
refDir="/work/00001/mattcowp/Hs_reference_datasets"
hgReference="$refDir/Homo_sapiens.GRCh37.72.dna.fa"
queue="normal"

# 
cd /scratch/00001/mattcowp/dakota
echo $pfx
cd $pfx
echo "$pfx:" `date` > ../$pfx.log

# we have to convert a bam into a fastq
# might have to change path to picard; check this in idev
echo "Running SamToFastQ on $pfx...."  `date` >> ../$pfx.log
(java -d64 -Xmx4g -jar /opt/apps/picard/1.98/SamToFastq.jar I=$pfx.bam F=$pfx.1.fastq F2=$pfx.2.fastq 2>&1) >> $pfx\_SamToFastQ.log
echo "SamToFastQ done: " `date` >> ../$pfx.log

exit

# the next step is to split and clean up the fastq files
# which is done by the script split_fastq_threaded.py
rm split.*
echo "" >> ../$pfx.log
echo "python $HOME/GBM_genomics/MAKE_VCF/split_fastq_threaded.py $pfx $splitSize" > split.script
echo "Waiting while splitting and cleaning the fastq file."  `date` >> ../$pfx.log
python $HOME/GBM_genoics/MAKE_VCF/launcher.py split.script 2
echo "Fastq split finished :" `date` >> ../$pfx.log
rm -rf pylauncher_tmp*

# not sure why we're making this note in the log here
echo "-----"
echo "$pfx:" `date` > $pfx.log

# I don't think we need this, but I will check
# we need to make sure we are in the right working directory
# which will depend if we ran unpack or not.
if [[ "$PWD" != *"$pfx"* ]];
then
    echo "Changing directory to /scratch/00001/mattcowp/dakota/$pfx"  `date` >> $pfx.log
    cd $SCRATCH/TCGA_NextGen_Analysis/$pfx
fi;

