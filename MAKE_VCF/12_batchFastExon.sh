#!/bin/bash

module load bwa/0.7.4
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

# this is all the prelim stuff that I'm not following
# big question is, what would the directory look like? Is there nothing but BAMs here? Because that is not at all what my directory looks like...
# little question is, where the hell does $pfx.log come from, anyway?
cd /scratch/00001/mattcowp/dakota
echo $pfx
cd $pfx
echo "$pfx:" `date` > $pfx.log

# we have to convert a bam into a fastq; possibly two, I'm not sure; these will be commands to do that...
echo "Running SamToFastQ on $pfx...."  `date` >> ../$pfx.log
(java -d64 -Xmx4g -jar /opt/apps/picard/1.92/SamToFastq.jar I=$pfx.bam F=$pfx.1.fastq F2=$pfx.2.fastq 2>&1) >> $pfx\_SamToFastQ.log
echo "SamToFastQ done: " `date` >> ../$pfx.log

# the next step is to split and clean up the fastq files, which will be here...
rm split.*
echo "" >> ../$pfx.log
echo "python $HOME/bin/split_fastq_threaded.py $pfx $splitSize"  `date` > split.script
echo "Waiting while splitting and cleaning the fastq file."  `date` >> ../$pfx.log
python $HOME/bin/launcher.py split.script 2
echo "Fastq split finished :" `date` >> ../$pfx.log
rm -rf pylauncher_tmp*