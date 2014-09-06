# start this script in $SCRATCH/DOWNSAMPLING_DATA
# this takes the downsampling data and makes a .bam, NOT A VCF

module load bwa/0.7.4
module load samtools
module load python
module load pylauncher

# pfx is the sample name (everything before .bam or .fastq), in this case GSAF barcodes
# split_size is 2,000,000
pfx=$1
splitSize=$2
batchSize=$3
coresPerRun=$4
refDir="/work/00001/mattcowp/Hs_reference_datasets"
hgReference="$refDir/Homo_sapiens.GRCh37.72.dna.fa"
queue="normal"

# make a directory into which the .bam file will go, and move to it
mkdir -p $pfx
cd $pfx
# split, and clean up the fastq files
rm split.*

echo "" >> ../$pfx.log
echo "python $SCRATCH/DOWNSAMPLING/split_fastq_threaded.py $pfx $splitSize"  `date` > split.script
echo "Waiting while splitting and cleaning the fastq file."  `date` >> ../$pfx.log
python $HOME/DOWNSAMPLING/launcher.py split.script 2
echo "Fastq split finished :" `date` >> ../$pfx.log
rm -rf pylauncher_tmp*

echo "---fastq split and cleaned---"
echo "$pfx:" `date` > $pfx.log

#  check to ensure that we go an equal number of splits, if there                                              
# are no r1 or r2 files, then try to re-run split script                                                       
numR1Files=`ls r1.* | wc -l`
echo "Number of r1 files: " $numR1Files  `date` >> ../$pfx.log
numR2Files=`ls r2.* | wc -l`
echo "Number of r2 files: " $numR2Files  `date` >> ../$pfx.log
if [ $numR1Files -eq 0 -o $numR2Files -eq 0 ]
then
    echo "python $HOME/DOWNSAMPLING/split_fastq_threaded.py $pfx $splitSize"  `date` > split.script
    echo "Initial split seems suspect:  r1=$numR1Files, r2=$numR2Files"  `date` >> $pfx.log
    echo "Waiting while splitting and cleaning the fastq file."  `date` >> ../$pfx.log
    python $HOME/DOWNSAMPLING/launcher.py split.script 2
    echo "Fastq split finished :" `date` >> ../$pfx.log
    rm -rf pylauncher_tmp*
    numR1Files=`ls r1.* | wc -l`
    echo "Number of r1 files: " $numR1Files  `date` >> ../$pfx.log
    numR2Files=`ls r2.* | wc -l`
    echo "Number of r2 files: " $numR2Files  `date` >> ../$pfx.log
fi

# check for equal numbers of files                                                                             
if [ $numR1Files -ne $numR2Files ]
then
    echo "Something wrong with split or fastq files.  Unqueal number of splits."  `date` >> ../$pfx.log
    exit 1;
fi

