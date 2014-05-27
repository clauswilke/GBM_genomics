#!/bin/bash

module load bwa/0.7.4
module load samtools
module load picard
module load python
module load pylauncher
module load gatk/2.5.2

pfx=$1
splitSize=$2
batchSize=$3
coresPerRun=$4
refDir="$WORK/Hs_reference_datasets"
hgReference="$refDir/Homo_sapiens.GRCh37.72.dna.fa"
queue="normal"
pkFile="$HOME/src/pk"

function decrypt {
      echo "gpg --no-tty --batch --yes --passphrase-file $3 --output $1 --decrypt $2/$1.gpg"  `date` 
      gpg --no-tty --batch --passphrase-file $3 --output $1 --decrypt $2/$1.gpg
      sleep 10
}

function encrypt {
    md5sum $1 > $1.md5
    gpg --no-tty --batch --cipher-algo AES256 --passphrase-file $2 --output $1.gpg --symmetric $1
    if [ -f $1.gpg ]
    then
        echo $1 " done."
        # rm $1
    else
        echo "Encryption of $1 failed!"  `date` 
    fi
}

# 0. Decrypt the bam file and extract, split and filter the reads (remove unpaired reads)
# Note:   I added a optional initial flag, "UNPACK_ONLY", that will run the decryption and SamToFastQ steps.
#         There are very low memory steps, but are also very slow.  Therefore, it's much more efficient
#         to unpack everything with 8 samples per node and then process everything else at 4 cores per
#         node.  Unpack is roughly 3/4 of the total processing time.
echo $1;
if [[ "$1" == *"UNPACK"* ]];
then

    # shift to right command-lines args                                                                                                                                
    pfx=$2
    splitSize=$3
    batchSize=$4
    coresPerRun=$5
    echo "$pfx:" `date` > $pfx.log
    echo "UNPACK:  Sample: " $pfx  `date` >> $pfx.log

    # make the dir and                                                                                                                                                  
    mkdir $pfx
    cd $pfx
    # delete old bams
    rm $pfx.bam
    echo "Waiting for $pfx to decrypt from ../tcga_bams"  `date` >> ../$pfx.log
    decrypt $pfx.bam ../tcga_bams $pkFile
    echo "Running SamToFastQ on $pfx...."  `date` >> ../$pfx.log
    (java -d64 -Xmx4g -jar /opt/apps/picard/1.92/SamToFastq.jar I=$pfx.bam F=$pfx.1.fastq F2=$pfx.2.fastq 2>&1) >> $pfx\_SamToFastQ.log
    echo "SamToFastQ done: " `date` >> ../$pfx.log

    # split and clean up the fastq files                                                                                                                                
    rm split.*
    echo "" >> ../$pfx.log
    echo "python $HOME/bin/split_fastq_threaded.py $pfx $splitSize"  `date` > split.script
    echo "Waiting while splitting and cleaning the fastq file."  `date` >> ../$pfx.log
    python $HOME/bin/launcher.py split.script 2
    echo "Fastq split finished :" `date` >> ../$pfx.log
    rm -rf pylauncher_tmp*

    if [ $1 == "UNPACK_ONLY" ];
    then
	   exit;
    fi;

fi;

echo "-----"
echo "$pfx:" `date` > $pfx.log

# we need to make sure we are in the right working directory
# which will depend if we ran unpack or not.
if [[ "$PWD" != *"$pfx"* ]];
then
    echo "Changing directory to $SCRATCH/TCGA_NextGen_Analysis/$pfx"  `date` >> $pfx.log
    cd $SCRATCH/TCGA_NextGen_Analysis/$pfx
fi;

#  check to ensure that we go an equal number of splits, if there
# are no r1 or r2 files, then try to re-run split script
numR1Files=`ls r1.* | wc -l`
echo "Number of r1 files: " $numR1Files  `date` >> ../$pfx.log
numR2Files=`ls r2.* | wc -l`
echo "Number of r2 files: " $numR2Files  `date` >> ../$pfx.log
if [ $numR1Files -eq 0 -o $numR2Files -eq 0 ]
then
    echo "python $HOME/bin/split_fastq_threaded.py $pfx $splitSize"  `date` > split.script
    echo "Initial split seems suspect:  r1=$numR1Files, r2=$numR2Files"  `date` >> $pfx.log
    echo "Waiting while splitting and cleaning the fastq file."  `date` >> ../$pfx.log
    python $HOME/bin/launcher.py split.script 2
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

# 2. Move a sets of splits into sub-directories to be processed independently
echo "Dividing splits into $batchSize per batch: " `date` >> ../$pfx.log
batchNum=0
subdir=""
for i in `seq 0 $((numR1Files-1))`; do
    if [ `expr $i % $batchSize` -eq 0 ]
    then
        subdir="b.$batchNum"
        mkdir $subdir
        batchNum=$((batchNum + 1))
    fi
    mv r2.$i $subdir
    mv r1.$i $subdir
done

# And the residual set, if any:
for file in $( ls r1.* ) ; do
    subdir="b.$batchNum"
    mkdir $subdir
    mv r1.* b.$i
    mv r2.* b.$i
done

# 3. Launch exome_step1.sh on each split within it's own directory; store job numbers; launch exome_step2.sh to combine chr files
subdirList=$( ls -d b.* )
rm map.script
for subdir in $subdirList ; do
  echo "Creating launcher for all files in $subdir: " `date` >> ../$pfx.log
  cd $subdir
  for file in $( ls r1.* ); do
     fileExt="${file##*.}"
     echo "Run exome_step1.sh on r1.$fileExt and r2.$fileExt" >> ../$pfx.log
     echo "$HOME/bin/exome_step1.bash $pfx $fileExt $hgReference $subdir >& $subdir/mapped.$fileExt.log" >> ../map.script
  done
  cd ..
done

echo "Waiting for set of mapping runs to finish: " `date` >> ../$pfx.log
python $HOME/bin/launcher.py map.script $coresPerRun
rm -rf pylauncher_tmp*
echo "Mapping Finished: " `date` >> ../$pfx.log

# 4. Launch job to combine final chr files across all directories
# split bam by chromosome
rm -f merge.script
chrList=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
for c in ${chrList[@]};
do
   echo "Merging chr$c "
   echo "samtools merge -f chr$c.bam b.*/chr$c.mapped.*.sorted.bam; samtools sort chr$c.bam chr$c.sorted; samtools index chr$c.sorted.bam" >> ./merge.script
done
echo "Waiting for per-chromosome merge to finish: " `date` >> ../$pfx.log
python $HOME/bin/launcher.py merge.script $coresPerRun
rm -rf pylauncher_tmp*
echo "Chromosome merge Finished: " `date` >> ../$pfx.log

# 5. Launch GATK on each chromosome, sorted bam file to realign around indels and do BQSR
# rm -f variants.script
# for c in ${chrList[@]};
# do
#     echo "GATK via exome_step2.bash on chr$c.sorted.bam"
#     echo "$HOME/bin/exome_step2.bash chr$c.sorted chr$c >& variants.chr$c.log" >> variants.script
# done
# echo "Waiting for variant calling to finish"
# python $HOME/bin/launcher.py variants.script $coresPerRun
# rm -rf pylauncher_tmp*
# echo "Variant Calling Finished: `date` "

# 6. Merge all bam files into single sample bam
echo "Waiting for final whole-sample merger to finish: " `date` >> ../$pfx.log
samtools merge -f $pfx.sort.bam *.sorted.bam
echo "Whole-sample merge Finished: " `date` >> ../$pfx.log
# samtools index $pfx.sort.bam

# 7.  Reduce the bam file size with GATK reduce reads
# echo "Waiting for GATK ReduceReads to finish"
# java -d64 -Xmx16g -jar $TACC_GATK_DIR/GenomeAnalysisTK.jar -R $hgReference -T ReduceReads -I $pfx.sort.bam -O $pfx.sort.reduce.bam
# echo "ReduceReads Finished: `date`.  Indexing... "

# 7.  Clean up and encrypt everything
# encrypt $pfx.sort.bam $pkFile
# encrypt $pfx.sort.reduce.bam $pkFile
# rm $pfx.bam
rm chr*
# rm -rf b.*
rm *.bai  # clean up any bam indices
# rm *.fastq
rm *.log
rm *.script

echo "Processing completed at: " `date` >> ../$pfx.log
exit;





