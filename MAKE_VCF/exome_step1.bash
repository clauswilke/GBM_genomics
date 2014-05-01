#!/bin/bash

module load bwa
module load samtools

date

pfx=$1
grpNum=$2
f3fastqfile="r1.$grpNum"
f5fastqfile="r2.$grpNum"
refindex=$3
outprefix="mapped.$grpNum"
workingdir=$4

f3fastqfile_prefix=$outprefix."r1"
f5fastqfile_prefix=$outprefix."r2"

echo "exome_step1.bash $f3fastqfile $f5fastqfile $refindex $outprefix $workingdir"
cd $workingdir
echo "bwa aln $refindex $f3fastqfile > $f3fastqfile_prefix.sai 2>$f3fastqfile_prefix.bwa.log "
bwa aln -q 30 -t 4 $refindex $f3fastqfile > $f3fastqfile_prefix.sai 2>$f3fastqfile_prefix.bwa.log
echo "bwa aln for $f3fastqfile done"
date
bwa aln -q 30 -t 4 $refindex $f5fastqfile > $f5fastqfile_prefix.sai 2>$f5fastqfile_prefix.bwa.log
echo "bwa aln for $f5fastqfile done"
date


echo "### BWA Generate Realigned bam file ###" >>$outprefix.sampe.log
# get the library information from the TCGA bam file
LB=`samtools view -H ../$pfx.bam | grep -m 1 '@RG' | awk '{ print substr($5,4) }'`
if [ ${#LB} -eq 0 ]  #some bams have no library information
then
	LB="Catch-1"
fi
RG="@RG\tID:1\tPL:ILLUMINA\tSM:$pfx\tLB:$LB\tDS:ref=hg19,pfx="$refindex
(bwa sampe -a 600 -P -r "$RG" $refindex $f3fastqfile_prefix.sai $f5fastqfile_prefix.sai $f3fastqfile $f5fastqfile | samtools view -bSh -o $outprefix.bam -) 2>$outprefix.sampe.log

echo "### Samtools sort bam file " >>$outprefix.sampe.log
samtools sort $outprefix.bam $outprefix.sorted >>$outprefix.sampe.log
echo "### Samtools index sorted bam file " >>$outprefix.sampe.log
samtools index $outprefix.sorted.bam >>$outprefix.sampe.log

# split bam by chromosome
for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y;
do
	echo "samtools view -bh $outprefix.sorted.bam "chr$c" > chr$c.$outprefix.sorted.bam">>$outprefix.sampe.log
	samtools view -bh $outprefix.sorted.bam -o chr$c.$outprefix.sorted.bam "chr$c"
	touch chr$c.$outprefix.sorted.bam # sometimes the files have initial size 0
	# as a last check remove any files that are size 0.  these arise when
    # no reads map to a particular chromosome this is important to keep the merge from crashing
	if [ ! -s chr$c.$outprefix.sorted.bam ]
	then
	    rm chr$c.$outprefix.sorted.bam;
	fi
done;