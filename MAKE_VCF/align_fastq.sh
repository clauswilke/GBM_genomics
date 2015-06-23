#!/bin/bash

# the purpose of this script is to align the two fastq files of a TCGA sample (paired ends) to hg19

module load bwa
module load samtools

# pfx is the sample name
pfx=$1
tissue=$2
outprefix="$pfx.out"
queue="normal"
# reference hg19
refDir="/work/00001/mattcowp/Hs_reference_datasets"
hgReference="$refDir/Homo_sapiens.GRCh37.72.dna.fa"
# what data we're running on...
data_dir="/scratch/00001/mattcowp/dakota"
if [ $tissue = "blood" ]
then
	data_dir="/scratch/00001/mattcowp/dakota/normal"
fi
fastq1="$data_dir/$pfx/$pfx.1.fastq"
fastq2="$data_dir/$pfx/$pfx.2.fastq"

# move to working directory
cd $data_dir/$pfx

# the bwa manual says this is important
# but Matt says take it out...
# bwa index -p RefSeqbwaidx -a bwtsw $hgReference

# bwa aln first fastq
bwa aln -q 30 -t 8 $hgReference $fastq1 > $fastq1.aln.sai 2>$fastq1.aln.bwa.log &
# bwa aln second fastq
bwa aln -q 30 -t 8 $hgReference $fastq2 > $fastq2.aln.sai 2>$fastq2.aln.bwa.log &
wait
echo "bwa aln for $pfx done"

# need to specify read groups (RG), which is the argument of  -r in bwa sempe
# first get the library information from the TCGA bam:
LB=`samtools view -H ./$pfx.bam | grep -m 1 '@RG' | awk '{print substr($5,4) }'`
if [ ${#LB} -eq 0 ] # this is the case if the bam has no library information at all
then
	LB="Catch-1"
fi
# then set the @RG header
RG="@RG\tID:1\tPL:ILLUMINA\tSM:$pfx\tLB:$LB\tDS:ref=hg19"

# bwa sampe | samtools view
(bwa sampe -a 600 -P -r "$RG" $hgReference $fastq1.aln.sai $fastq2.aln.sai $fastq1 $fastq2 | samtools view -bSh -o $outprefix.bam -) 2>$outprefix.log
# samtools sort
samtools sort -@ 16 $outprefix.bam $outprefix.sorted 2>>$outprefix.log
# samtools index
samtools index $outprefix.sorted.bam 2>>$outprefix.log
