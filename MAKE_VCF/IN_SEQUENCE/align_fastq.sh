#!/bin/bash

# the purpose of this script is to align the two fastq files of a TCGA sample (paired ends) to hg19

module load bwa
module load samtools

# pfx is the sample name
pfx=$1
refDir="/work/00001/mattcowp/Hs_reference_datasets"
hgReference="$refDir/Homo_sapiens.GRCh37.72.dna.fa"
queue="normal"

# move to working directory
cd /scratch/00001/mattcowp/dakota/$pfx

# bwa aln first fastq
bwa aln -q 30 -t 4 $refindex $f3fastqfile > $f3fastqfile_prefix.sai 2>$f3fastqfile_prefix.bwa.log
# bwa aln second fastq
bwa aln -q 30 -t 4 $refindex $f5fastqfile > $f5fastqfile_prefix.sai 2>$f5fastqfile_prefix.bwa.log
echo "bwa aln for $f5fastqfile done"

# bwa sampe | samtools view
(bwa sampe -a 600 -P -r "$RG" $refindex $f3fastqfile_prefix.sai $f5fastqfile_prefix.sai $f3fastqfile $f5fastqfile | samtools view -bSh -o $outprefix.bam -) 2>$outprefix.sampe.log
# samtools sort
samtools sort $outprefix.bam $outprefix.sorted >>$outprefix.sampe.log
# samtools index
samtools index $outprefix.sorted.bam >>$outprefix.sampe.log



* rename bams so aligned to hg19