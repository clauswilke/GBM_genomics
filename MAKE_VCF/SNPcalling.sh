#!/bin/bash

## command line variables (this is taken from the previous script, GATK_SNPs.sh)
pfx=$1
tumor=$2
blood=$3

## make paths to working directory, data directory, the data, and other important places and files
dataDir="/scratch/00001/mattcowp/dakota"
bloodDir="/scratch/00001/mattcowp/dakota/normal"
workDir="$dataDir/$pfx$tumor.WORK"
tumorbam="$workDir/$pfx$tumor.realn.recal.bam"
normalbam="$workDir/$pfx$blood.realn.recal.bam"

## make paths to files you will need, including hg19, 
refDir="/work/00001/mattcowp/Hs_reference_datasets"
hgReference="$refDir/Homo_sapiens.GRCh37.72.dna.fa"

## run SomaticSniper
cd $workDir
echo "run SomaticSniper: " `date` >> ./$pfx.SNPcalls.log
$somaticsniper -q 40 -Q 40 -J -s 0.001 -F vcf -f $hgReference $tumorbam $normalbam $pfx$tumor.SS.vcf >> ./$pfx.SNPcalls.log
echo "SomaticSniper finished: " `date` >> ./$pfx.SNPcalls.log
