#!/bin/bash

## command line variables (this is taken from the previous script, GATK_SNPs.sh)
tumor_pfx=$1
blood_pfx=$2

## working directory
tumor_dir="$tumor_pfx.WORK"
cd /scratch/00001/mattcowp/dakota/$tumor_dir

## files you may need...
ref_dir="/work/00001/mattcowp/Hs_reference_datasets"
hgReference="$ref_dir/Homo_sapiens.GRCh37.72.dna.fa"
tumor_bam="/scratch/00001/mattcowp/dakota/$tumor_dir/$tumor_pfx.realn.recal.bam"
normal_bam="/scratch/00001/mattcowp/dakota/$tumor_dir/$blood_pfx.realn.recal.bam"

## run SomaticSniper on a set of T-N pairs
echo "SomaticSniper: " `date` >> ./$tumor_pfx.ss.log
/scratch/00001/mattcowp/dakota/bam-somaticsniper -q 40 -Q 40 -J -s 0.001 -F vcf -f $hgReference $tumor_bam $normal_bam $tumor_pfx.SS.vcf &>> ./$tumor_pfx.ss.log
echo "SomaticSniper finished: " `date` >> ./$tumor_pfx.ss.log

exit;