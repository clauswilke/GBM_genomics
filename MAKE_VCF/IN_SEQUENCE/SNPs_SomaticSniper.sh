#!/bin/bash

## load modules

## command line variables (this is taken from the previous script, GATK_SNPs.sh)
tumor_pfx=$1
blood_pfx=$2
tumor_sample=${tumor_pfx:0:17}
blood_sample=${blood_pfx:0:17}
echo $tumor_pfx $tumor_sample
echo $blood_pfx $blood_sample
if [ $tumor_sample != $blood_sample ]; then
	echo "Non-matching sample names:"
	echo "	$tumor_sample"
	echo "	$blood_sample"
	exit 1;
fi

## files you may need...
$output_dir=""
ref_dir="/work/00001/mattcowp/Hs_reference_datasets"
hgReference="$ref_dir/Homo_sapiens.GRCh37.72.dna.fa"

## run SomaticSniper on a set of T-N pairs
echo "SomaticSniper: " 'date' >> ../$tumor_sample.log
/scratch/00001/mattcowp/dakota/bam-somaticsniper -q 40 -Q 40 -J -s 0.001 -F vcf -f $hgReference $tumor_pfx.realn.recal.bam $blood_pfx.realn.recal.bam $tumor_pfx.SS.vcf

## make sure the file exists, was made, then exit
if [ -f $tumor_pfx.SS.vcf ]; then
	cd ..
	mv tumor_sample/$tumor_pfx.SS.vcf $output_dir
	mv tumor_sample.log $output_dir/$tumor_pfx.log
	echo "SomaticSniper run: " 'date' >> $output_dir/$tumor_pfx.log
	exit 0;
fi;

exit1;