#!/bin/bash

## command line variables (this is taken from the previous script, GATK_SNPs.sh)
pfx=$1
tumor=$2
blood=$3
SS=$4
MT=$5
ST=$6
VS=$7

## make paths to working directory, data directory, the data, and other important places and files
dataDir="/scratch/00001/mattcowp/dakota"
bloodDir="/scratch/00001/mattcowp/dakota/normal"
workDir="$dataDir/$pfx$tumor.WORK"
tumorbam="$workDir/$pfx$tumor.realn.recal.bam"
normalbam="$workDir/$pfx$blood.realn.recal.bam"

## make paths to files you will need, including hg19, 
ref_dir="/work/00001/mattcowp/Hs_reference_datasets"
hgReference="$ref_dir/Homo_sapiens.GRCh37.72.dna.fa"

#####################################
## First SNP-Caller: SomaticSniper ##
#####################################

## locate or build the SomaticSniper executable
function LoB_SS {
	echo "Load or build SomaticSniper: " `date` >> ./$pfx.SNPcalls.log
	somaticsniper="/scratch/00001/mattcowp/dakota/bam-somaticsniper"
}

## run SomaticSniper
function runSS {
	cd $workDir
	echo "run SomaticSniper: " `date` >> ./$pfx.SNPcalls.log
	$somaticsniper -q 40 -Q 40 -J -s 0.001 -F vcf -f $hgReference $tumorbam $normalbam $pfx$tumor.SS.vcf >> ./$pfx.SNPcalls.log
	echo "SomaticSniper finished: " `date` >> ./$pfx.SNPcalls.log
}

###############################
## Second SNP-Caller: MuTect ##
###############################

## locate or build the MuTect executable
function LoB_MT {
	exit
}

## run MuTect
function runMT {
	exit
}

###############################
## Third SNP-Caller: Strelka ##
###############################

## locate or build the Strelka executable
function LoB_ST {
	exit
}

## run Strelka
function runST {
	exit
}

################################
## Fourth SNP-Caller: VarScan ##
################################

## locate or build the VarScan executable
function LoB_VS {
	exit
}

## run VarScan
function runVS {
	exit
}

###################
## MAIN FUNCTION ##
###################

if [ "$4" == "yes" ]; then
	LoB_SS
	runSS
fi

exit;