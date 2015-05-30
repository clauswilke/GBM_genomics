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
refDir="/work/00001/mattcowp/Hs_reference_datasets"
hgReference="$refDir/Homo_sapiens.GRCh37.72.dna.fa"

## additional references for Strelka

## additional references for VarScan

#####################################
## First SNP-Caller: SomaticSniper ##
#####################################

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

## run MuTect
function runMT {
	
	# go to the working directory
	cd $workDir
	
	# set variables
	muTectDir="/work/01839/dakotaz/MuTect"
	muTectJar="$muTectDir/muTect-1.1.1.jar"
	myRefDir="/work/01839/dakotaz/referenceDS"
	dbSNPgz="$myRefDir/dbsnp_132_b37.leftAligned.vcf.gz"
	COSMIC="$myRefDir/Cosmic.hg19.vcf"
	
	# unzip dbSNP in scratch (because it is huge)
	gunzip -c $dbSNPgz > $workDir/dbsnp_132_b37.leftAligned.vcf
	dbSNP="$workDir/dbsnp_132_b37.leftAligned.vcf"
	
	# and finally the run command...
	java -Xmx12g -jar $muTectJar --analysis_type MuTect --reference_sequence $hgReference --cosmic $COSMIC --dbsnp $dbSNP --input_file:normal $normalbam --input_file:tumor $tumorbam --out $pfx$tumor.MT.out --coverage_file $pfx$tumor.MTcoverage.wig.txt
	
	# clean up 
	rm $dbSNP
	
}

###############################
## Third SNP-Caller: Strelka ##
###############################

## run Strelka
function runST {
	exit
}

################################
## Fourth SNP-Caller: VarScan ##
################################

## run VarScan

function runVS {
	exit
}

###################
## MAIN FUNCTION ##
###################

if [ "$4" == "yes" ]; then
	runSS
fi

if [ "$5" == "yes" ]; then
	runMT
fi

if [ "$6" == "yes" ]; then
	runST
fi

if [ "$7" == "yes" ]; then
	runVS
fi

exit;