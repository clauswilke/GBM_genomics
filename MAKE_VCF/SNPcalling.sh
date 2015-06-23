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
	muTectDir="/work/01839/dakotaz/local/bin/MuTect"
	muTectJar="$muTectDir/muTect-1.1.1.jar"
	myRefDir="/work/01839/dakotaz/referenceDS"
	COSMIC="$myRefDir/Cosmic.hg19.vcf"
	dbSNP="$refDir/dbsnp_137.b37.vcf"
	
	# unzip dbSNP in scratch (because it is huge)
	# gunzip -c $dbSNPgz > $workDir/dbSNP.human_9606_b142_GRCh37p13.All.vcf
	# dbSNP="$workDir/dbSNP.human_9606_b142_GRCh37p13.All.vcf"
	
	# and finally the run command...
	java -Xmx12g -jar $muTectJar --analysis_type MuTect --reference_sequence $hgReference --cosmic $COSMIC --dbsnp $dbSNP --input_file:normal $normalbam --input_file:tumor $tumorbam --out $pfx$tumor.MT.out --coverage_file $pfx$tumor.MTcoverage.wig.txt
	
	# clean up 
	# rm $dbSNP
	
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
## note that this is the VarScan protocol to call Somatic Variants in CANCER; that is not the same thing as lots of other potential VarScan protocols. See their website for all their protocols.

function runVS {

	## set up environment
	ml samtools

	## important variables
	VarScanJar="/home1/01839/dakotaz/June2015/SNPcalling0169/VarScan/VarScan.v2.3.8.jar"
	outPfx="$pfx.VarScan"
	
	## samtools the important input files
	puNormal="$workDir/$pfx$blood.pileup"
	puTumor="$workDir/$pfx$tumor.pileup"
	samtools mpileup -q 1 -f $hgReference $normalbam > $puNormal
	samtools mpileup -q 1 -f $hgReference $tumorbam > $puTumor

	## and run VarScan
	java -jar $VarScanJar somatic $puNormal $puTumor $outPfx
}

###################
## MAIN FUNCTION ##
###################

if [ "$SS" == "yes" ]; then
	runSS
fi

if [ "$MT" == "yes" ]; then
	runMT
fi

if [ "$ST" == "yes" ]; then
	runST
fi

if [ "$VS" == "yes" ]; then
	runVS
fi

exit;