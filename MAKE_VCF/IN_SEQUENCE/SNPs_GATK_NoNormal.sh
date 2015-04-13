#!/bin/bash

# load modules (gatk current default is 2.7.2)
module load jdk64
module load gatk
module load samtools/0.1.19
module load picard 

## this script does the GATK processing and variant calibration on all of the sorted .bam and .bam.bai files

dataDir="/scratch/00001/mattcowp/dakota"
cd $dataDir

## get sample names and file names, and make sure the two files have the same sample name
tumorPfx1=$1
tumorPfx2=$2
sampPfx=$3
tumorSamp1=${tumorPfx1:0:17}
tumorSamp2=${tumorPfx2:0:17}
tumorDir="$sampPfx.NoNormal.WORK"
# echo $tumor_dir
echo $tumorPfx1 $tumorSamp1
echo $tumorPfx2 $tumorSamp2
if [ $tumorSamp1 != $tumorSamp2 ]; then
	echo "Non-matching sample names:"
	echo "	$tumorSamp1"
	echo "	$tumorSamp2"
	exit 1;
fi

## locations of relevant and needed files
gatkJar="$TACC_GATK_DIR/GenomeAnalysisTK.jar"
refDir="/work/00001/mattcowp/Hs_reference_datasets"
hgReference="$refDir/Homo_sapiens.GRCh37.72.dna.fa"
dbSNP="$refDir/dbsnp_137.b37.vcf"
hapMap="$refDir/hapmap_3.3.b37.vcf"
G1000HiConfSNPs="$refDir/1000G_phase1.snps.high_confidence.b37.vcf"
G1000Phase1Indels="$refDir/1000G_phase1.indels.b37.vcf"
G1000OmniSNPs="$refDir/1000G_omni2.5.b37.vcf"
G1000Mills="$refDir/Mills_and_1000G_gold_standard.indels.b37.vcf"

echo "Generating somatic calls for $sampPfx" > $sampPfx.gatk.log

# make a temporary directory in which this will actually run
mkdir $tumorDir
cd $tumorDir

# remove duplicate reads from bamfiles, leaving only the read with the highest map quality
samtools rmdup ../$tumorPfx1/$tumorPfx1.out.sorted.matefixed.bam $tumorPfx1.dedup.bam &
samtools rmdup ../$tumorPfx2/$tumorPfx2.out.sorted.matefixed.bam $tumorPfx2.dedup.bam &
wait;

# reindex the de-duplicated bamfiles
samtools index $tumorPfx1.dedup.bam &
samtools index $tumorPfx2.dedup.bam &
wait;

## indel re-alignment
# generate a common interval file for the blood and tumor data
echo "Start RealignerTargetCreator: " >> ../$sampPfx.gatk.log
java -d64 -Xmx20g -jar $gatkJar -R $hgReference -nt 12 -T RealignerTargetCreator -rf BadCigar -known $G1000Mills -known $G1000Phase1Indels -o $sampPfx.intervals -I $tumorPfx1.dedup.bam -I $tumorPfx2.dedup.bam &>>../$sampPfx.gatk.log
echo "$0: `date`"

# (local realignment around indels) tumor
echo "IndelRealigner: " `date` >> ../$sampPfx.gatk.log
java -d64 -Xmx10g -jar $gatkJar -R $hgReference -T IndelRealigner -rf BadCigar -I $tumorPfx1.dedup.bam -known $G1000Mills -known $G1000Phase1Indels -targetIntervals $sampPfx.intervals -o $tumorPfx1.realn.bam &>>../$sampPfx.gatk.log &
# (local realignment around indels) blood
java -d64 -Xmx10g -jar $gatkJar -R $hgReference -T IndelRealigner -rf BadCigar -I $tumorPfx2.dedup.bam -known $G1000Mills -known $G1000Phase1Indels -targetIntervals $sampPfx.intervals -o $tumorPfx2.realn.bam &>>../$sampPfx.gatk.log &
wait;

# index the re-aligned bams
samtools index $tumorPfx1.realn.bam
samtools index $tumorPfx2.realn.bam

## base recalibration
# "BaseRecalibrator"
echo "BaseRecalibrator on $tumorPfx1: " `date`  >> ../$sampPfx.gatk.log
java -d64 -Xmx10g -jar $gatkJar -nct 8 -T BaseRecalibrator -rf BadCigar -I $tumorPfx1.realn.bam -R $hgReference -knownSites $dbSNP -o $tumorPfx1.recal.grp &>>../$sampPfx.gatk.log &
echo "BaseRecalibrator on $tumorPfx2: " `date`  >> ../$tumorPfx.gatk.log
java -d64 -Xmx10g -jar $gatkJar -nct 8 -T BaseRecalibrator -rf BadCigar -I $tumorPfx2.realn.bam -R $hgReference -knownSites $dbSNP -o $tumorPfx2.recal.grp &>>../$sampPfx.gatk.log &
echo "Waiting for BQSR to finish... "`date` >> ../$sampPfx.gatk.log
wait;
 
# "PrintReads"   
echo "PrintReads on $tumorPfx1: " `date` >> ../$sampPfx.gatk.log
java -d64 -Xmx10g -jar $gatkJar -nct 8 -T PrintReads -rf BadCigar -R $hgReference -I $tumorPfx1.realn.bam -BQSR $tumorPfx1.recal.grp -o $tumorPfx1.realn.recal.bam 2>>../$sampPfx.gatk.log &
echo "PrintReads on $tumorPfx2: " `date` >> ../$tumorPfx.gatk.log
java -d64 -Xmx10g -jar $gatkJar -nct 8 -T PrintReads -rf BadCigar -R $hgReference -I $tumorPfx2.realn.bam -BQSR $tumorPfx2.recal.grp -o $tumorPfx2.realn.recal.bam 2>>../$sampPfx.gatk.log &
echo "Waiting for PrintReads to finish... "`date` >> ../$sampPfx.gatk.log
wait;

# index the re-aligned and re-calibrated bams
samtools index $tumorPfx1.realn.recal.bam &
samtools index $tumorPfx2.realn.recal.bam &
echo "Indexing recaled and realigned bam files... "`date` >> ../$sampPfx.gatk.log
wait;

exit;

