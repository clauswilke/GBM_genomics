#!/bin/bash

# load modules (gatk current default is 2.7.2)
module load gatk
module load samtools
module load picard 

## this script does the GATK processing and variant calibration on all of the sorted .bam and .bam.bai files

dataDir="/scratch/00001/mattcowp/dakota"
cd $dataDir

## get sample names and file names, and make sure the two files have the same sample name
tumorPfx=$1
bloodPfx=$2
tumorSamp=${tumorPfx:0:17}
bloodSamp=${bloodPfx:0:17}
tumorDir="$tumorPfx.WORK"
# echo $tumor_dir
echo $tumorPfx $tumorSamp
echo $bloodPfx $bloodSamp
if [ $tumorSamp != $bloodSamp ]; then
	echo "Non-matching sample names:"
	echo "	$tumorSamp"
	echo "	$bloodSamp"
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

echo "Generating somatic calls for $tumorPfx" > $tumorPfx.gatk.log

# make a temporary directory in which this will actually run
mkdir $tumorDir
cd $tumorDir

# remove duplicate reads from bamfiles, leaving only the read with the highest map quality
samtools rmdup ../$tumorPfx/$tumorPfx.out.sorted.matefixed.bam $tumorPfx.dedup.bam &
samtools rmdup ../normal/$bloodPfx/$bloodPfx.out.sorted.matefixed.bam $bloodPfx.dedup.bam &
wait;

# reindex the de-duplicated bamfiles
samtools index $tumorPfx.dedup.bam &
samtools index $bloodPfx.dedup.bam &
wait;

## indel re-alignment
# generate a common interval file for the blood and tumor data
echo "Start RealignerTargetCreator: " >> ../$tumorPfx.gatk.log
java -d64 -jar $gatkJar -R $hgReference -nt 12 -T RealignerTargetCreator -rf BadCigar -known $G1000Mills -known $G1000Phase1Indels -o $tumorPfx.intervals -I $tumorPfx.dedup.bam -I $bloodPfx.dedup.bam &>>../$tumorPfx.gatk.log
echo "$0: `date`"

# (local realignment around indels) tumor
echo "IndelRealigner: " `date` >> ../$tumor_pfx.gatk.log
java -d64 -jar $gatkJar -R $hgReference -T IndelRealigner -rf BadCigar -I $tumorPfx.dedup.bam -known $G1000Mills -known $G1000Phase1Indels -targetIntervals $tumorPfx.intervals -o $tumorPfx.realn.bam &>>../$tumorPfx.gatk.log &
# (local realignment around indels) blood
java -d64 -jar $gatkJar -R $hgReference -T IndelRealigner -rf BadCigar -I $bloodPfx.dedup.bam -known $G1000Mills -known $G1000Phase1Indels -targetIntervals $tumorPfx.intervals -o $bloodPfx.realn.bam &>>../$tumorPfx.gatk.log &
wait;

# index the re-aligned bams
samtools index $tumorPfx.realn.bam
samtools index $bloodPfx.realn.bam

## base recalibration
# "BaseRecalibrator"
echo "BaseRecalibrator on $tumorPfx: " `date`  >> ../$tumorPfx.gatk.log
java -d64 -jar $gatkJar -nct 8 -T BaseRecalibrator -rf BadCigar -I $tumorPfx.realn.bam -R $hgReference -knownSites $dbSNP -o $tumorPfx.recal.grp &>>../$tumorPfx.gatk.log &
echo "BaseRecalibrator on $bloodPfx: " `date`  >> ../$tumorPfx.gatk.log
java -d64 -jar $gatkJar -nct 8 -T BaseRecalibrator -rf BadCigar -I $bloodPfx.realn.bam -R $hgReference -knownSites $dbSNP -o $bloodPfx.recal.grp &>>../$tumorPfx.gatk.log &
echo "Waiting for BQSR to finish... "`date` >> ../$tumorPfx.gatk.log
wait;
 
# "PrintReads"   
echo "PrintReads on $tumorPfx: " `date` >> ../$tumorPfx.gatk.log
java -d64 -jar $gatkJar -nct 8 -T PrintReads -rf BadCigar -R $hgReference -I $tumorPfx.realn.bam -BQSR $tumorPfx.recal.grp -o $tumorPfx.realn.recal.bam 2>>../$tumorPfx.gatk.log &
echo "PrintReads on $bloodPfx: " `date` >> ../$tumorPfx.gatk.log
java -d64 -jar $gatkJar -nct 8 -T PrintReads -rf BadCigar -R $hgReference -I $bloodPfx.realn.bam -BQSR $bloodPfx.recal.grp -o $bloodPfx.realn.recal.bam 2>>../$tumorPfx.gatk.log &
echo "Waiting for PrintReads to finish... "`date` >> ../$tumorPfx.gatk.log
wait;

# index the re-aligned and re-calibrated bams
samtools index $tumorPfx.realn.recal.bam &
samtools index $bloodPfx.realn.recal.bam &
echo "Indexing recaled and realigned bam files... "`date` >> ../$tumorPfx.gatk.log
wait;

exit;

