#!/bin/bash

# load modules (gatk current default is 2.7.2)
# module load jdk64 # this line added just for 190 run
module load gatk
module load samtools/0.1.19 # used older version of samtools only on the 190 run
module load picard 

## this script does the GATK processing and variant calibration on all of the out.sorted.matefixed.bam (and corresponding .bam.bai) files

dataDir="/scratch/00001/mattcowp/dakota"
cd $dataDir

echo "Begin $0: `date`">>$tumorPfx.gatk.log
echo "-----------------------------------------------------------------------------">>$tumorPfx.gatk.log

## get sample names and file names, and make sure the two files have the same sample name
tumorPfx=$1
bloodPfx=$2
tumorSamp=${tumorPfx:0:17}
bloodSamp=${bloodPfx:0:17}
tumorDir="$tumorPfx.WORK"
# echo $tumor_dir
echo $tumorPfx $tumorSamp &>>$tumorPfx.gatk.log
echo $bloodPfx $bloodSamp &>>$tumorPfx.gatk.log
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

echo "Generating somatic calls for $tumorPfx" &>>$tumorPfx.gatk.log

# make a temporary directory in which this will actually run
mkdir $tumorDir
echo "Moving to working directory $tumorDir" &>>$tumorPfx.gatk.log
cd $tumorDir

# remove duplicate reads from bamfiles, leaving only the read with the highest map quality
# echo "Now running samtools rmdup..." >> $tumorPfx.gatk.log
# samtools rmdup ../$tumorPfx/$tumorPfx.out.sorted.matefixed.bam $tumorPfx.dedup.bam &>>../$tumorPfx.gatk.log &
# samtools rmdup ../normal/$bloodPfx/$bloodPfx.out.sorted.matefixed.bam $bloodPfx.dedup.bam &>>../$tumorPfx.gatk.log &
# wait;
# echo "`date`">>../$tumorPfx.gatk.log
# echo "">>../$tumorPfx.gatk.log

# how about mark duplicates with picard, instead of removing with samtools?
echo "-----------------------------------------------------------------------------" &>>../$tumorPfx.gatk.log
echo "Marking duplicates with picard..." &>>../$tumorPfx.gatk.log
echo "" &>>../$tumorPfx.gatk.log
java -d64 -Xmx22g -jar $TACC_PICARD_DIR/MarkDuplicates.jar MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=10000 INPUT="../$tumorPfx/$tumorPfx.out.sorted.matefixed.bam" OUTPUT="$tumorPfx.dupmark.bam" METRICS_FILE="$tumorPfx.dupmark.metrics" REMOVE_DUPLICATES="FALSE" &>>../$tumorPfx.gatk.log
java -d64 -Xmx22g -jar $TACC_PICARD_DIR/MarkDuplicates.jar MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=10000 INPUT="../normal/$bloodPfx/$bloodPfx.out.sorted.matefixed.bam" OUTPUT="$bloodPfx.dupmark.bam" METRICS_FILE="$bloodPfx.dupmark.metrics" REMOVE_DUPLICATES="FALSE" &>>../$tumorPfx.gatk.log

echo ""&>>../$tumorPfx.gatk.log
echo "-----------------------------------------------------------------------------" &>>../$tumorPfx.gatk.log
echo "Re-indexing the dupmarked bams..." &>>../$tumorPfx.gatk.log
echo ""&>>../$tumorPfx.gatk.log

# reindex the de-duplicated bamfiles
# this step is not necessary when we don't do samtools rmdup
# echo "Now running samtools index..." >>../$tumorPfx.gatk.log
# I don't know if we need to re-index dup marked stuff, but just in case I'm doing it...
samtools index $tumorPfx.dupmark.bam &>>../$tumorPfx.gatk.log 
samtools index $bloodPfx.dupmark.bam &>>../$tumorPfx.gatk.log 
echo "`date`">>../$tumorPfx.gatk.log
echo "">>../$tumorPfx.gatk.log

echo ""&>>../$tumorPfx.gatk.log
echo "-----------------------------------------------------------------------------" &>>../$tumorPfx.gatk.log
echo "Beginning Indel Realignment..." &>>../$tumorPfx.gatk.log
echo ""&>>../$tumorPfx.gatk.log

## indel re-alignment
# generate a common interval file for the blood and tumor data
echo "Start RealignerTargetCreator: " &>> ../$tumorPfx.gatk.log
java -d64 -jar $gatkJar -R $hgReference -nt 12 -T RealignerTargetCreator -rf BadCigar -known $G1000Mills -known $G1000Phase1Indels -o $tumorPfx.intervals -I $tumorPfx.dupmark.bam -I $bloodPfx.dupmark.bam &>>../$tumorPfx.gatk.log
echo "$0: `date`">>../$tumorPfx.gatk.log
echo "">>../$tumorPfx.gatk.log

# (local realignment around indels) tumor
echo "IndelRealigner: tumor " `date` >> ../$tumorPfx.gatk.log
java -d64 -jar $gatkJar -R $hgReference -T IndelRealigner -rf BadCigar -I $tumorPfx.dupmark.bam -known $G1000Mills -known $G1000Phase1Indels -targetIntervals $tumorPfx.intervals -o $tumorPfx.realn.bam &>>../$tumorPfx.gatk.log
echo "">>../$tumorPfx.gatk.log
# (local realignment around indels) blood
echo "IndelRealigner: blood " `date` >> ../$tumorPfx.gatk.log
java -d64 -jar $gatkJar -R $hgReference -T IndelRealigner -rf BadCigar -I $bloodPfx.dupmark.bam -known $G1000Mills -known $G1000Phase1Indels -targetIntervals $tumorPfx.intervals -o $bloodPfx.realn.bam &>>../$tumorPfx.gatk.log
echo "">>../$tumorPfx.gatk.log

echo ""&>>../$tumorPfx.gatk.log
echo "-----------------------------------------------------------------------------" &>>../$tumorPfx.gatk.log
echo "Re-indexing the realigned bams..." &>>../$tumorPfx.gatk.log
echo ""&>>../$tumorPfx.gatk.log

# index the re-aligned bams
samtools index $tumorPfx.realn.bam
samtools index $bloodPfx.realn.bam

echo ""&>>../$tumorPfx.gatk.log
echo "-----------------------------------------------------------------------------" &>>../$tumorPfx.gatk.log
echo "Beginning base recalibration..." &>>../$tumorPfx.gatk.log
echo ""&>>../$tumorPfx.gatk.log

## base recalibration
# "BaseRecalibrator"
echo "BaseRecalibrator on $tumorPfx: " `date`  >> ../$tumorPfx.gatk.log
java -d64 -jar $gatkJar -nct 8 -T BaseRecalibrator -rf BadCigar -I $tumorPfx.realn.bam -R $hgReference -knownSites $dbSNP -o $tumorPfx.recal.grp &>>../$tumorPfx.gatk.log
echo ""&>>../$tumorPfx.gatk.log
echo "BaseRecalibrator on $bloodPfx: " `date`  >> ../$tumorPfx.gatk.log
java -d64 -jar $gatkJar -nct 8 -T BaseRecalibrator -rf BadCigar -I $bloodPfx.realn.bam -R $hgReference -knownSites $dbSNP -o $bloodPfx.recal.grp &>>../$tumorPfx.gatk.log
echo ""&>>../$tumorPfx.gatk.log
echo "Waiting for BQSR to finish... "`date` >> ../$tumorPfx.gatk.log

echo ""&>>../$tumorPfx.gatk.log
echo "-----------------------------------------------------------------------------" &>>../$tumorPfx.gatk.log
echo "GATK PrintReads..." &>>../$tumorPfx.gatk.log
echo ""&>>../$tumorPfx.gatk.log

# "PrintReads"   
echo "PrintReads on $tumorPfx: " `date` >> ../$tumorPfx.gatk.log
java -d64 -jar $gatkJar -nct 8 -T PrintReads -rf BadCigar -R $hgReference -I $tumorPfx.realn.bam -BQSR $tumorPfx.recal.grp -o $tumorPfx.realn.recal.bam 2>>../$tumorPfx.gatk.log
echo ""&>>../$tumorPfx.gatk.log
echo "PrintReads on $bloodPfx: " `date` >> ../$tumorPfx.gatk.log
java -d64 -jar $gatkJar -nct 8 -T PrintReads -rf BadCigar -R $hgReference -I $bloodPfx.realn.bam -BQSR $bloodPfx.recal.grp -o $bloodPfx.realn.recal.bam 2>>../$tumorPfx.gatk.log
echo ""&>>../$tumorPfx.gatk.log
echo "Waiting for PrintReads to finish... "`date` >> ../$tumorPfx.gatk.log

echo ""&>>../$tumorPfx.gatk.log
echo "-----------------------------------------------------------------------------" &>>../$tumorPfx.gatk.log
echo "Re-indexing the recalibrated bams..." &>>../$tumorPfx.gatk.log
echo ""&>>../$tumorPfx.gatk.log

# index the re-aligned and re-calibrated bams
samtools index $tumorPfx.realn.recal.bam
samtools index $bloodPfx.realn.recal.bam
echo "Indexing recaled and realigned bam files... "`date` >> ../$tumorPfx.gatk.log

echo ""&>>../$tumorPfx.gatk.log
echo "-----------------------------------------------------------------------------" &>>../$tumorPfx.gatk.log
echo "Finally validating the realigned and recalibrated (but not de-duped) bams..." &>>../$tumorPfx.gatk.log
echo ""&>>../$tumorPfx.gatk.log

echo "Tumor bam:" &>>../$tumorPfx.gatk.log
java -d64 -Xmx4g -jar $TACC_PICARD_DIR/ValidateSamFile.jar MODE="SUMMARY" INPUT=$tumorPfx.realn.recal.bam
echo "" &>>../$tumorPfx.gatk.log
echo "Blood bam: " &>>../$tumorPfx.gatk.log
java -d64 -Xmx4g -jar $TACC_PICARD_DIR/ValidateSamFile.jar MODE="SUMMARY" INPUT=$bloodPfx.realn.recal.bam

echo "End $0: `date`">>../$tumorPfx.gatk.log
echo "-----------------------------------------------------------------------------">>../$tumorPfx.gatk.log

exit;

