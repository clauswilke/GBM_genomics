#!/bin/bash

# load modules (gatk current default is 2.7.2)
module load gatk/2.7.2
module load samtools
module load picard 

## this script does the GATK processing and variant calibration on all of the sorted .bam and .bam.bai files

cd /scratch/00001/mattcowp/dakota

## get sample names and file names, and make sure the two files have the same sample name
tumor_pfx=$1
blood_pfx=$2
tumor_sample=${tumor_pfx:0:17}
blood_sample=${blood_pfx:0:17}
tumor_dir="$tumor_pfx.WORK"
# echo $tumor_dir
echo $tumor_pfx $tumor_sample
echo $blood_pfx $blood_sample
if [ $tumor_sample != $blood_sample ]; then
	echo "Non-matching sample names:"
	echo "	$tumor_sample"
	echo "	$blood_sample"
	exit 1;
fi

## locations of relevant and needed files
gatk_jar="$TACC_GATK_DIR/GenomeAnalysisTK.jar"
ref_dir="/work/00001/mattcowp/Hs_reference_datasets"
hgReference="$ref_dir/Homo_sapiens.GRCh37.72.dna.fa"
dbSNP="$ref_dir/dbsnp_137.b37.vcf"
hapMap="$ref_dir/hapmap_3.3.b37.vcf"
G1000_HiConf_SNPs="$ref_dir/1000G_phase1.snps.high_confidence.b37.vcf"
G1000_Phase1_Indels="$ref_dir/1000G_phase1.indels.b37.vcf"
G1000_Omni_SNPs="$ref_dir/1000G_omni2.5.b37.vcf"
G1000_Mills="$ref_dir/Mills_and_1000G_gold_standard.indels.b37.vcf"

echo "Generating somatic calls for $tumor_pfx" > $tumor_pfx.gatk.log

# make a temporary directory in which this will actually run
mkdir $tumor_dir
cd $tumor_dir

# remove duplicate reads from bamfiles, leaving only the read with the highest map quality
samtools rmdup ../$tumor_pfx/$tumor_pfx\_out.sorted.bam $tumor_pfx.dedup.bam &
samtools rmdup ../normal/$blood_pfx/$blood_pfx\_out.sorted.bam $blood_pfx.dedup.bam &
wait;

# reindex the de-duplicated bamfiles
samtools index $tumor_pfx.dedup.bam &
samtools index $blood_pfx.dedup.bam &
wait;

## indel re-alignment
# generate a common interval file for the blood and tumor data
echo "Start RealignerTargetCreator: " >> ../$tumor_pfx.gatk.log
java -d64 -jar $gatk_jar -R $hgReference -nt 12 -T RealignerTargetCreator -rf BadCigar -known $G1000_Mills -known $G1000_Phase1_Indels -o $tumor_pfx.intervals -I $tumor_pfx.dedup.bam -I $blood_pfx.dedup.bam &>>../$tumor_pfx.gatk.log
echo "$0: `date`"

# (local realignment around indels) tumor
echo "IndelRealigner: " `date` >> ../$tumor_pfx.gatk.log
java -d64 -jar $gatk_jar -R $hgReference -T IndelRealigner -rf BadCigar -I $tumor_pfx.dedup.bam -known $G1000_Mills -known $G1000_Phase1_Indels -targetIntervals $tumor_pfx.intervals -o $tumor_pfx.realn.bam &>>../$tumor_pfx.gatk.log &

# (local realignment around indels) blood
java -d64 -jar $gatk_jar -R $hgReference -T IndelRealigner -rf BadCigar -I $blood_pfx.dedup.bam -known $G1000_Mills -known $G1000_Phase1_Indels -targetIntervals $tumor_pfx.intervals -o $blood_pfx.realn.bam &>>../$tumor_pfx.gatk.log &
wait;

# index the re-aligned bams
samtools index $tumor_pfx.realn.bam
samtools index $blood_pfx.realn.bam

## base recalibration
# "BaseRecalibrator"
echo "BaseRecalibrator on $tumor_pfx: " `date`  >> ../$tumor_pfx.gatk.log
java -d64 -jar $gatk_jar -nct 8 -T BaseRecalibrator -rf BadCigar -I $tumor_pfx.realn.bam -R $hgReference -knownSites $dbSNP -o $tumor_pfx.recal.grp &>>../$tumor_pfx.gatk.log &

echo "BaseRecalibrator on $blood_pfx: " `date`  >> ../$tumor_pfx.gatk.log
java -d64 -jar $gatk_jar -nct 8 -T BaseRecalibrator -rf BadCigar -I $blood_pfx.realn.bam -R $hgReference -knownSites $dbSNP -o $blood_pfx.recal.grp &>>../$tumorSample.log &

echo "Waiting for BQSR to finish... "`date` >> ../$tumor_pfx.gatk.log
wait;
 
# "PrintReads"   
echo "PrintReads on $tumor_pfx: " `date` >> ../$tumor_pfx.gatk.log
java -d64 -jar $gatk_jar -nct 8 -T PrintReads -rf BadCigar -R $hgReference -I $tumor_pfx.realn.bam -BQSR $tumor_pfx.recal.grp -o $tumor_pfx.realn.recal.bam 2>>../$tumor_pfx.gatk.log &

echo "PrintReads on $blood_pfx: " `date` >> ../$tumor_pfx.gatk.log
java -d64 -jar $gatk_jar -nct 8 -T PrintReads -rf BadCigar -R $hgReference -I $blood_pfx.realn.bam -BQSR $blood_pfx.recal.grp -o $blood_pfx.realn.recal.bam 2>>../$tumor_pfx.gatk.log &

echo "Waiting for PrintReads to finish... "`date` >> ../$tumor_pfx.gatk.log
wait;

# index the re-aligned and re-calibrated bams
samtools index $tumor_pfx.realn.recal.bam &
samtools index $blood_pfx.realn.recal.bam &
echo "Indexing recaled and realigned bam files... "`date` >> ../$tumor_pfx.gatk.log
wait;

exit;

