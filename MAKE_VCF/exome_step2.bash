#!/bin/bash

# module load gatk/2.5.2
module load samtools
module load picard/1.92

tumorPfx=$1
tumorSample=${tumorPfx:0:17}
echo $tumorPfx $tumorSample
bloodPfx=$2
bloodSample=${bloodPfx:0:17}
echo $bloodPfx $bloodSample
if [ $tumorSample != $bloodSample ]; then
    echo "Non-matching sample names:"
    echo "    $tumorSample"
    echo "    $bloodSample"
    exit 1;
fi

gatkJar="$HOME/bin/gatk-2.7.1/GenomeAnalysisTK.jar"
refDir="$WORK/Hs_reference_datasets"
hgReference="$refDir/Homo_sapiens.GRCh37.72.dna.fa"
dbsnp="$refDir/dbsnp_137.b37.vcf"
hapMap="$refDir/hapmap_3.3.b37.vcf"
G1000_HiConf_SNPs="$refDir/1000G_phase1.snps.high_confidence.b37.vcf"
G1000_Phase1_Indels="$refDir/1000G_phase1.indels.b37.vcf"
G1000_Omni_SNPs="$refDir/1000G_omni2.5.b37.vcf"
G1000_Mills="$refDir/Mills_and_1000G_gold_standard.indels.b37.vcf"

echo "Generating somatic calls for $tumorSample" > $tumorSample.log 

# make temporary director
mkdir $tumorSample
cd $tumorSample
samtools rmdup ../$tumorPfx/$tumorPfx.sort.bam $tumorPfx.bam &
samtools rmdup ../$bloodPfx/$bloodPfx.sort.bam $bloodPfx.bam &
wait;

# re-index the de-duped bamfile
samtools index $tumorPfx.bam &
samtools index $bloodPfx.bam &
wait;

# local realignment around indels -- generate common interval file for blood and tumor
echo "Start RealignerTargetCreator: " >> ../$tumorSample.log
java -d64 -jar $gatkJar -R $hgReference -nt 12 -T RealignerTargetCreator -rf BadCigar  \
     -known $G1000_Mills -known $G1000_Phase1_Indels -o $tumorSample.intervals \
     -I $tumorPfx.bam -I $bloodPfx.bam &>>../$tumorSample.log
echo "$0: `date`"

echo "IndelRealigner: " `date` >> ../$tumorSample.log
java -d64 -jar $gatkJar -R $hgReference -T IndelRealigner -rf BadCigar -I $tumorPfx.bam \
    -known $G1000_Mills -known $G1000_Phase1_Indels -targetIntervals $tumorSample.intervals \
    -o $tumorPfx.realn.bam &>>../$tumorSample.log &

java -d64 -jar $gatkJar -R $hgReference -T IndelRealigner -rf BadCigar -I $bloodPfx.bam \
    -known $G1000_Mills -known $G1000_Phase1_Indels -targetIntervals $tumorSample.intervals \
    -o $bloodPfx.realn.bam &>>../$tumorSample.log &
wait;

#index the realigned BAMS
samtools index $tumorPfx.realn.bam
samtools index $bloodPfx.realn.bam

echo "BaseRecalibrator on $tumorPfx: " `date`  >> ../$tumorSample.log
java -d64 -jar $gatkJar -nct 8 -T BaseRecalibrator -rf BadCigar -I $tumorPfx.realn.bam \
    -R $hgReference -knownSites $dbsnp -o $tumorPfx.recal.grp &>>../$tumorSample.log &

echo "BaseRecalibrator on $bloodPfx: " `date`  >> ../$tumorSample.log
java -d64 -jar $gatkJar -nct 8 -T BaseRecalibrator -rf BadCigar -I $bloodPfx.realn.bam \
    -R $hgReference -knownSites $dbsnp -o $bloodPfx.recal.grp &>>../$tumorSample.log &

echo "Waiting for BQSR to finish... "`date` >> ../$tumorSample.log
wait;
    
echo "PrintReads on $tumorPfx: " `date` >> ../$tumorSample.log
java -d64 -jar $gatkJar -nct 8 -T PrintReads -rf BadCigar -R $hgReference -I $tumorPfx.realn.bam \
    -BQSR $tumorPfx.recal.grp -o $tumorPfx.realn.recal.bam 2>>../$tumorSample.log &

echo "PrintReads on $bloodPfx: " `date` >> ../$tumorSample.log
java -d64 -jar $gatkJar -nct 8 -T PrintReads -rf BadCigar -R $hgReference -I $bloodPfx.realn.bam \
    -BQSR $bloodPfx.recal.grp -o $bloodPfx.realn.recal.bam 2>>../$tumorSample.log &

echo "Waiting for PrintReads to finish... "`date` >> ../$tumorSample.log
wait;

samtools index $tumorPfx.realn.recal.bam &
samtools index $bloodPfx.realn.recal.bam &
echo "Indexing recaled and realigned bam files... "`date` >> ../$tumorSample.log
wait;

# run somatic sniper on the T-N pairs
echo "SomaticSniper: " `date` >> ../$tumorSample.log
~/bin/bam-somaticsniper -q 40 -Q 40 -J -s 0.001 -F vcf -f $hgReference \
    $tumorPfx.realn.recal.bam $bloodPfx.realn.recal.bam $tumorSample.SS.vcf

if [ -f $tumorSample.SS.vcf ];
then
    cd ..
    mv $tumorSample hg19_bams
    mv $tumorSample.log hg19_bams/$tumorSample/
    echo "Step2 Completed: " `date` >> hg19_bams/$tumorSample/$tumorSample.log
    exit 0;
fi;

exit 1;
