#!bin/bash

## the purpose of this script is to gather the summary stats needed to do the last figures and the summary table in the GBM_Genomics paper revision

## environment
# module load samtools

## first you need to loop over each filename
for typeDir in /scratch/01839/dakotaz/alignments/2015-05.TCGA-GBM-53techreps/*/; do
	
	# get the type variable
	type1=${typeDir#/scratch/01839/dakotaz/alignments/2015-05.TCGA-GBM-53techreps/}
	type=${type1%/}
	# echo $type
	
	# get each sample
	for sampDir in $typeDir/*/; do
	
		#get the sampName variable
		sampName1=${sampDir#$typeDir/}
		sampName=${sampName1%/}
		# echo "$type, $sampName"
		
		# now go to the bamfile in the sampDir
		sampbam="$sampDir/$sampName.realn.recal.bam"
		if [ -e $sampbam ]
		then
			# number mapped v. unmapped reads for each reference
			samtools idxstats $sampbam > $sampDir/$sampName.realn.recal.idxstats.txt
			
			## and now the depth of coverage/percent of genome covered stuff
			echo "Genome size (number of bps in the genome, denominator of coverage): " > $sampDir/$sampName.realn.recal.depth.txt
			samtools view -H $sampbam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}' >> $sampDir/$sampName.realn.recal.depth.txt
			samtools view -H $sampbam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}' >> $sampDir/$sampName.realn.recal.depth.txt
			echo "" >> $sampDir/$sampName.realn.recal.depth.txt

			echo "Number of total bases mapped (numerator of coverage, where genome size is denominator):" >> $sampDir/$sampName.realn.recal.depth.txt
			samtools depth $sampbam |  awk '{sum+=$3} END {print sum}' >> $sampDir/$sampName.realn.recal.depth.txt
			echo "" >> $sampDir/$sampName.realn.recal.depth.txt

			echo "The incorrect coverage depth (doesn't include uncovered bases, so is high) with stdev:" >> $sampDir/$sampName.realn.recal.depth.txt
			samtools depth $sampbam | awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $sampDir/$sampName.realn.recal.depth.txt
			samtools depth $sampbam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' >> $sampDir/$sampName.realn.recal.depth.txt
			echo "" >> $sampDir/$sampName.realn.recal.depth.txt

			echo "Number of positions in the genome covered (numerator in % genome covered):" >> $sampDir/$sampName.realn.recal.depth.txt
			samtools depth C484.TCGA-14-0790-10A-01D-1494-08.5.realn.recal.bam | wc -l >> $sampDir/$sampName.realn.recal.depth.txt
			echo "" >> $sampDir/$sampName.realn.recal.depth.txt			
		fi 
	
	done
	
done