GBM_genomics
============

###############
## MAKE_VCFs ##
###############

-- these are the scripts used to make the SomaticSniper VCFs from the bamfiles
-- they were run on stampede

1. gtdownload : this is the script used to download the TCGA bamfiles from cghub
2. bam_to_fastq.sh : this is the script used to run picard to convert bamfiles to fastqs
3. align_fastq.sh : this is the script used to align the re-generated fastqs to hg19 using bwa
4. SNPs_GATK.sh : this is the script used to run GATK base recalibration and indel alignmnet
5. SNPcalling.sh : this is the script used to run SomaticSniper on all samples/home1/01839/dakotaz/GBM_genomics/MAKE_VCF

#################
## JOB_SCRIPTS ##
#################

-- these are the slurm scripts used to submit jobs to stampede

#################
## PROCESS_VCF ##
#################

-- These are the scripts that were used to filter the VCFs
-- they were run on phylocluster

##############
## ANALYSIS ##
##############
        
-- these are the custom python scripts used to analyze the VCF data
-- they were run on phylocluster

1. compare_overlap_to_filtered.py : 
2. compare_size_of_overlap.py :
3. double_compare_number_mutations.py : 
4. filtered_by_filter.py : 
5. germline_snps.py : 
6. how_many_muts.py : 
7. overlap_by_chrom.py : 
8. overlap_only_filter_by_filter.py : 
9. overlap.py : 
10. sort_doubles.py : 
11. VCF_list.py : 

#############
## FIGURES ##
#############

-- these are the custom scripts, run in R, used to make the figures in the paper
-- they were run on a Macbook

###########
## PAPER ##
###########

-- these are the files used to generate the paper and tables in the paper
-- they were all created and run on a Macbook

############
## POSTER ##
############

-- these are the files uesd to generate a poster on this work presented at the Big Data in Biology symposium
-- they were all created and run on a Macbook

