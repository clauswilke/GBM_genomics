# make mutation calls for each sample

from sets import Set
import sys
import os

data_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/"
vcf_dir = "/share/WilkeLab/work/MattCC/TCGA_NextGen_Data/SomaticSniper_Data_FINAL/VCF/"
indel_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/INDELS"
dbsnp_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/DBSNP"
# filter_dir and mutation_dir made in the script

## some VCF information
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=IGT,Number=1,Type=String,Description="Genotype when called independently (only filled if called in joint prior mode)">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
##FORMAT=<ID=BCOUNT,Number=4,Type=Integer,Description="Occurrence count for each base at this site (A,C,G,T)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=JGQ,Number=1,Type=Integer,Description="Joint genotype quality (only filled if called in join prior mode)">
##FORMAT=<ID=VAQ,Number=1,Type=Integer,Description="Variant allele quality">
##FORMAT=<ID=BQ,Number=.,Type=Integer,Description="Average base quality">
##FORMAT=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality across all reads">
##FORMAT=<ID=AMQ,Number=.,Type=Integer,Description="Average mapping quality for each allele present in the genotype">
##FORMAT=<ID=SS,Number=1,Type=Integer,Description="Variant status relative to non-adjacent Normal, 0=wildtype,1=germline,2=somatic,3=LOH,4=unknown">
##FORMAT=<ID=SSC,Number=1,Type=Integer,Description="Somatic Score">

## import filters
sys.path.insert( 0, './FILTERS' )
import filter_2_somatic_score_40
import filter_3_indels
import filter_4_10bp_window
import filter_5_VAQ
import filter_6_dbSNP_coverage
import filter_7_LOH
import filter_8_alt_coverage 

# open the VCF, and put the information in a list called "data"
def make_sample_file_list( vcf_dir ):

    sample_list = []
    # for testing used this: [ "C484.TCGA-32-2616.SS.vcf" ]
    vcf_list = os.path.join( data_dir, "VCFs.txt" )
    vcf_fh = open( vcf_list, 'r' )
    for line in vcf_fh:
        line = line.strip( "\r\n" )
        fields = line.split ( "\t" )
        sample = fields[ 1 ]
        sample_list.append( sample )

    sample_file_list = []
    for sample in sample_list:
        # sample_ext = sample + ".SS.vcf"
        sample_file = os.path.join( vcf_dir, sample )
        sample_file_list.append( sample_file )

    return sample_file_list

def make_data( sample_file ):
    
    sample_fh = open( sample_file, 'r' )
    data = sample_fh.readlines()

    return data

# make the homes where the data will live
def make_data_directories( sample_file ):

    VCF_name = sample_file[ -24:-7 ]
    
    mutation_dir = os.path.join( data_dir, "MUTATION_CALLS", VCF_name )
    if not os.path.exists( mutation_dir ):
        os.makedirs( mutation_dir )
        
    filter_dir = os.path.join( data_dir, "MUTATION_CALLS", VCF_name, 'FILT_LISTS' )
    if not os.path.exists( filter_dir ):
        os.makedirs( filter_dir )
        
    return VCF_name, mutation_dir, filter_dir

# make it a more work-with-able dictionary
def read_data( data ):

    VCF_dict = {}

    for line in data:
        if "##" in line:
            #print line
            continue
        else:
            fields = line.split( '\t' )
            chrom = fields[ 0 ]
            loc = fields[ 1 ]
            # these are ID, QUAL, FILTER, and INFO
            # they pretty much appear to all be blank
            # it might be worth time to ensure that they are all blank (print something, or break if they're not...)
            # blank is a "."
            # print fields[ 2 ], fields[ 5 ], fields[ 6 ], fields[ 7 ]
            ref_allele = fields[ 3 ]
            alt_allele = fields[ 4 ]
            fields_8 = fields[ 8 ].split( ":" )

            fields_normal = fields[ 9 ].split( ":" )
            fields_tumor = fields[ 10 ].split( ":" )

            # putting it all in the dictionary
            if chrom == "#CHROM" or len( fields_8 )<=8:
                continue
            if chrom not in VCF_dict:
                VCF_dict[ chrom ] = {}
                VCF_dict[ chrom ][ loc ] = {}
                VCF_dict[ chrom ][ loc ][ "SS_tumor" ] = int( fields_tumor[ 12 ].strip( '\n' ) )
                VCF_dict[ chrom ][ loc ][ "SS_normal" ] = fields_normal[ 12 ]
                VCF_dict[ chrom ][ loc ][ "VAQ_tumor" ] = fields_tumor[ 7 ]
                VCF_dict[ chrom ][ loc ][ "VAQ_normal" ] = fields_normal[ 7 ]
                VCF_dict[ chrom ][ loc ][ "Variant_status" ] = fields_tumor[ 11 ]
                VCF_dict[ chrom ][ loc ][ "read_depth_tumor" ] = fields_tumor[ 2 ]
                VCF_dict[ chrom ][ loc ][ "read_depth_tumor_detail" ] = fields_tumor[ 3 ]
                VCF_dict[ chrom ][ loc ][ "read_depth_normal" ] = fields_normal[ 2 ]
                # VCF_dict[ chrom ][ loc ][ "read_depth_normal_detail" ] = fields_normal[ 3 ]
                VCF_dict[ chrom ][ loc ][ "genotype_tumor" ] = fields_tumor[ 0 ]
                # VCF_dict[ chrom ][ loc ][ "genotype_normal" ] = fields_normal[ 0 ]
                VCF_dict[ chrom ][ loc ][ "REF" ] = ref_allele
                VCF_dict[ chrom ][ loc ][ "ALT" ] = alt_allele
            else:
                VCF_dict[ chrom ][ loc ] = {}
                VCF_dict[ chrom ][ loc ][ "SS_tumor" ] = int( fields_tumor[ 12 ].strip( '\n' ) )
                VCF_dict[ chrom ][ loc ][ "SS_normal" ] = fields_normal[ 12 ]
                VCF_dict[ chrom ][ loc ][ "VAQ_tumor" ] = fields_tumor[7 ]
                VCF_dict[ chrom ][ loc ][ "VAQ_normal" ] = fields_normal[ 7 ]
                VCF_dict[ chrom ][ loc ][ "Variant_status" ] = fields_tumor[ 11]
                VCF_dict[ chrom ][ loc ][ "read_depth_tumor" ] = fields_tumor[ 2 ]
                VCF_dict[ chrom ][ loc ][ "read_depth_tumor_detail" ] = fields_tumor[ 3 ]
                VCF_dict[ chrom ][ loc ][ "read_depth_normal" ] = fields_normal[ 2 ]
                # VCF_dict[ chrom ][ loc ][ "read_depth_normal_detail" ] = fields_normal[ 3 ]
                VCF_dict[ chrom ][ loc ][ "genotype_tumor" ] = fields_tumor[ 0 ]
                # VCF_dict[ chrom ][ loc ][ "genotype_normal" ] =fields_normal[ 0 ]
                VCF_dict[ chrom ][ loc ][ "REF" ] = ref_allele
                VCF_dict[ chrom ][ loc ][ "ALT" ] = alt_allele

    return VCF_dict

# write results to file
#################********************* this function has not yet been tested because all the filters aren't finished, so right now it doesn't work
def mutation_calls( VCF_dict, VCF_name, summary ):

    # all the mutations in the VCF dict
    total_unfiltered = 0
    total_filtered = 0

    unfiltered_file = os.path.join ( mutation_dir, "unfiltered.txt" )
    filtered_file = os.path.join( mutation_dir, "filtered.txt" )

    unfiltered_fh = open( unfiltered_file, 'w' )
    filtered_fh = open( filtered_file, 'w' )

    for chrom in VCF_dict:
        for loc in VCF_dict[ chrom ]:
            total_unfiltered = total_unfiltered + 1
            unfiltered_fh.write( "%s\t%s\n" % ( chrom, loc ) )
            
            if VCF_dict[ chrom ][ loc ][ "FILT_2" ] == "PASS":
                if VCF_dict[ chrom ][ loc ][ "FILT_3" ] == "PASS":
                    if VCF_dict[ chrom ][ loc ][ "FILT_4" ] == "PASS":
                        if VCF_dict[ chrom ][ loc ][ "FILT_5" ] == "PASS":
                            if VCF_dict[ chrom ][ loc ][ "FILT_6" ] == "PASS":
                                if VCF_dict[ chrom ][ loc ][ "FILT_7" ] == "PASS":
                                    if VCF_dict[ chrom ][ loc ][ "FILT_8" ] == "PASS":
                                        total_filtered = total_filtered + 1
                                        filtered_fh.write( "%s\t%s\n" % ( chrom, loc ) )

    unfiltered_fh.close()
    filtered_fh.close()

    summary[ VCF_name ] = {}
    summary[ VCF_name ][ "UNFILT" ] = total_unfiltered
    summary[ VCF_name ][ "FILT" ] = total_filtered
    # summary
    # when filters done, add to the summary how many passed each filter, also as columns
    # summary_file = os.path.join( mutation_dir, "summary.txt" )
    # summary_fh = open( summary_file, 'w' )
    # summary_fh.write( "SAMPLE\tUNFILTERED\tFILTERED\n" )
    # summary_fh.write( "%s\t%s\t%s\n" % ( VCF_sample, total_unfiltered, total_filtered ) )
    # summary_fh.close()

    return

###################
## MAIN FUNCTION ##
###################

summary = {}

sample_file_list = make_sample_file_list( vcf_dir )  ## makes the list of all the VCFs that will be analyzed

for sample_file in sample_file_list:

    if not os.path.exists( sample_file ):
        print "does not exist; try again: " + sample_file
        continue
    
    VCF_name, mutation_dir, filter_dir = make_data_directories( sample_file ) ## makes the directories in ../DATA/ where the output of this script will be stored
    print "begin processing %s..." % VCF_name
    data = make_data( sample_file ) ## reads the VCF into a file_header
    VCF_dict = read_data( data ) ## puts the relevant data into dictionary form

    ## run the filters of VCF dict
    # filter one, removes all reads with mapping quality <40, is part of makeing the VCF and does not need to be repeated here
    VCF_dict = filter_2_somatic_score_40.filter_two( VCF_dict, VCF_name, filter_dir )  ## runs filter two, no SNP with somatic score <50 in blood or normal tissue
    VCF_dict = filter_3_indels.filter_three( VCF_dict, VCF_name, indel_dir, mutation_dir, filter_dir ) ## this runs filter three, no SNP w/in 10 bp of a predicted indel of quality >=50
    VCF_dict = filter_4_10bp_window.filter_four( VCF_dict, VCF_name, filter_dir ) ## runs filter four, no SNP w/in 10 bp of another SNP, and saves output
    VCF_dict = filter_5_VAQ.filter_five( VCF_dict, VCF_name, filter_dir ) ## runs filter five, no SNP with VAQ <= 20
    VCF_dict = filter_6_dbSNP_coverage.filter_six( VCF_dict, VCF_name, filter_dir, dbsnp_dir ) ## this runs filter six, removing sites with rsIDs if coverage is very low
    VCF_dict = filter_7_LOH.filter_seven( VCF_dict, VCF_name, filter_dir ) ## runs filter seven, not counting any LOH 'cause they're probably artifacts
    VCF_dict = filter_8_alt_coverage.filter_eight( VCF_dict, VCF_name, filter_dir ) ## runs filter eight gets rid of all heterozygous sites where the coverage of the alternate allele is less than 10% of the major allele

    # print VCF_dict    
    mutation_calls( VCF_dict, VCF_name, summary )

    print "have now processed " +  VCF_name

print summary
summary_file = os.path.join( data_dir, "MUTATION_CALLS", "summary.txt" )
summary_fh = open( summary_file, 'w' )
summary_fh.write( "SAMPLE\tUNFILTERED\tFILTERED\n" )
for sample in summary:
    summary_fh.write( "%s\t%s\t%s\n" % ( sample, summary[ sample ][ "UNFILT" ], summary[ sample ][ "FILT" ] ) )
summary_fh.close()


################################################################################
## Somatic Mutation Filtering                                                 ##
##                                                                            ##
## 1. Filter reads with mapping quality < 40 (Matt did this)                  ##
## 2. Filter reads with somatic quality score < 40                            ##
## 3. mutations less than 10 bp from an indel with quality score >50          ##
## 4. mutations with 3 or more SNV calls in a 10 bp window around site        ##
## 5. Read quality score                                                      ##
##      a. SSnipter: VAQ >= 20                                                ##
##      b. GATK Quality score >= THRESHOLD                                    ##
## 6. dbSNP filter                                                            ##
##      a. filter out anything that has an rsID and anything with coverage<10 ##
##      b. filter out anything that has coverage<20 (nothing with rsID)       ##
## 7. Extra quality                                                           ##
##      a. filter 0/1 blood, when tumor is 1/1 or 0/0 (LOH)                   ##
##      b. skip 8                                                             ##
## 8. Alternate allele coverage                                               ##
##      a. Skip                                                               ##
##      b. alternate allele must have >= 10% of total reads at site           ##
################################################################################
