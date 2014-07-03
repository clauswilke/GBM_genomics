## this script does the same thing as filter_by_filter.py,
## which is to see which filters are removing which mutations,
## but it does so for only the putative mutations in the 55 double samples that are in the unfiltered overlap

import os
import sys
from sets import Set
import VCF_list


data_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/MUTATION_CALLS"
overlap_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/OVERLAP"
fig_data_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/GBM_genomics/FIGURES/FIGURE_DATA"

# this function is taken directly from overlap.py, in this file
# it checks to make sure the files I want to analyze actually exist
def check_file_existence( filter_level, double_samples ):

    useable_double_samples = {}
    num_doub = 0

    for sample in double_samples:
        existence = []
        # print double_samples
        for dir in double_samples[ sample ]:
            file = os.path.join( dir, filter_level )
            if os.path.exists( file ):
                existence.append( "TRUE" )
        if existence == [ "TRUE", "TRUE" ]:
            num_doub += 1
            useable_double_samples[ sample ] = double_samples[ sample ]

    print "there are %g useable samples..." % num_doub
    
    # sort so all the 282s (A) come before all the 484s (B)
    # and in the next set all the XXX (A) come before the XXX (B) 
    for sample in useable_double_samples:
        useable_double_samples[ sample ] = sorted( useable_double_samples[ sample ] )

    return useable_double_samples

def make_SNP_dict( dir, filter_level ):

	SNP_dict = {}
	#SNP_set = Set([])

	# read each line in the file
	file = os.path.join( dir, filter_level )
	file_fh = open( file, 'r' )
	for line in file_fh:
		line = line.strip( '\r\n' )
		fields = line.split( '\t' )
		chrom = fields[ 0 ]
		loc = fields[ 1 ]
		
		# save the info in a dict
		if chrom in SNP_dict:
			SNP_dict[ chrom ].add( loc )
		else:
			SNP_dict[ chrom ] = Set([ loc ])
	
	file_fh.close()

	return SNP_dict

# find the overlap between the two samples, using intersection of sets, and return it in dictionary form
def overlap( SNP_dict_1, SNP_dict_2 ):

	overlap_dict = {}
	for chrom in SNP_dict_1:
		if chrom in SNP_dict_2:
			overlap_dict[ chrom ] = ( SNP_dict_1[ chrom ] & SNP_dict_2[ chrom ] )
		else:
			overlap_dict[ chrom ] = Set([])
	for chrom in SNP_dict_2:
		if chrom not in SNP_dict_1:
			overlap_dict[ chrom ] = Set([])
			
	overlap_length = 0
	for chrom in overlap_dict:
		length = len( overlap_dict[ chrom ] )
		overlap_length = overlap_length + length
			
	return overlap_dict, overlap_length

# find which overlap SNPs were filtered by which filters; this function will run once each for filters two through 8
# note that there may be redundancy; as in one SNP may have been filtered by more than one filter, and so the total number of filtered SNPs may be less than the sum of the filter counts
def overlap_filter_count( i, overlap_dict, sample, overlap_filter_dict ):

	# make the file name for the filter
	filt_file = os.path.join( data_dir, "C282.%s" % sample, "FILT_LISTS", "%s_excluded.txt" % i )
	
	# open the filter file, and make a dictionary of everything filtered
	filt_dict = {}
	filt_fh = open( filt_file, 'r' )
	for line in filt_fh:
		line = line.strip( '\r\n' )
		fields = line.split( '\t' )
		if len( fields ) != 2:
			# print fields
			continue
		chrom = fields[ 0 ]
		loc = fields[ 1 ]
		if chrom not in filt_dict:
			filt_dict[ chrom ] = Set([ loc ])
		else:
			filt_dict[ chrom ].add( loc )
	#print filt_dict
	
	# compare and store
	over_count = 0
	over_filt_count = 0
	for chrom in overlap_dict:
		if chrom not in filt_dict:
			continue
		over_filt_count = over_filt_count + len( overlap_dict[ chrom ] & filt_dict[ chrom ] )
		over_count =  over_count + len( overlap_dict[ chrom ] )
	
	if over_filt_count != 0 and over_count != 0:
		percent_filtered = float( over_filt_count )/ float( over_count )
	else:
		percent_filtered = 0 
	tuple = ( over_filt_count, percent_filtered )
	overlap_filter_dict[ sample ].append( tuple )

	return overlap_filter_dict

def write_overlap( overlap_dict, sample ):

	overlap_file = os.path.join( overlap_dir, "%s.txt" % sample )
	overlap_fh = open( overlap_file, 'w' )
	overlap_fh.write( "%s\tWGS/WGA\n\n" % sample )
	overlap_fh.write( "CHROM\tLOC\n" )
	for chrom in overlap_dict:
		for loc in overlap_dict[ chrom ]:
			overlap_fh.write( "%s\t%s\n" % ( chrom, loc ) )
	overlap_fh.close()

	return

def save_to_figure_data( overlap_filter_dict ):

	new_file = os.path.join( fig_data_dir, "overlap_numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt" )

	new_fh = open( new_file, 'w' )
	new_fh.write( "SAMPLE\tFILT_2_COUNT\tFILT_2_PERCENT\tFILT_3_COUNT\tFILT_3_PERCENT\tFILT_4_COUNT\tFILT_4_PERCENT\tFILT_5_COUNT\tFILT_5_PERCENT\tFILT_6_COUNT\tFILT_6_PERCENT\tFILT_7_COUNT\tFILT_7_PERCENT\tFILT_8_COUNT\tFILT_8_PERCENT\n" )
	for sample in overlap_filter_dict:
		new_fh.write( ( "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" ) % ( sample, overlap_filter_dict[ sample ][ 0 ][ 0 ], overlap_filter_dict[ sample ][ 0 ][ 1 ], overlap_filter_dict[ sample ][ 1 ][ 0 ], overlap_filter_dict[ sample ][ 1 ][ 1 ], overlap_filter_dict[ sample ][ 2 ][ 0 ], overlap_filter_dict[ sample ][ 2 ][ 1 ], overlap_filter_dict[ sample ][ 3 ][ 0 ], overlap_filter_dict[ sample ][ 3 ][ 1 ], overlap_filter_dict[ sample ][ 4 ][ 0 ], overlap_filter_dict[ sample ][ 4 ][ 1 ], overlap_filter_dict[ sample ][ 5 ][ 0 ], overlap_filter_dict[ sample ][ 5 ][ 1 ], overlap_filter_dict[ sample ][ 6 ][ 0 ], overlap_filter_dict[ sample ][ 6 ][ 1 ] ) )	
	new_fh.close()

	return

###############################
######## MAIN FUNCTION ########
###############################

## the output of this file should be, for each double sample, 
## a file called <sample_name>_overlap.txt in the overlap_dir
## and 
## (since by fiat and for ease, this script is doing the counting also)
## the appropriate summaries in "FIGURES/FIGURE_SCRIPTS"

## first get the list of doubles samples and make sure they are all there and working
## this used the functions imported from VCF_list

data_set = "MUTATION_CALLS"
double_samples, doub_count = VCF_list.analyzed_doubles( data_set )

filter_level = "unfiltered.txt"
useable_double_samples = check_file_existence( filter_level, double_samples )

## then what this script is for:
## which is creating a list of the unfiltered_overlap for each double sample

overlap_filter_dict = {}

# print useable_double_samples
for sample in useable_double_samples:
	dir1 = useable_double_samples[ sample ][ 0 ]
	dir2 = useable_double_samples[ sample ][ 1 ]
	SNP_dict_1 = make_SNP_dict( dir1, filter_level )
	#print SNP_dict_1
	SNP_dict_2 = make_SNP_dict( dir2, filter_level )
	#print SNP_dict_2
	overlap_dict, overlap_length = overlap( SNP_dict_1, SNP_dict_2 )
	# print overlap_dict
	# now save the overlap SNPs to a file, just so that you have them
	write_overlap( overlap_dict, sample )
	overlap_filter_dict[ sample ] = [  ]
	for i in range( 2, 9 ):
		overlap_filter_dict = overlap_filter_count( i, overlap_dict, sample, overlap_filter_dict )
	
print overlap_filter_dict
save_to_figure_data( overlap_filter_dict )