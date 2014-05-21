## this script takes mutation files of doubles and finds the overlap
## then stores the results in ...

# for now this is only working on the WGS/WGA double files
# for data_set and filter_level, see command line argument explanations in "MAIN FUNCTION"

import os
import VCF_list
import sys

data_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/"
mutation_dir = os.path.join( data_dir, "MUTATION_CALLS" )

# next we make sure that all the files that we need are there
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

    return useable_double_samples

## this function makes a dictionary of all the mutations in a file, which can then be fed into a comparison function
def data_locs( file ):

    dict = {}

    file_fh = open( file, 'r' )
    for line in file_fh:
        line = line.strip( '\r\n' )
        fields = line.split( '\t' )
        chrom = fields[ 0 ]
        loc = int( fields[ 1 ] )

        if chrom not in dict:
            dict[ chrom ] = [ loc ]
        else:
            dict[ chrom ].append( loc )

    for chrom in dict:
        dict[ chrom ] = sorted( dict[ chrom ] )
        # print dict[ chrom ]

    return dict

## this function
## the hypothesis states that a substantial portion of the mutations in dictA (time_0 and WGS) should also be in dictB (time_1 and WGA)
def compare_locs( dictA, dictB ):

    all = {}
    overlap = {}

    # set up chrom keys in overlap and all
    chroms = []
    for chrom in dictA:
        if chrom not in chroms:
            chroms.append( chrom )
    for chrom in dictB:
        if chrom not in chroms:
            chroms.append( chrom )
    for chrom in chroms:
        all[ chrom ] = []
        if chrom in dictA and chrom in dictB:
            overlap[ chrom ] = []

    # now do the comparison, from the point of view of dictA
    for chrom in dictA:
    	if chrom not in dictB:
    		all[ chrom ] = dictA[ chrom ]
    		continue
    	else:
        	for loc in dictA[ chrom ]:
        		if loc in dictB[ chrom ]:
        			overlap[ chrom ].append( loc )
        			all[ chrom ].append( loc )
        		else:
        			all[ chrom ].append( loc )
    # and from the point of view of dictB, given all overlap with dictA has already been found
    for chrom in dictB:
    	for loc in dictB[ chrom ]:
    		if loc not in all[ chrom ]:
    			all[ chrom ].append( loc )
    			
    # print overlap

    # sort all and overlap for my peace of mind
    for chrom in all:
        all[ chrom ] = sorted( all[ chrom ] )
    for chrom in overlap:
        overlap[ chrom ] = sorted( overlap[ chrom ] )

    return all, overlap

###################
## MAIN FUNCTION ##
###################

# this function runs with several command line options, as follows:
# sys.argv1 = data_set, and should refer to which filter settings were used, and the file that houses that dataset
## it should be "MUTATION_CALLS" (original), "MUTATION_CALLS_STD40 (the std40 settings), ...
# sys.argv[2] = filter_level, and should be either "filtered.txt" or "unfiltered.txt"
## this gives you the filtered data set or the unfiltered data set for the run that you are analyzing

## pretty sure C282 is WGA and C484 WGS. Will double check that with Matt and update this note, however. 

## debugging...
# print sys.argv[ 0 ]
# print sys.argv[ 1 ], type( sys.argv[ 1 ] )
# print sys.argv[ 2 ], type( sys.argv[ 2 ] )

# first thing is to load the doubles files
data_set = sys.argv[ 1 ]
double_samples, doub_count = VCF_list.analyzed_doubles( data_set )
# print doub_count

# and check to make sure that all the files you need are there
filter_level = sys.argv[ 2 ]
useable_double_samples = check_file_existence( filter_level, double_samples )

# now run the comparison
# currently assuming that C282 is WGA and C484 is WGS
for sample in useable_double_samples:
    C282_dict = {}
    C484_dict = {}
    dirs = useable_double_samples[ sample ]
    
    # get each of the two file_names
    for dir in dirs:
        file = os.path.join( dir, filter_level )
        
        # analyze the C282 file
        if "C282" in file:
        	C282_dict = data_locs( file )
        	print "made C282 dict..."
        elif "C484" in file:
        	C484_dict = data_locs( file )
        	print "made C484 dict..."
        else:
        	print "error... not a C282 or a C484... what?"
    # print C282_dict
    # print C484_dict

	# and finally get around to comparing the two...
	all, overlap = compare_locs( C282_dict, C484_dict )
	print all, overlap

"""


# and save the data in /GBM_genomics/FIGURES/FIGURE_DATA, so that you can generate figures...
print 

"""