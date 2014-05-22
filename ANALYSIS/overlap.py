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

    difference = {}
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
        difference[ chrom ] = []
        if chrom in dictA and chrom in dictB:
            overlap[ chrom ] = []

    # now do the comparison, from the point of view of dictA
    for chrom in dictA:
    	if chrom not in dictB:
    		difference[ chrom ] = dictA[ chrom ]
    		continue
    	else:
        	for loc in dictA[ chrom ]:
        		if loc in dictB[ chrom ]:
        			overlap[ chrom ].append( loc )
        		else:
        			difference[ chrom ].append( loc )
    # and from the point of view of dictB, given all overlap with dictA has already been found
    for chrom in dictB:
    	if chrom not in dictA:
    		difference[ chrom ] = dictB[ chrom ]
    	else:
	    	for loc in dictB[ chrom ]:
    			if loc not in overlap[ chrom ]:
    				difference[ chrom ].append( loc )
    			
    # print overlap

    # sort all and overlap for my peace of mind
    for chrom in difference:
        difference[ chrom ] = sorted( difference[ chrom ] )
    for chrom in overlap:
        overlap[ chrom ] = sorted( overlap[ chrom ] )

    return difference, overlap

def write_files( sample, data_set, filter_level, overlap, difference ):

	data_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/DOUBLES/WGA-S"
	dataset_dir = os.path.join( data_dir, data_set )
	if not os.path.exists( dataset_dir ):
		os.makedirs( dataset_dir )
	samp_dir = os.path.join( dataset_dir, sample )
	if not os.path.exists( samp_dir ):
		os.makedirs( samp_dir )
	
	diff_file = os.path.join( samp_dir, "difference_" + filter_level )
	overlap_file = os.path.join( samp_dir, "overlap_" + filter_level )
	
	diff_fh = open( diff_file, 'w' )
	for chrom in difference:
		for loc in difference[ chrom ]:
			diff_fh.write( "%s\t%s\n" % ( chrom, loc ) )
	diff_fh.close()
	
	overlap_fh = open( overlap_file, 'w' )
	for chrom in overlap:
		for loc in overlap[ chrom ]: 
			overlap_fh.write( "%s\t%s\n" % ( chrom, loc ) )
	overlap_fh.close()

	return 

def summarize_for_figs( sample, data_set, filter_level, overlap, difference, summary ):

	diff_muts = 0
	over_muts = 0
	
	for chrom in overlap:
		for loc in overlap[ chrom ]:
			over_muts += 1
			
	for chrom in difference:
		for loc in difference[ chrom ]:
			diff_muts += 1
	
	all_muts = diff_muts + over_muts
	
	p_o = float( over_muts )/float( over_muts + diff_muts )
	p_d = float( diff_muts )/float( over_muts + diff_muts )
	
	print sample, over_muts, diff_muts, all_muts, p_o, p_d
	summary[ sample ] = [ over_muts, diff_muts, all_muts ]

	return

###################
## MAIN FUNCTION ##
###################

## run command: python overlap.py <data_set> <filter_level>
## example run command: python overlap.py MUTATION_CALLS_STD40 unfiltered.txt

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

summary = {}
# now run the comparison
for sample in useable_double_samples:
	dictA = data_locs( os.path.join( useable_double_samples[ sample ][ 0 ], sys.argv[ 2 ] ) )
	dictB = data_locs( os.path.join( useable_double_samples[ sample ][ 1 ], sys.argv[ 2 ] ) )
	difference, overlap = compare_locs( dictA, dictB )
	write_files( sample, data_set, filter_level, overlap, difference )
	summarize_for_figs( sample, data_set, filter_level, overlap, difference, summary )
	
# write summary to file
figure_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/GBM_genomics/FIGURES/FIGURE_DATA"
figure_data_file = os.path.join( figure_dir, "STD40_" + "unfiltered_" + "overlap.txt" )
figure_fh = open( figure_data_file, 'w' )
figure_fh.write( "SAMPLE\tOVERLAP\tDIFFERENCE\tTOTAL\tPERCENT_OVERLAP\tPERCENT_DIFFERENCE\n" )
for sample in summary:
	over = summary[ sample ][ 0 ]
	diff = summary[ sample ][ 1 ]
	all = summary[ sample ][ 2 ]
	percent_over = str( float( over )/float( all ) )
	percent_diff = str( float( diff )/float( all ) )
	figure_fh.write( "%s\t%s\t%s\t%s\t%s\t%s\n" % ( sample, over, diff, all, percent_over, percent_diff ) )
figure_fh.close()



"""


# and save the data in /GBM_genomics/FIGURES/FIGURE_DATA, so that you can generate figures...
print 

"""