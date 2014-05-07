## this script takes mutation files of doubles and finds the overlap

# for now this is only working on the WGS/WGA double files
# and is also only working on the unfiltered data
# eventually (when I have data) it should work on different filter levels and different doubles samples

import os
import VCF_list

data_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/"
mutation_dir = os.path.join( data_dir, "MUTATION_CALLS" )

# next we make sure that all the files that we need are there
# filter_level is basically file_name (e.g. "unfiltered.txt" or "filtered_04.txt")
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
    for chrom in dictB:
        if chrom not in chroms:
            chroms.append( chrom )
    for chrom in dictA:
        if chrom not in chroms:
            chroms.append( chrom )
    for chrom in chroms:
        all[ chrom ] = []
        if chrom in dictA and chrom in dictB:
            overlap[ chrom ] = []

    # now do the compare from the point of view of dictB
    for chrom in dictB:
        for loc in dictB:
            if chrom in dictA:
                if loc in dictA[ chrom ]:
                    overlap[ chrom ].append( loc )
                    all[ chrom ].append( loc )
                else:
                    all[ chrom ].append( loc )
            else:
                all[ chrom ].append( loc )

    # and from the point of view of dictA, given all overlap with dictB has already been found
    for chrom in dictA:
        for loc in dictA[ chrom ]:
            if loc not in all[ chrom ]:
                all[ chrom ].append( loc )

    # sort all and overlap for my peace of mind
    for chrom in all:
        all[ chrom ] = sorted( all[ chrom ] )
    for chrom in overlap:
        overlap[ chrom ] = sorted( overlap[ chrom ] )

    return all, overlap

###################
## MAIN FUNCTION ##
###################

## pretty sure C282 is WGA and C484 WGS. Will double check that with Matt and update this note, however. 

# first thing is to load the doubles files
double_samples, doub_count = VCF_list.analyzed_doubles()

# and check to make sure that all the files you need are there
filter_level = "unfiltered.txt"
useable_double_samples = check_file_existence( filter_level, double_samples )

# now run the comparison
# currently assuming that C282 is WGA and C484 is WGS
for sample in useable_double_samples:
    C282_dict = {}
    C484_dict = {}
    for dir in useable_double_samples[ sample ]:
        file = os.path.join( dir, filter_level )
        if "C282" in file:
            C282_dict = data_locs( file )
            # print "made WGA daict..."
        elif "C484" in file:
            C484_dict = data_locs( file )
            # print "made WGS dict..."
    # print C484_dict
    all, overlap = compare_locs( C484_dict, C282_dict )
    print overlap

# and save the data in /GBM_genomics/FIGURES/FIGURE_DATA, so that you can generate figures...
print 
