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
        print double_samples[ sample ]
        for dir in double_samples[ sample ]:
            file = os.path.join( dir, filter_level )
            if os.path.exists( file ):
                existence.append( "TRUE" )
        if existence == [ "TRUE", "TRUE" ]:
            num_doub += 1
            useable_double_samples[ sample ] = double_samples[ sample ]

    print "there are" + num_doub + "useable samples"

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
def compare_locs( filt_dict, unfilt_dict ):

    overlap = {}

    

    return overlap

###################
## MAIN FUNCTION ##
###################

## pretty sure C282 is WGA and C484 WGS. Will double check that with Matt and update this note, however. 

# first thing is to load the doubles files
double_samples = VCF_list.analyzed_doubles()

# and check to make sure that all the files you need are there
filter_level = "unfiltered.txt"
useable_double_samples = check_file_existence( filter_level, double_samples )

# now run the comparison



# and save the data in /ANALYSIS/FIGURES, so that you can generate figures...
