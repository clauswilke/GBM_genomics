# this script takes mutation files of doubles and finds the overlap
## should take different filter levels, will work out how when I have all the levels...

import os

data_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/"
doubles_file = os.path.join( data_dir, "doubles.txt" )
mutation_dir = os.path.join( data_dir, "MUTATION_CALLS" )

## this function checks to make sure that all the amplification v. non-amplification files that should be there, are there
## it returns a dictionary of sample names as keys with two paths (list, sorted) as values
def check_for_amp_doubles():

    doubles_fh = open( doubles_file, 'r' )
    for line in doubles_fh:
        line = line.strip( "\r\n" )
        fields = line.split( "\t" )
        sample_name = 
        file_one = fields[ 1 ][ :17 ]
        file_two = fields[ 2 ][ :17 ]
        path_one = os.path.join( mutation_dir, file_one )
        path_two = os.path.join( mutation_dir, file_two )
        if os.path.exists( path_one ):
            one = "has one"
        if os.path.exists( path_two ):
            two = "has two"
        if one == "has one" and two == "has two":
            double = "TRUE"
        else:
            double = "FALSE"
        print double, one, two

    return 

## this function
def calc_overlaps
