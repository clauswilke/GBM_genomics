# this script takes mutation files of doubles and finds the overlap
## should take different filter levels, will work out how when I have all the levels...

import os

data_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/"
doubles_file = os.path.join( data_dir, "doubles.txt" )
mutation_dir = os.path.join( data_dir, "MUTATION_CALLS" )

## this function checks to make sure that all the amplification v. non-amplification files that should be there, are there
## it returns a dictionary of sample names as keys with two paths (list, sorted) as values
def check_for_amp_doubles():

    amp_doubles = {}
    doubles_fh = open( doubles_file, 'r' )
    for line in doubles_fh:
        line = line.strip( "\r\n" )
        fields = line.split( "\t" )
        sample_name = fields[ 0 ]
        file_one = fields[ 1 ][ :17 ]
        file_two = fields[ 2 ][ :17 ]
        path_one = os.path.join( mutation_dir, file_one )
        path_two = os.path.join( mutation_dir, file_two )

        # if os.path.exists( path_one ):
        #     one = "has one"
        # if os.path.exists( path_two ):
        #     two = "has two"
        # if one == "has one" and two == "has two":
        #     double = "TRUE"
        # else:
        #     double = "FALSE"
        # print double, one, two
        
        amp_doubles[ sample_name ] = sorted( [ path_one, path_two ] )

    return amp_doubles

## this function
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

amp_doubles = check_for_amp_doubles()
for double in amp_doubles:

    # these are the constructs you'll be working with; this will change as you have filt300 and filt100 samples to add to these
    WGA_filt_dict = {}
    WGA_unfilt_dict = {}
    WGS_filt_dict = {}
    WGS_unfilt_dict = {}

    # for each sample for which you have a WGS and a WGA
    # make the WGA (C282) or WGS(C484) dictionaries
    for dir in double:
        if "C282" in dir:
            files = [ os.path.join( dir, "filtered.txt" ), os.path.join( dir, "unfiltered.txt" ) ]
            for file in files:
                if "un" in file and os.path.exists( file ):
                    WGA_unfilt_dict = data_locs( file )
                elif os.path.exists( file ):
                    WGA_filt_dict = data_locs( file )
        elif "C484" in dir:
            files = [ os.path.join( dir, "filtered.txt" ), os.path.join( dir, "unfiltered.txt" ) ]
            for file in files:
                if "un" in file and os.path.exists( file ):
                    WGS_unfilt_dict = data_locs( file )
                elif os.path.exists( file ):
                    WGS_filt_dict = data_locs( file )

    # now run the comparison
    if WGA_unfilt_dict and WGS_unfilt_dict:
        print "yay!"


        #filt_dict, unfilt_dict = data_locs( type, dir )
        #overlap_55 = compare_locs( filt_dict, unfilt_dict )
