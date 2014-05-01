# this is filter four
# it finds all mutations that are within a 10 bp window of another mutation

# you must import Set from sets to run this function

def filter_four( VCF_dict, VCF_name, filter_dir ):

    from sets import Set
    import os.path

    # make a dictionary all_locs[ chrom ] : [ list of integer locations on the chromosome with a SNP ]
    all_locs = {}
    filtered_4_10bp_window = {}
    for chrom in VCF_dict:
        filtered_4_10bp_window[ chrom ] = []
        all_locs[ chrom ] = []

        for loc in VCF_dict[ chrom ]:
            int_loc = int( loc )
            all_locs[ chrom ].append( int_loc )

        # figure out which ones don't pass the filter
        sort_all_locs = sorted( all_locs[ chrom ] )
        filt_set = Set([])
        for i in range( 1, len( sort_all_locs )-1 ):
            loc = sort_all_locs[ i ]
            lower_bound = loc-10
            if sort_all_locs[ i-1 ] >= lower_bound:
                filt_set.add( sort_all_locs[ i-1 ] )
                filt_set.add( sort_all_locs[ i ] )

        # now store those values in the dictionary
        for loc in filt_set:
            filtered_4_10bp_window[ chrom ].append( loc )

    # save the output in a file in DATA/sample/filters/filter_4_excluded.txt
    new_file = os.path.join( filter_dir, "4_excluded.txt" )
    new_fh = open( new_file, 'w' )
    for chrom in filtered_4_10bp_window:
        for loc in filtered_4_10bp_window[ chrom ]:
            new_fh.write( "%s\t%s\n" % ( chrom, loc ) )
            
    # and add a column to the VCF_dict
    for chrom in VCF_dict:
        for loc in VCF_dict[ chrom ]:
            if int( loc ) in filtered_4_10bp_window[ chrom ]:
                VCF_dict[ chrom ][ loc ][ "FILT_4" ] = "FAIL"
            else:
                VCF_dict[ chrom ][ loc ][ "FILT_4" ] = "PASS"

    print "filter four has run"

    return VCF_dict
