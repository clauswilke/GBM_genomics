# this is filter 2
# it removes any SNP that has a somatic score (somatic quality score) <50 for either the tumor or the normal sample

# we originally said 40, but I changed it to 50
# this was a completely arbitrary decision based on the fact that I wasn't catching anything at 40
# and also on the fact that Matt wanted to me get the number of mutations per file somewhat below 300, closer to 100, so I needed to catch more

def filter_two( VCF_dict, VCF_name, filter_dir ):

    import os.path

    filtered_2_SSC_50 = {}

    for chrom in VCF_dict:
        filtered_2_SSC_50[ chrom ] = []
        for loc in VCF_dict[ chrom ]:
            SS_tumor = int( VCF_dict[ chrom ][ loc ][ "SS_tumor" ] ) 
            ## SS_normal = VCF_dict[ chrom ][ loc ][ "SS_normal" ]  ## this was removed because it winds up being an empty field

            if SS_tumor < 50:
                filtered_2_SSC_50[ chrom ].append( loc )
            #print SS_tumor, type( SS_tumor )

    # save the output in DATA/mutation_calls/sample/filt_lists/2_excluded.txt  
    new_file = os.path.join( filter_dir, "2_excluded.txt" )
    new_fh = open( new_file, 'w' )
    for chrom in filtered_2_SSC_50:
        for loc in filtered_2_SSC_50[ chrom ]:
            new_fh.write( "%s\t%s\n" % ( chrom, loc ) )

    # and add a value to the dictionary VCF_dict for passing or failing the filter
    for chrom in VCF_dict:
        for loc in VCF_dict[ chrom ]:
            if loc in filtered_2_SSC_50[ chrom ]:
                VCF_dict[ chrom ][ loc ][ "FILT_2" ] = "FAIL"
            else:
                VCF_dict[ chrom][ loc ][ "FILT_2" ] = "PASS"
    
    print "filter two has run"
            
    return VCF_dict
