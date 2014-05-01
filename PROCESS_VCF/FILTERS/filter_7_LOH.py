## the purpose of filter 7 is to remove anything that is LOH (loss of heterozygosity),
## I believe on the assumption that these variants are more likely to be artifacts than they are to be real
## functionally, this filters out any SNP whose tumor genotype is 0/0 or 1/1, if the blood is 0/1   

## the alternative to this filter is filter 8, as we talked about it with Max, but I can do whatever I want

def filter_seven( VCF_dict, VCF_name, filter_dir ):

    import os.path

    new_file = os.path.join( filter_dir, "7_excluded.txt" )
    new_fh = open( new_file, 'w' )

    for chrom in VCF_dict:
        for loc in VCF_dict[ chrom ]:
            if VCF_dict[ chrom ][ loc ][ "Variant_status" ] == "3":
                VCF_dict[ chrom ][ loc ][ "FILT_7" ] = "FAIL"
                new_fh.write( "%s\t%s\n" % ( chrom, loc ) )
            else:
                VCF_dict[ chrom ][ loc ][ "FILT_7" ] = "PASS"

    print "filter seven has run"

    return VCF_dict
