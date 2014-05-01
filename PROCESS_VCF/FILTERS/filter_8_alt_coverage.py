## tis filter removes all SNPs where the tumor is heterozygous (0/1) and
## coverage of the alternate allele is less than 10% of the coverage of the major allele
## it is used in place of 7, according to the talk we had with Max that one time
## but we can always just use them all, also

## especially becuase I'm not convinced this one actually does anything at all. I think SS does this step, inclusively

##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">

def filter_eight( VCF_dict, VCF_name, filter_dir ):

    import os.path

    existence = []

    for chrom in VCF_dict:
        for loc in VCF_dict[ chrom ]:
            tumor_detail = VCF_dict[ chrom ][ loc ][ "read_depth_tumor_detail" ].split( ',' )
            t_ref_f = tumor_detail[ 0 ]
            t_ref_r = tumor_detail[ 1 ]
            t_alt_f = tumor_detail[ 2 ]
            t_alt_r = tumor_detail[ 3 ]
            t_gt = VCF_dict[ chrom ][ loc ][ "genotype_tumor" ]
            REF = VCF_dict[ chrom ][ loc ][ "REF" ]
            ALT = VCF_dict[ chrom ][ loc ][ "ALT" ]

            if t_gt == "0/1" or t_gt == "0/2":
                ref_cov = float( t_ref_r ) + float( t_ref_f )
                alt_cov = float( t_alt_r ) + float( t_alt_f )
                # print t_gt, ref_cov, alt_cov
                if alt_cov <= ref_cov/10.0:
                    existence.append( [ chrom, loc ] )
                    VCF_dict[ chrom ][ loc ][ "FILT_8" ] = "FAIL"
                else:
                    VCF_dict[ chrom ][ loc ][ "FILT_8" ] = "PASS"
            else:
                VCF_dict[ chrom ][ loc ][ "FILT_8" ] = "PASS"

    # writing it all to file...
    new_file = os.path.join( filter_dir, "8_excluded.txt" )
    new_fh = open( new_file, 'w' )
    if len( existence ) >= 1:
        for SNP in existence:
            new_fh.write( "%s\t%s\n" % ( SNP[ 0 ], SNP[ 1 ] ) )
    else:
        new_fh.write( "NONE" )
    new_fh.close

                    
    print "filter eight has run"

    return VCF_dict
