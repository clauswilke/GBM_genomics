## the purpose of this filter is to remove any SNPs with an rsID in dbSNP and coverage less than 10 reads
## this threshold can be increased for finer filtering

## the alternative we discussed with Max, not included here, is to just remove any SNP with coverage less than 20 or some other higher threshold

def filter_six( VCF_dict, VCF_name, filter_dir, dbsnp_dir ):

    import os.path

    # the test is twofold, for rs_ID and coverage
    for chrom in VCF_dict:

        # first make a list of all the SNPs in dbSNP for that chrom
        dbsnp_snps = []
        filename = "%s.txt" % chrom
        chrom_file = os.path.join( dbsnp_dir, filename )
        chrom_fh = open( chrom_file, 'r' )
        for line in chrom_fh:
            loc = int( line.strip( "\r\n" ) )
            dbsnp_snps.append( loc )

        # now make a list of all the locations for this chrom in VCF_dict
        # with a sorted list of a tuple of n_dp and t_dp
        vcf_snps = []
        vcf_coverage = []
        for loc in VCF_dict[ chrom ]:
            # print chrom
            # print VCF_dict[ chrom ]
            vcf_snps.append( int( loc ) )
            n_dp = VCF_dict[ chrom ][ loc ][ "read_depth_normal" ]
            t_dp = VCF_dict[ chrom ][ loc ][ "read_depth_tumor" ]
            vcf_coverage.append( [ int( n_dp ), int( t_dp ) ] )
            
        # then compare them... if ( loc in SNP_dict and reads<20 ) or (reads < 15 ): exclude
        # this is super messy, I know
        # print len( vcf_snps )
        # print len( dbsnp_snps )
        i=0
        j=0
        while i < len( dbsnp_snps ):
            # print i, j 
            while j < len( vcf_snps ):
                # print i,j
                if vcf_snps[ j ] > dbsnp_snps[ i ]:
                    i = i+1
                    continue
                if vcf_snps[ j ] < dbsnp_snps[ i ]:
                    # then vcf_snps[ j ] is not a site with a dbSNP entry, so...
                    coverage = vcf_coverage[ j ]
                    if coverage[ 0 ] >=15 and coverage[ 1 ] >= 15:
                        VCF_dict[ chrom ][ loc ][ "FILT_6" ] = "PASS"
                    else:
                        VCF_dict[ chrom ][ loc ][ "FILT_6" ] = "FAIL"
                    j = j+1
                    continue
                if vcf_snps[ j ] == dbsnp_snps[ i ]:
                    # then vcf_snps[ j ] is a site with a dbSNP entry, so...
                    coverage = vcf_coverage[ j ]
                    if coverage[ 0 ] >= 20 and coverage[ 1 ] >= 20:
                        VCF_dict[ chrom ][ loc ][ "FILT_6" ] = "PASS"
                    else:
                        VCF_dict[ chrom ][ loc ][ "FILT_6" ] = "FAIL"
                    j = j+1
                    continue
            # if we reach the end of j before we reach the end of i,
            # what do we do to stop the infinite loop?
            i=i+1

    new_file = os.path.join( filter_dir, "6_excluded.txt" )
    new_fh = open( new_file, 'w' )
    for chrom in VCF_dict:
        for loc in VCF_dict[ chrom ]:
            if VCF_dict[ chrom ][ loc ][ "FILT_6" ] == "FAIL":
                new_fh.write( "%s\t%s\n" % ( chrom, loc ) )


    print "where is filter six? finally running..."

    return VCF_dict
