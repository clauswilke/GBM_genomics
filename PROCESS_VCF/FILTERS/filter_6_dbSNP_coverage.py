## the purpose of this filter is to remove any SNPs with an rsID in dbSNP and coverage less than 10 reads
## this threshold can be increased for finer filtering

## the alternative we discussed with Max, not included here, is to just remove any SNP with coverage less than 20 or some other higher threshold

def filter_six( VCF_dict, VCF_name, filter_dir ):

    import os.path

    # get and parse the list of SNPs, and make a dictionary of all the SNP locations by chromosome


    # the test is twofold, for rs_ID and coverage
    for chrom in VCF_dict:
        for loc in VCF_dict[ chrom ]:
            n_dp = VCF_dict[ chrom ][ loc ][ "read_depth_normal" ]
            t_dp = VCF_dict[ chrom ][ loc ][ "read_depth_tumor" ]
            
    # then it will be for loc on chrom in dict, if ( loc in SNP_dict and reads<=20 ) or (reads <= 15 ): exclude
    # this is super messy, I know


    new_file = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/MUTATION_CALLS/%s/FILT_LISTS/6_excluded.txt" % VCF_name

    print "where is filter six? finally running..."

    return VCF_dict
