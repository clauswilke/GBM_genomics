## this is filter five, which gets rid of any SNP that has a VAQ <= 20
## where VAQ is the Variant Allele Quality, some parameter calculated by SomaticSniper
## this is the same as the SomaticSniper where we got all the filtering information from

## tumor and blood appear to have different VAQs for each variant;
## for this filter, they must both exceed 20

## when we discussed this with Max, he proposed an alternative that I have declined to code, partly becuase I don't know how
## his alternative is to filter everything with a GATK quality score less than a certain THRESHOLD (which is not specified)

def filter_five( VCF_dict, VCF_name ):

    new_file = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/MUTATION_CALLS/%s/FILT_LISTS/5_excluded.txt" % VCF_name
    new_fh = open( new_file, 'w' )

    # the filter and adding it to the VCF_dict are practically the same thing, so I've combined them into one step here
    # also with writing the file
    for chrom in VCF_dict:
        for loc in VCF_dict[ chrom ]:
            if int( VCF_dict[ chrom ][ loc ][ "VAQ_tumor" ] ) >=20 and int( VCF_dict[ chrom ][ loc ][ "VAQ_normal"] ) >= 20:
                VCF_dict[ chrom ][ loc ][ "FILT_5" ] = "PASS"
            else:
                VCF_dict[ chrom][ loc ][ "FILT_5" ] = "FAIL"
                new_fh.write( "%s\t%s\n" % ( chrom, loc ) )

    new_fh.close()
    
    print "filter five has run"

    return VCF_dict
