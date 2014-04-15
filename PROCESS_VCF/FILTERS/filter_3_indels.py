# this is filter three
# it finds all mutations within 10 bp of an indel of quality 50 or greater
# it also parses the idel file, cause that had to happe sometime

def filter_three( VCF_dict, VCF_name ):

    normal_indel_file = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/INDELS/%s_blood.txt" % VCF_name
    tumor_indel_file = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/INDELS/%s_tumor.txt" % VCF_name

    # make a chrom dictionary of loc lists of indel locations, the lists should be sorted
    indel_locs = {}
    tumor_fh = open( tumor_indel_file, 'r' )
    normal_fh = open( normal_indel_file, 'r' )
    headers = [ tumor_fh, normal_fh ]
    for header in headers:
        data = header.readlines()
        for line in data:
            line = line.strip( '\n' )
            fields = line.split( '\t' )
            chrom = fields[ 0 ]
            loc = int( fields[ 1 ] )
            if chrom not in indel_locs:
                indel_locs[ chrom ] = []
                indel_locs[ chrom ].append( loc )
            else:
                indel_locs[ chrom ].append( loc )

    for chrom in indel_locs:
        indel_locs[ chrom ].sort()
    
    # make a chrom dictionary loc list of SNP mutations, the lists should be sorted
    SNP_locs = {}
    for chrom in VCF_dict:
        for loc in VCF_dict[ chrom ]:
            int_loc = int( loc )
            if chrom not in SNP_locs:
                SNP_locs[ chrom ] = []
                SNP_locs[ chrom ].append( int_loc )
            else:
                SNP_locs[ chrom ].append( int_loc )

    for chrom in SNP_locs:
        SNP_locs[ chrom ].sort()

    # compare them to find which SNPs to filter
    filtered_3_indels = {}
    
    for chrom in indel_locs:
        j = 0
        k = 0
        filtered_3_indels[ chrom ] = []

        while j < len( indel_locs[ chrom ] ):
            while k < len( SNP_locs[ chrom ] ) and SNP_locs[ chrom ][ k ] <= indel_locs[ chrom ][ j ] + 10:
                if SNP_locs[ chrom ][ k ] < indel_locs[ chrom ][ j ] - 10:
                    k = k + 1
                    continue
                if indel_locs[ chrom ][ j ] - 10 <= SNP_locs[ chrom ][ k ] <= indel_locs[ chrom ][ j ] + 10:
                    filtered_3_indels[ chrom ].append( SNP_locs[ chrom ][ k ] )
                    k = k + 1
                    continue
            j = j + 1

    # save to file
    new_file = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/MUTATION_CALLS/%s/FILT_LISTS/3_excluded.txt" % VCF_name
    new_fh = open( new_file, 'w' )
    for chrom in filtered_3_indels:
        for loc in filtered_3_indels[ chrom ]:
            new_fh.write( "%s\t%s\n" % ( chrom, loc ) )

    # add column to dictionary VCF_dict
    for chrom in VCF_dict:
        for loc in VCF_dict[ chrom ]:
            if int( loc ) in filtered_3_indels[ chrom ]:
                VCF_dict[ chrom ][ loc ][ "FILT_3" ] = "FAIL"
                # print "got 'em"
            else:
                VCF_dict[ chrom][ loc ][ "FILT_3" ] = "PASS"

    print "filter three has run"

    return VCF_dict

 
