# this script counts the incidence of each of things filtered by filter 6
# answers the question, which, if any, germline SNPs appear in this sample a lot?

import os

mutations_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/MUTATION_CALLS"

sample_list = [ "C484.TCGA-19-4065" ]

germline_SNP_dict = {}
still_to_analyze = 0

for sample in sample_list:
    filter6_file = os.path.join( mutations_dir, sample, "FILT_LISTS", "6_excluded.txt" )

    if not os.path.exists( filter6_file ):
        still_to_analyze = still_to_analyze + 1
    
    6_fh = open( filter6_file, 'r' )
    for line in 6_fh:
        line = line.strip( "\n\r" )
        fields = line.split( "\t" )
        chrom = fields[ 0 ]
        loc = fields[ 1 ]
        if chrom not in germline_SNP_dict:
            germline_SNP_dict[ chrom ] = {}
        if loc not in germline_SNP_dict[ chrom ]:
            germline_SNP_dict[ chrom ][ loc ] = 1
        else:
            germline_SNP_dict[ chrom ][ loc ] = germline_SNP_dict[ chrom ][ loc ] + 1

print germline_SNP_dict
print still_to_analyze


