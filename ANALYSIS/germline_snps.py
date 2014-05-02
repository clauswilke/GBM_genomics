# this script counts the incidence of each of things filtered by filter 6
# answers the question, which, if any, germline SNPs appear in this sample a lot?

import os
import VCF_list

sample_paths = VCF_list.sample_dir_paths()

# print sample_paths

germline_SNP_dict = {}
still_to_analyze = 0
analyzed = 0

for path in sample_paths:
    filter6_file = os.path.join( path, "FILT_LISTS", "6_excluded.txt" )

    if not os.path.exists( filter6_file ):
        still_to_analyze = still_to_analyze + 1
        # print "NOT HERE"
        continue

    # print "HERE"
    analyzed = analyzed + 1
    filt_fh = open( filter6_file, 'r' )
    for line in filt_fh:
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

# print germline_SNP_dict
print still_to_analyze
print analyzed

repeat_germline_snps = 0
for chrom in germline_SNP_dict:
    for loc in germline_SNP_dict[ chrom ]:
        if germline_SNP_dict[ chrom ][ loc ] >1:
            print chrom, loc, germline_SNP_dict[ chrom ][ loc ]
            repeat_germline_snps = repeat_germline_snps + 1

print repeat_germline_snps
