# this is the file that reads dbSNP build 137's VCF and makes it into something that I can use...

import os

data_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA"
dbsnp_vcf = os.path.join( data_dir, "dbsnp_137.b37.vcf" )

dbsnp_dict = {}

dbsnp_fh = open( dbsnp_vcf, 'r' )
for i, line in enumerate( dbsnp_fh ):

    if i == 143:
        header = line.strip( "\r\n" )
        print header

    if i >= 144:
        line = line.strip( "\r\n" )
        fields = line.split( "\t" )
        chrom = fields[ 0 ]
        loc = int( fields[ 1 ] )

        if chrom in dbsnp_dict:
            dbsnp_dict[ chrom ].append( loc )
        else:
            dbsnp_dict[ chrom ] = []

# make dbSNP_dir in DATA/
dbsnp_dir = os.path.join( data_dir, "DBSNP" )
if not os.path.exists( dbsnp_dir ):
    os.makedirs( dbsnp_dir )

# and make a file for each chromosome (by the name of the chromosome.txt) with the locations of all the snps on that chromosome
for chrom in dbsnp_dict:
    dbsnp_dict[ chrom ] = sorted( dbsnp_dict[ chrom ] )
    chrom_file_name = "%s.txt" % chrom
    chrom_file = os.path.join( dbsnp_dir, chrom_file_name )
    chrom_fh = open( chrom_file, 'w' )
    for loc in dbsnp_dict[ chrom ]:
        chrom_fh.write( "%s\n" % str( loc ) )
    
