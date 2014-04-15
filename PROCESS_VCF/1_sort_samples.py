## 

import os
import string
from subprocess import call
import shutil
from sets import Set

## this is where unencrypted files live
path = "/share/WilkeLab/work/MattCC/TCGA_NextGen_Data/SomaticSniper_Data_FINAL/"
path_dir = os.listdir( path )
# make a list of all unencrypted sample files in that directory
all_sample_files = []
for file in path_dir:
    if ".vcf" in file:
        if ".gpg" not in file:
            all_sample_files.append( file )

# make a dictionary with all the relevant information for each sample name
sample_files = {}
for file in all_sample_files:
    sample_name = file[ 5:-7 ]
    if sample_name not in sample_files:
        sample_files[ sample_name ] = []
        sample_files[ sample_name ].append( file )
    else:
        sample_files[ sample_name ].append( file )

# pair each sample with a TCGA barcode, in a dictionary barcodes
barcodes = {}
barcode_file = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/TCGA_barcodes.txt"
barcode_fh = open( barcode_file, 'r' )
for line in  barcode_fh.readlines():
    line = line.strip( '\r\n' )
    fields = line.split( '/' )
    file_name = fields[ 0 ]
    barcode = fields[ 1 ]
    barcodes[ file_name ] =  barcode

# write the double samples file
doubles_file = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/doubles.txt"
doubles_fh = open( doubles_file, 'w' )
for sample in sample_files:
    if len( sample_files[ sample ] ) != 1:
        # print sample_files[ sample ]
        doubles_fh.write( "%s\t%s\t%s\n" % ( sample, barcodes[ sample_files[ sample ][ 0 ][ :-7 ] ], barcodes[ sample_files[ sample ][ 1 ][ :-7 ] ] ) )
        # print sample, sample_files[ sample ]
doubles_fh.close()

# also write a file name with all the sample_name, file combinations, TCGA barcode
# call this file VCFs
VCFs_file = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/VCFs.txt"
VCFs_fh = open( VCFs_file, 'w' )
for sample in sample_files:
    for file in sample_files[ sample ]:
        VCFs_fh.write( "%s\t%s\t%s\n" % ( sample, file, barcodes[ file[ :-7 ] ] ) )
VCFs_fh.close()

#for file in all_sample_files:
    #print file
