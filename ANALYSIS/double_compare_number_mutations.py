# this is the script that copares the number of mutations in each of the pairs of doubles files

## this script is used to make the histogram "distribution of frequency of mutations"
## and the "directionalty" figure (plot of condition one against condition2, to see if directionalty hypothesis is supported)

## for now, it is running only on the unfiltered data of the WGS and WGA data set (55 samples)
## but this can be changes using command line arguments later

import os
import VCF_list
import how_many_muts

# get the paths to the relevant directories
sample_dirs, double_count = VCF_list.analyzed_doubles()

# next let's make three lists: (i) the number of mutaitons per file, (ii) the number of mutations per C282 file, and (iii) the number of mutaitons per C484 file
# these are for the histogram
C282_files = []
C484_files = []
all_files = []

# then there will be a list of tuples (C282, C484)
paired_counts = []

# this is the function that propogates these data structures
for sample in sample_dirs:

    C282_dir = sample_dirs[ sample ][ 0 ]
    C282_file = os.path.join( C282_dir, "unfiltered.txt" )
    C282_count = how_many_muts.count_muts( C282_file )

    C484_dir = sample_dirs[ sample ][ 1 ]
    C484_file = os.path.join( C484_dir, "unfiltered.txt" )
    C484_count = how_many_muts.count_muts( C484_file )

    all_files.append( C282_count )
    all_files.append( C484_count )
    C282_files.append( C282_count )
    C484_files.append( C484_count )
    pair = ( sample, C282_count, C484_count )
    paired_counts.append( pair )

# and finally save them in a file /GBM_genomics/FIGURES, where I can make them into figures using R...
figure_data_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/GBM_genomics/FIGURES/FIGURE_DATA"
file = os.path.join( figure_data_dir, "unfiltered_WGS-WGA_counts.txt" )
file_fh = open( file, 'w' )
file_fh.write( "SAMPLE\tC282_count\tC484_count\n" )
for tuple in paired_counts:
    file_fh.write( "%s\t%s\t%s\n" % ( tuple[ 0 ], tuple[ 1 ], tuple[ 2 ] ) )
file_fh.close()
