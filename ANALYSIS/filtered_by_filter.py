import os
import sys
import glob


# the purpose of this script is to look at what the individual filters are taking out and make comparisons

GBM_genomics = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/GBM_genomics"
data_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA"
data_home = os.path.join( data_dir, sys.argv[ 1 ] )
figure_data = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/GBM_genomics/FIGURES/FIGURE_DATA"

## MAIN FUNCTION ##
##
## which filtered data set to run this on is sys.argv[ 1 ]
## but the default is MUTATION_CALLS_STD40

# the directory that will hold all the data is here:
filter_counts = {}

# for every directory in data_home
sampdirs = [ name for name in os.listdir( data_home ) ]

# for each sample
for samp in sampdirs:
	samp_filters = os.path.join( data_home, samp, "FILT_LISTS" )
	filter_counts[ samp ] = []
	# now for each filter:
	for i in range( 2, 9 ):
		file_name = "%s_excluded.txt" % i
		#print samp, file_name (checked for existence in all of MUTATION_CALLS_STD40)
		file = os.path.join( samp_filters, file_name )
		line_count = 0
		file_fh = open( file, 'r' )
		for line in file_fh:
			line_count += 1
		file_fh.close()
		filter_counts[ samp ].append( line_count )
		
#print filter_counts

# now make a file with this information in /share/WilkeLab/work/dzd58/TCGA_Reanalysis/GBM_genomics/FIGURES/FIGURE_DATA
data_file_name =  "numbers_excluded_by_filter_%s.txt" % sys.argv[ 1 ]
data_file = os.path.join( figure_data, data_file_name )
data_fh = open( data_file, 'w' )
data_fh.write( "SAMPLE\tFILTER_2\tFILTER_3\tFILTER_4\tFILTER_5\tFILTER_6\tFILTER_7\tFILTER_8\n" )
for samp in filter_counts:
	data_fh.write( "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ( samp, filter_counts[ samp ][ 0 ], filter_counts[ samp ][ 1 ], filter_counts[ samp ][ 2 ], filter_counts[ samp ][ 3 ], filter_counts[ samp ][ 4 ], filter_counts[ samp ][ 5 ], filter_counts[ samp ][ 6 ] ) )
data_fh.close()
		
		
		