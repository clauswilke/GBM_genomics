# the purpose of this script is to make data sets of filtered data (STD_40) that also includes the LOH, the VAQ, and the LOH/VAQ, to see if there is better agreement with or whithout these filters

import os

## these are the directories where data lives
data_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA"
STD_40_dir = os.path.join( data_dir, "MUTATION_CALLS_STD40" )
LOH_dir = os.path.join( data_dir, "STD_40_LOH" )
VAQ_dir = os.path.join( data_dir, "STD_40_VAQ" )
LOH_VAQ_dir = os.path.join( data_dir, "STD_40_LOH_VAQ" )

def read_file( filepath ):

	mutation_list = []
	
	file_fh = open( filepath, 'r' )
	for line in file_fh:
		mutation_list.append( line )
	
	return mutation_list

def write_newfile( newpath, inputs ):

	newfile_fh = open( newpath, 'w' )
	for input in inputs:
		for line in input:
			newfile_fh.write( line )
			
	newfile_fh.close()

	return

#####################
### MAIN FUNCTION ###
#####################

samples = [ name for name in os.listdir( STD_40_dir ) ]

## for each sample we are doing this for, make the appropriate paths...
for sample in samples:
	samp_dir = os.path.join( STD_40_dir, sample )
	filter_dir = os.path.join( samp_dir, "FILT_LISTS" )
	filtered_filepath = os.path.join( samp_dir, "filtered.txt" )
	LOH_filepath = os.path.join( filter_dir, "7_excluded.txt" )
	VAQ_filepath = os.path.join( filter_dir, "5_excluded.txt" )
	
	# now make a list of each of the three file_contents
	filtered_muts = read_file( filtered_filepath )
	LOH_muts = read_file( LOH_filepath )
	VAQ_muts = read_file( VAQ_filepath )
	#print VAQ_muts
	
	# and save the appropriate lists in the right new spots
	sample_file = "%s.txt" % sample
	
	LOH_newpath = os.path.join( LOH_dir, sample_file )
	LOH_inputs = [ filtered_muts, LOH_muts ]
	#write_newfile( LOH_newpath, LOH_inputs, LOH_muts, VAQ_muts, filtered_muts )
	write_newfile( LOH_newpath, LOH_inputs )
	
	VAQ_newpath = os.path.join( VAQ_dir, sample_file )
	VAQ_inputs = [ filtered_muts, VAQ_muts ]
	write_newfile( VAQ_newpath, VAQ_inputs )
	
	LOH_VAQ_newpath = os.path.join( LOH_VAQ_dir, sample_file )
	LOH_VAQ_inputs = [ filtered_muts, LOH_muts ]
	write_newfile( LOH_VAQ_newpath, LOH_VAQ_inputs )