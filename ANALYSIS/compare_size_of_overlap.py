# the purpose of this script is to make a summary of the size of the overlap for each sample,
# the filtered data,
# the unfiltered data,
# the filtered data plus LOH,
# the filtered data plus VAQ,
# and the filtered data plus LOH and VAQ

import os
from sets import Set

# first, relevant data directories
data_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA"
STD_40_dir = os.path.join( data_dir, "MUTATION_CALLS_STD40" )
LOH_dir = os.path.join( data_dir, "STD_40_LOH" )
VAQ_dir = os.path.join( data_dir, "STD_40_VAQ" )
LOH_VAQ_dir = os.path.join( data_dir, "STD_40_LOH_VAQ" )

# to be a sample worth analyzing, 10 files must exist
# this function spits out a list of the samples worth analyzing (in tuples, since each sample will actually be a pair of samples)
def find_samples():
	
	useable_samples = []
	all_samples = Set([])
	
	# first find all the sample names (without the C282 or the C484)
	STD_40_samps = [ name for name in os.listdir( STD_40_dir ) ]
	for sample in STD_40_samps:
		sample_name = sample[ 5: ]
		all_samples.add( sample_name )
	
	# now test for the 10 files, and if they all exist, add them to useable_samples
	for sample in all_samples:
	
		# the unfiltered data
		C282_dict = "C282.%s" % sample
		C484_dict = "C484.%s" % sample
		C282_STD40_unfilt = os.path.exists( os.path.join( STD_40_dir, C282_dict, "unfiltered.txt" ) )
		C484_STD40_unfilt = os.path.exists( os.path.join( STD_40_dir, C484_dict, "unfiltered.txt" ) )
		# print C282_STD40_unfilt, C484_STD40_unfilt
		
		# the filtered data
		C282_STD40_filt = os.path.exists( os.path.join( STD_40_dir, C282_dict, "filtered.txt" ) )
		C484_STD40_filt = os.path.exists( os.path.join( STD_40_dir, C484_dict, "filtered.txt" ) )
		# print C282_STD40_filt, C484_STD40_filt
		
		# the LOH data
		C282_file = "C282.%s.txt" % sample
		C484_file = "C484.%s.txt" % sample
		C282_LOH = os.path.exists( os.path.join( LOH_dir, C282_file ) )
		C484_LOH = os.path.exists( os.path.join( LOH_dir, C484_file ) )
		#print C282_LOH, C484_LOH
		
		# the VAQ data
		C282_VAQ = os.path.exists( os.path.join( VAQ_dir, C282_file ) )
		C484_VAQ = os.path.exists( os.path.join( VAQ_dir, C484_file ) )
		#print C282_VAQ, C484_VAQ
		
		# the LOH/VAQ data
		C282_LOH_VAQ = os.path.exists( os.path.join( LOH_VAQ_dir, C282_file ) )
		C484_LOH_VAQ = os.path.exists( os.path.join( LOH_VAQ_dir, C484_file ) )
		# print C282_LOH_VAQ, C484_LOH_VAQ
		
		# if they're all there, add it to the list...
		if C282_STD40_unfilt == True:
			#print "yes"
			if C484_STD40_unfilt == True:
				if C282_STD40_filt == True:
					if C484_STD40_filt == True:
						if C282_LOH == True:
							if C484_LOH == True:
								if C282_VAQ == True:
									if C484_VAQ == True:
										if C282_LOH_VAQ == True:
											if C484_LOH_VAQ == True:
												#print "FOUND"
												files = [ os.path.join( STD_40_dir, C282_dict, "unfiltered.txt" ), os.path.join( STD_40_dir, C484_dict, "unfiltered.txt" ), os.path.join( STD_40_dir, C282_dict, "filtered.txt" ), os.path.join( STD_40_dir, C484_dict, "filtered.txt" ), os.path.join( LOH_dir, C282_file ), os.path.join( LOH_dir, C484_file ), os.path.join( VAQ_dir, C282_file ), os.path.join( VAQ_dir, C484_file ), os.path.join( LOH_VAQ_dir, C282_file ), os.path.join( LOH_VAQ_dir, C484_file ) ]
												useable_samples.append( files )
	return useable_samples
	
def find_overlap_length( C282_file, C484_file ):

	# make a set of the C282_file
	C282_muts = Set([])
	C282_fh = open( C282_file, 'r' )
	for line in C282_fh:
		line = line.strip( '\r\n' )
		fields = line.split( '\t' )
		tuple = ( fields[ 0 ], fields[ 1 ] )
		C282_muts.add( tuple )
	
	# make a set of the C484_file
	C484_muts = Set([])
	C484_fh = open( C484_file, 'r' )
	for line in C484_fh:
		line = line.strip( '\r\n' )
		fields = line.split( '\t' )
		tuple = ( fields[ 0 ], fields[ 1 ] )
		C484_muts.add( tuple )
	
	# and find the overlap between the two sets
	overlap = ( C282_muts & C484_muts )
	overlap_length = len( overlap )

	return overlap_length



###########################
###### MAIN FUNCTION ######
###########################

samples = find_samples()

# for each sample
for sample in samples:

	# find the lengths of the overlaps between them
	unfilt_overlap = find_overlap_length( sample[ 0 ], sample[ 1 ] )
	#print unfilt_overlap
	filt_overlap = find_overlap_length( sample[ 2 ], sample[ 3 ] )
	#print filt_overlap
	LOH_overlap = find_overlap_length( sample[ 4 ], sample[ 5 ] )
	#print LOH_overlap
	VAQ_overlap = find_overlap_length( sample[ 6 ], sample[ 7 ] )
	#print VAQ_overlap
	LOH_VAQ_overlap = find_overlap_length( sample[ 8 ], sample[ 9 ] )
	print filt_overlap, VAQ_overlap, LOH_overlap, LOH_VAQ_overlap, unfilt_overlap
	


## THERE IS STILL A BUG IN THIS, BECAUSE THE LOH AND THE VAQ_LOH ARE COMING OUT IDENTICALLY, WHICH THEY ARE NOT ALLOWED TO DO BECAUSE THAT IS NOT RIGHT!!!!!!!!!!!!!!!!!!