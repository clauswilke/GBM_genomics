# the purpose of this script is to calculate the overlap by chromosome of a sample
# it is possible that I am wasting my time... but if one chromosome is bigger than the rest in terms ov NEVER OVERLAPPING that would be cool

## COMMAND LINE ARGUMENTS ##
## sys.argv[1] = filtered or not; it is entered as "unfiltered.txt" or "filtered.txt"

import sys
import os

all_overlaps = {}

# first load the data...
overlap_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/DOUBLES/WGA-S/MUTATION_CALLS_STD40"
sample_dirs = [ x[ 0 ] for x in os.walk( overlap_dir ) ]
for dir in sample_dirs:
	overlap_file = os.path.join( dir, "overlap_" + sys.argv[ 1 ] )
	
	# make sure the file exists...
	if not os.path.exists( overlap_file ):
		print "This one, real is not: " + overlap_file
		continue
	else:
		#print "GOOD ONE"
		overlap_fh = open( overlap_file, 'r' )
		for line in overlap_fh:
			line = line.strip( '\r\n' )
			fields = line.split( '\t' )
			chrom = fields[ 0 ]
			loc = fields[ 1 ]
			if chrom not in all_overlaps:
				all_overlaps[ chrom ] = [ loc ]
			else:
				all_overlaps[ chrom ].append( loc )
				
# now start to do some calculations...
chroms = []
counts = []
for chrom in all_overlaps:
	chroms.append( chrom )
	counts.append( len( all_overlaps[ chrom ] ) )

print chroms
print counts

## for a rough estimate, I do still need to normalize this by chromosome length...
