## runs on stampede
## extracts information about how mapping quality
## outputs table in form <sample ID> <patient ID> <%mapped> <%covered> <avg depth> <avg depth of covered>
## outputs this data in /home1/01839/dakotaz/GBM_genomics/

## environment
import os

## set up, variables, structures
print "now I will analyze the sequences..."
dataDict={}

## commands
dataDir="/scratch/01839/dakotaz/alignments/2015-05.TCGA-GBM-53techreps/"
# for type in dataDir...
types = os.listdir( dataDir )
# print types
for type in types:
	# for sample in type (not whole genome)
	typeDir=os.path.join( dataDir, type )
	samples = os.listdir( typeDir )
	# print type, samples
	for sample in samples:
		# print type, sample
		sampDir=os.path.join( typeDir, sample )
		depth=sample+".realn.recal.depth.txt"
		idx=sample+".realn.recal.idxstats.txt"
		depthFile=os.path.join( sampDir, depth )
		idxFile=os.path.join( sampDir, idx )
		if os.path.isfile( depthFile ) and os.path.isfile( idxFile ):
			# print sample, " is a real file"
			# put the sample names into the dictionary
			
			# read idxstats file to get %mapped
			depthFh=open( depthFile, 'r' )
			linenumber=1
			coverage="zero"
			for line in depthFh:
				if linenumber==10:
					coverage=line.strip('\n')
					coverage=coverage.strip( 'Average: = ' )
				linenumber=linenumber+1
			if linenumber==16:
			
				dataDict[ sample ]={}
				dataDict[ sample ][ "type" ] = type
				patientID=sample[5:17]
				# print patientID
				dataDict[ sample ][ "patientID" ]= patientID
			
				dataDict[ sample ][ "cover" ]=coverage
				# print sample, coverag3
				# read idxstats file to get mapped/unmapped and calculate percent
				idxFh=open( idxFile, 'r' )
				mapped=0.00
				unmapped=0.00
				for line in idxFh:
					line = line.strip( '\n' )
					line = line.split( '\t' )
					mapped=mapped+float( line[ 2 ] )
					unmapped=unmapped+float( line[ 3 ] )
					perm=mapped/(mapped+unmapped)
				dataDict[ sample ][ "percentMapped" ] = perm
					
# now put all that stuff in a file you can read into R
outFile="/home1/01839/dakotaz/GBM_genomics/FIGURES/FIGURE_DATA/qualstats.txt"
outFh=open( outFile, 'w' )
outFh.write( 'SAMPLE\tPATIENT_ID\tCOVERAGE\tPERCENT_MAPPED\n' )
for sample in dataDict:
	outFh.write( sample )
	outFh.write( '\t' )
	outFh.write( dataDict[ sample ][ "patientID" ] )
	outFh.write( '\t' )
	outFh.write( dataDict[ sample ][ "cover" ] )
	outFh.write( '\t' )
	outFh.write( str( dataDict[ sample ][ "percentMapped" ] ) )
	outFh.write( '\n' )
	print sample, dataDict[ sample ][ "patientID" ], dataDict[ sample ][ "cover" ], dataDict[ sample ][ "percentMapped" ]
outFh.close() 
