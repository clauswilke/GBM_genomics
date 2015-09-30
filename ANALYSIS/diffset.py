## this is a new script
import os

## in GBM_Genomics live ANALYSIS (home of this script) and FIGURES
filtData="../FIGURES/FIGURE_DATA/numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt"
overData="../FIGURES/FIGURE_DATA/overlap_numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt"

## set up a dictionary, smapDict, with sample: { filter: { total, overlap, diffset } }
## also set up a dictionary, filtDict, with filter: { total:[], overlap:[], diffset:[] }
sampDict={}

filtDict={}
filtDict[ "F2" ]={}
filtDict[ "F3" ]={}
filtDict[ "F4" ]={}
filtDict[ "F5" ]={}
filtDict[ "F6" ]={}
filtDict[ "F7" ]={}
filtDict[ "F8" ]={}

for filter in filtDict:
	filtDict[ filter ][ "total" ] = []
	filtDict[ filter ][ "overlap" ] = []
	filtDict[ filter ][ "diffset" ] = []
	## some calculated properties in the dictionary...
	filtDict[ filter ][ "percentOverlap" ]=[]
	filtDict[ filter ][ "percentDiffset" ]=[]
	filtDict[ filter ][ "ratioDiff2Over" ]=[]
	filtDict[ filter ][ "ratioOver2Diff" ]=[]
	
## NUMBER IN OVERLAP
overFh=open( overData, 'r' )
for line in overFh:
	line=line.strip( '\n' )
	fields=line.split( '\t' )
	sample=fields[ 0 ]
	samp1="C484."+sample
	samp2="C282."+sample
	Ftwo=fields[ 1 ]
	Fthree=fields[ 3 ]
	Ffour=fields[ 5 ]
	Ffive=fields[ 7 ]
	Fsix=fields[ 9 ]
	Fseven=fields[ 11 ]
	Feight=fields[ 13 ]
	
	if sample=="SAMPLE":
		continue
	
	# samp1
	sampDict[ samp1 ]={}
	sampDict[ samp1 ][ "F2" ]={}
	sampDict[ samp1 ][ "F2" ][ "overlap" ]=Ftwo
		
	sampDict[ samp1 ][ "F3" ]={}
	sampDict[ samp1 ][ "F3" ][ "overlap" ]=Fthree

	sampDict[ samp1 ][ "F4" ]={}
	sampDict[ samp1 ][ "F4" ][ "overlap" ]=Ffour
			
	sampDict[ samp1 ][ "F5" ]={}
	sampDict[ samp1 ][ "F5" ][ "overlap" ]=Ffive
		
	sampDict[ samp1 ][ "F6" ]={}
	sampDict[ samp1 ][ "F6" ][ "overlap" ]=Fsix
		
	sampDict[ samp1 ][ "F7" ]={}
	sampDict[ samp1 ][ "F7" ][ "overlap" ]=Fseven
		
	sampDict[ samp1 ][ "F8" ]={}
	sampDict[ samp1 ][ "F8" ][ "overlap" ]=Feight		

	# samp2		
	sampDict[ samp2 ]={}
	sampDict[ samp2 ][ "F2" ]={}
	sampDict[ samp2 ][ "F2" ][ "overlap" ]=Ftwo
		
	sampDict[ samp2 ][ "F3" ]={}
	sampDict[ samp2 ][ "F3" ][ "overlap" ]=Fthree

	sampDict[ samp2 ][ "F4" ]={}
	sampDict[ samp2 ][ "F4" ][ "overlap" ]=Ffour
		
	sampDict[ samp2 ][ "F5" ]={}
	sampDict[ samp2 ][ "F5" ][ "overlap" ]=Ffive

	sampDict[ samp2 ][ "F6" ]={}
	sampDict[ samp2 ][ "F6" ][ "overlap" ]=Fsix
		
	sampDict[ samp2 ][ "F7" ]={}
	sampDict[ samp2 ][ "F7" ][ "overlap" ]=Fseven

	sampDict[ samp2 ][ "F8" ]={}
	sampDict[ samp2 ][ "F8" ][ "overlap" ]=Feight
	
overFh.close()
		
## FILTER TOTALS
filtFh=open( filtData, 'r' )
for line in filtFh:
	line=line.strip( '\n' )
	fields=line.split( '\t' )
	sample=fields[ 0 ]
	Ftwo=fields[ 1 ]
	Fthree=fields[ 2 ]
	Ffour=fields[ 3 ]
	Ffive=fields[ 4 ]
	Fsix=fields[ 5 ]
	Fseven=fields[ 6 ]
	Feight=fields[ 7 ]

	if sample in sampDict.keys():
		sampDict[ sample ][ "F2" ][ "total" ]=Ftwo
		sampDict[ sample ][ "F3" ][ "total" ]=Fthree
		sampDict[ sample ][ "F4" ][ "total" ]=Ffour
		sampDict[ sample ][ "F5" ][ "total" ]=Ffive
		sampDict[ sample ][ "F6" ][ "total" ]=Fsix
		sampDict[ sample ][ "F7" ][ "total" ]=Fseven
		sampDict[ sample ][ "F8" ][ "total" ]=Feight
	else:
		continue

filtFh.close()

## NUMBER IN DIFFSET
for sample in sampDict:
	for filter in sampDict[ sample ]:
		sampDict[ sample ][ filter ][ "diffset" ] = int( sampDict[ sample ][ filter ][ "total" ] ) - int( sampDict[ sample ][ filter ][ "overlap" ] )
		sampDict[ sample ][ filter ][ "total" ]=int( sampDict[ sample ][ filter ][ "total" ] )
		sampDict[ sample ][ filter ][ "overlap" ]=int( sampDict[ sample ][ filter ][ "overlap" ] )

## FILT_DICT
for sample in sampDict:
	for filter in sampDict[ sample ]:
		tot = sampDict[ sample ][ filter ][ "total" ]
		over = sampDict[ sample ][ filter ][ "overlap" ]
		diff = sampDict[ sample ][ filter ][ "diffset" ]
		
		if filter in filtDict.keys():
			filtDict[ filter ][ "total" ].append( float( tot ) )
			filtDict[ filter ][ "overlap" ].append( float( over ) )
			filtDict[ filter ][ "diffset" ].append( float( diff ) )
		else:
			print "MISSED", filter
			
## MAKE SOME MORE QUANTITIES IN FILT DICT
for filter in filtDict:
# check the lists are all the same length
#	if len( filtDict[ filter ][ "total" ] ) == len( filtDict[ filter ][ "overlap" ] ) == len( filtDict[ filter ][ "diffset" ] ):
#		print "YES"
#	else:
#		print "NO"
# they are.
	samLength = len( filtDict[ filter ][ "total" ] )
	i=0
	while i<samLength:
		## the quantities we'd like to calculate
		if filtDict[ filter ][ "total" ][ i ] != 0:
			OVER = filtDict[ filter ][ "overlap" ][ i ]/filtDict[ filter ][ "total" ][ i ]
		else:
			OVER = "INF"
		filtDict[ filter ][ "percentOverlap" ].append( OVER )
		
		if filtDict[ filter ][ "total" ][ i ] != 0:
			DIFF = filtDict[ filter ][ "diffset" ][ i ]/filtDict[ filter ][ "total" ][ i ]
		else:
			DIFF = "INF"
		filtDict[ filter ][ "percentDiffset" ].append( DIFF )
		
		if OVER != "INF" and OVER != 0:
			DIFFRATE = filtDict[ filter ][ "percentDiffset" ][ i ]/filtDict[ filter ][ "percentOverlap" ][ i ]
		elif OVER == 0:
			DIFFRATE = "INF"
		elif OVER == "INF":
			DIFFRATE = 0.00
		filtDict[ filter ][ "ratioDiff2Over" ].append( DIFFRATE )

		if DIFF != "INF" and DIFF != 0:		
			OVERRATE = filtDict[ filter ][ "percentOverlap" ][ i ]/filtDict[ filter ][ "percentDiffset" ][ i ]
		elif DIFF == 0:
			OVERRATE = "INF"
		elif DIFF == "INF":
			OVERRATE = 0.00
		filtDict[ filter ][ "ratioOver2Diff" ].append( OVERRATE )

		print i
		print filtDict[ filter ][ "ratioDiff2Over" ]

		i = i+1
	
## print some stuff for my own edification
## nm this will be easier in R...
# print "filt_7: " + 
# print "filt_5: " + 
	
## save all the info in FiltDict into FIGURE_DATA/DIFFSET.TXT, for purposes of R analysis
#diffsetData="../FIGURES/FIGURE_DATA/diffset.txt"
#diffFh.open( diffsetData, 'w' )
#header = "SAMPLE\tF2-tot\tF2-over\tF2-diff\tF3-tot\tF3-over\tF3-diff\tF4-tot\tF4-over\tF4\t\t\t\t\t\t\t\t\t\t\t\t\t"
#diffFh.write( header )

#samLength = len( filtDict[ filter ][ "total" ] )
#i=0
#while i<samLength:
	
	
	
#	i=i+1






