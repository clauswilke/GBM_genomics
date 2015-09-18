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
filtDict[ "F2" ][ "total" ]=[]
filtDict[ "F2" ][ "overlap" ]=[]
filtDict[ "F2" ][ "diffset" ]=[]
filtDict[ "F3" ]={}
filtDict[ "F3" ][ "total" ]=[]
filtDict[ "F3" ][ "overlap" ]=[]
filtDict[ "F3" ][ "diffset" ]=[]
filtDict[ "F4" ]={}
filtDict[ "F4" ][ "total" ]=[]
filtDict[ "F4" ][ "overlap" ]=[]
filtDict[ "F4" ][ "diffset" ]=[]
filtDict[ "F5" ]={}
filtDict[ "F5" ][ "total" ]=[]
filtDict[ "F5" ][ "overlap" ]=[]
filtDict[ "F5" ][ "diffset" ]=[]
filtDict[ "F6" ]={}
filtDict[ "F6" ][ "total" ]=[]
filtDict[ "F6" ][ "overlap" ]=[]
filtDict[ "F6" ][ "diffset" ]=[]
filtDict[ "F7" ]={}
filtDict[ "F7" ][ "total" ]=[]
filtDict[ "F7" ][ "overlap" ]=[]
filtDict[ "F7" ][ "diffset" ]=[]
filtDict[ "F8" ]={}
filtDict[ "F8" ][ "total" ]=[]
filtDict[ "F8" ][ "overlap" ]=[]
filtDict[ "F8" ][ "diffset" ]=[]

## FILTER TOTALS
filtFh=open( filtData, 'r' )
for line in filtFh:
	line=line.strip( '\n' )
	fileds=line.split( '\t' )
	sample=fileds[ 0 ]
	Ftwo=fileds[ 1 ]
	Fthree=fileds[ 2 ]
	Ffour=fileds[ 3 ]
	Ffive=fileds[ 4 ]
	Fsix=fileds[ 5 ]
	Fseven=fileds[ 6 ]
	Feight=fileds[ 7 ]

	sampDict[ sample ]={}
	sampDict[ sample ][ "F2" ]={}
	sampDict[ sample ][ "F2" ][ "total" ]=Ftwo
	
	sampDict[ sample ][ "F3" ]={}
	sampDict[ sample ][ "F3" ][ "total" ]=Fthree
	
	sampDict[ sample ][ "F4" ]={}
	sampDict[ sample ][ "F4" ][ "total" ]=Ffour
	
	sampDict[ sample ][ "F5" ]={}
	sampDict[ sample ][ "F5" ][ "total" ]=Ffive
	
	sampDict[ sample ][ "F6" ]={}
	sampDict[ sample ][ "F6" ][ "total" ]=Fsix
	
	sampDict[ sample ][ "F7" ]={}
	sampDict[ sample ][ "F7" ][ "total" ]=Fseven
	
	sampDict[ sample ][ "F8" ]={}
	sampDict[ sample ][ "F8" ][ "total" ]=Feight
	
## NUMBER IN OVERLAP VERSUS DIFFSET


## OVERLAP PERCENTAGES