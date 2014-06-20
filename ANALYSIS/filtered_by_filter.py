import path
import sys

# the purpose of this script is to look at what the individual filters are taking out and make comparisons

GBM_genomics = /share/WilkeLab/work/dzd58/TCGA_Reanalysis/GBM_genomics
data_dir = /share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA

## make sure you're using the correct data
if sys.argv[ 1 ] = "":
	data_level = "MUTATION_CALLS_STD40"
else:
	data_level = sys.argv[ 1 ] 

data = os.path.join( data_dir, data_level )

##



## MAIN FUNCTION ##
##
## which filtered data set to run this on is sys.argv[ 1 ]
## but the default is MUTATION_CALLS_STD40