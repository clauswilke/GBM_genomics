# this scripts counts the number of "mutations" (locations?) in a file
# it is a single function meant to be used by other functions
# functionally, it counts the lines in a file and nothing else

def count_muts( file ):

    count = 0

    file_fh = open( file, 'r' )
    for line in file_fh:
        count +=1

    return count

##
## TESTING
##

# file = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/MUTATION_CALLS/C282.TCGA-19-2624/unfiltered.txt"
# count = count_muts( file )
# print count
