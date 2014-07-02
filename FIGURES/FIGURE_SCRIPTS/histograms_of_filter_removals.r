# this script looks in different ways at which mutations are removed when I run the filters

# first look at the STD_40 data
# first, just everything removed, each file individually, how many by each filter 3-8
data = read.table( "../FIGURE_DATA/unfiltered_WGS-WGA_counts.txt", header=T )