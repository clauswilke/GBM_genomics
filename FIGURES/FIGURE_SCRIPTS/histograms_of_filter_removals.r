# this script looks in different ways at which mutations are removed when I run the filters

# first look at the STD_40 data
# first, just everything removed, each file individually, how many by each filter 3-8
data = read.table( "../FIGURE_DATA/numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt", header=T )

# make density of each filter removal
d3 = density( data$FILTER_3 )
d4 = density( data$FILTER_4 )
d5 = density( data$FILTER_5 )
d6 = density( data$FILTER_6 )
d7 = density( data$FILTER_7 )
d8 = density( data$FILTER_8 )

# make the box plot that will show everything
par(bty='n')
boxplot( data[,c(5,7,4,3,6,8)], log='y', ylim=c(1,100000), xaxt='n' )
axis( 1, at=1:6, labels=c( "Variant Allele\nQuality (VAQ)\n>20", "Removes\nLOH", "No putative\nSNPs within\n10 bp window\nof other SNPs", "No putative\nSNPs within\n10 bp window\nof indels", "No overlap\nwith dbSNP\ncoverage", "Alternate\nallele\ncoverage\n>=10%") )

### still working on two problems:
### how to make the figure horizontally longer so that all labels fit, and 
### how to move the x-axis lables down relative to the x-axis, again so they fit


#axis( 1, at=1:6 )
#mtext( side=1, text="Variant Allele\nQuality (VAQ)\n>20", line=4, cex=-20 )
#, "Removes\nLOH", "No putative\nSNPs within\n10 bp window\nof other SNPs", "No putative\nSNPs within\n10 bp window\nof indels", "No overlap\nwith dbSNP\ncoverage", "Alternate\nallele\ncoverage\n>=10%"), line=4 )

## this second part of the script makes the same figure, but for only the putative SNPs that are in the overlaps in the unfiltered data