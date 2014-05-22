# load the data
data = read.table( "../FIGURE_DATA/STD40_unfiltered_overlap.txt", header=T )

# density plot of the percent overlap for the 55 WGA-S samples
pdf( "../FIGURE_PDFS/unfiltered_overlap_density.pdf" )
par(bty = 'n')
overlap = density( data$PERCENT_OVERLAP_by_WGA.282, from=0, to=1 )
plot( overlap, main="Density of [WGS/WGA overlap]/[Total WGA mutations]", xlab="[WGS/WGA overlap]/[Total WGA mutations]" )
dev.off()

# density plot of the percent difference for the 55 WGA-S samples
pdf( "../FIGURE_PDFS/unfiltered_difference_density.pdf" )
par(bty='n')
difference = density( data$PERCENT_DIFFERENCE_by_WGA.282, from=0, to=1 )
plot( difference, main="Density of [WGS/WGA difference]/[Total WGA mutations]", xlab="[WGS/WGA difference]/[Total WGA mutations]" )
dev.off()

# plot of the size of the overlap v. total # unfiltered mutations
pdf( "../FIGURE_PDFS/unfiltered_total_muts_v_percent_overlap.pdf" )
par(bty = 'n')
plot( data$TOTAL, data$PERCENT_OVERLAP_by_WGA.282, log='x', pch=20, main="No. of putative SNPs does not correlate with % WGS/WGA overlap", xlab="log of number of putative SNPs (unfiltered)", ylab="percent overlap between WGS and WGA (unfiltered)" )
abline( 0,1 )
dev.off()

# would also like to look at the split of the difference -- how much from 282, how much from 484?