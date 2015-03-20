## This is for figure 4 and 5 in the paper
## that is the density plot (4)
## and the log figure that is random

# load the data
data = read.table( "../FIGURE_DATA/STD40_unfiltered_overlap.txt", header=T )

# density plot of the percent overlap for the 55 WGA-S samples
png( "../FIGURE_PDFS/unfiltered_overlap_density.png" )
# pdf( "../FIGURE_PDFS/unfiltered_overlap_density.pdf" )
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
#png( "../FIGURE_PDFS/unfiltered_total_muts_v_percent_overlap.png" )
pdf( "../FIGURE_PDFS/unfiltered_total_muts_v_percent_overlap.pdf" )
# when necessary, use the following title: No. of putative SNPs does not correlate with % WGS/WGA overlap
par(bty = 'n')
plot( data$TOTAL, data$PERCENT_OVERLAP_by_WGA.282, log='x', pch=20, main="", xlab="log of no. of putative SNPs (unfiltered)", ylab="percent overlap between WGS and WGA (unfiltered)" )
abline( 0,1 )
dev.off()

# for each of WGA(282) and WGS(484), is there a clear difference in percentage of mutations that are in the overlap?
# in other words, what is teh density of (a_only/a_all) - (b_only/b_all)?
#png( "../FIGURE_PDFS/unfiltered_overlap_WGS_WGA_together_densities.png" )
pdf( "../FIGURE_PDFS/unfiltered_overlap_WGS_WGA_together_densities.pdf" )
# figure title (should it ever be needed): Greater % of WGS sample found in overlap than WGA
amp_frac = density( data$OVERLAP/data$TOTAL_WGA.C282., from=0 )
seq_frac = density( data$OVERLAP/data$TOTAL_WGS.C484., from=0 )
par(bty='n')
plot( seq_frac$x, seq_frac$y, type='l', col='blue', main="", xlab="% of sample obtained in overlap", ylab="Density")
legend( "topright", c( "WGS", "WGA" ), col=c( "blue", "red" ), lty=1 )
lines( amp_frac$x, amp_frac$y, col='red' )
dev.off()

