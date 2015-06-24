## This is for figure 3 in the paper
## and hte plot with the trendline (2)

# setwd( "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/GBM_genomics/FIGURES/" )

# this script makes histograms of mutation frequency for the different data sets

# first load the data
data = read.table( "../FIGURE_DATA/unfiltered_WGS-WGA_counts.txt", header=T )

# density plot of C282
pdf( "../FIGURE_PDFS/unfiltered_C282_density.pdf" )
d1 = density( data$C282_count, from=0, to=10000 )
plot(d1)
dev.off()

# density plot of C484
pdf( "../FIGURE_PDFS/unfiltered_C484_density.pdf" )
d2 = density( data$C484_count, from=0, to=10000 )
plot(d2)
dev.off()

# plot of (paired) C282 against C484, with a line at 1 to see how things fall
# anything below the line is in line with the hypothesis, above is bad
#png( "../FIGURE_PDFS/C282_v_C484.png" )
pdf( "../FIGURE_PDFS/C282_v_C484.pdf" )
par(bty='n')
plot(data$C282_count, data$C484_count, log='xy', pch=20, xlab="No. of putative SNVs, WGA", ylab="No. of putative SNVs, WGS" )
abline(0,1)
dev.off()

# the two density plots on top of each other, to see what is going on
pdf( "../FIGURE_PDFS/unfiltered_C282_C484_density.pdf" )
par(bty='n')
plot( d2$x, d2$y, type='l', col='red', xlab="Density", ylab="% of sample contained in overlap")
lines( d1$x, d1$y, col='blue' )
dev.off()

# and finally, the differene: how much bigger is the one than the other?
pdf( "../FIGURE_PDFS/unfiltered_C282-C484.pdf" ) 
d3 = density( data$C282_count-data$C484_count )
plot( d3$x, d3$y, type='l', col='blue')
dev.off()





