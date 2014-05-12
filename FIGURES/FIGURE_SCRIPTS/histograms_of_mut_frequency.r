setwd( "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/GBM_genomics/FIGURES/" )

# this script makes histograms of mutation frequency for the different data sets

# first load the data
data = read.table( "./FIGURE_DATA/unfiltered_WGS-WGA_counts.txt", header=T )


# histogram of C282
pdf( "./FIGURE_PDFS/C282_count_histogram.pdf" )
hist( data$C282_count )
dev.off()

# histogram of C484
pdf( "./FIGURE_PDFS/C484_count_histogram.pdf" )
hist( data$C484_count )
dev.off

# plot of C282 against C484


