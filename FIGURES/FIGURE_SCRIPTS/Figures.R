## This is the R script (run on a Macbook) to create the figures from the data in FIGURE_DATA

###############
## LOAD DATA ##
###############

data23 = read.table( "../FIGURE_DATA/unfiltered_WGS-WGA_counts.txt", header=T )

#####################################################
## FIGURE 2: Number of putative SNVs in WGS v. WGA ##
#####################################################

# pdf( "../FIGURE_PDFS/C282_v_C484.pdf" )
par(bty='n')
plot(data23$C282_count, data23$C484_count, log='xy', pch=20, xlab="No. of putative SNVs, WGA", ylab="No. of putative SNVs, WGS" )
abline(0,1)
# dev.off()