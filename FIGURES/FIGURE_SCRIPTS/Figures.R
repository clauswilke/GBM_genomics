## This is the R script (run on a Macbook) to create the figures from the data in FIGURE_DATA

###############
## LOAD DATA ##
###############

data2 = read.table( "../FIGURE_DATA/unfiltered_WGS-WGA_counts.txt", header=T )
data34 = read.table( "../FIGURE_DATA/STD40_unfiltered_overlap.txt", header=T )

#########################################################
## FIGURE 2: Number of putative SNVs in WGS versus WGA ##
#########################################################

# pdf( "../FIGURE_PDFS/Figure2.pdf" )
par(bty='n')
plot(data2$C282_count, data2$C484_count, log='xy', pch=20, xlab="Number of putative SNVs, WGA", ylab="Number of putative SNVs, WGS" )
abline(0,1)
# dev.off()

############################
## FIGURE 2 CORERELATIONS ##
############################

# First the Pearson correlation on the raw data
# cor.test(data2$C282_count, data2$C484_count)
# OUTPUT: t = 4.2622, df = 53, p-value = 8.35e-05, 95 percent confidence interval: 0.2770757 0.6794588, cor=0.5052344

# Next the Pearson correlation on the log-transformed data
# cor.test( log(data2$C282_count), log(data2$C484_count) )
# OUTPUT: t = 3.6456, df = 53, p-value = 0.0006086, 95 percent confidence interval: 0.2070558 0.6373449, cor=0.4477571

# Finally the Spearman correlation
# cor.test(data2$C282_count, data2$C484_count, method="spearman")
# OUTPUT: S = 16142.79, p-value = 0.001511, rho=0.4176482

#########################################################################
## FIGURE 3: Density of percentage overlap between WGS and WGA samples ##
#########################################################################

# pdf( "../FIGURE_PDFS/Figure3.pdf" )
# figure title (should it ever be needed): Greater % of WGS sample found in overlap than WGA
amp_frac = density( data34$OVERLAP/data34$TOTAL_WGA.C282., from=0 )
seq_frac = density( data34$OVERLAP/data34$TOTAL_WGS.C484., from=0 )
par(bty='n')
plot( seq_frac$x, seq_frac$y, type='l', col='blue', main="", xlab="Percentage of putative SNVs in a sample also recovered in the samples's technical replicate", ylab="Density")
legend( "topright", c( "WGS", "WGA" ), col=c( "blue", "red" ), lty=1 )
lines( amp_frac$x, amp_frac$y, col='red' )
# dev.off()

#####################################################################################
## FIGURE 4: Percentage of overlap versus (log) number of putative SNVs per sample ##
#####################################################################################

pdf( "../FIGURE_PDFS/Figure4.pdf" )
# when necessary, use the following title: Number of putative SNPs does not correlate with percentage WGS/WGA overlap
par(bty = 'n')
plot( data34$TOTAL, data34$PERCENT_OVERLAP_by_WGA.282, log='x', pch=20, main="", xlab="Number of putative SNVs in unfiltered data (log)", ylab="Percentage overlap between WGS and WGA in unfiltered data" )
abline( 0,1 )
dev.off()

############################
## FIGURE 4 CORERELATIONS ##
############################

# First the Pearson correlation on the raw data
# cor.test(data34$TOTAL, data34$PERCENT_OVERLAP_by_WGA.282)
# OUTPUT: t = -0.2989, df = 53, p-value = 0.7662, 95 percent confidence interval: -0.3030232  0.2267410, cor=-0.04102397

# Next the Pearson correlation on the log-transformed data
# cor.test( log(data34$TOTAL), data34$PERCENT_OVERLAP_by_WGA.282 )
# OUTPUT: t = -0.043, df = 53, p-value = 0.9659, 95 percent confidence interval: -0.2707813  0.2597958, cor=-0.005908626

# Finally the Spearman correlation
# cor.test(data34$TOTAL, data34$PERCENT_OVERLAP_by_WGA.282, method="spearman")
# OUTPUT: S = 29268, p-value = 0.6847, rho=-0.05584416

## FIGURE 5: ##

## FIGURE 6: ##

## FIGURE 7: ##

## FIGURE 8: ##








