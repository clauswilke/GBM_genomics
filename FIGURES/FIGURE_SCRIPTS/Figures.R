## This is the R script (run on a Macbook) to create the figures from the data in FIGURE_DATA

###############
## LOAD DATA ##
###############

data23 = read.table( "../FIGURE_DATA/unfiltered_WGS-WGA_counts.txt", header=T )

#####################################################
## FIGURE 2: Number of putative SNVs in WGS v. WGA ##
#####################################################

pdf( "../FIGURE_PDFS/Figure2.pdf" )
par(bty='n')
plot(data23$C282_count, data23$C484_count, log='xy', pch=20, xlab="Number of putative SNVs, WGA", ylab="Number of putative SNVs, WGS" )
abline(0,1)
dev.off()

############################
## FIGURE 2 CORERELATIONS ##
############################

# First the Pearson correlation on the raw data
cor.test(data23$C282_count, data23$C484_count)
# OUTPUT: t = 4.2622, df = 53, p-value = 8.35e-05, 95 percent confidence interval: 0.2770757 0.6794588, cor=0.5052344

# Next the Pearson correlation on the log-transformed data
cor.test( log(data23$C282_count), log(data23$C484_count) )
# OUTPUT: t = 3.6456, df = 53, p-value = 0.0006086, 95 percent confidence interval: 0.2070558 0.6373449, cor=0.4477571

# Finally the Spearman correlation
cor.test(data23$C282_count, data23$C484_count, method="spearman")
# OUTPUT: S = 16142.79, p-value = 0.001511, rho=0.4176482

## FIGURE 3: ##