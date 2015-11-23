## This is the R script (run on a Macbook) to create the figures from the data in FIGURE_DATA

###############
## LOAD DATA ##
###############

data2 = read.table( "../FIGURE_DATA/unfiltered_WGS-WGA_counts.txt", header=T )
data34 = read.table( "../FIGURE_DATA/STD40_unfiltered_overlap.txt", header=T )
data5 = read.table( "../FIGURE_DATA/numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt", header=T )
data6 = read.table( "../FIGURE_DATA/overlap_numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt", header=T )
data78a = read.table( "../FIGURE_DATA/STD40_unfiltered_overlap.txt", header=T )
data78a[ 3:14 ] <- list( NULL )
data78b = read.table( "../FIGURE_DATA/overlap_numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt", header=T )
data78b[ 2:7 ] <- list( NULL )
data78b[ 4:5 ] <- list( NULL )
data78b[ 6:7 ] <- list( NULL )

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

##################################################################################
## FIGURE 3 ORIGINAL: Density of percentage overlap between WGS and WGA samples ##
## THIS IS NOW FIGURE S1														##
##################################################################################

#pdf( "../FIGURE_PDFS/FigureS1.pdf" )
# figure title (should it ever be needed): Greater % of WGS sample found in overlap than WGA
amp_frac = density( data34$OVERLAP/data34$TOTAL_WGA.C282., from=0 )
seq_frac = density( data34$OVERLAP/data34$TOTAL_WGS.C484., from=0 )
par(bty='n')
plot( seq_frac$x, seq_frac$y, type='l', col='orange', main="", xlab="Percentage of putative SNVs in a sample also recovered in its technical replicate", ylab="Density")
legend( "topright", c( "WGS", "WGA" ), col=c( "green", "orange" ), lty=1 )
lines( amp_frac$x, amp_frac$y, col='green' )
#dev.off()

#############################################################################################################
## FIGURE 3 MOD (ORIGINAL FIGURE 4): Percentage of overlap versus (log) number of putative SNVs per sample ##
#############################################################################################################

#pdf( "../FIGURE_PDFS/Figure3.pdf" )
# when necessary, use the following title: Number of putative SNPs does not correlate with percentage WGS/WGA overlap
par(bty = 'n')
plot( data34$TOTAL, data34$PERCENT_OVERLAP_by_WGA.282, log='x', pch=20, col="orange", main="", xlab="Number of putative SNVs in unfiltered data (log)", ylab="Percentage overlap between WGS and WGA in unfiltered data" )
points ( data34$TOTAL, data34$PERCENT_OVERLAP_by_WGS.484, pch="+", col="green" )
abline( 0,1 )
legend( 'topleft', c( "overlap/WGS total", "overlap/WGA total" ), col=c( 'green', 'orange' ) )
legend( 'topleft', c( "overlap/WGS total", "overlap/WGA total" ), col=c( 'green', 'orange' ), pch=c(3,20) )
#dev.off()

############################
## FIGURE 3 CORERELATIONS ##
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

############################################################################################
## FIGURE 4 (ORIG FIGURE 5 and 6): Number of putative SNPs removed by each of six filters ##
############################################################################################

# pdf( "../FIGURE_PDFS/Figure4.pdf", width=6, height=5 )
par(bty='n', mfrow=c(2,1), mar=c(5.1,4.1,1.1,0))
boxplot( data5[,c(5,7,4,3,6,8)], log='y', ylim=c(1,100000), xaxt='n', cex.lab=0.7, pch=20, ylab="Number putative SNVs filtered" )
axis( 1, at=1:6, padj=0.8, cex.axis=0.7, labels=c( "VAQ", "LOH", "10bp-SNV", "10bp-INDEL", "dbSNP", "<10%" ) )
boxplot( data6[,c(9,13,7,5,11,15)], ylim=c(0,1), xaxt='n', cex.lab=0.7, ylab="Percentage overlap filtered", pch=20 )
axis( 1, at=1:6, padj = 0.8, cex.axis=0.7, labels=c( "VAQ", "LOH", "10bp-SNV", "10bp-INDEL", "dbSNP", "<10%" ) )
# dev.off()

##############################################################
## FIGURE 5 IS MADE IN THE SCRIPT Jacquard.R IN THIS FOLDER ##
##############################################################

##############################################################################################################################################
## FIGURE 6 (ORIG FIGURES 7 AND 8): Percentage of replicate SNVs filtered out by LOH, as a function of the total number of overlapping SNVs ##
##############################################################################################################################################

LOH_VAQ_data <- merge( data78a, data78b, by="SAMPLE" )
names( LOH_VAQ_data )[ names( LOH_VAQ_data )=="FILT_5_COUNT" ] <- "VAQ_COUNT"
names( LOH_VAQ_data )[ names( LOH_VAQ_data )=="FILT_7_COUNT" ] <- "LOH_COUNT"
names( LOH_VAQ_data )[ names( LOH_VAQ_data )=="FILT_5_PERCENT" ] <- "VAQ_PERCENT"
names( LOH_VAQ_data )[ names( LOH_VAQ_data )=="FILT_7_PERCENT" ] <- "LOH_PERCENT"

# note that this is pretty close to 1.00 every time, hence the inversion...
LOH_VAQ_data$SUM=LOH_VAQ_data$LOH_PERCENT + LOH_VAQ_data$VAQ_PERCENT

#pdf( "../FIGURE_PDFS/Figure6.pdf", width=10, height=5 )
par(bty='n', mfrow=c(1,2), mar=c(5,7,4,2) )
plot( LOH_VAQ_data$OVERLAP, LOH_VAQ_data$LOH_PERCENT, xlab="Number of SNVs in both replicates", ylab="Percentage of SNVs in both replicates\nfiltered out by LOH", pch=20 )
plot( LOH_VAQ_data$OVERLAP, LOH_VAQ_data$VAQ_PERCENT, xlab="Number of SNVs in both replicates", ylab="Percentage of SNVs in both replicates\nfiltered out by VAQ", pch=20 )
#dev.off()


##############
## OLD CODE ##
##############

# make density of each filter removal

#pdf( "../FIGURE_PDFS/Figure4.pdf", width=6, height=5 )
par(bty='n', mar=c(5.1,4.1,1.1,0))
boxplot( data5[,c(5,7,4,3,6,8)], log='y', ylim=c(1,100000), xaxt='n', pch=20, ylab="Number of putative SNVs filtered, per sample" )
axis( 1, at=1:6, padj=0.8, cex.axis=0.7, labels=c( "VAQ", "LOH", "10bp-SNV", "10bp-INDEL", "dbSNP", "<10%" ) )
#dev.off()


#pdf( "../FIGURE_PDFS/Figure5.pdf", width=6, height=5 )
par(bty='n', mar=c(5.1,4.1,1.1,0))
boxplot( data6[,c(9,13,7,5,11,15)], ylim=c(0,1), xaxt='n', ylab="Percentage of overlap filtered, per sample", pch=20 )
axis( 1, at=1:6, padj = 0.8, cex.axis=0.7, labels=c( "VAQ", "LOH", "10bp-SNV", "10bp-INDEL", "dbSNP", "<10%" ) )
#dev.off()



