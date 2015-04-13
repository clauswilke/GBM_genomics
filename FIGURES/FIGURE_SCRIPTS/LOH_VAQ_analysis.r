# we are redoing the LOH/VAQ analysis, because the figures are screwed...

#######################################
# STEP 1: make the correct data frame #
#######################################
 
overlap_data = read.table( "../FIGURE_DATA/STD40_unfiltered_overlap.txt", header=T )
overlap_data[ 3:14 ] <- list( NULL )

filter_data = read.table( "../FIGURE_DATA/overlap_numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt", header=T )
filter_data[ 2:7 ] <- list( NULL )
filter_data[ 4:5 ] <- list( NULL )
filter_data[ 6:7 ] <- list( NULL )

LOH_VAQ_data <- merge( overlap_data, filter_data, by="SAMPLE" )
names( LOH_VAQ_data )[ names( LOH_VAQ_data )=="FILT_5_COUNT" ] <- "VAQ_COUNT"
names( LOH_VAQ_data )[ names( LOH_VAQ_data )=="FILT_7_COUNT" ] <- "LOH_COUNT"
names( LOH_VAQ_data )[ names( LOH_VAQ_data )=="FILT_5_PERCENT" ] <- "VAQ_PERCENT"
names( LOH_VAQ_data )[ names( LOH_VAQ_data )=="FILT_7_PERCENT" ] <- "LOH_PERCENT"

# note that this is pretty close to 1.00 every time, hence the inversion...
LOH_VAQ_data$SUM=LOH_VAQ_data$LOH_PERCENT + LOH_VAQ_data$VAQ_PERCENT

#############################################################
# STEP 2: make figure 8, x-axis=OVERLAP, y-axis=LOH_PERCENT #
#############################################################

par(bty='n')
plot( LOH_VAQ_data$OVERLAP, LOH_VAQ_data$LOH_PERCENT, xlab="Length of overlap", ylab="% of overlap filtered out by LOH", pch=20 )

#############################################################
# STEP 3: make figure 9, x-axis=OVERLAP, y-axis=VAQ_PERCENT #
#############################################################

par(bty='n')
plot( LOH_VAQ_data$OVERLAP, LOH_VAQ_data$VAQ_PERCENT, xlab="Length of overlap", ylab="% of overlap filtered out by VAQ", pch=20 )