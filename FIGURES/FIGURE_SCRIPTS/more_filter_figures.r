# the purpose of this script is to make some figures...

## 1. For LOH, VAQ, plot size of overlap (x-axis) against percent of overlap removed by relevant filter (y-axis)

overlap_data = read.table( "../FIGURE_DATA/STD40_unfiltered_overlap.txt", header=T )
filter_data = read.table( "../FIGURE_DATA/overlap_numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt", header=T )

all_data <- merge( overlap_data, filter_data, by="SAMPLE" )

# LOH (filter 7)
LOH_vars = c( "OVERLAP", "FILT_7_COUNT" )
LOH_data = all_data[ LOH_vars ]
FILT_7_PERCENT_new = LOH_data$FILT_7_COUNT / LOH_data$OVERLAP
LOH_data$FILT_7_PERCENT_new = FILT_7_PERCENT_new

pdf( "/Users/dakota/github/GBM_genomics/FIGURES/FIGURE_PDFS/LOH_VAQ/LOH_all.pdf" )
par(bty='n')
plot( LOH_data$OVERLAP, LOH_data$FILT_7_PERCENT_new, xlab="Length of overlap", ylab="% of overlap filtered out by LOH", pch=20 )
dev.off()

#png( "/Users/dakota/Desktop/FIG_RD/LOH_partial.png" )
#par(bty='n')
#plot( LOH_data$OVERLAP, LOH_data$FILT_7_PERCENT_new, xlim=c( 0, 500 ), xlab="Length of overlap", ylab="% of overlap filtered out by LOH", pch=20 )
#dev.off()

# VAQ (filter 5)
VAQ_vars = c( "OVERLAP", "FILT_5_PERCENT" )
VAQ_data = all_data[ VAQ_vars ]

pdf( "/Users/dakota/github/GBM_genomics/FIGURES/FIGURE_PDFS/LOH_VAQ/VAQ_all.pdf" )
par(bty='n')
plot( VAQ_data$OVERLAP, VAQ_data$FILT_5_PERCENT, xlab="Length of overlap", ylab="% of overlap filtered out by VAQ", pch=20 )
dev.off()

#png( "/Users/dakota/Desktop/FIG_RD/VAQ_partial.png" )
#par(bty='n')
#plot( VAQ_data$OVERLAP, VAQ_data$FILT_5_PERCENT, xlim=c( 0, 1000 ), xlab="Length of overlap", ylab="% of overlap filtered out by VAQ", pch=20 )
#dev.off()

## 2. Plot LOH against VAQ: are they removing the same thing? 

compare_vars = c( "FILT_5_COUNT", "FILT_7_COUNT" )
comparison_data = all_data[ compare_vars ]

#png( "/Users/dakota/Desktop/FIG_RD/LOHvVAQ_all.png" )
#par(bty='n')
#plot( comparison_data$FILT_5_COUNT , comparison_data$FILT_7_COUNT, xlab="% of overlap filtered out by VAQ", ylab="% of overlap filtered out by LOH", pch=20 )
#dev.off()

#png( "/Users/dakota/Desktop/FIG_RD/LOHvVAQ_partial.png" )
#par(bty='n')
#plot( comparison_data$FILT_5_COUNT , comparison_data$FILT_7_COUNT, xlim=c( 0, 150 ), ylim=c( 0, 1200 ), xlab="% of overlap filtered out by VAQ", ylab="% of overlap filtered out by LOH", pch=20 )
#dev.off()

## 2.5 Instead of just looking at quantity, can we look at whether or not they are the same mutations being filtered by each filter? (this is a python script that I will feed into this point... when everyone in lab SHUTS THE FUCK UP. ALL OF YOU. SHUT THE FUCK UP.

## 3. What I was working on, on phylocluster, figure this out tonight...