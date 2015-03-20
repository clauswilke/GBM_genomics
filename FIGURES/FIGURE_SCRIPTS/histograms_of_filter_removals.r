## FIGURES IN PAPER: 6 and 7
## it's the two box plots

# this script looks in different ways at which mutations are removed when I run the filters

# first look at the STD_40 data
# first, just everything removed, each file individually, how many by each filter 3-8
data = read.table( "../FIGURE_DATA/numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt", header=T )

# make density of each filter removal
d3 = density( data$FILTER_3 )
d4 = density( data$FILTER_4 )
d5 = density( data$FILTER_5 )
d6 = density( data$FILTER_6 )
d7 = density( data$FILTER_7 )
d8 = density( data$FILTER_8 )

# make the box plot that will show everything
png( "../FIGURE_PDFS/boxplot_number_filtered.png", width = 480, height = 360, res=55 )
par(bty='n')
boxplot( data[,c(5,7,4,3,6,8)], log='y', ylim=c(1,100000), xaxt='n', pch=20, ylab="# of putative SNPs filtered, per sample" )
axis( 1, at=1:6, padj=0.8, cex.axis=0.7, labels=c( "Removes Variant\nAllele Quality\n(VAQ)<20", "Removes\nLOH", "Removes putative\nSNPs within\n10 bp window\nof other SNPs", "Removes putative\nSNPs within\n10 bp window\nof indels", "Removes\noverlap\nwith dbSNP\ncoverage", "Removes when\nalt. allele\ncoverage<=10%" ) )
dev.off()

pdf( "../FIGURE_PDFS/boxplot_number_filtered.pdf", width=6, height=5 )
par(bty='n', mar=c(5.1,4.1,1.1,0))
boxplot( data[,c(5,7,4,3,6,8)], log='y', ylim=c(1,100000), xaxt='n', pch=20, ylab="# of putative SNPs filtered, per sample" )
axis( 1, at=1:6, padj=0.8, cex.axis=0.7, labels=c( "VAQ", "LOH", "10bp-SNP", "10bp-INDEL", "dbSNP", "<10%" ) )
dev.off()

### still working on two problems:
### how to make the figure horizontally longer so that all labels fit, and 
### how to move the x-axis lables down relative to the x-axis, again so they fit


#axis( 1, at=1:6 )
#mtext( side=1, text="Variant Allele\nQuality (VAQ)\n>20", line=4, cex=-20 )
#, "Removes\nLOH", "No putative\nSNPs within\n10 bp window\nof other SNPs", "No putative\nSNPs within\n10 bp window\nof indels", "No overlap\nwith dbSNP\ncoverage", "Alternate\nallele\ncoverage\n>=10%"), line=4 )

## this second part of the script makes the same figure, but for only the putative SNPs that are in the overlaps in the unfiltered data

overlap_data = read.table( "../FIGURE_DATA/overlap_numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt", header=T )

# make the same plot for the overlaps
png( "../FIGURE_PDFS/boxplot_percent_overlap_filtered.png", width = 480, height = 360, res=55 )
par( bty='n' )
boxplot( overlap_data[,c(9,13,7,5,11,15)], ylim=c(0,1), xaxt='n', ylab="% of unfiltered overlap filtered, per sample", pch=20 )
axis( 1, at=1:6, padj = 0.8, cex.axis=0.7, labels=c( "Removes Variant\nAllele Quality\n(VAQ)<20", "Removes\nLOH", "Removes putative\nSNPs within\n10 bp window\nof other SNPs", "Removes putative\nSNPs within\n10 bp window\nof indels", "Removes\noverlap\nwith dbSNP\ncoverage", "Removes when\nalt. allele\ncoverage<=10%" ) )
dev.off()

# make the same plot for the overlaps
pdf( "../FIGURE_PDFS/boxplot_percent_overlap_filtered.pdf", width=6, height=5 )
par(bty='n', mar=c(5.1,4.1,1.1,0))
boxplot( overlap_data[,c(9,13,7,5,11,15)], ylim=c(0,1), xaxt='n', ylab="% of unfiltered overlap filtered, per sample", pch=20 )
axis( 1, at=1:6, padj = 0.8, cex.axis=0.7, labels=c( "VAQ", "LOH", "10bp-SNP", "10bp-INDEL", "dbSNP", "<10%" ) )
dev.off()
