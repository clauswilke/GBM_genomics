###################################
## first do the Jaccard analysis ##
###################################

## this table has the size of the overlap and the size of the total, everything for the "before" Jaccard, for each sample
data1 = read.table( "../FIGURE_DATA/STD40_unfiltered_overlap.txt", header=T, stringsAsFactors=FALSE )
## this table has the number of mutations excluded from the total for each sample for each filter
data2 = read.table( "../FIGURE_DATA/numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt", header=T, stringsAsFactors=FALSE )
## this table has the number of mutations for each sample/filter that were excluded FOR THE OVERLAP
data3 = read.table( "../FIGURE_DATA/overlap_numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt", header=T, stringsAsFactors=FALSE )

## from data1, we want the following columns: SAMPLE, OVERLAP, TOTAL
## (Jaccard^before( WGS, WGA )=OVERLAP/TOTAL
data1$DIFFERENCE<-NULL
data1$ONLY_WGA.C282.<-NULL
data1$TOTAL_WGA.C282.<-NULL
data1$ONLY_WGS.C484.<-NULL
data1$TOTAL_WGS.C484.<-NULL
data1$PERCENT_OVERLAP_by_ALL<-NULL
data1$PERCENT_OVERLAP_by_WGA.282<-NULL
data1$PERCENT_OVERLAP_by_WGS.484<-NULL
data1$PERCENT_DIFFERENCE_by_ALL<-NULL
data1$PERCENT_DIFFERENCE_by_WGA.282<-NULL
data1$PERCENT_DIFFERENCE_by_WGS.484<-NULL

data1$Jaccard.before<-data1$OVERLAP/data1$TOTAL

## from data2 (excluded by filter from SAMPLE), we want SAMPLE, FILT_3, FILT_4, FILT_5, FILT_6, FILT_7 (renamed tot.)
# change C282 and C484 to a new column, type
data2$newtype<-"NA"
data2$newname<-"NA"
for (i in 1:length(data2$SAMPLE)){
	oldname<-data2$SAMPLE[ i ]
	newname<-substr( oldname, 6, 18 )
	data2$newname[ i ]<-newname
}
for (i in 1:length(data2$SAMPLE)){
	oldname<-data2$SAMPLE[ i ]
	newtype<-substr( oldname, 1, 5 )
	data2$newtype[ i ]<-newtype
}
# divide data2 into two data frames, by type
WGA<-data2[data2$newtype=="C282.",]
WGS<-data2[data2$newtype=="C484.",]
# WGA change the name of the columns in each new dataframe, and get rid of the type column
for (i in 1:length(WGA$SAMPLE)){
	WGA$SAMPLE[ i ]<-WGA$newname[ i ]
}
WGA$newname<-NULL
WGA$newtype<-NULL
names( WGA )[ 2 ] <- "WGA.tot.F2"
names( WGA )[ 3 ] <- "WGA.tot.F3"
names( WGA )[ 4 ] <- "WGA.tot.F4"
names( WGA )[ 5 ] <- "WGA.tot.F5"
names( WGA )[ 6 ] <- "WGA.tot.F6"
names( WGA )[ 7 ] <- "WGA.tot.F7"
names( WGA )[ 8 ] <- "WGA.tot.F8"
# WGS change the name of the columns in each new dataframe, and get rid of the type column
for (i in 1:length(WGS$SAMPLE)){
	WGS$SAMPLE[ i ]<-WGS$newname[ i ]
}
WGS$newname<-NULL
WGS$newtype<-NULL
names( WGS )[ 2 ] <- "WGS.tot.F2"
names( WGS )[ 3 ] <- "WGS.tot.F3"
names( WGS )[ 4 ] <- "WGS.tot.F4"
names( WGS )[ 5 ] <- "WGS.tot.F5"
names( WGS )[ 6 ] <- "WGS.tot.F6"
names( WGS )[ 7 ] <- "WGS.tot.F7"
names( WGS )[ 8 ] <- "WGS.tot.F8"
# merge the two data frames with data1
merged1<-merge( data1, WGA )
merged2<-merge( merged1, WGS )

## from data3 (excluded by sample from OVERLAP), we want SAMPLE, FILT_3, FILT_4, FILT_5, FILT_6, FILT_7 (renamed over.)
names( data3 )[ 2 ] <- "over.F2"
names( data3 )[ 4 ] <- "over.F3"
names( data3 )[ 6 ] <- "over.F4"
names( data3 )[ 8 ] <- "over.F5"
names( data3 )[ 10 ] <- "over.F6"
names( data3 )[ 12 ] <- "over.F7"
names( data3 )[ 14 ] <- "over.F8"
data3$FILT_2_PERCENT<-NULL
data3$FILT_3_PERCENT<-NULL
data3$FILT_4_PERCENT<-NULL
data3$FILT_5_PERCENT<-NULL
data3$FILT_6_PERCENT<-NULL
data3$FILT_7_PERCENT<-NULL
data3$FILT_8_PERCENT<-NULL
mergedData<-merge( merged2, data3 )

## per filter, 
## (Jaccard^after( WGS, WGA )=(OVERLAP - over.F3)/(TOTAL - WGS.tot.F3 - WGA.tot.F3 + over.F3)
mergedData$F2.overlap<-mergedData$OVERLAP-mergedData$over.F2
mergedData$F3.overlap<-mergedData$OVERLAP-mergedData$over.F3
mergedData$F4.overlap<-mergedData$OVERLAP-mergedData$over.F4
mergedData$F5.overlap<-mergedData$OVERLAP-mergedData$over.F5
mergedData$F6.overlap<-mergedData$OVERLAP-mergedData$over.F6
mergedData$F7.overlap<-mergedData$OVERLAP-mergedData$over.F7
mergedData$F8.overlap<-mergedData$OVERLAP-mergedData$over.F8

mergedData$F2.total<-mergedData$TOTAL-mergedData$WGS.tot.F2-mergedData$WGA.tot.F2+mergedData$over.F2
mergedData$F3.total<-mergedData$TOTAL-mergedData$WGS.tot.F3-mergedData$WGA.tot.F3+mergedData$over.F3
mergedData$F4.total<-mergedData$TOTAL-mergedData$WGS.tot.F4-mergedData$WGA.tot.F4+mergedData$over.F4
mergedData$F5.total<-mergedData$TOTAL-mergedData$WGS.tot.F5-mergedData$WGA.tot.F5+mergedData$over.F5
mergedData$F6.total<-mergedData$TOTAL-mergedData$WGS.tot.F6-mergedData$WGA.tot.F6+mergedData$over.F6
mergedData$F7.total<-mergedData$TOTAL-mergedData$WGS.tot.F7-mergedData$WGA.tot.F7+mergedData$over.F7
mergedData$F8.total<-mergedData$TOTAL-mergedData$WGS.tot.F8-mergedData$WGA.tot.F8+mergedData$over.F8

mergedData$Jaccard.F2<-mergedData$F2.overlap/mergedData$F2.total
mergedData$Jaccard.F3<-mergedData$F2.overlap/mergedData$F3.total
mergedData$Jaccard.F4<-mergedData$F2.overlap/mergedData$F4.total
mergedData$Jaccard.F5<-mergedData$F2.overlap/mergedData$F5.total
mergedData$Jaccard.F6<-mergedData$F2.overlap/mergedData$F6.total
mergedData$Jaccard.F7<-mergedData$F2.overlap/mergedData$F7.total
mergedData$Jaccard.F8<-mergedData$F2.overlap/mergedData$F8.total

## finally, calculate delta( Jaccard )

mergedData$delta.Jaccard.F2<-((mergedData$Jaccard.F2-mergedData$Jaccard.before)/mergedData$Jaccard.before)*100
mergedData$delta.Jaccard.F3<-((mergedData$Jaccard.F3-mergedData$Jaccard.before)/mergedData$Jaccard.before)*100
mergedData$delta.Jaccard.F4<-((mergedData$Jaccard.F4-mergedData$Jaccard.before)/mergedData$Jaccard.before)*100
mergedData$delta.Jaccard.F5<-((mergedData$Jaccard.F5-mergedData$Jaccard.before)/mergedData$Jaccard.before)*100
mergedData$delta.Jaccard.F6<-((mergedData$Jaccard.F6-mergedData$Jaccard.before)/mergedData$Jaccard.before)*100
mergedData$delta.Jaccard.F7<-((mergedData$Jaccard.F7-mergedData$Jaccard.before)/mergedData$Jaccard.before)*100
mergedData$delta.Jaccard.F8<-((mergedData$Jaccard.F8-mergedData$Jaccard.before)/mergedData$Jaccard.before)*100

##################################################
## THEN DO THE DIFFSET ANALYSIS WITH THE RATIOS ##
##################################################

## LOAD DATA

data4=read.table( "../FIGURE_DATA/overlap_numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt", header=T, stringsAsFactors=FALSE )
data5=read.table( "../FIGURE_DATA/numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt", header=T, stringsAsFactors=FALSE )

## now re-name and delete columns
names( data4 )[ 2 ] <- "F2overlap"
names( data4 )[ 4 ] <- "F3overlap"
names( data4 )[ 6 ] <- "F4overlap"
names( data4 )[ 8 ] <- "F5overlap"
names( data4 )[ 10 ] <- "F6overlap"
names( data4 )[ 12 ] <- "F7overlap"
names( data4 )[ 14 ] <- "F8overlap"

data4$FILT_2_PERCENT<-NULL
data4$FILT_3_PERCENT<-NULL
data4$FILT_4_PERCENT<-NULL
data4$FILT_5_PERCENT<-NULL
data4$FILT_6_PERCENT<-NULL
data4$FILT_7_PERCENT<-NULL
data4$FILT_8_PERCENT<-NULL

## now merge the two data frames
for (i in 1:length(data4$SAMPLE)){
	short<-data4$SAMPLE[ i ]
	newname<-paste( "C282", short, sep="." )
	data4$SAMPLE[ i ]<-newname
}

# duplicate the data.frame
newdata <- data4

# change C282 to C484
for (i in 1:length(newdata$SAMPLE)){
	short <- newdata$SAMPLE[ i ]
	old <- substring( short, 6, )
	newname <- paste( "C484", old, sep="." )
	newdata$SAMPLE[ i ]<-newname
}

# add the two together
overlapdata <- merge(data4, newdata, all=TRUE )

# merge with data5...
names( data5 )[ 2 ] <- "F2all"
names( data5 )[ 3 ] <- "F3all"
names( data5 )[ 4 ] <- "F4all"
names( data5 )[ 5 ] <- "F5all"
names( data5 )[ 6 ] <- "F6all"
names( data5 )[ 7 ] <- "F7all"
names( data5 )[ 8 ] <- "F8all"

all <- merge( overlapdata, data5 )

## then calculate the diffset size
all$F2diff <- all$F2all - all$F2overlap
all$F3diff <- all$F3all - all$F3overlap
all$F4diff <- all$F4all - all$F4overlap
all$F5diff <- all$F5all - all$F5overlap
all$F6diff <- all$F6all - all$F6overlap
all$F7diff <- all$F7all - all$F7overlap
all$F8diff <- all$F8all - all$F8overlap

## a ratio ...
all$F2ratio <- all$F2diff / all$F2overlap
all$F3ratio <- all$F3diff / all$F3overlap
all$F4ratio <- all$F4diff / all$F4overlap
all$F5ratio <- all$F5diff / all$F5overlap
all$F6ratio <- all$F6diff / all$F6overlap
all$F7ratio <- all$F7diff / all$F7overlap
all$F8ratio <- all$F8diff / all$F8overlap

## add normalized ratio...
#all$F2ratioNORM <- (all$F2diff/all$F2all) / (all$F2overlap/all$F2all)
#all$F3ratioNORM <- (all$F3diff/all$F2all) / (all$F3overlap/all$F2all)
#all$F4ratioNORM <- (all$F4diff/all$F2all) / (all$F4overlap/all$F2all)
#all$F5ratioNORM <- (all$F5diff/all$F2all) / (all$F5overlap/all$F2all)
#all$F6ratioNORM <- (all$F6diff/all$F2all) / (all$F6overlap/all$F2all)
#all$F7ratioNORM <- (all$F7diff/all$F2all) / (all$F7overlap/all$F2all)
#all$F8ratioNORM <- (all$F8diff/all$F2all) / (all$F8overlap/all$F2all)

## then make some averages...

## ratio averages...
all$F3ratio[all$F3ratio=="Inf"] <- NA
all$F3ratio[all$F3ratio=="NaN"] <- NA
mean(all$F3ratio, na.rm=TRUE)
median(all$F3ratio, na.rm=TRUE)

all$F4ratio[all$F4ratio=="Inf"] <- NA
all$F4ratio[all$F4ratio=="NaN"] <- NA
mean(all$F4ratio, na.rm=TRUE)
median(all$F4ratio, na.rm=TRUE)

all$F5ratio[all$F5ratio=="Inf"] <- NA
all$F5ratio[all$F5ratio=="NaN"] <- NA
mean(all$F5ratio, na.rm=TRUE)
median(all$F5ratio, na.rm=TRUE)

all$F6ratio[all$F6ratio=="Inf"] <- NA
all$F6ratio[all$F6ratio=="NaN"] <- NA
mean(all$F6ratio, na.rm=TRUE)
median(all$F6ratio, na.rm=TRUE)

all$F7ratio[all$F7ratio=="Inf"] <- NA
all$F7ratio[all$F7ratio=="NaN"] <- NA
mean(all$F7ratio, na.rm=TRUE)
median(all$F7ratio, na.rm=TRUE)

########################################################
## THEN MAKE FIGURE 6 (COMPOUND FIGURE) FOR THE PAPER ##
########################################################

#pdf( "../FIGURE_PDFS/Figure5.pdf", width=6, height=5 )

par(bty='n', mfrow=c(2,2), mar=c(3,4.1,2,1) )

plot( all$F7ratio, ylim=c(0,80), pch=20, col='black', ylab="Ratio of difference to overlap", xaxt='n', xlab="Samples" )
points( all$F5ratio, pch=20, col='navy' )
points( all$F3ratio, pch=20, col='red' )
points( all$F4ratio, pch=20, col='yellow' )
points( all$F6ratio, pch=20, col='orange' )
legend( 'topright', c( "VAQ", "LOH", "10bp-SNV", "10bp-INDEL", "dbSNP" ), col=c( 'black', 'navy', 'yellow', 'red', 'orange' ), pch=20 )

boxplot( all[,c(26,28,25,24,27)], ylim=c(1,50), xaxt='n', pch=20, ylab="" )
axis( 1, at=1:5, padj=0.8, cex.axis=0.7, labels=c( "VAQ", "LOH", "10bp-SNV", "10bp-INDEL", "dbSNP" ) )

boxplot( mergedData[,c(52,50)], ylim=c(0,450), xaxt='n', ylab=expression(paste( Delta, "Jaccard( WGS, WGA )")), pch=20 )
axis( 1, at=1:2, padj = 0.8, cex.axis=0.7, labels=c( "LOH", "VAQ" ) )

boxplot( mergedData[,c(49,48,51,53)], ylim=c(0,15), xaxt='n', ylab="", pch=20 )
axis( 1, at=1:4, padj = 0.8, cex.axis=0.7, labels=c( "10bp-SNV", "10bp-INDEL", "dbSNP", "<10%" ) )

#dev.off()




##############
## OLD CODE ##
##############

#pdf( "../FIGURE_PDFS/Jaccard1.pdf", width=6, height=5 )
par(bty='n', mar=c(5.1,4.1,1.1,0) )
boxplot( mergedData[,c(50,52)], ylim=c(0,450), xaxt='n', ylab=expression(paste( "Percentage ", Delta, "Jaccard( WGS, WGA )")), pch=20 )
axis( 1, at=1:2, padj = 0.8, cex.axis=0.7, labels=c( "VAQ", "LOH" ) )
#dev.off()

#pdf( "../FIGURE_PDFS/Jaccard2.pdf", width=6, height=5 )
par(bty='n', mar=c(5.1,4.1,1.1,0))
boxplot( mergedData[,c(49,48,51,53)], ylim=c(0,15), xaxt='n', ylab=expression(paste( "Percentage ", Delta, "Jaccard( WGS, WGA )")), pch=20 )
axis( 1, at=1:4, padj = 0.8, cex.axis=0.7, labels=c( "10bp-SNV", "10bp-INDEL", "dbSNP", "<10%" ) )
#dev.off()

#pdf( "../FIGURE_PDFS/Diffset_bypoints.pdf" )
par(bty='n')
plot( all$F7ratio, ylim=c(0,80), pch=20, col='black' )
points( all$F5ratio, pch=20, col='navy' )
points( all$F3ratio, pch=20, col='red' )
points( all$F4ratio, pch=20, col='yellow' )
points( all$F6ratio, pch=20, col='orange' )
#dev.off()

# and here is a boxplot...
#24-28
#pdf( "../FIGURE_PDFS/Diffset_bybars.pdf" )
par(bty='n')
boxplot( all[,c(26,28,25,24,27)], ylim=c(1,50), xaxt='n', pch=20, ylab="Ratio of the diffset to the overlap, per filter" )
axis( 1, at=1:5, padj=0.8, cex.axis=0.7, labels=c( "VAQ", "LOH", "10bp-SNV", "10bp-INDEL", "dbSNP" ) )
#dev.off()



