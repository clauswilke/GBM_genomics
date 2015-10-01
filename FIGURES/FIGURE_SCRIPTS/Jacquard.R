## first get data

## this table has the size of the overlap and the size of the total, everything for the "before" Jacquard, for each sample
data1 = read.table( "../FIGURE_DATA/STD40_unfiltered_overlap.txt", header=T, stringsAsFactors=FALSE )
## this table has the number of mutations excluded from the total for each sample for each filter
data2 = read.table( "../FIGURE_DATA/numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt", header=T, stringsAsFactors=FALSE )
## this table has the number of mutations for each sample/filter that were excluded FOR THE OVERLAP
data3 = read.table( "../FIGURE_DATA/overlap_numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt", header=T, stringsAsFactors=FALSE )

## from data1, we want the following columns: SAMPLE, OVERLAP, TOTAL
## (Jacquard^before( WGS, WGA )=OVERLAP/TOTAL
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

data1$Jacquard.before<-data1$OVERLAP/data1$TOTAL

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
## (Jacquard^after( WGS, WGA )=(OVERLAP - over.F3)/(TOTAL - WGS.tot.F3 - WGA.tot.F3 + over.F3)
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

mergedData$Jacquard.F2<-mergedData$F2.overlap/mergedData$F2.total
mergedData$Jacquard.F3<-mergedData$F2.overlap/mergedData$F3.total
mergedData$Jacquard.F4<-mergedData$F2.overlap/mergedData$F4.total
mergedData$Jacquard.F5<-mergedData$F2.overlap/mergedData$F5.total
mergedData$Jacquard.F6<-mergedData$F2.overlap/mergedData$F6.total
mergedData$Jacquard.F7<-mergedData$F2.overlap/mergedData$F7.total
mergedData$Jacquard.F8<-mergedData$F2.overlap/mergedData$F8.total

## finally, calculate delta( Jacquard )

mergedData$delta.Jacquard.F2<-((mergedData$Jacquard.F2-mergedData$Jacquard.before)/mergedData$Jacquard.before)*100
mergedData$delta.Jacquard.F3<-((mergedData$Jacquard.F3-mergedData$Jacquard.before)/mergedData$Jacquard.before)*100
mergedData$delta.Jacquard.F4<-((mergedData$Jacquard.F4-mergedData$Jacquard.before)/mergedData$Jacquard.before)*100
mergedData$delta.Jacquard.F5<-((mergedData$Jacquard.F5-mergedData$Jacquard.before)/mergedData$Jacquard.before)*100
mergedData$delta.Jacquard.F6<-((mergedData$Jacquard.F6-mergedData$Jacquard.before)/mergedData$Jacquard.before)*100
mergedData$delta.Jacquard.F7<-((mergedData$Jacquard.F7-mergedData$Jacquard.before)/mergedData$Jacquard.before)*100
mergedData$delta.Jacquard.F8<-((mergedData$Jacquard.F8-mergedData$Jacquard.before)/mergedData$Jacquard.before)*100


## then average everything for the paper






