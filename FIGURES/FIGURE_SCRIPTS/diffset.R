## LOAD DATA

data=read.table( "../FIGURE_DATA/overlap_numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt", header=T, stringsAsFactors=FALSE )
data2=read.table( "../FIGURE_DATA/numbers_excluded_by_filter_MUTATION_CALLS_STD40.txt", header=T, stringsAsFactors=FALSE )

## now re-name and delete columns
names( data )[ 2 ] <- "F2overlap"
names( data )[ 4 ] <- "F3overlap"
names( data )[ 6 ] <- "F4overlap"
names( data )[ 8 ] <- "F5overlap"
names( data )[ 10 ] <- "F6overlap"
names( data )[ 12 ] <- "F7overlap"
names( data )[ 14 ] <- "F8overlap"

data$FILT_2_PERCENT<-NULL
data$FILT_3_PERCENT<-NULL
data$FILT_4_PERCENT<-NULL
data$FILT_5_PERCENT<-NULL
data$FILT_6_PERCENT<-NULL
data$FILT_7_PERCENT<-NULL
data$FILT_8_PERCENT<-NULL

## now merge the two data frames
for (i in 1:length(data$SAMPLE)){
	short<-data$SAMPLE[ i ]
	newname<-paste( "C282", short, sep="." )
	data$SAMPLE[ i ]<-newname
}

# duplicate the data.frame
newdata <- data

# change C282 to C484
for (i in 1:length(newdata$SAMPLE)){
	short <- newdata$SAMPLE[ i ]
	old <- substring( short, 6, )
	newname <- paste( "C484", old, sep="." )
	newdata$SAMPLE[ i ]<-newname
}

# add the two together
overlapdata <- merge(data, newdata, all=TRUE )

# merge with data2...
names( data2 )[ 2 ] <- "F2all"
names( data2 )[ 3 ] <- "F3all"
names( data2 )[ 4 ] <- "F4all"
names( data2 )[ 5 ] <- "F5all"
names( data2 )[ 6 ] <- "F6all"
names( data2 )[ 7 ] <- "F7all"
names( data2 )[ 8 ] <- "F8all"

all <- merge( overlapdata, data2 )

## then calculate the diffset size
all$F2diff <- all$F2all - all$F2overlap
all$F3diff <- all$F3all - all$F3overlap
all$F4diff <- all$F4all - all$F4overlap
all$F5diff <- all$F5all - all$F5overlap
all$F6diff <- all$F6all - all$F6overlap
all$F7diff <- all$F7all - all$F7overlap
all$F8diff <- all$F8all - all$F8overlap

## add ratio...
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

## finally make some (supplementary) figures

# here is a plot of points...
par(bty='n')
plot( all$F7ratio, ylim=c(0,80), pch=20, col='black' )
points( all$F5ratio, pch=20, col='navy' )
points( all$F3ratio, pch=20, col='red' )
points( all$F4ratio, pch=20, col='yellow' )
points( all$F6ratio, pch=20, col='orange' )

# and here is a boxplot...
#24-28
par(bty='n')
boxplot( all[,c(26,28,25,24,27)], ylim=c(1,50), xaxt='n', pch=20, ylab="Ratio of the diffset to the overlap, per filter" )
axis( 1, at=1:6, padj=0.8, cex.axis=0.7, labels=c( "VAQ", "LOH", "10bp-SNV", "10bp-INDEL", "dbSNP" ) )

