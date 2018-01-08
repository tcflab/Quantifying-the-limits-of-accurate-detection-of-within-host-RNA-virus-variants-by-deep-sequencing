require(ggplot2)
require(grid)
require(stringr)
require(reshape2)
require(format)
require(plyr)

# blue/red pallette
color1 <- rgb(0,148,194, maxColorValue=255)
color2 <- rgb(74,204,245, maxColorValue=255)
color3 <- rgb(0,90,117, maxColorValue=255)
color5 <- rgb(97,7,4, maxColorValue=255)
color4 <- rgb(207,48,31, maxColorValue=255)
color6 <- rgb(255,128,114, maxColorValue=255)
color7 <- rgb(255,255,255, maxColorValue=255)


###############################################################################
# Figure 6: Normalization does not impact SNP detection
###############################################################################
data_1000x=read.table("~/Documents/Graduate_School/CA04_SNP_detection/revisions/data_reanalysis_with_illumina_pipeline/coverage_depth_normalization/combined_vcfs_amplicon_1000x_171221.txt", header=T, stringsAsFactors=F, na.string=c(NA,"NA","na"), sep="\t", col.names=c("site","reference_allele","variant_allele","frequency","input","replicate"))

data_1000x$subsampled = "1000x"

data=read.table("~/Documents/Graduate_School/CA04_SNP_detection/revisions/data_reanalysis_with_illumina_pipeline/combined_vcfs_amplicon_171221.txt", header=T, stringsAsFactors=F, na.string=c(NA,"NA","na"), sep="\t", col.names=c("site","reference_allele","variant_allele","frequency","input","replicate"))

data$subsampled = "no"

# combine/merge dataframes; all=TRUE will make it so that all rows are merged, and if a SNP was only detected in either 1000x or the original data, the field that it was not detected in will be filled in with an NA and not just omitted from the dataframe
combined_df = merge(data, data_1000x, by=c("site","reference_allele","variant_allele","input","replicate"), all=TRUE)




### TO MAKE FIGURE 6B, SUBSET DATA TO INCLUDE ONLY SNPs < 10% in frequency, unhash the next 3 lines and run  ###
#low_freq_df = data[data$frequency <= 10, ]
#low_freq_1000x_df = data_1000x[data_1000x$frequency <= 10, ]
#combined_df = merge(low_freq_df, low_freq_1000x_df, by=c("site","reference_allele","variant_allele","input","replicate"), all=TRUE)



# fill in columns to remove NAs and replace with values or 0s
combined_df$subsampled.x = "no"
combined_df$subsampled.y = "1000x"
combined_df$frequency.x[is.na(combined_df$frequency.x)] <- 0
combined_df$frequency.y[is.na(combined_df$frequency.y)] <- 0

combined_df$input <- as.character(as.numeric(combined_df$input))
combined_df$input <- gsub("10000","1e+04",combined_df$input)
combined_df$input <- gsub("1000","1e+03",combined_df$input)
combined_df$input <- gsub("100","1e+02",combined_df$input)
combined_df$input <- gsub("10","1e+01",combined_df$input)

combined_df$input <- factor(combined_df$input, levels=c("1e+01","1e+02","1e+03","1e+04","1e+05","1e+06"))


# calculate regression between input and standard deviation
r2.lm = lm(frequency.x~frequency.y, data=combined_df)
r2.lm$residuals #get residualsa
summary(r2.lm)

# compare subsampled vs. not as a scatter plot (correlation)
ggplot(combined_df, aes(x=frequency.x,y=frequency.y,color=input))+geom_point(size=2)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+theme(plot.title=element_text(size=13))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title=element_text(size=13, vjust=0.2))+theme(axis.text=element_text(size=11, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=11))+theme(legend.title=element_text(size=13, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+labs(x="SNP frequency (%)",y="SNP frequency (%) after subsampling")+scale_x_continuous(breaks=seq(0,100,10), limits=c(0,100))+scale_y_continuous(breaks=seq(0,100,10), limits=c(0,100))+scale_colour_manual(values=c(color1,color2,color3,color4,color5,color6))


# calculate the average difference between subsampled vs. not SNP frequencies

combined_df$difference <- abs(combined_df$frequency.x - combined_df$frequency.y)
mean(combined_df$difference)

df_1e2 = combined_df[combined_df$input == "1e+02",]
df_1e3 = combined_df[combined_df$input == "1e+03",]
df_1e4 = combined_df[combined_df$input == "1e+04",]
df_1e5 = combined_df[combined_df$input == 1e5,]
df_1e6 = combined_df[combined_df$input == 1e6,]
df_1e7 = combined_df[combined_df$input == 1e7,]

# perform paired t-tests comparing frequencies detectedd before and after coverage normalization
t.test(combined_df$frequency.x, combined_df$frequency.y, paired=TRUE)
t.test(df_1e7$frequency.x,df_1e7$frequency.y,paired=TRUE)
t.test(df_1e6$frequency.x,df_1e6$frequency.y,paired=TRUE)
t.test(df_1e5$frequency.x,df_1e5$frequency.y,paired=TRUE)
t.test(df_1e4$frequency.x,df_1e4$frequency.y,paired=TRUE)
t.test(df_1e3$frequency.x,df_1e3$frequency.y,paired=TRUE)
t.test(df_1e2$frequency.x,df_1e2$frequency.y,paired=TRUE)
