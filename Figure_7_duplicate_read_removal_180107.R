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
# Figure 7: Duplicate read removal does not improve reproducibility
###############################################################################

### FIGURE 7E: DUPLICATE READ REMOVAL SNP DETECTION REPRODUCIBILITY ##

# read in dataframes and combine
bbtools=read.table("~/Documents/Graduate_School/CA04_SNP_detection/revisions/data_reanalysis_with_illumina_pipeline/duplicate_read_removal/combined_vcfs_amplicon_dedupe_180107.txt", header=T, stringsAsFactors=F, na.string=c(NA,"NA","na"), sep="\t", col.names=c("reference_position","reference","allele","frequency","input","replicate"))
bbtools$method = "bbtools"

picard=read.table("~/Documents/Graduate_School/CA04_SNP_detection/revisions/data_reanalysis_with_illumina_pipeline/duplicate_read_removal/combined_vcfs_amplicon_picard_180107.txt", header=T, stringsAsFactors=F, na.string=c(NA,"NA","na"), sep="\t", col.names=c("reference_position","reference","allele","frequency","input","replicate"))
picard$method = "picard"

original=read.table("~/Documents/Graduate_School/CA04_SNP_detection/revisions/data_reanalysis_with_illumina_pipeline/combined_vcfs_amplicon_171221.txt", header=T, stringsAsFactors=F, na.string=c(NA,"NA","na"), sep="\t", col.names=c("reference_position","reference","allele","frequency","input","replicate"))
original$method = "duplicate_reads_kept"

# combined dataframes
data = rbind(original, bbtools, picard)

# calculate the number of times a specific SNP is detected in each dilution among replicates (1,2, or 3 times) and see if that differs depending on duplicate removal method. 

method_list = unique(data$method)
dilution_list = unique(data$input)
count = 0

for (i in method_list)
{
  for (d in dilution_list)
  {
    df <- data[data$method == i & data$input == d,]
    table <- as.data.frame(table(df$reference_position))
    table$input = d
    table$method = i
    
    number_sites = nrow(table)
    fraction_in_all_3 = (sum(table$Freq == 3)/number_sites)*100
    fraction_in_2 = (sum(table$Freq == 2)/number_sites)*100
    fraction_in_1 = (sum(table$Freq == 1)/number_sites)*100
    fraction = c(fraction_in_all_3, fraction_in_2, fraction_in_1)
    times_detected = c(3,2,1)
    input = c(d,d,d)
    method = c(i,i,i)
    df2 = data.frame(input,method,fraction,times_detected)
    
    if (count == 0){
      counts_df = table
      percentage_df = df2
      count = count + 1
    } else {
      counts_df = rbind(counts_df,table)
      percentage_df = rbind(df2,percentage_df)
    }
  }
}

# Plot the fraction of SNPs detected 1, 2 or 3 times for each method
percentage_df$times_detected = as.character(as.numeric(percentage_df$times_detected))
percentage_df$input = as.character(as.numeric(percentage_df$input))
percentage_df$method = factor(percentage_df$method, levels=c("duplicate_reads_kept","bbtools","picard"))
percentage_df$d <- factor(percentage_df$input, levels = c(10, 100, 1000, 10000, 1e+05, 1e+06))

ggplot(percentage_df, aes(x=input, y=fraction, color=times_detected, fill=times_detected))+geom_bar(position="stack", stat = "identity")+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+labs(x="input RNA copies",y="fraction of SNPs")+theme(plot.title=element_text(size=16))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title.y=element_text(size=16, vjust=0.5))+theme(axis.title.x=element_text(size=16, vjust=0.5))+theme(axis.text=element_text(size=16, colour="black"))+theme(axis.text.x=element_text(hjust=1, angle=45))+theme(legend.text=element_text(size=16))+theme(legend.title=element_text(size=16, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(legend.key.size=unit(0.7, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+facet_wrap(~method, scales="free")+scale_y_continuous(breaks=c(0,20,40,60,80,100), limits=c(0,100))+scale_colour_manual(values=c("grey80","black","grey40"))+scale_fill_manual(values=c("grey80","black","grey40"))


#### FIGURE 7A:-D CORRELATION PLOTS  ########################################

# read in dataframes and combine
bbtools=read.table("~/Documents/Graduate_School/CA04_SNP_detection/revisions/data_reanalysis_with_illumina_pipeline/duplicate_read_removal/combined_vcfs_amplicon_dedupe_180107.txt", header=T, stringsAsFactors=F, na.string=c(NA,"NA","na"), sep="\t", col.names=c("reference_position","reference","allele","frequency","input","replicate"))
bbtools$method = "bbtools"

picard=read.table("~/Documents/Graduate_School/CA04_SNP_detection/revisions/data_reanalysis_with_illumina_pipeline/duplicate_read_removal/combined_vcfs_amplicon_picard_180107.txt", header=T, stringsAsFactors=F, na.string=c(NA,"NA","na"), sep="\t", col.names=c("reference_position","reference","allele","frequency","input","replicate"))
picard$method = "picard"

original=read.table("~/Documents/Graduate_School/CA04_SNP_detection/revisions/data_reanalysis_with_illumina_pipeline/combined_vcfs_amplicon_171221.txt", header=T, stringsAsFactors=F, na.string=c(NA,"NA","na"), sep="\t", col.names=c("reference_position","reference","allele","frequency","input","replicate"))
original$method = "duplicate_reads_kept"

# combined dataframes
data = rbind(original, bbtools, picard)
data <- data[c(1,4:7)]

# convert to long format again
df1 = data[data$method=="duplicate_reads_kept",]
df2 = data[data$method=="bbtools",]
df3 = data[data$method=="picard",]

df = merge(df1, df2, by=c("reference_position","input","replicate"), all=TRUE)
df$method.x = "duplicate_reads_kept"
df$method.y = "bbtools"
df[is.na(df)] <- 0

df = merge(df, df3, by=c("reference_position","input","replicate"), all=TRUE)
df$method.x = "duplicate_reads_kept"
df$method.y = "bbtools"
df$method = "picard"
df[is.na(df)] <- 0

# calculate regression between not subsampled and bbtools
r2.lm = lm(frequency.x~frequency.y, data=df)
r2.lm$residuals #get residuals
summary(r2.lm)

# calculate regression between not subsampled and picard
r2.lm = lm(frequency.x~frequency, data=df)
r2.lm$residuals #get residuals
summary(r2.lm)

df$input <- as.character(as.numeric(df$input))
df$input <- gsub("10000","1e+04",df$input)
df$input <- gsub("1000","1e+03",df$input)
df$input <- gsub("100","1e+02",df$input)
df$input <- gsub("10","1e+01",df$input)
df$input <- factor(df$input, levels=c("1e+01","1e+02","1e+03","1e+04","1e+05","1e+06","1e+07"))

# compare subsampled vs. not as a scatter plot (correlation)
ggplot(df, aes(x=frequency.x,y=frequency.y,color=input))+geom_point(size=2)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+theme(plot.title=element_text(size=15))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=15))+theme(axis.title=element_text(size=15, vjust=0.2))+theme(axis.text=element_text(size=15, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=15))+theme(legend.title=element_text(size=15, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+labs(x="SNP frequency (%)",y="SNP frequency (%) after bbtools")+scale_x_continuous(breaks=seq(0,100,10), limits=c(0,100))+scale_y_continuous(breaks=seq(0,100,10), limits=c(0,100))+scale_colour_manual(values=c(color1,color2,color3,color4,color5,color6))

# original vs picard (7b)
ggplot(df, aes(x=frequency.x,y=frequency,color=input))+geom_point(size=2)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+theme(plot.title=element_text(size=15))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=15))+theme(axis.title=element_text(size=15, vjust=0.2))+theme(axis.text=element_text(size=15, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=15))+theme(legend.title=element_text(size=15, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+labs(x="SNP frequency (%)",y="SNP frequency (%) after picard")+scale_x_continuous(breaks=seq(0,100,10), limits=c(0,100))+scale_y_continuous(breaks=seq(0,100,10), limits=c(0,100))+scale_colour_manual(values=c(color1,color2,color3,color4,color5,color6))

# only SNPs 1-10%
df2 = df[df$frequency <= 10,]
df2 = df[df$frequency.x <= 10,]
df2 = df[df$frequency.y <= 10,]

# calculate regression between not subsampled and bbtools
r2.lm = lm(frequency.x~frequency.y, data=df2)
r2.lm$residuals #get residuals
summary(r2.lm)

# calculate regression between not subsampled and picard
r2.lm = lm(frequency.x~frequency, data=df2)
r2.lm$residuals #get residuals
summary(r2.lm)

# compare subsampled vs. not as a scatter plot (correlation)
ggplot(df2, aes(x=frequency.x,y=frequency.y,color=input))+geom_point(size=2)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+theme(plot.title=element_text(size=15))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=15))+theme(axis.title=element_text(size=15, vjust=0.2))+theme(axis.text=element_text(size=15, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=15))+theme(legend.title=element_text(size=15, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+labs(x="SNP frequency (%)",y="SNP frequency (%) after bbtools")+scale_x_continuous(breaks=seq(0,10,1), limits=c(0,10))+scale_y_continuous(breaks=seq(0,10,1), limits=c(0,10))+scale_colour_manual(values=c(color1,color2,color3,color4,color5,color6))

# original vs picard
ggplot(df2, aes(x=frequency.x,y=frequency,color=input))+geom_point(size=2)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+theme(plot.title=element_text(size=15))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=15))+theme(axis.title=element_text(size=15, vjust=0.2))+theme(axis.text=element_text(size=15, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=15))+theme(legend.title=element_text(size=15, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+labs(x="SNP frequency (%)",y="SNP frequency (%) after picard")+scale_x_continuous(breaks=seq(0,10,1), limits=c(0,10))+scale_y_continuous(breaks=seq(0,10,1), limits=c(0,10))+scale_colour_manual(values=c(color1,color2,color3,color4,color5,color6))

