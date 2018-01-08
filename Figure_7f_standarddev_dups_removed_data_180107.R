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


#### FIGURE 7F: STANDARD DEVIATION VS INPUT DNA COPIES ################

# read in the combined dataframe. 
data=read.table("~/Documents/Graduate_School/CA04_SNP_detection/revisions/data_reanalysis_with_illumina_pipeline/duplicate_read_removal/combined_vcfs_amplicon_dedupe_and_picard_180107.txt", header=T, stringsAsFactors=F, na.string=c(NA,"NA","na"), sep="\t", col.names=c("reference_position","reference","allele","frequency","input","replicate", "method"))

original=read.table("~/Documents/Graduate_School/CA04_SNP_detection/revisions/data_reanalysis_with_illumina_pipeline/combined_vcfs_amplicon_171221.txt", header=T, stringsAsFactors=F, na.string=c(NA,"NA","na"), sep="\t", col.names=c("reference_position","reference","allele","frequency","input","replicate"))
original$method = "duplicate_reads_kept"

data = rbind(data,original)

# for this, need to convert the input to a character
data$input = as.character(as.numeric(data$input))
data$input <- gsub("10000","1e+04",data$input)
data$input <- gsub("1000","1e+03",data$input)
data$input <- gsub("100","1e+02",data$input)
data$input <- gsub("10","1e+01",data$input)


# Subset the dataframe to remove nonduplicated sites. 
# 1st step creates a list of sites that are duplicated; 2nd step subsets data to include only rows with those nonunique sites
nonunique <- subset(data$reference_position, duplicated(data$reference_position), )
nonunique = sort(nonunique)
nonunique = unique(nonunique)
data_nonunique <- data[data$reference_position %in% nonunique, ]

# Now calculate the mean and standard deviation for the nonunique data
dilution_list = unique(data_nonunique$input)
method_list = unique(data_nonunique$method)
count = 0

data_nonunique$frequency = as.numeric(as.character(data_nonunique$frequency))

for (m in method_list)
{
  for (d in dilution_list)
  {
    sites <- data_nonunique[data_nonunique$input == d & data_nonunique$method == m,]
    site_list = unique(sites$reference_position)
    
    for (i in site_list)
    {
      df <- data_nonunique[data_nonunique$reference_position == i & data_nonunique$input == d & data_nonunique$method == m,]
      if (nrow(df) == 1){
        all_values = list(df$frequency, 0, 0)
        times_detected = 1
        non_mean_frequency = df$frequency
      } else if (nrow(df) == 2) {
        all_values = list(df$frequency[1], df$frequency[2], 0)
        times_detected=2
        non_mean_frequency = NA
      }
      else if (nrow(df) == 3) {
        all_values = list(df$frequency[1], df$frequency[2], df$frequency[3])
        times_detected = 3
        non_mean_frequency = NA
      }
      # now calculate the mean, variance, and predicted variance
      mean = mean(as.numeric(all_values))
      st_dev = sd(as.numeric(all_values))
      replicate1 = all_values[[1]]
      replicate2 = all_values[[2]]
      replicate3 = all_values[[3]]
      df2 = data.frame(i, replicate1,replicate2,replicate3,mean, st_dev, times_detected, non_mean_frequency, d, m)
      
      if (count == 0){
        mean_variance_df = df2
        count = count + 1
      } else {
        mean_variance_df = rbind(mean_variance_df, df2)
      }
    }
  }
}

ggplot(mean_variance_df, aes(x=d,y=log10(st_dev),color=d, shape=m))+geom_point(size=2, position=position_dodge(width=0.5))+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+theme(plot.title=element_text(size=13))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title=element_text(size=13, vjust=0.2))+theme(axis.text=element_text(size=11, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=11))+theme(legend.title=element_text(size=13, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+labs(x="log10(dilution group)",y="log10(standard deviation)")+scale_y_continuous(breaks=seq(-2,2,1), limits=c(-2,2), minor_breaks=c(-2,-1,0,1,2))+scale_colour_manual(values=c(color1,color2,color3,color4,color5,color6))+scale_shape_manual(values=c(4,1,19))

