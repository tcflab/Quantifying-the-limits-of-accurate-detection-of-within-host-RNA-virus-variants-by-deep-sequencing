###############################################################################
# This includes the code for generating figures for the CA04 illumina manuscript
###############################################################################
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
color4 <- rgb(117,11,0, maxColorValue=255)
color5 <- rgb(207,19,0, maxColorValue=255)
color6 <- rgb(245,89,74, maxColorValue=255)
color7 <- rgb(255,255,255, maxColorValue=255)

###############################################################################
# Figure 2: coverage plots and quality score plots for random hexamer dataset
###############################################################################

####### FIGURE 2B: QSCORES PLOT FOR RANDOM HEXAMER DATA ###########################

data = read.table("~/Documents/Graduate_School/CA04_SNP_detection/revisions/data_reanalysis_with_illumina_pipeline/coverage_depth_and_quality_scores/combined.aqhist.randomhex_170831.txt", header=T, sep="\t", col.names=c("id","quality","count","fraction"))

total_reads = sum(data$count)
data$total_fraction = (data$count/total_reads)*100

data = data[data$id != "1e7_A_S1" ,] 
data = data[data$id != "1e7_B_S2" ,] 
data = data[data$id != "1e7_C_S3" ,] 


# cast the data into long format
data = dcast(data, quality~id, value.var="count", fun.aggregate=sum)

data$total = rowSums(data)
data$corrected = data$total - data$quality
data$corrected_fraction = (data$corrected/total_reads)*100

# make a bar plot
ggplot(data, aes(x=quality, y=corrected_fraction))+geom_bar(position="dodge", stat = "identity")+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+labs(x="average q-score",y="fraction of reads")+theme(plot.title=element_text(size=16))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title.y=element_text(size=16, vjust=0.5))+theme(axis.title.x=element_text(size=16, vjust=0.5))+theme(axis.text=element_text(size=16, colour="black"))+theme(axis.text.x=element_text(hjust=1, angle=45))+theme(legend.text=element_text(size=16))+theme(legend.title=element_text(size=16, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(legend.key.size=unit(0.7, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+scale_x_continuous(limits=c(25,40), breaks=seq(25,40,1))+scale_y_continuous(breaks=seq(0,35,5), limits=c(0,35))



####### FIGURE 2D: COVERAGE PLOT FOR RANDOM HEXAMER DATA ##########################

# read in combined pileup files
data=read.table("~/Documents/Graduate_School/CA04_SNP_detection/revisions/data_reanalysis_with_illumina_pipeline/coverage_depth_and_quality_scores/combined2_randomhex_170830.pileup", header=F, sep="\t", col.names=c("dilution","replicate","reference","site","ref_base","coverage"))

# the sites where there is 0 coverage are just omitted in the pileup file and therefore will give skewed values for the average coverage. So I need to loop through each site, replicate, and dilution and add lines with 0 coverage for the missing ones
## NOTE THIS TAKES FOREVER BUT IT DOES WORK

data = data[data$dilution != 1e+07 ,] 

sites_list = (1:1701)
dilution_list = unique(data$dilution)
replicate_list = unique(data$replicate)
count = 0

for (dilution in dilution_list)
{
  for (replicate in replicate_list)
  {
    for (site in sites_list)
    {
      df <- data[data$site == site & data$dilution == dilution & data$replicate == replicate,]
      
      if (nrow(df) == 0){
        coverage = 0
        row_to_append = data.frame(dilution,replicate,"CA04_HA_GQ117044",site,"N",coverage)
      } else {
        coverage = df$coverage
        row_to_append = data.frame(dilution,replicate,"CA04_HA_GQ117044",site,"N",coverage)
      }
      
      if (count == 0){
        df2 = row_to_append
        count = count + 1
      } else {
        df2 = rbind(df2,row_to_append)
      }
    }
  }
}


mean(df2$coverage)
sd(df2$coverage)

#df_1e7 = df2[df2$dilution == 1e7,]
df_1e6 = df2[df2$dilution == 1e6,]
df_1e5 = df2[df2$dilution == 1e5,]
df_1e4 = df2[df2$dilution == 1e4,]
df_1e3 = df2[df2$dilution == 1000,]
df_1e2 = df2[df2$dilution == 100,]

mean(df_1e7$coverage)
sd(df_1e7$coverage)
mean(df_1e6$coverage)
sd(df_1e6$coverage)
mean(df_1e5$coverage)
sd(df_1e5$coverage)
mean(df_1e4$coverage)
sd(df_1e4$coverage)
mean(df_1e3$coverage)
sd(df_1e3$coverage)
mean(df_1e2$coverage)
sd(df_1e2$coverage)


## plot the mean and standard deviation for each dilution on 1 plot (averaging across replicates)

data = df2
site_list = unique(data$site)
dilution_list = unique(data$dilution)
count = 0

for (d in dilution_list)
{
  for (s in site_list)
  {
    df <- data[data$site == s & data$dilution == d,]
    mean_coverage = mean(df$coverage)
    sd_coverage = sd(df$coverage)
    df2 = data.frame(d,s,mean_coverage,sd_coverage)
    
    if (count == 0){
      coverage_df = df2
      count = count + 1
    } else {
      coverage_df = rbind(coverage_df,df2)
    }
  }
}

coverage_df["lower"] = coverage_df$mean_coverage - coverage_df$sd_coverage
coverage_df["upper"] = coverage_df$mean_coverage + coverage_df$sd_coverage
coverage_df[is.na(coverage_df)] <- 0

ggplot(coverage_df, aes(x=s, y=mean_coverage))+
  facet_wrap(~d, scales="free")+
  theme(panel.grid.major=element_line(colour=NA,size=NA))+
  theme(panel.grid.minor=element_line(colour=NA,size=NA))+
  labs(x="nucleotide site",y="coverage depth")+
  theme(plot.title=element_text(size=13))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(axis.line.x=element_line(colour="black"))+
  theme(axis.line.y=element_line(colour="black"))+
  theme(strip.text.x=element_text(size=11))+
  theme(axis.title=element_text(size=13, vjust=0.2))+
  theme(axis.text=element_text(size=11, colour="black"))+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(legend.text=element_text(size=11))+
  theme(legend.title=element_text(size=13, face="plain"))+
  theme(panel.margin=unit(0.9, "lines"))+
  theme(legend.key.size=unit(0.5, "cm"))+
  theme(panel.background=element_rect(fill=NA))+
  theme(legend.key=element_rect(fill=NA))+
  scale_y_log10(limits=c(1,10000), breaks=c(1,10,100,1000,10000))+
  scale_x_continuous(breaks=seq(0,1800,300), limits=c(0,1800))+
  geom_ribbon(aes(x=s, ymin=lower, ymax=upper), fill="grey80", linetype=0, alpha=0.6)+
  geom_line(size=0.7)
