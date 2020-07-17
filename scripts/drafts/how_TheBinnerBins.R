

# binning: a visual example 
# this script gives a visual example of 
# how the differential abundance of contigs through time
# determines which contigs fall within the same bin 
# not to forget: 
# metabat2 also uses coverage and GC content (or k-mer frequency) info 
# to bin! 

library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(ggforce)

setwd("~/Desktop/metapigs_dry/drafts_out")
depth_clean <- read_csv("~/Desktop/bins_clustering_parsing_dataframes/out_new_test/14159/metabat/depth_clean.csv")
bins_to_contigs <- read_csv("~/Desktop/bins_clustering_parsing_dataframes/out_new_test/14159/metabat/bins_to_contigs.csv")
colnames(bins_to_contigs) <- c("pig","bin","contigName")
bins_to_contigs <- bins_to_contigs[,2:3]

z <- inner_join(bins_to_contigs,depth_clean)

z <- z[, -grep(".var", colnames(z))]

z <- z %>%
  pivot_longer(cols=ends_with(".bam"), names_to = "date", values_to = "depth") %>%
  group_by(date) %>%
  mutate(depth=depth/sum(depth))

# # correction for contig length
# z <- z %>%
#   mutate(depth=depth/contigLen)

z0 <- z[1:20000,]


all_20k_contigs <- ggplot(z0, aes(x=date,y=depth, group=contigName, fill=bin,color=bin))+
  geom_point()+
  geom_line()+
  theme(legend.position="top",
        axis.text.x=element_text(angle=30, hjust=1))

facet_20k_contigs <- ggplot(z0, aes(x=date,y=depth, group=contigName, fill=bin,color=bin))+
  geom_point()+
  geom_line()+
  theme(legend.position="top",
        axis.text.x=element_blank())+
  facet_wrap(~bin,scales = "free_y")

both <- ggarrange(all_20k_contigs,facet_20k_contigs,nrow=2,common.legend = TRUE)

pdf("Desktop/metapigs_dry/how_binning_works.pdf")
both 
dev.off()


# looking at the one contig that does something different 
z1 <- z
z2 <- z1 %>% filter(bin=="bins.100.fa")
ggplot(z2, aes(x=date,y=depth, group=contigName, fill=bin,color=bin))+
  geom_point()+
  geom_line()+
  theme(legend.position="top",
        axis.text.x=element_blank())+
  facet_wrap(~bin,scales = "free_y")


# split by bin
multiple_DFs <- split( z1 , f = z1$bin ,drop = TRUE)

pdf("bins_contigs.pdf")
for (single_DF in multiple_DFs) {
  
  single <- as.data.frame(single_DF)
  
  print(ggplot(single, aes(x=date,y=depth, group=contigName, fill=bin,color=bin))+
    geom_point()+
    geom_line()+
    theme(legend.position="top",
          axis.text.x=element_blank())+
    facet_wrap(~bin,scales = "free_y"))
}
dev.off()