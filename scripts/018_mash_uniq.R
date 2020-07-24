library(dplyr)
library(splitstackshape)
library(readr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(pheatmap)
library(treemap)
library(data.table)
library(tidyverse)
library(magrittr)
library(reshape)
library(EnvStats)



basedir = "~/Desktop/metapigs_dry/"
setwd("~/Desktop/metapigs_dry/mash")

######################################################################

# counts data 

no_reps_all <- read.csv(paste0(basedir,"no_reps_all.csv"), 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)


# remove .fa extension to match bins in checkm df 
no_reps_all$bin <- gsub(".fa","", no_reps_all$bin)
head(no_reps_all)
NROW(no_reps_all)

colnames(no_reps_all)[colnames(no_reps_all)=="pig"] <- "isolation_source"

# metadata to get cohort info only 
cohorts <- no_reps_all %>%
  dplyr::select(isolation_source,cohort) %>%
  distinct()

######################################################################


# counts data  - pos controls

no_reps_pos_controls <- read.csv(paste0(basedir,"pos_controls/no_reps_pos_controls.csv"), 
                                 na.strings=c("","NA"),
                                 check.names = FALSE,
                                 header = TRUE)

# remove .fa extension to match bins in checkm df 
no_reps_pos_controls$bin <- gsub(".fa","", no_reps_pos_controls$bin)


######################################################################

# mash data 

merged_uniq <- read_table2(paste0(basedir,"merged_uniq.txt"), col_names = FALSE)
NROW(merged_uniq)
head(merged_uniq)

mash_out <- merged_uniq %>%
  dplyr::filter(!X4==1)
NROW(mash_out)
head(mash_out)

sink(file = "merged_uniq.txt", 
     append = FALSE, type = c("output"))
paste0("Total number of mash comparisons : ",NROW(merged_uniq))
paste0("Total number of mash comparisons where p-value is not 1 : ",NROW(mash_out))
sink()

mash_out <- cSplit(mash_out,"X1","/") 
mash_out <- cSplit(mash_out,"X2","_")
mash_out <- cSplit(mash_out,"X5","/")

# remove .fa extension to match bins in checkm df 
mash_out$X1_7 <- gsub(".fa","", mash_out$X1_7)
mash_out$X2_3 <- gsub("fastq","", mash_out$X2_3)

mash_out <- cSplit(mash_out,"X1_7","_")

mash_out <- mash_out %>%
  dplyr::mutate(match_percentage=(X5_1/X5_2)*100 ) %>%
  dplyr::select(X2_3,X1_7_1,X1_7_2,X3,X4,match_percentage) 

colnames(mash_out) <- c("reference","isolation_source", "bin", "mash_dist","p_value","match_percentage")

summary(mash_out$match_percentage)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1000  0.1000  0.1000  0.9994  0.1000 71.3000 

summary(mash_out$mash_dist)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.008732 0.295981 0.295981 0.272260 0.295981 0.295981 

summary(mash_out$p_value)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000e+00 4.087e-05 1.903e-04 1.613e-04 2.557e-04 4.880e-04

mash_out1 <- mash_out %>%
  group_by(isolation_source,reference) %>%
  filter(p_value<1e-300) 


######################################################################

# lib size normalisation

df1 <- no_reps_pos_controls %>%
  dplyr::select(cohort,isolation_source,bin,variable,value)

#################################
# STEP 1.

# normalization for library size 
df2 <- df1 %>%
  dplyr::group_by(isolation_source,variable) %>% # variable is like date , a pig sample is: pig+date
  dplyr::mutate(norm_value=value/sum(value))    # a positive control sample is: poscontrol+replicate
head(df2)

#################################
# STEP 2.

# sum all the norm values that fall within same isolation_source,variable,bin

# not necessary because this was already done (dereplicated)
NROW(unique(paste0(df2$isolation_source,df2$variable,df2$bin)))==NROW(df2)

# ready to merge to mash output! 
##################################################################

# MERGE MASH with no_reps

# merge
df3 <- merge(mash_out1,df2) 
NROW(df3) # why 378 instead of 432? because we removed some rows in the mash output that had p_value 1  
head(df3)


pdf("mash_uniq_pos_controls.pdf")
df3 %>% 
  filter(!variable=="R9") %>%
  filter(!isolation_source=="MockCommunity") %>%
  ggplot(., aes(x=variable,y=norm_value,fill=reference))+
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~isolation_source)
df3 %>% 
  filter(!variable=="R9") %>%
  filter(!isolation_source=="MockCommunity") %>%
  ggplot(., aes(x=variable,y=norm_value,fill=isolation_source))+
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~reference)
dev.off()



######################################################################
######################################################################
######################################################################


# further filtering for piglets



# cohorts


######################################################################

# lib size normalisation

df1 <- no_reps_all %>%
  dplyr::select(cohort,isolation_source,bin,date,value)

#################################
# STEP 1.

# normalization for library size 
df2 <- df1 %>%
  dplyr::group_by(isolation_source,date) %>% # variable is like date , a pig sample is: pig+date
  dplyr::mutate(norm_value=value/sum(value))    # a positive control sample is: poscontrol+replicate
head(df2)

#################################
# STEP 2.

# sum all the norm values that fall within same isolation_source,variable,bin

# not necessary because this was already done (dereplicated)
NROW(unique(paste0(df2$isolation_source,df2$date,df2$bin)))==NROW(df2)

# ready to merge to mash output! 
##################################################################

# MERGE MASH with no_reps

unique(mash_out1$reference)

# merge
df3 <- merge(mash_out1,df2) 
NROW(df3) # why 378 instead of 432? because we removed some rows in the mash output that had p_value 1  
head(df3)

unique(df3$reference)

df3 <- df3 %>% 
  filter(!cohort=="Mothers") %>% 
  filter(date=="t0"|date=="t2"|
           date=="t4"|date=="t6"|
           date=="t8")


gg <- df3 %>% 
  ggplot(., aes(x=date,y=norm_value,fill=reference))+
  geom_bar(position="dodge", stat="identity")+
  stat_n_text(size = 2, angle = 90, y.pos = 0.037) +
  facet_wrap(~cohort)

# manually picking colours to match colors in the plot of the pos controls
cols <- c("LD." = "#F25E5A", "LS." = "#B95EFF")
gg2 <- gg + scale_fill_manual(values = cols)

pdf("mash_uniq_cohorts.pdf")
gg2
dev.off()


######################################################################
######################################################################
######################################################################