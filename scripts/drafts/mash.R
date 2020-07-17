library(dplyr)
library(splitstackshape)
library(readr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(pheatmap)


basedir = "~/Desktop/metapigs_dry/"

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

# mash data 

mash_output_all <- read_table2(paste0(basedir,"mash_output_all.csv"), col_names = FALSE)
NROW(mash_output_all)
head(mash_output_all)

mash_out <- mash_output_all %>%
  dplyr::filter(!X4==1)
NROW(mash_out)
head(mash_out)

sink(file = "mash_output.txt", 
     append = FALSE, type = c("output"))
paste0("Total number of mash comparisons : ",NROW(mash_output_all))
paste0("Total number of mash comparisons where p-value is not 1 : ",NROW(mash_out))
sink()

mash_out <- cSplit(mash_out,"X1","/") 
mash_out <- cSplit(mash_out,"X2","/")
mash_out <- cSplit(mash_out,"X5","/")

# remove .fa extension to match bins in checkm df 
mash_out$X2_7 <- gsub(".fa","", mash_out$X2_7)
mash_out$X1_7 <- gsub(".contigs.fasta","", mash_out$X1_7)

mash_out <- cSplit(mash_out,"X2_7","_")

mash_out <- mash_out %>%
  dplyr::mutate(match_percentage=(X5_1/X5_2)*100 ) %>%
  dplyr::select(X1_7,X2_7_1,X2_7_2,X3,X4,match_percentage)

colnames(mash_out) <- c("reference","isolation_source", "bin", "mash_dist","p_value","match_percentage")

summary(mash_out$match_percentage)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.100   0.100   0.100   1.095   0.100  95.100 

summary(mash_out$mash_dist)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001211 0.295981 0.295981 0.274610 0.295981 0.295981

######################################################################


# lib size normalisation

df1 <- no_reps_all %>%
  dplyr::select(cohort,isolation_source,bin,date,value)

#################################
# STEP 1.

# normalization for library size 
df2 <- df1 %>%
  dplyr::group_by(isolation_source,date) %>% # a pig sample is: pig+date
  dplyr::mutate(norm_value=value/sum(value))    
head(df2)

# # test:
# test <- df2 %>%
#   filter(isolation_source=="29912") %>%
#   filter(date=="t0") %>%
#   dplyr::mutate(norm_value=value/sum(value))
# head(test)
# sum(test$norm_value)

#################################
# STEP 2.

# sum all the norm values that fall within same isolation_source,variable,bin

# not necessary because this was already done (dereplicated)
NROW(unique(paste0(df2$isolation_source,df2$date,df2$bin)))==NROW(df2)

# ready to merge to mash output! 
##################################################################


# MERGE MASH with no_reps

unique(mash_out$isolation_source)
unique(df2$isolation_source)
NROW(mash_out)
NROW(df2)

# merge
df3 <- merge(mash_out,df2) %>%
  dplyr::filter(!date=="tM") %>%
  dplyr::filter(!date=="t10")
NROW(df3) 
head(df3)

# add cohorts info 
df4 <- inner_join(df3,cohorts)
NROW(df4) 

unique(df4$reference)


pdf("mash_cohorts.pdf")
df4 %>% 
  filter(reference=="Lactobacillus_salivarius_140318") %>%
  ggplot(data=., mapping=aes(x=date, y=sqrt(norm_value))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, size=0.1)+
  facet_wrap(~cohort) +
  ggtitle("Lactobacillus_salivarius_140318")
df4 %>% 
  filter(reference=="Lactobacillus_salivarius_140318") %>%
  ggplot(data=., mapping=aes(x=date, y=log(norm_value))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, size=0.1)+
  facet_wrap(~cohort) +
  ggtitle("Lactobacillus_salivarius_140318")
df4 %>% 
  filter(reference=="Lactobacillus_rhamnosus_NCIMB8010_ATCC7469") %>%
  ggplot(data=., mapping=aes(x=date, y=sqrt(norm_value))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, size=0.1)+
  facet_wrap(~cohort) +
  ggtitle("Lactobacillus_rhamnosus_NCIMB8010_ATCC7469")
df4 %>% 
  filter(reference=="Lactobacillus_rhamnosus_NCIMB8010_ATCC7469") %>%
  ggplot(data=., mapping=aes(x=date, y=log(norm_value))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, size=0.1)+
  facet_wrap(~cohort) +
  ggtitle("Lactobacillus_rhamnosus_NCIMB8010_ATCC7469")
df4 %>% 
  filter(reference=="Lactobacillus_plantarum_140318") %>%
  ggplot(data=., mapping=aes(x=date, y=sqrt(norm_value))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, size=0.1)+
  facet_wrap(~cohort) +
  ggtitle("Lactobacillus_plantarum_140318")
df4 %>% 
  filter(reference=="Lactobacillus_plantarum_140318") %>%
  ggplot(data=., mapping=aes(x=date, y=log(norm_value))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, size=0.1)+
  facet_wrap(~cohort) +
  ggtitle("Lactobacillus_plantarum_140318")
df4 %>% 
  filter(reference=="Lactobacillus_plantarum_NCIMB11974_ATCC14917") %>%
  ggplot(data=., mapping=aes(x=date, y=sqrt(norm_value))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, size=0.1)+
  facet_wrap(~cohort) +
  ggtitle("Lactobacillus_plantarum_NCIMB11974_ATCC14917")
df4 %>% 
  filter(reference=="Lactobacillus_plantarum_NCIMB11974_ATCC14917") %>%
  ggplot(data=., mapping=aes(x=date, y=log(norm_value))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, size=0.1)+
  facet_wrap(~cohort) +
  ggtitle("Lactobacillus_plantarum_NCIMB11974_ATCC14917")
df4 %>% 
  filter(reference=="Lactobacillus_acidophilus_NCIMB701748_ATCC4356") %>%
  ggplot(data=., mapping=aes(x=date, y=sqrt(norm_value))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, size=0.1)+
  facet_wrap(~cohort) +
  ggtitle("Lactobacillus_acidophilus_NCIMB701748_ATCC4356")
df4 %>% 
  filter(reference=="Lactobacillus_acidophilus_NCIMB701748_ATCC4356") %>%
  ggplot(data=., mapping=aes(x=date, y=log(norm_value))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, size=0.1)+
  facet_wrap(~cohort) +
  ggtitle("Lactobacillus_acidophilus_NCIMB701748_ATCC4356")
dev.off()
