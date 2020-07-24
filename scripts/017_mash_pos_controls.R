

library(data.table)
library(tidyverse)
library(magrittr)
library(reshape)
library(splitstackshape)
library(dplyr)



basedir = "/Users/12705859/Desktop/metapigs_dry/"
setwd("/Users/12705859/Desktop/metapigs_dry/mash/")

######################################################################

# counts data 

no_reps_pos_controls <- read.csv(paste0(basedir,"pos_controls/no_reps_pos_controls.csv"), 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)

# remove .fa extension to match bins in checkm df 
no_reps_pos_controls$bin <- gsub(".fa","", no_reps_pos_controls$bin)

######################################################################

# mash data 


mash_output_all <- read_table2(paste0(basedir,"mash_output_all.csv"), col_names = FALSE)
NROW(mash_output_all)
head(mash_output_all)

mash_out <- mash_output_all %>%
  dplyr::filter(!X4==1)
NROW(mash_out)
head(mash_out)

sink(file = "....", 
     append = TRUE, type = c("output"))
paste0("Total number of mash comparisons : ",NROW(mash_output_all))
paste0("Total number of mash comparisons where p-value is not 1 : ",NROW(mash_out))

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

# # test:
# test <- df2 %>%
#   filter(isolation_source=="MockCommunity") %>%
#   filter(variable=="R3") %>%
#   dplyr::mutate(norm_value=value/sum(value))
# head(test)
# sum(test$norm_value)

#################################
# STEP 2.

# sum all the norm values that fall within same isolation_source,variable,bin

# not necessary because this was already done (dereplicated)
NROW(unique(paste0(df2$isolation_source,df2$variable,df2$bin)))==NROW(df2)

# ready to merge to mash output! 
##################################################################


# MERGE MASH with no_reps

unique(mash_out$isolation_source)
unique(df2$isolation_source)
NROW(mash_out)
NROW(df2)

# merge
df3 <- merge(mash_out,df2) 
NROW(df3) # why 378 instead of 432? because we removed some rows in the mash output that had p_value 1  
head(df3)




# positive control coliguard 
# has got 2 main strains: 
# L. salivarius (%90)
# L. plantarum (10%)

# lessons learned: 

# the less complex a sample is, the better it will show matches against specific strains, because
# their relative abundance will be higher
# reason why rarefaction is a thing ...

# matches that are below 10% are non specific 

# L. salivarius is most descriptive of ColiGuard, but also the best Lactobacillus genome we have to describe Protexin 
# L. rhamnosus is discriminative to describe Protexin alone 


pdf("mash_pos_controls.pdf")
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

