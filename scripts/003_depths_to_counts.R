
# old name: 7.R script                                    #
# bins without a cluster assigned will be "no_cluster"    #
# depth to counts                                         #
# removes samples replicates                              #
# removal of SD piglets                                   #
# saves no replicated bins                                #

library(data.table)
library(tidyverse)
library(magrittr)
library(reshape)
library(splitstackshape)
library(dplyr)

# PATHS
middle_dir = "/Users/12705859/metapigs_dry/middle_dir/"

# input files: 
# merged_all_clustered_wa_bins_with_cohorts.csv

# output files: 
# no_reps_all.csv

# upload input file
df <- read.csv(paste0(middle_dir,"merged_all_clustered_wa_bins_with_cohorts.csv"),
                                                      na.strings=c("","NA"),
                                                      check.names = FALSE,
                                                      header = TRUE)

# cohort names edits: 
df$cohort <- plyr::mapvalues(as.character(df$cohort), 
                                  from = c("Neomycin+D-scour","Neomycin+ColiGuard","D-scour"), 
                                  to = c("NeoD","NeoC","D-Scour"))

# bins that don't have secondary cluster assigned are filled with "no_cluster"
df <- df %<>% 
  dplyr::mutate(secondary_cluster = fct_explicit_na(secondary_cluster, na_level = "no_cluster"))

# get the column names before entering function
original_colnames <- colnames(df)

# transform depth values into counts
# where {counts = depth*median bin length/300} (300 is the read pair length) 
A <- function(x) x * median(df[,2]) / 300    # A <- function(x) x * df[,2] / 300 we used to multiply by the bin length but that would un-normalize data 
counts <- data.frame(df[1:5], apply(df[6:ncol(df)],2, A) )
counts
#return the original colnames to new dataframe
colnames(counts) <- original_colnames
#now we don't need binLen anymore we can remove this column
counts <- counts[,-2]

NROW(counts)

#make long (one value per row)

counts_long <- reshape::melt(counts, id=c("cohort", "pig", "bin", "secondary_cluster"))
NROW(counts_long)

# remove NAs (NAs were automatically created when pivoting large)
counts_long <- na.omit(counts_long)
NROW(counts_long)
head(counts_long)

# Merge the replicates (some samples are duplicate, same day and pig, so bins are also duplicate)
#split column "variable" using dot separator
#install.packages("splitstackshape")

NL2 <- cSplit(counts_long, "variable", ".")
NROW(NL2)
head(NL2)

# rename new columns
colnames(NL2)[colnames(NL2)=="variable_1"] <- "collection_date"
colnames(NL2)[colnames(NL2)=="variable_2"] <- "replicate"
#reformat to be recognized as dates
NL2$collection_date <- as.Date(NL2$collection_date, format = "%y-%m-%d")

# checking how many replicates we have
#return rows that have replicate==1, ==2, ==3, ==4
NLrep1 <- NL2[which(NL2$replicate == 1), ]
NLrep2 <- NL2[which(NL2$replicate == 2), ]
NLrep3 <- NL2[which(NL2$replicate == 3), ]
NLrep4 <- NL2[which(NL2$replicate == 4), ]
NROW(NLrep1)
NROW(NLrep2)
NROW(NLrep3)
NROW(NLrep4)
#concatenate the above dataframes
total <- rbind(NLrep1, NLrep2, NLrep3, NLrep4)
NROW(total)

# 1st de-replication (take mean out of replicates): 
NL2$replicate <- NULL
NROW(NL2)
no_reps <- stats::aggregate(NL2$value,by=list(cohort=NL2$cohort,
                                       pig=NL2$pig,bin=NL2$bin, 
                                       collection_date=NL2$collection_date,
                                       secondary_cluster=NL2$secondary_cluster),data=NL2,FUN=mean)
NROW(no_reps)

#replace column name x with "value"
names(no_reps)[names(no_reps) == 'x'] <- 'value'


# checking if mean of replicates is working fine
a <- dplyr::filter(NLrep1, pig == "29951", bin == "bins.1.fa", collection_date == "2017-01-31")
a
b <- dplyr::filter(NLrep2, pig == "29951", bin == "bins.1.fa", collection_date == "2017-01-31")
b
c <- dplyr::filter(no_reps, pig == "29951", bin == "bins.1.fa", collection_date == "2017-01-31")
c

# is the value of c the mean of a and b ?
(a$value+b$value)/2==c$value
# yes

########################################################

# Important change 20191108: because some piglets could not be sampled on given sampling date
# (or low quality sample was taken, so sample was taken again the next day),
# they were re-sampled the morning after (before treatment was administered)
# for this reason we are going to make groups of sampling dates as follows: 

# duplicate date column 
head(no_reps)
no_reps$date <- no_reps$collection_date 
head(no_reps)

from = c("2017-01-30",
         "2017-01-31","2017-02-01",
         "2017-02-03",
         "2017-02-06","2017-02-07","2017-02-08",
         "2017-02-10",
         "2017-02-14",
         "2017-02-16","2017-02-17",
         "2017-02-21",
         "2017-02-24", 
         "2017-02-28",
         "2017-03-03",
         "2017-03-06","2017-03-07","2017-03-08","2017-03-09","2017-03-10")

to = c("tM",
       "t0","t0",
       "t1",
       "t2","t2","t2",
       "t3",
       "t4",
       "t5","t5",
       "t6",
       "t7", 
       "t8",
       "t9",
       "t10","t10","t10","t10","t10")

# replace collection dates (date format) with groups of collection dates (character format)
no_reps$date <- plyr::mapvalues(as.character(no_reps$date), from, to)
unique(no_reps$date)

# select columns to keep
no_reps <- no_reps %>%
  dplyr::select(cohort,pig,bin,date,secondary_cluster,value)

######################

# 2nd de-replication: now there should be again replicates because we joined some dates. 
# (take mean out of replicates)
NROW(NL2)
no_reps$sampleID <- paste0(no_reps$date,"_",no_reps$cohort,"_",no_reps$pig,"_",no_reps$secondary_cluster,no_reps$bin)
NROW(no_reps) == length(unique(no_reps$sampleID))
NROW(no_reps) # 362863
length(unique(no_reps$sampleID)) # 361054
# here above we can see the amount of replicates (same pig, same date)
no_reps2 <- stats::aggregate(no_reps$value,by=list(cohort=no_reps$cohort,
                                            pig=no_reps$pig,
                                            bin=no_reps$bin,
                                            secondary_cluster=no_reps$secondary_cluster,
                                            date=no_reps$date),
                      data=no_reps,FUN=mean)

#replace column name x with "value"
names(no_reps2)[names(no_reps2) == 'x'] <- 'value'

# sanity check: 
# number of unique pigs
length(unique(paste0(no_reps2$pig)))
# 168
# number of bins in total across all pigs
length(unique(paste0(no_reps2$pig,no_reps2$bin)))
# 50787
# number of unique pig,date,bin
length(unique(paste0(no_reps2$pig,no_reps2$date, no_reps2$bin)))
# 361054


# final number of samples per date and cohort : 
az <- no_reps2 %>%
  dplyr::select(date, pig, cohort)
bz <- unique(az)
cz2 <- data.frame(table(bz$cohort, bz$date))
# subset to non zero containing groups 
dz2 <- cz2 %>% dplyr::filter(Freq != 0 )
head(dz2)

# Write out to play with it (this is with secondary_clusters assigned only)
# fwrite(x = tot_counts_dereplicated, file = "~/Desktop/bins_clustering_parsing_dataframes/tot_counts_dereplicated.csv")

# remove piglets with SD (swine dysentery): 
no_reps3 <- no_reps2 %>%
  filter(!pig=="29665"&
           !pig=="29865"&
           !pig=="29702")

# Write out to play with it (this is everything)
fwrite(x = no_reps3, file = paste0(middle_dir,"no_reps_all.csv"))



