#upload all libraries
library(base)
library(data.table)
library(dplyr)
library(stringr)
library(utils)
library(tidyverse)
library(plyr)
library(splitstackshape)
library(reshape)

#library(magrittr)



out.dir = "/shared/homes/12705859/contig_abundances" # on HPC
#out.dir = "/Users/12705859/Desktop/bins_clustering_parsing_DFs/contig_abundances" # on local


df0 <- read.table(file = file.path(out.dir,"merged_all_wa_contigs_with_cohorts.csv"), header = TRUE, sep = ",", row.names = NULL)

df1 <- df0

# cohort names edits: 
df1$cohort <- plyr::mapvalues(as.character(df1$cohort), 
                             from = c("Neomycin+D-scour","Neomycin+ColiGuard","D-scour"), 
                             to = c("NeoD","NeoC","D-Scour"))

# remove totAvgDepth column
df1 <- df1 %>% 
  dplyr::select(!totalAvgDepth)

# move bin column as first (might be useful at later stages)
df1 <- df1 %>%
  dplyr::select(bin, everything())

# # get the column names before entering function
# original_colnames <- colnames(df1)
# 
# # artificially inflate the contig lengths (it will help get counts away from zeros)
# df1 <- df1 %>%
#   dplyr::mutate(contigLen=contigLen*1000000)

# transform depth values into counts
# where {counts = depth*median contig length/300} (300 is the read pair length) 
A <- function(x) x * median(df1[,3]) / 300 
counts <- data.frame(df1[1:5], apply(df1[6:ncol(df1)],2, A) )

# now we don't need contigLen anymore we can remove this column
counts <- counts %>% 
  dplyr::select(!contigLen)

# make long (one value per row)
counts_long <- reshape::melt(counts, id=c("bin", "cohort", "pig","contigName"))

# remove NAs (NAs were automatically created when pivoting wide)
counts_long <- na.omit(counts_long)



########################################################

# First de-replication : 

# some info on contigs
# hist(counts_long$value, xlim = c(0,1000000), breaks=1000)
# summary(counts_long$value)

# Merge the replicates (some samples are duplicate, same day and pig)
# split column "variable" using dot separator
NL2 <- cSplit(counts_long, "variable", ".")
NL2$variable_1 <- as.numeric(gsub( "X", "", NL2$variable_1))
NL2$collection_date <- paste0(NL2$variable_1,"-",NL2$variable_2,"-",NL2$variable_3)
colnames(NL2)[colnames(NL2)=="variable_4"] <- "replicate"

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

# remove replicate column
NL2$replicate <- NULL

no_reps_contigs <- stats::aggregate(NL2$value,by=list(cohort=NL2$cohort,
                                              pig=NL2$pig,bin=NL2$bin, contigName=NL2$contigName,
                                              collection_date=NL2$collection_date),data=NL2,FUN=mean)

#replace column name x with "value"
names(no_reps_contigs)[names(no_reps_contigs) == 'x'] <- 'value'

##### TEST: 
# checking if mean of replicates is working fine
a <- dplyr::filter(NLrep1, pig == "29951", contigName == "k141_1000014", collection_date == "2017-01-31")
a
b <- dplyr::filter(NLrep2, pig == "29951", contigName == "k141_1000014", collection_date == "2017-01-31")
b
c <- dplyr::filter(no_reps_contigs, pig == "29951", contigName == "k141_1000014", collection_date == "2017-01-31")
c

# is the value of c the mean of a and b ?
c$value==(a$value+b$value)/2
# yes
##### END OF TEST. 

########################################################

# Second de-replication : 

# Important change 20191108: because some piglets could not be sampled on given sampling date
# (or low quality sample was taken, so sample was taken again the next day),
# they were re-sampled the morning after (before treatment was administered)
# for this reason we are going to make groups of sampling dates as follows: 

# duplicate date column 
no_reps_contigs$date <- no_reps_contigs$collection_date 

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
no_reps_contigs$date <- plyr::mapvalues(as.character(no_reps_contigs$date), from, to)
unique(no_reps_contigs$date)

# select columns to keep
no_reps_contigs <- no_reps_contigs %>%
  dplyr::select(cohort,pig,bin,contigName,date,value)

######################

# 2nd de-replication: now there should be again replicates because we joined some dates. 
# (take mean out of replicates)
no_reps_contigs$sampleID <- paste0(no_reps_contigs$date,"_",
                                   no_reps_contigs$cohort,"_",
                                   no_reps_contigs$pig,"_",
                                   no_reps_contigs$contigName,
                                   no_reps_contigs$bin)
NROW(no_reps_contigs) == length(unique(no_reps_contigs$sampleID))
NROW(no_reps_contigs) 
length(unique(no_reps_contigs$sampleID)) 

# here above we can see the amount of replicates (same pig, same date)
no_reps_contigs2 <- stats::aggregate(no_reps_contigs$value,
                                     by=list(cohort=no_reps_contigs$cohort,
                                             pig=no_reps_contigs$pig,
                                             bin=no_reps_contigs$bin,
                                             contigName=no_reps_contigs$contigName,
                                             date=no_reps_contigs$date),
                                     data=no_reps_contigs,FUN=mean)

#replace column name x with "value"
names(no_reps_contigs2)[names(no_reps_contigs2) == 'x'] <- 'value'

# sanity check: 
# number of unique pigs
length(unique(paste0(no_reps_contigs2$pig)))

# number of bins in total across all pigs
length(unique(paste0(no_reps_contigs2$pig,no_reps_contigs2$bin)))

# number of bins in total across all pigs
length(unique(paste0(no_reps_contigs2$pig,no_reps_contigs2$contigName)))

# # remove piglets with SD (swine dysentery): 
# no_reps3 <- no_reps2 %>%
#   filter(!pig=="29665"&
#            !pig=="29865"&
#            !pig=="29702")

# Write out to play with it (this is everything)
fwrite(x = no_reps_contigs2,
       file = file.path(out.dir,"no_reps_contigs.csv"))


