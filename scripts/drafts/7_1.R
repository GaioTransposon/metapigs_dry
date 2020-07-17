# 7_1.R 
# takes all dereplicated weighted bins (6.R output)       #
# merges with taxa info (kraken and gtdb)                 #
# merges with dRep clustering data                        #
# saves                                                   #
# for each classification of bins (gtdb, kraken or dRep), #
# it makes wide and aggregates                            #


library(dplyr)
library(data.table)
no_reps_all <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/no_reps_all.csv", 
                                                      na.strings=c("","NA"),
                                                      check.names = FALSE,
                                                      header = TRUE)
NROW(no_reps_all)

######################################################################################################

# load taxonomy from 200gb gtdb_kraken
gtdb <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/all_concatenated_220_essential.kraken", 
                 na.strings=c("","NA"),
                 check.names = FALSE,
                 header = FALSE,
                 sep = "\t")
NROW(gtdb)
head(gtdb)
View(gtdb)

# separate V2 column into pig and bin
library(splitstackshape)
gtdb2 <- cSplit(gtdb, "V2", "_")
gtdb2$V1 <- NULL
gtdb2
length(unique(gtdb2$V2_1)) # yes, 172 pig IDs, which is 126 piglets + 42 mothers + 3 pos ctrl + 1 neg ctrl
length(unique(gtdb2$V3)) # mapping all bins against the gtdb krakeb db (200GB) created 1582 unique taxa 

######################################################################################################

# load taxonomy from 8gb kraken
kraken <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/all_concatenated_8_essential.kraken", 
                 na.strings=c("","NA"),
                 check.names = FALSE,
                 header = FALSE,
                 sep = "\t")
NROW(kraken)
head(kraken)
View(kraken)

# separate V2 column into pig and bin
library(splitstackshape)
kraken2 <- cSplit(kraken, "V2", "_")
kraken2$V1 <- NULL
kraken2
length(unique(kraken2$V2_1)) # yes, 172 pig IDs, which is 126 piglets + 42 mothers + 3 pos ctrl + 1 neg ctrl
length(unique(kraken2$V3)) # mapping all bins against the gtdb krakeb db (200GB) created 1582 unique taxa 

#remove depth cols (we only need one) and "C" cols from dfs
gtdb$V1 <- NULL
kraken$V1 <- NULL
gtdb$V4 <- NULL
kraken$V4 <- NULL
gtdb
kraken

# joins the two taxa datasets
merged_taxa <- merge(gtdb, kraken, by = "V2")
# split V2 column to pig and bins
library(splitstackshape)
merged2 <- cSplit(merged_taxa, "V2", "_")
merged2
#bring pig and bin cols to the front
merged3 <- merged2 %>% select(V2_2, everything())
merged_taxa <- merged3 %>% select(V2_1, everything())
merged_taxa
#rename cols to meaningful names
library(data.table)
setnames(merged_taxa, old = c('V2_1','V2_2', 'V3.x', 'V3.y'), new = c('pig','bin', 'gtdb_taxa', 'kraken_taxa'))
head(merged_taxa)
NROW(merged_taxa)
View(merged_taxa)

# save merged taxa
fwrite(x = merged_taxa, file = "~/Desktop/bins_clustering_parsing_dataframes/merged_taxa.csv")


######################################################################################################

# necessary a semi_join before making wide
tot_counts_all_taxa <-  merge(no_reps_all, merged_taxa, by = c("pig", "bin"))
NROW(tot_counts_all_taxa)
head(tot_counts_all_taxa)

# checking for duplicate rows 
a <- tot_counts_all_taxa
# check how many uniq pig, bin, cohort, date. 
b <- paste(a$date, a$cohort, a$pig, a$bin)
NROW(b)
length(unique(b))

# save bins with gtdb and kraken taxa
fwrite(x = tot_counts_all_taxa, file = "~/Desktop/bins_clustering_parsing_dataframes/gtdb_kraken_bins.csv")

################################

# keep only kraken taxa
kraken_bins <- tot_counts_all_taxa
kraken_bins$gtdb_taxa <- NULL
NROW(kraken_bins)
head(kraken_bins)

# move taxa to front
kraken_bins <- kraken_bins %>% 
  select(kraken_taxa, everything())

# save bins with gtdb taxa
fwrite(x = kraken_bins, file = "~/Desktop/bins_clustering_parsing_dataframes/kraken_bins.csv")

################################

# keep only gtdb taxa
gtdb_bins <- tot_counts_all_taxa
gtdb_bins$kraken_taxa <- NULL
NROW(gtdb_bins)
head(gtdb_bins)

# move taxa to front
gtdb_bins <- gtdb_bins %>% 
  select(gtdb_taxa, everything())

# save bins with kraken taxa
fwrite(x = gtdb_bins, file = "~/Desktop/bins_clustering_parsing_dataframes/gtdb_bins.csv")

############################################################################################################

# load secondary_clusters info 

Cdb = read.csv(
  file="/Users/12705859/Desktop/Cdb.csv"
)
length(unique(Cdb$genome))
# split genome column into pig id and bin\
library(stringr)
Cdb_split <- data.frame(Cdb, str_split_fixed(Cdb$genome, "_", 2))

# rename new columns
colnames(Cdb_split)[colnames(Cdb_split)=="X1"] <- "pig"
colnames(Cdb_split)[colnames(Cdb_split)=="X2"] <- "bin"

#subset df: pig, bin, secondary_cluster
Cdb_essential_cols <- Cdb_split[c("pig", "bin", "secondary_cluster")]
NROW(Cdb_essential_cols)

#change to character pig and bin items in dataframes
Cdb_essential_cols <- Cdb_essential_cols %>%  
  mutate(pig = as.character(pig))
Cdb_essential_cols <- Cdb_essential_cols %>%  
  mutate(bin = as.character(bin))

# merge with Cdb3
tot_counts_all_dRep <- dplyr::left_join(no_reps_all, Cdb_essential_cols, by=c("pig","bin"))

#move secondary_cluster column to first position of dataframe
tot_counts_all_dRep <- tot_counts_all_dRep %>% 
  select(secondary_cluster, everything())
NROW(tot_counts_all_dRep)

View(tot_counts_all_dRep)

# exclude NAs (bins that are not assigned a cluster by dREp)
tot_counts_all_dRep <- na.omit(tot_counts_all_dRep)
NROW(tot_counts_all_dRep)

# save 
fwrite(x = tot_counts_all_dRep, file = "~/Desktop/bins_clustering_parsing_dataframes/dRep_bins.csv")


##############################

# Your aggregation depends on the taxa assignment method (dRep, gtdb or kraken). 
# In this case we are aggregating based on gtdb_taxa assignment

# make widest, one row per gtdb_taxa + bin
library(tidyr)
widest <- gtdb_bins %>% 
  pivot_wider(names_from = c("date", "cohort", "pig"), values_from = "value")
NROW(widest)

bbb <- widest
colsToAggregate <- colnames(bbb[,3:ncol(bbb)])
aggregateBy <- c("gtdb_taxa")
dummyaggfun <- function(v, na.rm = TRUE) {
  c(mean = mean(v, na.rm = TRUE))
}
ccc <- aggregate(bbb[colsToAggregate], by = bbb[aggregateBy], FUN = dummyaggfun)

#replace missing values with zeros (NA values are automatically assigned when long -> wide)
ccc[is.na(ccc)] <- 0
NROW(ccc) # counts of 1578 uniquely gtdb-identified taxa 
NCOL(ccc) # 801 samples

colSums(ccc[,-1])
rowSums(ccc[,-1])

# save 
fwrite(x = ccc, file = "~/Desktop/bins_clustering_parsing_dataframes/gtdb_agg_bins.csv")

##############################

# make widest, one row per kraken_taxa + bin
library(tidyr)
widest <- kraken_bins %>% 
  pivot_wider(names_from = c("date", "cohort", "pig"), values_from = "value")
NROW(widest)

bbb <- widest
colsToAggregate <- colnames(bbb[,3:ncol(bbb)])
aggregateBy <- c("kraken_taxa")
dummyaggfun <- function(v, na.rm = TRUE) {
  c(mean = mean(v, na.rm = TRUE))
}
ccc <- aggregate(bbb[colsToAggregate], by = bbb[aggregateBy], FUN = dummyaggfun)

#replace missing values with zeros (NA values are automatically assigned when long -> wide)
ccc[is.na(ccc)] <- 0
NROW(ccc) # counts of 1578 uniquely gtdb-identified taxa 
NCOL(ccc) # 801 samples

colSums(ccc[,-1])
rowSums(ccc[,-1])

# save 
fwrite(x = ccc, file = "~/Desktop/bins_clustering_parsing_dataframes/kraken_agg_bins.csv")

##############################

# make widest, one row per kraken_taxa + bin
library(tidyr)
widest <- tot_counts_all_dRep %>% 
  pivot_wider(names_from = c("date", "cohort", "pig"), values_from = "value")
NROW(widest)

bbb <- widest
colsToAggregate <- colnames(bbb[,3:ncol(bbb)])
aggregateBy <- c("secondary_cluster")
dummyaggfun <- function(v, na.rm = TRUE) {
  c(mean = mean(v, na.rm = TRUE))
}
ccc <- aggregate(bbb[colsToAggregate], by = bbb[aggregateBy], FUN = dummyaggfun)
View(ccc)
#replace missing values with zeros (NA values are automatically assigned when long -> wide)
ccc[is.na(ccc)] <- 0
NROW(ccc) # counts of 1578 uniquely gtdb-identified taxa 
NCOL(ccc) # 801 samples

colSums(ccc[,-1])
rowSums(ccc[,-1])

# save 
fwrite(x = ccc, file = "~/Desktop/bins_clustering_parsing_dataframes/dRep_agg_bins.csv")

##############################



dRep <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/dRep_bins.csv")
gtdb <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/gtdb_bins.csv")
kraken <- read.csv("~/Desktop/bins_clustering_parsing_dataframes/kraken_bins.csv")
NROW(dRep)
NROW(gtdb)
NROW(kraken)
merge1 <- merge(kraken, gtdb, by = c("pig","bin"))
View(merge1)
merge2 <- merge(merge1, dRep, by = c("pig","bin"))
colnames(merge2)

all <- merge2 %>%
  select(kraken_taxa, gtdb_taxa, secondary_cluster, pig, bin, date, cohort)

# save 
fwrite(x = all, file = "~/Desktop/bins_clustering_parsing_dataframes/bins_gtdb_dRep_kraken.csv")




######
# checking correspondance of ids (e.g. what gtdb taxa is "x" sec cluster)
NROW(all)
length(unique(all$secondary_cluster))
length(unique(all$kraken_taxa))
length(unique(all$gtdb_taxa))
selection <- all %>% filter(
  secondary_cluster == "909_6",
  date == "t2"
)
View(selection)
######
