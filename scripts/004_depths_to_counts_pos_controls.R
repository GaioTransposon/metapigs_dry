
# no_reps_pos_controls <- counts_long.R script                             #
# bins without a cluster assigned will be "no_cluster"                     #
# depth to counts                                                          #
# removes samples replicates  (none in this case, for pos controls)        #
# saves no replicated bins                                                 #

library(data.table)
library(tidyverse)
library(magrittr)
library(reshape)
library(splitstackshape)
library(dplyr)


middle_dir = "/Users/12705859/Desktop/metapigs_dry/middle_dir/"


# upload input file
MockCommunity <- read.csv(paste0(middle_dir,"MockCommunity_clustered_wa_bins.csv"),
               na.strings=c("","NA"),
               check.names = FALSE,
               header = TRUE)
colnames(MockCommunity) = c("secondary_cluster","bin","binLen","pig","R1","R2","R3","R4","R5","R6","R7","R8","R9")
MockCommunity$isolation_source = "MockCommunity"
MockCommunity <- MockCommunity %>%
  dplyr::select(binLen,isolation_source,bin,secondary_cluster,R1,R2,R3,R4,R5,R6,R7,R8,R9)

# upload input file
ColiGuard <- read.csv(paste0(middle_dir,"ColiGuard_clustered_wa_bins.csv"),
                          na.strings=c("","NA"),
                          check.names = FALSE,
                          header = TRUE)
colnames(ColiGuard) = c("secondary_cluster","bin","binLen","pig","R1","R2","R3","R4","R5","R6","R7","R8")
ColiGuard$isolation_source = "ColiGuard"
ColiGuard$R9 = 0
ColiGuard <- ColiGuard %>%
  dplyr::select(binLen,isolation_source,bin,secondary_cluster,R1,R2,R3,R4,R5,R6,R7,R8,R9)

# upload input file
Protexin <- read.csv(paste0(middle_dir,"Protexin_clustered_wa_bins.csv"),
                          na.strings=c("","NA"),
                          check.names = FALSE,
                          header = TRUE)
colnames(Protexin) = c("secondary_cluster","bin","binLen","pig","R1","R2","R3","R4","R5","R6","R7","R8")
Protexin$isolation_source = "Protexin"
Protexin$R9 = 0
Protexin <- Protexin %>%
  dplyr::select(binLen,isolation_source,bin,secondary_cluster,R1,R2,R3,R4,R5,R6,R7,R8,R9)


# ready all to be concatenated ! 

df <- rbind(MockCommunity, Protexin, ColiGuard)
df$cohort = "PosCtrl"
df <- df %>%
  dplyr::select(cohort,binLen,isolation_source,bin,secondary_cluster,R1,R2,R3,R4,R5,R6,R7,R8,R9)

# bins that don't have secondary cluster assigned are filled with "no_cluster"
df <- df %<>% dplyr::mutate(secondary_cluster = fct_explicit_na(secondary_cluster, na_level = "no_cluster"))


################


# get the column names before entering function
original_colnames <- colnames(df)

# transform depth values into counts
# where {counts = depth*binLen/300} (300 is the read pair length) #2 is binLen 1:5 is the non-depths data
A <- function(x) x * df[,2] / 300
counts <- data.frame(df[1:5], apply(df[6:ncol(df)],2, A) )
counts
#return the original colnames to new dataframe
colnames(counts) <- original_colnames
#now we don't need binLen anymore we can remove this column
counts <- counts[,-2]

NROW(counts)


#make long (one value per row)

counts_long <- reshape::melt(counts, id=c("cohort", "isolation_source", "bin", "secondary_cluster"))
NROW(counts_long)

# remove NAs (NAs were automatically created when pivoting large)
counts_long <- na.omit(counts_long)
NROW(counts_long)
head(counts_long)



# no de-replication needed for the positive controls (differently from 7.R)
# "no_reps" as a prefix just for consistency

no_reps_pos_controls <- counts_long


# ready! 



# Write out to play with it (this is everything)
fwrite(x = no_reps_pos_controls, file = paste0(middle_dir,"no_reps_pos_controls.csv"))

# ########################################################
