

library(readr)
library(dplyr)
library(tidyr)
library(data.table)
library(splitstackshape)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(treemapify)
library(ggrepel)
library(factoextra)
library(stringr)
library(zip)


source_dir = "/Users/12705859/metapigs_dry/source_data/" # git 
middle_dir = "/Users/12705859/metapigs_dry/middle_dir/" # git 
out_dir_git = "/Users/12705859/metapigs_dry/out/" # git 
out_dir = "/Users/12705859/Desktop/metapigs_dry/dbcan/"  # local

######################################################################

# load dbcan_hmmer data (output of mapping our-bins-AAs against CAZy db)

## First extract
tmp2 <- tempfile()
unzip(zipfile=paste0(middle_dir,"hmmer.out.zip"), exdir = tmp2)

hmmer <- read_table2(paste0(tmp2,"/hmmer.out"), col_names = FALSE)


##########################################################

# HMMER data cleaning 

# remove .hmm from enzyme ID
hmmer <- cSplit(hmmer,"X1",".")

# separate by _
hmmer <- cSplit(hmmer,"X1_1","_")

# copy enzyme ID into new column 
hmmer$enzymeID <- hmmer$X1_1_1

# extra enzyme ID - suffixes together 
hmmer <- hmmer  %>%
  dplyr::mutate_at(c('X1_1_2', 'X1_1_3','X1_1_4','X1_1_5'), ~str_replace_na(., "")) %>%
  dplyr::mutate(combo_var = paste0(X1_1_2,".", X1_1_3,".",X1_1_4,".",X1_1_5))

# subject/bin/contig
hmmer <- cSplit(hmmer,"X3","_")

# columns selection
hmmer <- hmmer %>%
  dplyr::select(enzymeID,combo_var,
                X3_1,X3_2, X3_3, # pig, bin, predicted_protein_ID
                X2,X4, # length enzyme, length subject protein
                X5,X6,X7,X8,X9,X10) # evalue, start query, end query, start sample, end sample, coverage  

# get enzyme class
hmmer$getclass <- hmmer$enzymeID
hmmer <- hmmer %>%
  separate(getclass, 
           into = c("text", "num"), 
           sep = "(?<=[A-Za-z])(?=[0-9])"
  )

hmmer$num <- NULL

# cols renaming
colnames(hmmer) <- c("enzymeID","enzymeID_suffixes",
                     "pig","bin","predicted_protein_ID",
                     "HMM_length","query_lenght",
                     "evalue",
                     "start HMM","end HMM",
                     "start sample","end sample",
                     "coverage", "enzymeNAME")

hmmer$bin <- gsub(".fa","", hmmer$bin)

# QUESTIONABLE: 
# when there is two or more hits for a certain predicted protein, 
# where you get two hits matching two enzymeIDs (with different e-values), 
# should you or should you not take only the best evalue? 
# problem is that if the match is just as good (99% coverage identity), 
# but the enzyme within the database is longer, 
# then whenever the choice is between a longer vs a shorter enzyme ID, 
# you would always get a better e-value for the longer enzyme.  
# up to now I opted not to filter. 

####################################################################
################# TEST start
#################
# let's try filtering; but let's first run a test

# reduce dataset for test (columns)
hmmer1 <- hmmer %>%
  dplyr::select(enzymeID,enzymeID_suffixes,pig,bin,predicted_protein_ID,evalue,enzymeNAME,coverage)

# reduce dataset for test (rows)
hmmer2 <- hmmer1[1:50,]

# for group selected, arrange evalue (you will see that for some predicted proteins there are two hits; spot them)
hmmer3 <- hmmer2 %>%
  group_by(pig,bin,predicted_protein_ID,coverage) %>%
  arrange(evalue) 
#View(hmmer3)

# now, for group selected, slice to obtain only the lowest evalue-hit. Compare this df with the one above to see if it works (it works)
hmmer4 <- hmmer3 %>%
  group_by(pig,bin,predicted_protein_ID) %>%
  arrange(evalue) %>%
  slice(1)
#View(hmmer4)
#################
################# TEST end
####################################################################

# Let's now run on whole dataset

NROW(hmmer)
hmmerClean <- hmmer %>%
  dplyr::select(enzymeID,enzymeID_suffixes,pig,bin,predicted_protein_ID,evalue,enzymeNAME,coverage) %>%
  group_by(pig,bin,predicted_protein_ID) %>%
  arrange(evalue) %>%
  slice(1)
NROW(hmmerClean)
head(hmmerClean)

#################

# number of lines is slightly reduced from original as we only took the best hit for each predicted protein 
# we also now have less columns 

# save hmmer CLEAN output
write.csv(hmmerClean, gzfile(paste0(middle_dir, "hmmer.out_Clean.csv.gz")), row.names = FALSE)

##########################################################
##########################################################

# get enzyme count: when the same enzyme is present twice within the same bin (and host), count should be 2. 
hmmerCleanCounts <- hmmerClean %>%
  dplyr::ungroup() %>%
  dplyr::select(enzymeID,pig,bin,enzymeNAME) %>% 
  group_by(enzymeNAME,enzymeID,pig,bin) %>%
  dplyr::summarise(enz_count=n()) %>%
  dplyr::select(enzymeID,pig,bin,enz_count,enzymeNAME) 

# save hmmer CleanCounts output
write.csv(hmmerCleanCounts, gzfile(paste0(middle_dir, "hmmer.out_CleanCounts.csv.gz")), row.names = FALSE)

##########################################################
##########################################################
