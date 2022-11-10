
# WHAT THIS SCRIPT DOES: 
# This script takes the contig depths 
# from each subject, and summarises them into bin depths. 

# STEPS:
## STEP 1. contigs depths to bins depths 
## STEP 2. Reads in the clustering info file (output of dRep) & 
# merges clustering info to bin depths. 


library(readr)
library(reshape)
library(tidyr)
library(splitstackshape)
library(readxl)
library(stringr)
library(dplyr)
library(data.table)


# input dir
pig.id.basedir = "/shared/homes/152324/out_new" # make sure you have permission to write in the "subjectDIR"/metabat folder
source_dir = "/shared/homes/152324/metapigs_function/source_data/" # should contain: Cdb.csv

#test here: 
#pig.id.basedir = "/Users/dgaio/cloudstor/Gaio/out_new_test" 
#source_dir = "/Users/dgaio/cloudstor/Gaio/github/metapigs_function/source_data/" # should contain: Cdb.csv




########################################################

# Read in clustering info:

# clustering input file
Cdb = read.csv(paste0(source_dir,"Cdb.csv"))

# split genome column into pig id and bin
Cdb2 <- data.frame(Cdb, str_split_fixed(Cdb$genome, "_", 2))

# rename new columns
colnames(Cdb2)[colnames(Cdb2)=="X1"] <- "pig"
colnames(Cdb2)[colnames(Cdb2)=="X2"] <- "bin"

#subset df: pig, bin, secondary_cluster
Cdb3 <- Cdb2[c("pig", "bin", "secondary_cluster")]


########################################################


# contigs depths to bins depths 
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  
  f <- paste(pig.id.dir, "contig_depths.csv", sep="/")

  # if the file exists, open it
  if (file_test("-f", f)) {
    print(paste0(f,"   file exists"))
    
    df <- read.table(
      file = f, 
      header = TRUE, 
      sep = ",", 
      row.names = NULL
    )
    
    
    ######### STEP 1 #######################################
    
    
    
    # take the median (more statistically robust, not affected by outliers)
    bins_depths <- df %>%
      group_by(pig,bin,date) %>% 
      summarise(value=median(value),
                binLen=sum(contigLen),
                binLen_sd=sd(contigLen))
    
    # # little sanity check (this is tested on a sow subject - because only one time point)
    # check <- df %>%
    #   dplyr::filter(bin=="bins.1.fa") %>%
    #   dplyr::filter(date=="tM")
    # out <- bins_depths %>%
    #   dplyr::filter(bin=="bins.1.fa") %>%
    #   dplyr::filter(date=="tM")
    # out$value==median(check$value)
    # out$binLen==sum(check$contigLen)
    
    
    ######### STEP 2 #######################################
    
    
    
    bins_depths$pig <- as.character(bins_depths$pig)
    
    # merge
    clu_bins_depths <- dplyr::left_join(bins_depths, Cdb3, by=c("pig","bin"))
    
    #move secondary_cluster column to first position of dataframe
    clu_bins_depths <- clu_bins_depths %>% 
      dplyr::select(secondary_cluster, everything())
    

    fwrite(
      x = clu_bins_depths,
      file = file.path(pig.id.dir, "bin_depths.csv"),
      row.names=FALSE)
    
  }
  
  else { print(paste0(f,"   file does not exist"))}
  
}

View(clu_bins_depths)


# # Metabat's peeps wrote this script to aggregate contig depths into bins depths. 
# # The script below produces the same bins depths produced above. 
# library(foreach)
# library(matrixStats)
# 
# depth <- read.table('/Users/dgaio/cloudstor/Gaio/out_new_test/Y09733/metabat/depth.txt', as.is=T, header=T)
# depth <- depth[,-which(grepl('var$',colnames(depth)))][,-c(2,3)]
# bins <- list.files('/Users/dgaio/cloudstor/Gaio/out_new_test/Y09733/metabat/', pattern='.fa', full.name=T)
# mems <- foreach(f = bins) %do% {
#   sub('>','',system(sprintf("grep '^>' %s", f), intern=T))
# }
# names(mems) <- sub('.fa','',basename(bins)) 
# 
# bins.cv <- do.call(rbind, foreach(m=mems) %do% {
#   colMedians(as.matrix(depth[match(m, depth$contigName),-1]))
# })
# dimnames(bins.cv) <- list(names(mems), colnames(depth)[-1])
# View(bins.cv)
# # compare bins.cv with '/Users/dgaio/cloudstor/Gaio/out_new_test/Y09733/metabat/bin_depths.csv'






    
    






