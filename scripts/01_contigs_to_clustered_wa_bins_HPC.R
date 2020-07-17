#!/usr/bin/Rscript
# runs from the HPC 
# old name: 1_5_HPC.R:
# collects contigs from .fa and merges with depth.txt files
# language: R 
# version this script was developed in: "R version 3.6.0 (2019-04-26)"
# platform this script was developed in: "x86_64-apple-darwin15.6.0"

#This script requires the following packages:
install.packages("base", repos = "http://cran.us.r-project.org")
install.packages("data.table", repos = "http://cran.us.r-project.org", dependencies = TRUE)
install.packages("dplyr", repos = "http://cran.us.r-project.org")
install.packages("seqinr", repos = "http://cran.us.r-project.org")
install.packages("stringr", repos = "http://cran.us.r-project.org")
install.packages("utils", repos = "http://cran.us.r-project.org")

#upload all libraries
library(base)
library(data.table)
library(dplyr)
library(seqinr)
library(stringr)
library(utils)

# clustering input file
Cdb = read.csv(
  file="/shared/homes/12705859/Cdb.csv"
)

# input dir
pig.id.basedir = "/shared/homes/12705859/out_new"

# gather all contig names (>k141") from fasta files and store them in a file "bins_to_contigs.csv"
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  
  df = data.frame(
    pig = character(),
    bin = character(),
    contig = character()
  )
  
  fasta.files = list.files(pig.id.dir, pattern=".fa")
  for (fasta.file in fasta.files){
    ver <- read.fasta(
      file.path(pig.id.dir, fasta.file),
      as.string = TRUE,
      strip.desc = TRUE
    )
    annotation <- getAnnot(ver)
    #unlist to create a vector from a list (goes faster)
    annotation.vector = unlist(annotation)
    
    df = rbind(
      df,
      data.frame(
        pig = rep(pig.id, length(annotation.vector)),
        bin = rep(fasta.file, length(annotation.vector)),
        contig = annotation.vector
      )
    )
  }
  
  fwrite(
    x = df,
    file = file.path(pig.id.dir, "bins_to_contigs.csv"),
    row.names=FALSE
  )
}

#this loop creates a file inside each sample directory from depth.txt files, stripping off 
#the sample IDs from all headers, leaving only sampling dates and extensions
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  
  depth.files = list.files(pig.id.dir, pattern=".txt")
  for (depth.file in depth.files) {
    nub <- read.table(
      file.path(pig.id.dir, depth.file), 
      header=TRUE
    )
    colnames(nub)
    #function to strip the sample ID from the headers
    colClean <- function(nub){ colnames(nub) <- gsub("^[0-9a-zA-z]+_[0-9a-zA-z]+.", " ", colnames(nub)); nub }
    new <- colClean(nub)
  }
  fwrite(
    x = new,
    file = file.path(
      pig.id.dir, 
      "depth_clean.csv"
    ),
    row.names=FALSE
  )
}

#merging depth_clean.csv with bins_to_contigs.csv
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  
  infile1 <- paste(pig.id.dir, "depth_clean.csv", sep="/")
  infile2 <- paste(pig.id.dir, "bins_to_contigs.csv", sep="/")
  
  a <- read.table(
    infile1, 
    header=TRUE, 
    sep=",", 
    row.names=NULL
  )
  b <- read.table(
    infile2, 
    header=TRUE, 
    sep=",", 
    row.names=NULL
  )
  
  a_b <- merge(
    a, 
    b, 
    by.x="contigName", 
    by.y="contig"
  )
  
  fwrite(
    x = a_b,
    file = file.path(pig.id.dir, "depth_with_bins.csv"),
    row.names=FALSE
  )
}

#calculating weighted averages for bins and saving as wa_bins.csv
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  
  old <- paste(pig.id.dir, "depth_with_bins.csv", sep="/")
  
  df <- read.table(
    file = old, 
    header = TRUE, 
    sep = ",", 
    row.names = NULL
  )
  
  df_new <- df %>% 
    select(-ends_with(".var"))
  dt <- as.data.table(df_new)
  colsTOweight <- dt %>% 
    select(contains('.bam'))
  
  # create bin length from contig lengths
  bins_lengths <- setNames(aggregate(dt$contigLen, by=list(bin=dt$bin), FUN=sum), c("bin","binLen"))
  
  dt_new <- dt[,lapply(.SD,weighted.mean,w=contigLen),by=list(bin, pig), .SDcols = colnames(colsTOweight)]
  
  # merge the current dataframe (with weighted bins) with dataframe containing bins lengths
  dt_new2 <- merge(bins_lengths,dt_new, by="bin")
  
  fwrite(
    x = dt_new2,
    file = file.path(pig.id.dir, "wa_bins.csv"),
    row.names=FALSE
  )
}

# split genome column into pig id and bin
Cdb2 <- data.frame(Cdb, str_split_fixed(Cdb$genome, "_", 2))

# rename new columns
colnames(Cdb2)[colnames(Cdb2)=="X1"] <- "pig"
colnames(Cdb2)[colnames(Cdb2)=="X2"] <- "bin"

#subset df: pig, bin, secondary_cluster
Cdb3 <- Cdb2[c("pig", "bin", "secondary_cluster")]

#change to character pig and bin items in dataframes
Cdb4 <- Cdb3 %>%  
  mutate(pig = as.character(pig))
Cdb5 <- Cdb4 %>%  
  mutate(bin = as.character(bin))

# merging bins clustering info to depths, creating a file in each
# pig dir containing the weighted depths per bin and its clustering info
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  
  # read in wa_bins.csv
  wa_bins <- paste(pig.id.dir, "wa_bins.csv", sep="/")
  wa_bins_df <- read.table(
    file = wa_bins, 
    header = TRUE, 
    sep = ",", 
    row.names = NULL
  )
  wa_bins_df2 <- wa_bins_df %>%  
    mutate(pig = as.character(pig))
  wa_bins_df3 <- wa_bins_df2 %>%  
    mutate(bin = as.character(bin))
  
  # merge with Cdb3
  clu_wa_bins <- dplyr::left_join(wa_bins_df3, Cdb5, by=c("pig","bin"))
  
  #move secondary_cluster column to first position of dataframe
  clustered_wa_bins <- clu_wa_bins %>% 
    select(secondary_cluster, everything())
  
  fwrite(
    x = clustered_wa_bins,
    file = file.path(pig.id.dir, "clustered_wa_bins.csv"),
    row.names=FALSE
  )
}
