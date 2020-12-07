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
install.packages("tidyr", repos = "http://cran.us.r-project.org")

#upload all libraries
library(base)
library(data.table)
library(dplyr)
library(seqinr)
library(stringr)
library(utils)
library(tidyr)

# input dir
pig.id.basedir = "/shared/homes/12705859/out_new"


#calculating weighted averages for contigs and saving as wa_contigs.csv
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  
  old <- paste(pig.id.dir, "depth_with_bins.csv", sep="/")
  
  df <- read.table(
    file = old, 
    header = TRUE, 
    sep = ",", 
    row.names = NULL
  )
  
  dt <- as.data.table(df)
  
  # remove the var cols
  dt_new <- select(dt, -contains(".var"))
  
  # weight the contigs abundances by the contig length
  dt_new2 <- dt_new %>% 
    pivot_longer(cols = contains(".bam")) %>% 
    dplyr::mutate(value=value/contigLen) %>% 
    pivot_wider(names_from = name)
  
  fwrite(
    x = dt_new2,
    file = file.path(pig.id.dir, "wa_contigs.csv"),
    row.names=FALSE
  )
}
