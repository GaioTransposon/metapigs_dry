#!/usr/bin/Rscript
# runs from the HPC 
# old name: 6_HPC.R: 
# merges depths from pigs together (incl. cohorts info)
# language: R 
# version this script was developed in: "R version 3.6.0 (2019-04-26)"
# platform this script was developed in: "x86_64-apple-darwin15.6.0"

# this script : 
# 1. merges all the weighted and clustered bins depths from all pig   #
#    directories (excl. positive controls) into one dataframe.        #            
# 2. merges this dataframe with cohorts.xlsx                          #
# 3. rearranges column positions, reformats, NA-fills missing values  #


# required packages:
install.packages("data.table", repos = "http://cran.us.r-project.org", dependencies = TRUE)
install.packages("plyr", repos = "http://cran.us.r-project.org")
install.packages("seqinr", repos = "http://cran.us.r-project.org")
install.packages("readxl", repos = "http://cran.us.r-project.org")

# upload all libraries
library(seqinr)
library(data.table)
library(plyr)
library(dplyr)
library(readxl)

# create empty df with all possible headers (dates) and interate through each sample's df
# containing depths from clustered weighted bins 
# concatenate all to one df > output as: new_headers_clustered_wa_bins.csv
pig.id.basedir = "/shared/homes/12705859/out_new"
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
    
  df = data.frame(
    pig = character(),
    bin = character(),
    secondary_cluster = character(),
    binLen = character(),
    X170130.01.bam = character(),
    X170131.01.bam = character(),
    X170131.02.bam = character(),
    X170131.03.bam = character(),
    X170201.01.bam = character(),
    X170201.02.bam = character(),
    X170203.01.bam = character(),
    X170206.01.bam = character(),
    X170206.02.bam = character(),
    X170207.01.bam = character(),
    X170207.02.bam = character(),
    X170208.01.bam = character(),
    X170210.01.bam = character(),
    X170214.01.bam = character(),
    X170214.02.bam = character(),
    X170214.03.bam = character(),
    X170216.01.bam = character(),
    X170216.02.bam = character(),
    X170217.01.bam = character(),
    X170221.01.bam = character(),
    X170224.01.bam = character(),
    X170228.01.bam = character(),
    X170303.01.bam = character(),
    X170306.01.bam = character(),
    X170307.01.bam = character(),
    X170308.01.bam = character(),
    X170309.01.bam = character(),
    X170310.01.bam = character()
  )
  #if statement to exclude NegativeControl, MockCommunity, Protexin and Coliguard
  if (grepl("[[:digit:]]", pig.id) == TRUE) {
    old <- paste(pig.id.dir, "clustered_wa_bins.csv", sep="/")
    df2 <- read.table(file = old, header = TRUE, sep = ",", row.names = NULL)
    
    new_df <- rbind.fill(df,df2)
    
    fwrite(
      x = new_df,
      file = file.path(pig.id.dir, "new_headers_clustered_wa_bins.csv"),
      row.names=FALSE
    )
  }
}

pig.id.basedir = "/shared/homes/12705859/out_new"
df_total = data.frame()
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  #if statement to exclude NegativeControl, MockCommunity, Protexin and Coliguard
  if (grepl("[[:digit:]]", pig.id) == TRUE) {
    old <- paste(pig.id.dir, "new_headers_clustered_wa_bins.csv", sep="/")
    df2 <- read.table(file = old, header = TRUE, sep = ",", row.names = NULL)
    
    df3 <- data.frame(df2)
    df_total <- rbind(df_total,df3)
    
    fwrite(
      x = df_total,
      file = file.path("/shared/homes/12705859/merged_all_clustered_wa_bins.csv"),
      row.names=FALSE
      )
  }
}


# merge merged_all_clustered_wa_bins.csv with cohorts.xlsx (from metapigs/source_data)

#input files
cohorts <- read_excel("/shared/homes/12705859/cohorts.xlsx")
df_total <- read.csv("/shared/homes/12705859/merged_all_clustered_wa_bins.csv")

# reformat dates
colnames(df_total)[colnames(df_total)=="X170130.01.bam"] <- "17-01-30.1"
colnames(df_total)[colnames(df_total)=="X170131.01.bam"] <- "17-01-31.1"
colnames(df_total)[colnames(df_total)=="X170131.02.bam"] <- "17-01-31.2"
colnames(df_total)[colnames(df_total)=="X170131.03.bam"] <- "17-01-31.3"
colnames(df_total)[colnames(df_total)=="X170201.01.bam"] <- "17-02-01.1"
colnames(df_total)[colnames(df_total)=="X170201.02.bam"] <- "17-02-01.2"
colnames(df_total)[colnames(df_total)=="X170203.01.bam"] <- "17-02-03.1"
colnames(df_total)[colnames(df_total)=="X170206.01.bam"] <- "17-02-06.1"
colnames(df_total)[colnames(df_total)=="X170206.02.bam"] <- "17-02-06.2"
colnames(df_total)[colnames(df_total)=="X170207.01.bam"] <- "17-02-07.1"
colnames(df_total)[colnames(df_total)=="X170207.02.bam"] <- "17-02-07.2"
colnames(df_total)[colnames(df_total)=="X170208.01.bam"] <- "17-02-08.1"
colnames(df_total)[colnames(df_total)=="X170210.01.bam"] <- "17-02-10.1"
colnames(df_total)[colnames(df_total)=="X170214.01.bam"] <- "17-02-14.1"
colnames(df_total)[colnames(df_total)=="X170214.02.bam"] <- "17-02-14.2"
colnames(df_total)[colnames(df_total)=="X170214.03.bam"] <- "17-02-14.3"
colnames(df_total)[colnames(df_total)=="X170216.01.bam"] <- "17-02-16.1"
colnames(df_total)[colnames(df_total)=="X170216.02.bam"] <- "17-02-16.2"
colnames(df_total)[colnames(df_total)=="X170217.01.bam"] <- "17-02-17.1"
colnames(df_total)[colnames(df_total)=="X170221.01.bam"] <- "17-02-21.1"
colnames(df_total)[colnames(df_total)=="X170207.01.bam"] <- "17-02-7.1"
colnames(df_total)[colnames(df_total)=="X170224.01.bam"] <- "17-02-24.1"
colnames(df_total)[colnames(df_total)=="X170228.01.bam"] <- "17-02-28.1"
colnames(df_total)[colnames(df_total)=="X170303.01.bam"] <- "17-03-3.1"
colnames(df_total)[colnames(df_total)=="X170306.01.bam"] <- "17-03-6.1"
colnames(df_total)[colnames(df_total)=="X170307.01.bam"] <- "17-03-7.1"
colnames(df_total)[colnames(df_total)=="X170308.01.bam"] <- "17-03-8.1"
colnames(df_total)[colnames(df_total)=="X170309.01.bam"] <- "17-03-9.1"
colnames(df_total)[colnames(df_total)=="X170310.01.bam"] <- "17-03-10.1"

#merge cohorts
df_total_2 <- merge.data.frame(df_total, cohorts, by.x="pig", by.y = "Animal ID")

#move binLen column to first position of dataframe
df_total_3 <- df_total_2 %>% 
  select("binLen", everything())

#move Cohort Name column to first position of dataframe
df_total_4 <- df_total_3 %>% 
  select("Cohort Name", everything())

#rename to eliminate space
colnames(df_total_4)[colnames(df_total_4)=="Cohort Name"] <- "cohort"

#fill in column "secondary_cluster" -empty cells with NA
df_total_5 <- mutate_all(df_total_4, funs(na_if(.,"")))

fwrite(
  x = df_total_5,
  file = file.path("/shared/homes/12705859/merged_all_clustered_wa_bins_with_cohorts.csv"),
  row.names=FALSE
)

