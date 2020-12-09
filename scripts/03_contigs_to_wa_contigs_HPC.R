# runs from the HPC 
# language: R 

#This script requires the following packages:
install.packages("base", repos = "http://cran.us.r-project.org")
install.packages("data.table", repos = "http://cran.us.r-project.org", dependencies = TRUE)
install.packages("dplyr", repos = "http://cran.us.r-project.org")
install.packages("stringr", repos = "http://cran.us.r-project.org")
install.packages("utils", repos = "http://cran.us.r-project.org")
install.packages("tidyr", repos = "http://cran.us.r-project.org")
install.packages("plyr", repos = "http://cran.us.r-project.org")
install.packages("readxl", repos = "http://cran.us.r-project.org")

#upload all libraries
library(base)
library(data.table)
library(dplyr)
library(stringr)
library(utils)
library(tidyr)
library(plyr)
library(readxl)

# input dir
pig.id.basedir = "/shared/homes/12705859/out_new" # on HPC
out.dir = "/shared/homes/12705859/contig_abundances" # on HPC
#pig.id.basedir = "/Users/12705859/Desktop/bins_clustering_parsing_DFs/out_new_test" # on local
#out.dir = "/Users/12705859/Desktop/bins_clustering_parsing_DFs/contig_abundances" # on local

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

# create empty df with all possible headers (dates) and interate through each sample's df
# containing depths from weighted contigs 
# concatenate all to one df > output as: new_headers_wa_contigs.csv
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  
  df = data.frame(
    pig = character(),
    contigName = character(),
    contigLen = character(),
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
    old <- paste(pig.id.dir, "wa_contigs.csv", sep="/")
    df2 <- read.table(file = old, header = TRUE, sep = ",", row.names = NULL)
    
    new_df <- rbind.fill(df,df2)
    
    write.csv(new_df, row.names=FALSE, quote=FALSE, file = file.path(pig.id.dir, "new_headers_wa_contigs.csv"))
    
    # fwrite(
    #   x = new_df,
    #   file = file.path(pig.id.dir, "new_headers_wa_contigs.csv"),
    #   row.names=FALSE
    # )
  }
  else {
    print("not a pig subject")
  }
}



df_total = data.frame()
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  #if statement to exclude NegativeControl, MockCommunity, Protexin and Coliguard
  if (grepl("[[:digit:]]", pig.id) == TRUE) {
    old <- paste(pig.id.dir, "new_headers_wa_contigs.csv", sep="/")
    df2 <- read.table(file = old, header = TRUE, sep = ",", row.names = NULL)
    
    df3 <- data.frame(df2)
    df_total <- rbind(df_total,df3)
    
    fwrite(
      x = df_total,
      file = file.path(out.dir,"merged_all_wa_contigs.csv"),
      row.names=FALSE
    )
  }
}


# merge merged_all_clustered_wa_bins.csv with cohorts.xlsx (from metapigs/source_data)
#input files
cohorts <- read_excel("/shared/homes/12705859/cohorts.xlsx") # on HPC
#cohorts <- read_excel("/Users/12705859/metapigs_dry/source_data/cohorts.xlsx") # on local
df_total <- read.csv(file.path(out.dir,"merged_all_wa_contigs.csv")) 

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
  dplyr::select("contigLen", everything())

#move Cohort Name column to first position of dataframe
df_total_4 <- df_total_3 %>% 
  dplyr::select("Cohort Name", everything())

#rename to eliminate space
colnames(df_total_4)[colnames(df_total_4)=="Cohort Name"] <- "cohort"

fwrite(
  x = df_total_4,
  file = file.path(out.dir,"merged_all_wa_contigs_with_cohorts.csv"), 
  row.names=FALSE
)


