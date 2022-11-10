# Install packages:
# install.packages("base", repos = "http://cran.us.r-project.org")
# install.packages("data.table", repos = "http://cran.us.r-project.org", dependencies = TRUE)
# install.packages("dplyr", repos = "http://cran.us.r-project.org")
# install.packages("seqinr", repos = "http://cran.us.r-project.org")
# install.packages("stringr", repos = "http://cran.us.r-project.org")
# install.packages("utils", repos = "http://cran.us.r-project.org")

#upload all libraries
library(base)
library(data.table)
library(dplyr)
library(seqinr)
library(stringr)
library(utils)

# input dir
pig.id.basedir = "/shared/homes/152324/out_new" #make sure you have permission to write in the subdir /metabat
## test here: 
#pig.id.basedir = "/Users/dgaio/cloudstor/Gaio/out_new_test"


# 1. gather all contig names (>k141") from fasta files and store them in a file "bins_to_contigs.csv"
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  
  
  # if the dir is not empty:
  if (file.exists(pig.id.dir)) {
    
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
    
  } else { print(paste0(pig.id.dir," directory is empty")) }
  
}




# 2. this loop creates a file inside each sample directory from depth.txt files, stripping off 
#the sample IDs from all headers, leaving only sampling dates and extensions
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  
  # if the dir is not empty:
  if (file.exists(pig.id.dir)) {
    
    depth.files = list.files(pig.id.dir, pattern=".txt")
    
    print(paste0(depth.files," exists"))
    
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
    
  } else { print(paste0(pig.id.dir," directory is empty")) }
  
}



# 3. merging depth_clean.csv with bins_to_contigs.csv
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")

  # if the dir is not empty:
  if (file.exists(pig.id.dir)) {
    
    infile1 <- paste(pig.id.dir, "depth_clean.csv", sep="/")
    infile2 <- paste(pig.id.dir, "bins_to_contigs.csv", sep="/")
    
    print(paste0(infile1, "   file exists", "  and  ",infile2, "   file exists"))
    
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
    
    # need to change the column name to match 
    id <- b$pig[1]
    b$pig <- NULL
    colnames(a)[colnames(a)=="contigName"] <- "contig"
    
    # merge bin info to depths  
    a_b <- full_join(a, b)
    a_b$pig <- as.character(id)
    
    # if contig has not been binnen to any bin, explicitly attribute "no_bin" string
    a_b$bin <- a_b$bin %>% replace_na("no_bin")
    
    fwrite(
      x = a_b,
      file = file.path(pig.id.dir, "depth_with_bins.csv"),
      row.names=FALSE
    )
    
  } else { print(paste0(pig.id.dir," directory is empty")) }
  
}

