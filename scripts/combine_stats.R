
###########################################################################################

library(readr)
library(openxlsx)


out_dir_git = "/Users/12705859/metapigs_dry/out/"



## 1. create workbook 
wb <- createWorkbook()



## 2. open csv and save as sheets of the workbook
filenames <- list.files(out_dir_git, pattern="*.csv", full.names=TRUE)

for (each_filename in filenames) {
  
  ldf <- lapply(each_filename, read.csv)
  
  clean_name <- gsub(".*out//\\s*|.csv.*", "", as.character(each_filename))
  
  addWorksheet(wb, clean_name)
  writeData(wb, sheet = clean_name, ldf, rowNames = FALSE)
  
}


## 3. save workbook 
saveWorkbook(wb, paste0(out_dir_git,"stats.xlsx"), overwrite=TRUE)


###########################################################################################

