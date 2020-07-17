library(readr)
library(splitstackshape)
library(dplyr)
library(data.table)
library(taxize)
library(tidyr)
library(treemap)
library(ggplot2) 

setwd("~/Desktop/metapigs_dry/gtdbtk")

# load gtdbtk assignments of the bins
all_concatenated_220_essential <- read_delim("~/Desktop/metapigs_dry/gtdbtk/all_concatenated_220_essential.kraken", 
                                             "\t", escape_double = FALSE, col_names = FALSE, 
                                             trim_ws = TRUE)

all_concatenated_220_essential <- cSplit(all_concatenated_220_essential, "X3","(")
all_concatenated_220_essential <- cSplit(all_concatenated_220_essential, "X3_2",")")
all_concatenated_220_essential <- cSplit(all_concatenated_220_essential, "X2","_")
all_concatenated_220_essential$X2_2 <- gsub(".fa","", all_concatenated_220_essential$X2_2)
all_concatenated_220_essential$X3_2_1 <- gsub("taxid ","", all_concatenated_220_essential$X3_2_1)

# separate the sp[digits]
all_concatenated_220_essential <- cSplit(all_concatenated_220_essential, "X3_1"," ")

# all containing numbers 
all_concatenated_220_essential$contain_digits_first <- grepl("[[:digit:]]", all_concatenated_220_essential$X3_1_1)

# separate underscore
all_concatenated_220_essential <- cSplit(all_concatenated_220_essential, "X3_1_1","_")
all_concatenated_220_essential <- cSplit(all_concatenated_220_essential, "X3_1_2","_")

# useless 
all_concatenated_220_essential$X3_1_1_2 <- NULL
all_concatenated_220_essential$X3_1_2_2 <- NULL

# all containing numbers 
all_concatenated_220_essential$contain_digits_second <- grepl("[[:digit:]]", all_concatenated_220_essential$X3_1_2)

# subset : 
df_simple <- subset(all_concatenated_220_essential , contain_digits_first == FALSE)

# subset : 
df_other <- anti_join(all_concatenated_220_essential, df_simple)

# awesome, it sums up: 
NROW(df_simple)+NROW(df_other)

# subset : 
df_simple_sp <- subset(df_simple , contain_digits_second == TRUE)

df_simplest <- anti_join(df_simple, df_simple_sp)

unique(df_simplest$X3_1_1_1)
unique(df_simplest$X3_1_2_1)

unique(df_simple_sp$X3_1_1_1)
unique(df_simple_sp$X3_1_2_1)

unique(df_other$X3_1_1_1)
unique(df_other$X3_1_2_1)

# awesome, it sums up: 
NROW(df_simplest)+NROW(df_simple_sp)+NROW(df_other)

df_simplest$name_type = "simplest"
df_simplest$gtdbtk_name <- paste0(df_simplest$X3_1_1_1," ",df_simplest$X3_1_2_1)
df_simplest$sp_id <- NA

df_simple_sp$name_type = "simple_sp"
df_simple_sp$gtdbtk_name = paste0(df_simple_sp$X3_1_1_1)
df_simple_sp$sp_id <- paste0(df_simple_sp$X3_1_2_1)
  

df_other$name_type = "other"
df_other$gtdbtk_name <- paste0(df_other$X3_1_1_1," ",df_other$X3_1_2_1)
df_other$sp_id <- paste0(df_other$X3_1_2_1)

df <- rbind(df_simplest,df_simple_sp,df_other)

df <- df %>%
  mutate(taxid = X3_2_1,
         pig = X2_1,
         bin = X2_2) %>%
  select(taxid,pig,bin,gtdbtk_name,name_type,sp_id)

View(df)


# add ncbi full names 
# nb: this can only be done for "simplest" and "simple_sp", not for "other" as these 
# extra info can't be retrieved

# what proportion are we talking about? 
df1 <- setDT(df)[, .(Freq = .N), by = .(name_type)]
df1[order(df1$name_type)]

classifiable <- df %>%
  dplyr::filter(name_type=="simplest"|name_type=="simple_sp") %>%
  select(pig,bin,gtdbtk_name)
NROW(classifiable)
head(classifiable)

colnames(classifiable)[colnames(classifiable) == 'gtdbtk_name'] <- 'name'

# remove "NA" string (NA was assigned next to genus when speies name not available, it can be removed)
classifiable$name <- gsub(" NA","", classifiable$name)
NROW(classifiable)


# Taxize will be used to find full taxonomy for each genus or species we have

# get the unique names
taxa_to_search <- unique(classifiable$name)
taxa_to_search <- as.data.frame(taxa_to_search)

# subset into smaller dataframes (570 can be divided into 15 = 38) otherwise we'll get a "Argument list too long" error message
taxa_search_multiple <- split(taxa_to_search, (0:NROW(taxa_to_search) %/% 57))

# remove the last item from the list as it's empty and it would otherwise cause issues in the loop below
taxa_search_multiple[11] <- NULL

mylist <- as.list(as.character(taxa_to_search$taxa_to_search))

# construct an empty dataframe to build on 
fin <- data.frame(
  name = character(),
  rank = character(),
  id = character(),
  index = character(),
  stringsAsFactors = FALSE
)

for(i in mylist){
  
  print(i)
  uids <- get_uid(i)
  out <- classification(uids)
  out2 <- do.call(rbind.data.frame, out)
  out2$index <- rownames(out2)
  if (out2$index!=1) {
    fin <- rbind(fin,out2)
  }
}


fin2 <- cSplit(fin, "index",".")

fin2$id <- NULL
require(reshape2)
fin3 <- reshape2::dcast(fin2, index_1 ~ index_2, value.var="name")
fin3 <- as.data.frame(fin3)

unique(fin3$taxa_simple) # 535/570 found

head(fin3)
ind <- !is.na(fin3)
fin3$taxa_simple <- tapply(fin3[ind], row(fin3)[ind], tail, 1)
fin3$index_1 <- NULL
NROW(fin3)
fwrite(x = fin3, file = "ncbi_retrieved_taxa")

unique(fin3$taxa_simple)


