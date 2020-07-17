
# AIM: create a dictionary of gtdb species with all the upper levels taxonomy

# this dictionary has been created from: 
# https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_sp_labels.tree
# https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar122.sp_labels.tree

# previous steps: 
# 1. imported archaea and bacterial dtdb databases into Archaeopteryx 
# 2. saved both as phyloXML files

# import libs: 
library(readr)
library(data.table)
library(dplyr)
library(zoo)

setwd("~/Desktop/metapigs_dry/gtdbtk")

# input files: 
arc_phyloXML <- read_csv("arc_phyloXML.txt", col_names = FALSE)
bac_phyloXML <- read_csv("bac_phyloXML.txt", col_names = FALSE)

# extract only rows containing names 
arc <- arc_phyloXML$X1[arc_phyloXML$X1 %like% "<name>" ]
bac <- bac_phyloXML$X1[bac_phyloXML$X1 %like% "<name>" ]

# list to extract rows containing names that starts with: 
l <- c("d__|p__|c__|o__|f__|g__|s__")

arc <- arc[arc %like% l ]
bac <- bac[bac %like% l ]

arc2 <- as.data.frame(arc)
bac2 <- as.data.frame(bac)

# cleaning 
arc2$arc <- arc2$arc %<>%
  stringr::str_remove_all("<name>") %>%
  stringr::str_remove_all("</name>") %>%
  gsub('.*:', '', .)  # removes everything up to __ (keeping only the most specific)
bac2$bac <- bac2$bac %<>%
  stringr::str_remove_all("<name>") %>%
  stringr::str_remove_all("</name>") %>%
  gsub('.*:', '', .)  # removes everything up to __ (keeping only the most specific)


# some are collapsed, for example: one row: c__Halobacteria; o__Halobacteriales
arc2 <- separate_rows(arc2, arc, convert = TRUE, sep = ";")
bac2 <- separate_rows(bac2, bac, convert = TRUE, sep = ";")

# preparing new (empty) columns
arc2$domain <- NA
arc2$phylum <- NA
arc2$class <- NA
arc2$order <- NA
arc2$family <- NA
arc2$genus <- NA
arc2$species <- NA
bac2$domain <- NA
bac2$phylum <- NA
bac2$class <- NA
bac2$order <- NA
bac2$family <- NA
bac2$genus <- NA
bac2$species <- NA


# place all in the right place 
arc3 <- arc2 %>%
  mutate(., genus=ifelse(grepl('^g__',arc), as.character(arc), "NA"))
arc3 <- arc3 %>%
  mutate(., species=ifelse(grepl('^s__',arc), as.character(arc), "NA"))
arc3 <- arc3 %>%
  mutate(., domain=ifelse(grepl('^d__',arc), as.character(arc), "NA"))
arc3 <- arc3 %>%
  mutate(., phylum=ifelse(grepl('^p__',arc), as.character(arc), "NA"))
arc3 <- arc3 %>%
  mutate(., class=ifelse(grepl('^c__',arc), as.character(arc), "NA"))
arc3 <- arc3 %>%
  mutate(., order=ifelse(grepl('o__',arc), as.character(arc), "NA"))
arc3 <- arc3 %>%
  mutate(., family=ifelse(grepl('^f__',arc), as.character(arc), "NA"))
bac3 <- bac2 %>%
  mutate(., genus=ifelse(grepl('^g__',bac), as.character(bac), "NA"))
bac3 <- bac3 %>%
  mutate(., species=ifelse(grepl('^s__',bac), as.character(bac), "NA"))
bac3 <- bac3 %>%
  mutate(., domain=ifelse(grepl('^d__',bac), as.character(bac), "NA"))
bac3 <- bac3 %>%
  mutate(., phylum=ifelse(grepl('^p__',bac), as.character(bac), "NA"))
bac3 <- bac3 %>%
  mutate(., class=ifelse(grepl('^c__',bac), as.character(bac), "NA"))
bac3 <- bac3 %>%
  mutate(., order=ifelse(grepl('o__',bac), as.character(bac), "NA"))
bac3 <- bac3 %>%
  mutate(., family=ifelse(grepl('^f__',bac), as.character(bac), "NA"))

# replace literal NA with actual NA
arc3[] <- lapply(arc3, gsub, pattern = "NA", replacement = NA, fixed = TRUE)
bac3[] <- lapply(bac3, gsub, pattern = "NA", replacement = NA, fixed = TRUE)

arc4 <- arc3
bac4 <- bac3

# when NA in a row, fill with value from above 
arc4 <- arc4 %>% 
  mutate(phylum = na.locf(phylum, fromLast=FALSE, na.rm = FALSE))
arc4 <- arc4 %>% 
  mutate(class = na.locf(class, fromLast=FALSE, na.rm = FALSE))
arc4 <- arc4 %>% 
  mutate(order = na.locf(order, fromLast=FALSE, na.rm = FALSE))
arc4 <- arc4 %>% 
  mutate(family = na.locf(family, fromLast=FALSE, na.rm = FALSE))
arc4 <- arc4 %>% 
  mutate(genus = na.locf(genus, fromLast=FALSE, na.rm = FALSE))
arc4 <- arc4 %>% 
  mutate(domain = na.locf(domain, fromLast=FALSE, na.rm = FALSE))

bac4 <- bac4 %>% 
  mutate(phylum = na.locf(phylum, fromLast=FALSE, na.rm = FALSE))
bac4 <- bac4 %>% 
  mutate(class = na.locf(class, fromLast=FALSE, na.rm = FALSE))
bac4 <- bac4 %>% 
  mutate(order = na.locf(order, fromLast=FALSE, na.rm = FALSE))
bac4 <- bac4 %>% 
  mutate(family = na.locf(family, fromLast=FALSE, na.rm = FALSE))
bac4 <- bac4 %>% 
  mutate(genus = na.locf(genus, fromLast=FALSE, na.rm = FALSE))
bac4 <- bac4 %>% 
  mutate(domain = na.locf(domain, fromLast=FALSE, na.rm = FALSE))


# remove all the rows where species is NA 
archaea <- arc4[!is.na(arc4$species),]
bacteria <- bac4[!is.na(bac4$species),]

# remove now useless column 
archaea$arc <- NULL
bacteria$bac <- NULL

# join dataframes 
bac120_arc122_dictionary <- rbind(archaea, bacteria)
NROW(bac120_arc122_dictionary)
View(bac120_arc122_dictionary)


bac120_arc122_dictionary <- bac120_arc122_dictionary %>% 
  mutate(domain = gsub("d__", "", domain)) %>% 
  mutate(phylum = gsub("p__", "", phylum)) %>% 
  mutate(class = gsub("c__", "", class)) %>% 
  mutate(order = gsub("o__", "", order)) %>% 
  mutate(family = gsub("f__", "", family)) %>% 
  mutate(genus = gsub("g__", "", genus)) %>% 
  mutate(species = gsub("s__", "", species)) %>% 
  mutate(node = row_number()) # I am assigning a row number a tree node, to get a rough idea of how 
                                    # distant two nodes are when, for instance, species names won't tell 
                                          # (thinking of plots and networks)

fwrite(x = bac120_arc122_dictionary, file = "bac120_arc122_dictionary")


