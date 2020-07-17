
library(readr)
library(splitstackshape)
library(dplyr)
library(data.table)

source_dir = "/Users/12705859/Desktop/metapigs_dry/source_data/" # git 
middle_dir = "/Users/12705859/Desktop/metapigs_dry/middle_dir/" # git 
out_dir = "/Users/12705859/Desktop/metapigs_dry/gtdbtk/"   # local (because file is 14.7GB!)


# load gtdbtk assignments of the bins
gtdbtk_bins <- read_delim(paste0(out_dir,"all_concatenated_220_essential.kraken"), 
                                             "\t", escape_double = FALSE, col_names = FALSE, 
                                             trim_ws = TRUE)

gtdbtk_bins <- gtdbtk_bins %>% 
  dplyr::mutate(X3 = gsub("taxid.*", "", X3)) %>% 
  dplyr::mutate(X2 = gsub(".fa", "", X2)) 

gtdbtk_bins <- cSplit(gtdbtk_bins, "X3"," (")
gtdbtk_bins$X1 <- NULL
gtdbtk_bins$X4 <- NULL

gtdbtk_bins <- cSplit(gtdbtk_bins, "X2","_")
colnames(gtdbtk_bins) <- c("species","pig","bin")


# load gtdb dictionary
bac_arc <- read_csv(paste0(middle_dir,"bac120_arc122_dictionary"),col_names = TRUE)


# SPECIES

# retrieve upper level taxa of species that have been assigned to my bins
named_byspecies <- inner_join(gtdbtk_bins,bac_arc)
NROW(named_byspecies) # were succesfully found as species (out of 51175)
NROW(unique(named_byspecies$species)) # <- these are the numbers of unique bins that have been assigned upper taxonomic levels
NROW(unique(gtdbtk_bins$species)) # <- these are the numbers of unique bins 

NOTnamed_byspecies <- anti_join(gtdbtk_bins,bac_arc)
colnames(NOTnamed_byspecies) <- c("genus","pig","bin")
NROW(NOTnamed_byspecies) # total number of bins (out of 51175) not found as species as these were not classified as species to start with
head(NOTnamed_byspecies)

# GENUS
bac_arcG <- bac_arc %>% 
  group_by(genus) %>% 
  dplyr::filter(row_number()==1)

named_bygenus <- inner_join(NOTnamed_byspecies,bac_arcG)
named_bygenus$species <- NA
NROW(named_bygenus)
NROW(unique(named_bygenus$genus)) 

# FAMILY
bac_arcF <- bac_arc %>% 
  group_by(family) %>% 
  dplyr::filter(row_number()==1)
colnames(NOTnamed_byspecies) <- c("family","pig","bin")
named_byfamily <- inner_join(NOTnamed_byspecies,bac_arcF)
named_byfamily$species <- NA
named_byfamily$genus <- NA
NROW(named_byfamily)
NROW(unique(named_byfamily$family)) 


# ORDER
bac_arcO <- bac_arc %>% 
  group_by(order) %>% 
  dplyr::filter(row_number()==1)
colnames(NOTnamed_byspecies) <- c("order","pig","bin")
named_byorder <- inner_join(NOTnamed_byspecies,bac_arcO)
named_byorder$species <- NA
named_byorder$genus <- NA
named_byorder$family <- NA
NROW(named_byorder)
NROW(unique(named_byorder$order)) 

# CLASS
bac_arcC <- bac_arc %>% 
  group_by(class) %>% 
  dplyr::filter(row_number()==1)
colnames(NOTnamed_byspecies) <- c("class","pig","bin")
named_byclass <- inner_join(NOTnamed_byspecies,bac_arcC)
named_byclass$species <- NA
named_byclass$genus <- NA
named_byclass$family <- NA
named_byclass$order <- NA
NROW(named_byclass)
NROW(unique(named_byclass$class)) 


# PHYLUM
bac_arcP <- bac_arc %>% 
  group_by(phylum) %>% 
  dplyr::filter(row_number()==1)
colnames(NOTnamed_byspecies) <- c("phylum","pig","bin")
named_byphylum <- inner_join(NOTnamed_byspecies,bac_arcP)
named_byphylum$species <- NA
named_byphylum$genus <- NA
named_byphylum$family <- NA
named_byphylum$order <- NA
named_byphylum$class <- NA
NROW(named_byphylum)
NROW(unique(named_byphylum$phylum)) 

# domain
bac_arcD <- bac_arc %>% 
  group_by(domain) %>% 
  dplyr::filter(row_number()==1)
colnames(NOTnamed_byspecies) <- c("domain","pig","bin")
named_bydomain <- right_join(NOTnamed_byspecies,bac_arcD)
named_bydomain$species <- NA
named_bydomain$genus <- NA
named_bydomain$family <- NA
named_bydomain$order <- NA
named_bydomain$class <- NA
named_bydomain$phylum <- NA
NROW(named_bydomain)
NROW(unique(named_bydomain$domain)) 

all <- rbind(named_byspecies,
      named_bygenus,
      named_byfamily,
      named_byorder,
      named_byclass,
      named_byphylum,
      named_bydomain)

NROW(all)

NROW(gtdbtk_bins)-NROW(all)

NROW(unique(paste0(gtdbtk_bins$pig,gtdbtk_bins$bin)))
NROW(unique(paste0(all$pig,all$bin)))


fwrite(x = all, file = paste0(middle_dir,"gtdb_bins_completeTaxa"))

