

library(readr)
library(tidyverse)
library(dplyr)
library(robCompositions)
library(microbiome)
library(phyloseq)
library(ggplot2)
library(ape)
library(splitstackshape)
library(magrittr)
library(scales)
library(ggpubr)
library(pheatmap)

library(EnvStats)
library(readxl)
library(data.table)
library(FSA)
library(openxlsx)
library(purrr)
library(broom)

library(cowplot)
library(ggsci)

library(compositions) # for clr transform 

source_dir = "/Users/12705859/metapigs_dry/source_data/" # git 
middle_dir = "/Users/12705859/metapigs_dry/middle_dir/" # git 
out_dir_git = "/Users/12705859/metapigs_dry/out/" # git 
out_dir = "/Users/12705859/Desktop/metapigs_dry/gtdbtk/"  # local


######################################################################


# input files: 
# gtdbtk_bins_completeTaxa
# no_reps_all.csv (BINS COUNTS)
# Metagenome.environmental_20190308_2.xlsx (metadata, necessary for last part, co-housing)


# OUTPUTS:

# gt_phylo_barplot.pdf
# gt_phylo_heatmap.pdf
# gt_phylo_ordination.pdf
# gt_phylo_diversity.pdf
# gt_phylo_network.pdf
# gt_phylo_heatmap_ProbioticCheck.pdf
# gt_ProbioticCheck.pdf
# gt_cohousing.pdf


######################################################################

# counts data 

no_reps_all <- read.csv(paste0(middle_dir,"no_reps_all.csv"), 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)


# remove .fa extension to match bins in checkm df 
no_reps_all$bin <- gsub(".fa","", no_reps_all$bin)
head(no_reps_all)
NROW(no_reps_all)

######################################################################

# load gtdbtk assignments of the bins

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv(paste0(middle_dir,"gtdb_bins_completeTaxa"),
                        col_types = cols(node = col_character(),
                                         pig = col_character()))


head(gtdbtk_bins)

######################################################################
######################################################################


# upload metadata for pen info

mdat <- read_excel(paste0(source_dir,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)

mdat$`*collection_date` <- as.character(mdat$`*collection_date`)

# cohort names edits: 
mdat$Cohort <- plyr::mapvalues(as.character(mdat$Cohort), 
                               from = c("Neomycin+D-scour","Neomycin+ColiGuard","D-scour"), 
                               to = c("NeoD","NeoC","D-Scour"))


mdat2 <- mdat %>%
  dplyr::filter(!Cohort=="Mothers")  %>%
  dplyr::filter(!`*collection_date`=="2017-01-31"|
                  `*collection_date`=="2017-02-01"|
                  `*collection_date`=="2017-02-03") %>%
  dplyr::select(isolation_source,PigPen)

colnames(mdat2) <- c("pig","pen")
mdat2 <- as.data.frame(mdat2)

mdat2 <- mdat2 %>%
  group_by(pig) %>%
  distinct()
NROW(mdat2)

mdat2$pen <- gsub("nan",NA,mdat2$pen)
mdat2 <- na.omit(mdat2)

# we need to keep only record of pigs that were not relocated. 
mdat2 <- setDT(mdat2)[,if(.N ==1) .SD,by=pig]


######################################################################


# upload breed and bday info 

suppl_piglets_details_mothers_weight <- read_excel(paste0(source_dir,"suppl_piglets_details_mothers&weight.xlsx"))

# select cols of interest
breed_bday <- suppl_piglets_details_mothers_weight %>%
  dplyr::select(TATTOO,BIRTH_DAY,...8,`Nursing Dam`,STIGDAM)

# rename columns
colnames(breed_bday) <- c("pig","birth_day","breed","nurse_mother","mother")

breed_bday$birth_day <- as.character(breed_bday$birth_day)

# clean names
breed_bday$pig <- gsub("G","", breed_bday$pig)
breed_bday$pig <- gsub("T","", breed_bday$pig)

breed_bday <- as.data.frame(breed_bday)

######################################################################

# upload weight info 


weights <- read_csv(paste0(source_dir,"weights.csv"), 
                    col_types = cols(Pig = col_character(), 
                                     Room = col_character()))
colnames(weights) <- c("room","pen","pig","t0","t2","t4","t6","t8")


weights_final <- read_csv(paste0(source_dir,"weights_final.csv"), 
                          col_types = cols(Pig = col_character(), 
                                           Room = col_character()))
colnames(weights_final) <- c("room","pen","pig","date","weight")
weights_final$date <- gsub("6-Mar","t10",weights_final$date)
weights_final$date <- gsub("7-Mar","t10",weights_final$date)
weights_final$date <- gsub("8-Mar","t10",weights_final$date)
weights_final$date <- gsub("9-Mar","t10",weights_final$date)

weights_final <- weights_final %>%
  dplyr::select(pig,date,weight) %>%
  dplyr::filter(!date=="10-Mar") # as it's NA

weights <- weights %>%
  dplyr::select(pig,t0,t2,t4,t6,t8) %>%
  pivot_longer(., cols = c(t0,t2,t4,t6,t8), names_to = "date", values_to = "weight")

weights <- as.data.frame(weights)

weights <- rbind(weights,weights_final)
NROW(weights)



######################################################################

# merge bins info to gtdbtk assignment info :  

NROW(gtdbtk_bins)
NROW(no_reps_all)
head(gtdbtk_bins)
head(no_reps_all)
df0 <- merge(no_reps_all, gtdbtk_bins, by=c("pig","bin"))

# rename node as gOTU and place "gOTU_" in front of node number: a separate genomic OTU identifier for each different genome

colnames(df0)[colnames(df0) == 'node'] <- 'gOTU'
df0$gOTU <- paste0("gOTU_",df0$gOTU)

NROW(unique(df0$gOTU))
NROW(df0)

######################################################################

# merge all other info: 
# add pen info (mdat2), breed and bday info (breed_bday) and weight info (weights)

# add breed and bday info (breed_bday)
df0 <- left_join(df0,breed_bday)
NROW(df0)

# add pen info (mdat2)
df0 <- left_join(df0,mdat2)
NROW(df0)

# add weight info (weights)
df0 <- left_join(df0,weights)
NROW(df0)




























######################################################################

# gOTU_mat

# columns to be kept 
keep <- c("cohort","pig","bin","date","value","species","genus","family","class")
df1 <- df0[ , (names(df0) %in% keep)]

# NA to zeros 
df1$value[is.na(df1$value)] <- 0

# as dates with NA was giving problems, change to class character and swap date NAs with "no-t"
df1$date <- as.character(df1$date)
df1$date[is.na(df1$date)] <- "no-t"

NROW(df1)

# sum up all the counts from the same sample (pig and date) that belong to the same OTU
df2 <- df1 %>%
  group_by(pig,date,species,genus,family,class) %>%
  dplyr::summarise(all_bins_value = sum(value)) 

# clr transform per sample
df2 <- df2 %>%
  group_by(pig,date) %>%
  dplyr::mutate(all_bins_value = as.numeric(clr(all_bins_value)))


######################################################################

# SAMPLES 

# create a sample table for phyloseq 

sample_df <- df0

sample_df$sample <- paste0(sample_df$date,"_",sample_df$pig)
NROW(unique(sample_df$sample))

sample_df <- sample_df %>%
  dplyr::select(sample,pig,date,cohort,pen,birth_day,breed,weight) %>%
  group_by(sample) %>%
  slice(1)

unique(sample_df$breed)

sample_df$breed <- gsub("Landrace x Cross bred [(]LW x D[])]","LxLWD", sample_df$breed)
sample_df$breed <- gsub("Duroc x Landrace","DxL", sample_df$breed)
sample_df$breed <- gsub("Duroc x Large white","DxLW", sample_df$breed)
sample_df$breed <- gsub("Large white x Duroc","LWxD", sample_df$breed)

sample_df$gOTU <- NULL
sample_df <- as.data.frame(sample_df)

NROW(sample_df)
head(sample_df)

# reorder dates 
sample_df$date  = factor(sample_df$date, levels=c("t0",
                                                  "t1", 
                                                  "t2",
                                                  "t3",
                                                  "t4",
                                                  "t5",
                                                  "t6",
                                                  "t7",
                                                  "t8",
                                                  "t9",
                                                  "t10",
                                                  "no-t"))

# reorder cohorts 
sample_df$cohort  = factor(sample_df$cohort, levels=c("Control", 
                                                      "D-Scour",
                                                      "ColiGuard", 
                                                      "Neomycin",
                                                      "NeoD",
                                                      "NeoC",
                                                      "Mothers"))


# ready

df2
sample_df


all <- inner_join(sample_df,df2) %>% 
  dplyr::filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"| # select those with weight data 
                  date=="t8"|date=="t10") %>%
  dplyr::select(weight,pig,date,all_bins_value,species,genus,family,class) %>% drop_na()

# reorder dates 
all$date  = factor(all$date, levels=c("t0","t2","t4","t6","t8","t10"))

######################################################################

# # subset to frequent species
# all$n=1
# all_frequent <- all %>% 
#   group_by(species,date) %>% 
#   dplyr::mutate(num=sum(n)) %>% 
#   dplyr::filter(num>20) 


last <- data.frame() # empty dataframe
require(plyr)
func <- function(xx)
{
  if (NROW(xx) > 10) {
    f <- cor.test(xx$weight, xx$all_bins_value, method = "pearson")
    pval <- f$p.value
    cor <- f$estimate
    species <- unique(xx$species)
    out <- data.frame(species,cor,pval)
    last <- rbind(last,out)
    return(out)
  }
}

res <- ddply(all, .(date,species), func)
head(res)
NROW(res)

# filter out weak correlations 
ress <- res %>% 
  dplyr::filter(pval<0.01) %>%
  dplyr::filter(!between(cor, -0.3, 0.3)) %>%
  dplyr::arrange(pval)

NROW(ress)
View(ress)

# fish out from the ress output 
# selecting species and plot 
pdf(paste0(out_dir,"gt_corr_weight_species.pdf"))
for (row in 1:nrow(ress)) {
  
  speciess <- ress[row,2]
  all0 <- subset(all, species %in% speciess) 
  
  print(all0 %>%
          ggplot(., aes(x=weight,y=all_bins_value))+
          geom_point(size=0.5) +
          facet_wrap(~date, scale="free", shrink = TRUE)+
          stat_smooth(method="lm", se=TRUE) +
          stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
          ggtitle(paste0(unique(all0$species),"__",
                         unique(all0$family),"__",
                         unique(all0$class))))
  
}
dev.off()

############################################################

# abundance of which taxa correlate with final weight? 


end_weight <- all %>%
  dplyr::filter(date=="t10") %>%
  dplyr::select(pig,weight,species,genus,family,class)

start_abund <- all %>%
  dplyr::filter(date=="t0") %>%
  dplyr::select(pig,all_bins_value,species,genus,family,class)

start_end <- inner_join(start_abund,end_weight) 


start_end_res <- ddply(start_end, .(species), func)
head(start_end_res)
NROW(start_end_res)

# filter out weak correlations 
start_end_ress <- start_end_res %>% 
  dplyr::filter(pval<0.05) %>%
  dplyr::filter(!between(cor, -0.3, 0.3)) %>%
  dplyr::arrange(pval)

NROW(start_end_ress)
head(start_end_ress)

# fish out from the ress output 
# selecting species and plot 
pdf(paste0(out_dir,"gt_STARTcorr_ENDweight_species.pdf"))
for (row in 1:nrow(start_end_ress)) {
  
  speciess <- start_end_ress[row,1]
  all0 <- subset(start_end, species %in% speciess) 
  
  print(all0 %>%
          ggplot(., aes(x=weight,y=all_bins_value))+
          geom_point(size=0.5) +
          stat_smooth(method="lm", se=TRUE) +
          stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
          ggtitle(paste0(unique(all0$species),"__",
                         unique(all0$family),"__",
                         unique(all0$class))))
  
}
dev.off()