
library(readr)
library(splitstackshape)
library(tidyverse)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(EnvStats)
library(ggpubr)
library(readxl)
library(data.table)
library(ape)
library(scales)
library(FSA)
library(openxlsx)
library(broom)

setwd("~/Desktop/metapigs_dry/gtdbtk")
basedir = "~/Desktop/metapigs_dry/"


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

no_reps_all <- read.csv(paste0(basedir,"no_reps_all.csv"), 
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
gtdbtk_bins <- read_csv("gtdbtk_bins_completeTaxa",
                        col_types = cols(node = col_character(),
                                         pig = col_character()))


head(gtdbtk_bins)

######################################################################
######################################################################


# upload metadata for pen info

mdat <- read_excel(paste0(basedir,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)

mdat$`*collection_date` <- as.character(mdat$`*collection_date`)
mdat$Cohort <- gsub("Sows","Sows",mdat$Cohort)

mdat2 <- mdat %>%
  dplyr::filter(!Cohort=="Sows")  %>%
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

suppl_piglets_details_mothers_weight <- read_excel("~/Desktop/metapigs_dry/suppl_piglets_details_mothers&weight.xlsx")

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


weights <- read_csv(paste0(basedir,"weights.csv"), 
                    col_types = cols(Pig = col_character(), 
                                     Room = col_character()))
colnames(weights) <- c("room","pen","pig","t0","t2","t4","t6","t8")


weights_final <- read_csv(paste0(basedir,"weights_final.csv"), 
                          col_types = cols(Pig = col_character(), 
                                           Room = col_character()))
colnames(weights_final) <- c("room","pen","pig","date","weight")
weights_final$date <- gsub("6-Mar","t10",weights_final$date)
weights_final$date <- gsub("7-Mar","t10",weights_final$date)
weights_final$date <- gsub("8-Mar","t10",weights_final$date)
weights_final$date <- gsub("9-Mar","t10",weights_final$date)
weights_final <- weights_final %>%
  dplyr::select(pig,date,weight) %>%
  filter(!date=="10-Mar") # as it's NA

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

###########################################################################################

# create workbook to add stats 

wb <- createWorkbook()

###########################################################################################


# TAXA


taxa_mat <- df0 %>%
  dplyr::select(gOTU,species,genus,family,order,class,phylum,domain) %>%
  group_by(gOTU) %>%
  slice(1)

NROW(taxa_mat)
NROW(unique(taxa_mat$gOTU))

taxa_mat_df <- as.data.frame(taxa_mat)
# to matrix 
taxa_mat <- taxa_mat_df
rownames(taxa_mat) <- taxa_mat[,1]
taxa_mat[,1] <- NULL
taxa_mat <- as.matrix(taxa_mat)
# ready 

NROW(unique(rownames(taxa_mat)))
head(taxa_mat_df)


######################################################################

# gOTU_mat

# columns to be kept 
keep <- c("cohort","pig","bin","date","value","gOTU")
df1 <- df0[ , (names(df0) %in% keep)]

# NA to zeros 
df1$value[is.na(df1$value)] <- 0

# as dates with NA was giving problems, change to class character and swap date NAs with "no-t"
df1$date <- as.character(df1$date)
df1$date[is.na(df1$date)] <- "no-t"

NROW(df1)

# sum up all the counts from the same sample (pig and date) that belong to the same OTU
df2 <- df1 %>%
  group_by(pig,date,gOTU) %>%
  dplyr::summarise(all_bins_value = sum(value))

NROW(df2)
NROW(unique(paste0(df2$pig,df2$date)))


# assign a unique sample name 
df2$sample <- paste0(df2$date,"_",df2$pig)
# remove now pig and date (redundant)
df2$pig <- NULL
df2$date <- NULL

# long to wide 
df3 <- df2 %>%
  pivot_wider(names_from = sample, values_from = all_bins_value, values_fill = list(all_bins_value = 0)) 

# check whether this list is empty(no NAs)
check_DF <- df3[rowSums(is.na(df3)) > 0,]
NROW(check_DF)


# to matrix 
gOTU_mat <- as.data.frame(df3)
rownames(gOTU_mat) <- gOTU_mat[,1]
gOTU_mat[,1] <- NULL
gOTU_mat <- as.matrix(gOTU_mat)
mode(gOTU_mat) <- "integer"
# ready 

NROW(unique(rownames(gOTU_mat)))
NROW(unique(colnames(gOTU_mat)))


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
                                                      "DScour",
                                                      "ColiGuard", 
                                                      "Neomycin",
                                                      "NeoD",
                                                      "NeoC"))

rownames(sample_df) <- sample_df[,1]

# ready


######################################################################


# create phyloseq object

gOTU = otu_table(gOTU_mat, taxa_are_rows = TRUE)
TAX = tax_table(taxa_mat)
samples = sample_data(sample_df)


############################################################################################################

######################################################################
######################################################################

my.theme <- theme(axis.title.x = element_text(size=8),
                  axis.title.y = element_text(size=8),
                  axis.text.x = element_text(size=7),
                  axis.text.y = element_text(size=7),
                  legend.text = element_blank(),
                  legend.position = "top",
                  title=element_text(size=7))
  

# ORDINATION  - effect of weight 


# t0
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
mydf <- as.data.frame(carbom.ord$vectors)
mydf$sample <- rownames(mydf)
mydf <- mydf %>%
  pivot_longer(., cols=-sample, names_to = "component", values_to = "value")
mydf <- cSplit(indt = mydf, "sample", sep = "_", drop = NA)
colnames(mydf) <- c("sample","component","value","date","pig")
mydf$pig <- as.character(mydf$pig)
mydf1 <- left_join(mydf,weights)
mydf2 <- na.omit(mydf1) 
# Spearman! 
z <- mydf2 %>% 
  nest(-component) %>% 
  mutate(cor=map(data,~cor.test(.x$value, .x$weight, method = "sp"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>%
  arrange(desc(estimate)) %>%
  dplyr::select(component,estimate,p.value,method,alternative)
head(z)
mid<-median(mydf2$weight)
t0 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight", axes = c(20,77)) + 
  geom_point(size=2)+
  scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                        high="red", space ="Lab" ) +
  ggtitle(paste0(mydf2$date[1],
                 "  ",
                 "Rho: ","x=",round(z$estimate[1],3),
                 " ",
                 "y=",round(z$estimate[2],3),
                 "   ",
                 "p-value: ","x=",round(z$p.value[1],3),
                 " ",
                 "y=",round(z$p.value[2],3)))+
  theme_bw()+
  my.theme
t0

# t2
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t2")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
mydf <- as.data.frame(carbom.ord$vectors)
mydf$sample <- rownames(mydf)
mydf <- mydf %>%
  pivot_longer(., cols=-sample, names_to = "component", values_to = "value")
mydf <- cSplit(indt = mydf, "sample", sep = "_", drop = NA)
colnames(mydf) <- c("sample","component","value","date","pig")
mydf$pig <- as.character(mydf$pig)
mydf1 <- left_join(mydf,weights)
mydf2 <- na.omit(mydf1) 
# Spearman! 
z <- mydf2 %>% 
  nest(-component) %>% 
  mutate(cor=map(data,~cor.test(.x$value, .x$weight, method = "sp"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>%
  arrange(desc(estimate)) %>%
  dplyr::select(component,estimate,p.value,method,alternative)
head(z)
head(z) 
mid<-median(mydf2$weight)
t2 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight", axes = c(6,34)) + 
  geom_point(size=2)+
  scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                        high="red", space ="Lab" )+
  ggtitle(paste0(mydf2$date[1],
                 "  ",
                 "Rho: ","x=",round(z$estimate[1],3),
                 " ",
                 "y=",round(z$estimate[2],3),
                 "   ",
                 "p-value: ","x=",round(z$p.value[1],3),
                 " ",
                 "y=",round(z$p.value[2],3)))+
  theme_bw()+
  my.theme
t2

# t4
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t4")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
mydf <- as.data.frame(carbom.ord$vectors)
mydf$sample <- rownames(mydf)
mydf <- mydf %>%
  pivot_longer(., cols=-sample, names_to = "component", values_to = "value")
mydf <- cSplit(indt = mydf, "sample", sep = "_", drop = NA)
colnames(mydf) <- c("sample","component","value","date","pig")
mydf$pig <- as.character(mydf$pig)
mydf1 <- left_join(mydf,weights)
mydf2 <- na.omit(mydf1) 
# Spearman! 
z <- mydf2 %>% 
  nest(-component) %>% 
  mutate(cor=map(data,~cor.test(.x$value, .x$weight, method = "sp"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>%
  arrange(desc(estimate)) %>%
  dplyr::select(component,estimate,p.value,method,alternative)
head(z)
mid<-median(mydf2$weight)
t4 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight", axes = c(17,5)) + 
  geom_point(size=2)+
  scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                        high="red", space ="Lab" )+
  ggtitle(paste0(mydf2$date[1],
                 "  ",
                 "Rho: ","x=",round(z$estimate[1],3),
                 " ",
                 "y=",round(z$estimate[2],3),
                 "   ",
                 "p-value: ","x=",round(z$p.value[1],3),
                 " ",
                 "y=",round(z$p.value[2],3)))+
  theme_bw()+
  my.theme
t4


# t6
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t6")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
mydf <- as.data.frame(carbom.ord$vectors)
mydf$sample <- rownames(mydf)
mydf <- mydf %>%
  pivot_longer(., cols=-sample, names_to = "component", values_to = "value")
mydf <- cSplit(indt = mydf, "sample", sep = "_", drop = NA)
colnames(mydf) <- c("sample","component","value","date","pig")
mydf$pig <- as.character(mydf$pig)
mydf1 <- left_join(mydf,weights)
mydf2 <- na.omit(mydf1) 
# Spearman! 
z <- mydf2 %>% 
  nest(-component) %>% 
  mutate(cor=map(data,~cor.test(.x$value, .x$weight, method = "sp"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>%
  arrange(desc(estimate)) %>%
  dplyr::select(component,estimate,p.value,method,alternative)
head(z)
mid<-median(mydf2$weight)
t6 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight", axes = c(15,59)) + 
  geom_point(size=2)+
  scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                        high="red", space ="Lab" )+
  ggtitle(paste0(mydf2$date[1],
                 "  ",
                 "Rho: ","x=",round(z$estimate[1],3),
                 " ",
                 "y=",round(z$estimate[2],3),
                 "   ",
                 "p-value: ","x=",round(z$p.value[1],3),
                 " ",
                 "y=",round(z$p.value[2],3)))+
  theme_bw()+
  my.theme
t6


# t8
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
mydf <- as.data.frame(carbom.ord$vectors)
mydf$sample <- rownames(mydf)
mydf <- mydf %>%
  pivot_longer(., cols=-sample, names_to = "component", values_to = "value")
mydf <- cSplit(indt = mydf, "sample", sep = "_", drop = NA)
colnames(mydf) <- c("sample","component","value","date","pig")
mydf$pig <- as.character(mydf$pig)
mydf1 <- left_join(mydf,weights)
mydf2 <- na.omit(mydf1) 
# Spearman! 
z <- mydf2 %>% 
  nest(-component) %>% 
  mutate(cor=map(data,~cor.test(.x$value, .x$weight, method = "sp"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>%
  arrange(desc(estimate)) %>%
  dplyr::select(component,estimate,p.value,method,alternative)
head(z)
# plot it! 
mid<-median(sample_data(carbom)$weight)
t8 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight", axes = c(17,41)) + 
  geom_point(size=2)+
  scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                        high="red", space ="Lab" )+
  ggtitle(paste0(mydf2$date[1],
                 "  ",
                 "Rho: ","x=",round(z$estimate[1],3),
                 " ",
                 "y=",round(z$estimate[2],3),
                 "   ",
                 "p-value: ","x=",round(z$p.value[1],3),
                 " ",
                 "y=",round(z$p.value[2],3)))+
  theme_bw()+
  my.theme
t8

# t10
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t10")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
carbom.ord <- ordinate(carbom, "PCoA", "bray")
mydf <- as.data.frame(carbom.ord$vectors)
mydf$sample <- rownames(mydf)
mydf <- mydf %>%
  pivot_longer(., cols=-sample, names_to = "component", values_to = "value")
mydf <- cSplit(indt = mydf, "sample", sep = "_", drop = NA)
colnames(mydf) <- c("sample","component","value","date","pig")
mydf$pig <- as.character(mydf$pig)
mydf1 <- left_join(mydf,weights)
mydf2 <- na.omit(mydf1) 
# Spearman! 
z <- mydf2 %>% 
  nest(-component) %>% 
  mutate(cor=map(data,~cor.test(.x$value, .x$weight, method = "sp"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>%
  arrange(desc(estimate)) %>%
  dplyr::select(component,estimate,p.value,method,alternative)
head(z) 
# plot it! 
mid<-median(mydf2$weight)

t10 <- plot_ordination(carbom, carbom.ord, type="samples", color="weight", axes = c(15,48)) + 
  geom_point(size=2)+
  scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                        high="red", space ="Lab" )+
  ggtitle(paste0(mydf2$date[1],
                 "  ",
                 "Rho: ","x=",round(z$estimate[1],3),
                 " ",
                 "y=",round(z$estimate[2],3),
                 "   ",
                 "p-value: ","x=",round(z$p.value[1],3),
                 " ",
                 "y=",round(z$p.value[2],3)))+
  theme_bw()+
  my.theme
t10


all_timepoints <- ggarrange(t0,t2,t4,t6,t8,t10,ncol=2,nrow=3,common.legend = TRUE)

pdf("gt_phylo_ordination_weight.pdf")
all_timepoints
dev.off()





