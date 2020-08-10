
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

###########################################################################################


# TAXA


taxa_mat <- df0 %>%
  dplyr::select(gOTU,species,genus,family,order,class,phylum,domain) %>%
  group_by(gOTU) %>%
  dplyr::slice(1)

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
                                                      "D-Scour",
                                                      "ColiGuard", 
                                                      "Neomycin",
                                                      "NeoD",
                                                      "NeoC",
                                                      "Mothers"))

rownames(sample_df) <- sample_df[,1]

# ready


######################################################################


# create phyloseq object

gOTU = otu_table(gOTU_mat, taxa_are_rows = TRUE)
TAX = tax_table(taxa_mat)
samples = sample_data(sample_df)


############################################################################################################
############################################################################################################
############################################################################################################

# PLOT

######################


# ORDINATION 

# NORMALIZATION BY MEDIAN SEQUENCING DEPTH
carbom <- phyloseq(gOTU,TAX,samples)
# SUBSETTING phyloseq obejct
carbom <- subset_samples(carbom, (date %in% c("t0","t1","t2","t3","t4","t5","t6","t7","t8","t9")))
# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)


# Remove taxa not seen more than 10 times in at least a 30% of the samples 
carbom_abund <- filter_taxa(carbom, function(x) sum(x > 10) > (0.3*length(x)), TRUE)


carbom_abund.ord <- ordinate(carbom_abund, "NMDS", "bray")

gt_ordination_plot <- plot_ordination(carbom_abund, carbom_abund.ord, type="samples", color="date") + 
  geom_point(size=1) +
  theme_bw()+
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=7))+
  guides(colour = guide_legend(nrow = 1))+
  theme(legend.position="top")

pdf(paste0(out_dir,"gt_phylo_ordination.pdf"))
gt_ordination_plot #+
  #facet_wrap(~cohort)
dev.off()


########################################################################################

# NETWORK ANALYSIS 

ig = make_network(carbom_abund, type = "samples", distance = "bray", max.dist = 0.55)
gt_network_plot <- plot_network(ig, carbom_abund, color = "date", shape = "cohort", line_weight = 0.3, 
                                label = NULL, point_size = 1)+
  theme(legend.position = "bottom")+
  guides(shape = guide_legend(nrow = 1))+
  guides(size = "legend", colour = "none")

pdf(paste0(out_dir,"gt_phylo_network.pdf"))
gt_network_plot
dev.off()


##############################################


# BAR PLOT

# BAR GRAPH - by time point
pdf(paste0(out_dir,"gt_phylo_barplot_time.pdf"))
plot_bar(carbom_abund, fill = "phylum") +
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack") +
  facet_grid(~date,scales="free_x") +
  theme(axis.text.x = element_blank())
plot_bar(carbom_abund, fill = "class") +
  geom_bar(aes(color=class, fill=class), stat="identity", position="stack") +
  facet_grid(~date,scales="free_x") +
  theme(axis.text.x = element_blank())
dev.off()


##############################################


# HEATMAP

# NORMALIZATION BY RAREFACTION
carbom <- phyloseq(gOTU,TAX,samples)
# SUBSETTING phyloseq obejct
carbom <- subset_samples(carbom, (date %in% c("t0","t2","t4","t6","t8","t10")))

# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE,
                           rngseed = 42)

# Remove taxa not seen more than 5 times in at least a 20% of the samples 
carbom_abund <- filter_taxa(carbom, function(x) sum(x > 5) > (0.2*length(x)), TRUE)


random_tree = rtree(ntaxa(carbom_abund), rooted=TRUE, tip.label=taxa_names(carbom_abund))
physeq1 = merge_phyloseq(carbom_abund,random_tree)


# HEATMAP time - genus, family, order, etc ...
pdf(paste0(out_dir,"gt_phylo_heatmap.pdf"))
plot_heatmap(physeq1, 
             taxa.label = "species", 
             sample.order = "date") +
  facet_grid(~ date, switch = "x", scales = "free_x", space = "free_x")+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "Microbe Species Diversity - log transformed data")
plot_heatmap(physeq1, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "species", taxa.order="species",
             sample.order = "date", trans=identity_trans(),
             low="blue", high="red", na.value="white") +
  facet_grid(~ date, switch = "x", scales = "free_x", space = "free_x")+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "Microbe Species Diversity")
dev.off()


# HEATMAP time - cohorts
pdf(paste0(out_dir,"gt_phylo_heatmap_cohorts.pdf"))
plot_heatmap(physeq1, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "species", taxa.order="species",
             sample.order = "date", trans=identity_trans(),
             low="blue", high="red", na.value="white") +
  facet_wrap(~ cohort, switch = "x", scales = "free_x")+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "Microbe Species Diversity")
plot_heatmap(physeq1, 
             taxa.label = "species", 
             sample.order = "date") +
  facet_wrap(~ cohort, switch = "x", scales = "free_x")+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "Microbe Species Diversity")
dev.off()

##############################################

# this is a small break in between the phyloseq analysis: 

# DISPLAYING MOST ABUNDANT SPECIES TIME-TREND

# this is a zoom in on the species shown in phylo_heatmap because most abundant 
# (taking gOTUs filtering selection from above

# 1. raw data is normalized by lib size and 
# 2. the mean is taken from bins assigned the same species

# 3. at this point those species are now selected and plotted using the log10

physeq1
keep_these_taxa <- as.data.frame(tax_table(physeq1))$species
keep_these_taxa <- as.character(keep_these_taxa)

# normalization and sum of same species all in one 
z <- df0 %>%
  dplyr::filter(!cohort=="Mothers") %>%
  group_by(pig,date) %>%
  dplyr::mutate(norm_value=value/sum(value)) %>% # lib size normalization
  group_by(cohort,pig,date,species) %>%
  dplyr::summarize(z = mean(norm_value)) # mean of bins falling within same species

z <- subset(z, (species %in% keep_these_taxa))
NROW(unique(z$species))

z <- as.data.frame(z)

# reorder dates 
z$date  = factor(z$date, levels=c("t0",
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
z$cohort  = factor(z$cohort, levels=c("Control",
                                      "D-Scour", 
                                      "ColiGuard",
                                      "Neomycin",
                                      "NeoD",
                                      "NeoC",
                                      "Mothers"))


# split df by species
multiple_DFs <- split( z , f = z$species )

pdf(paste0(out_dir,"gt_zoomIN_on_phylo_heatmap.pdf"))
for (single_DF in multiple_DFs) {
  
  DF <- as.data.frame(single_DF)
  this_species <- DF$species
  
  p <- ggplot(DF, aes(x=date, y=log10(z), fill=cohort)) + 
    geom_boxplot(outlier.shape = NA) +
    #geom_jitter(size=0.2)+       # keep this if you want to see how many samples behind measurement
    #geom_line() + geom_point(size=0.8)+
    theme_bw()+
    theme(legend.position="right")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "collection date",
         y = "relative abundance (log10)",
         color = "Phylum")  +
    theme(legend.title = element_text()) +
    ggtitle(this_species)
  
  print(p)
  
}
dev.off()



######################


# LARGE HEATMAP

# Remove taxa not seen more than 1 times in at least a 20% of the samples 
carbom_abund <- filter_taxa(carbom, function(x) sum(x > 1) > (0.2*length(x)), TRUE)


random_tree = rtree(ntaxa(carbom_abund), rooted=TRUE, tip.label=taxa_names(carbom_abund))
physeq1 = merge_phyloseq(carbom_abund,random_tree)

# HEATMAP time - genus, family, order, etc ...
pdf(paste0(out_dir,"gt_phylo_heatmap_large.pdf"))
plot_heatmap(physeq1,
             taxa.label = "species", 
             sample.order = "date") +
  facet_grid(~ date, switch = "x", scales = "free_x", space = "free_x")+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = paste0("Microbe Species Diversity in the piglet population \n(log transformed data) \n (n=",
                         NROW(tax_table(carbom_abund)),
                         " GTDB species)"))
dev.off()




######################

# DIVERSITY 

# NORMALIZATION BY RAREFACTION
carbom <- phyloseq(gOTU,TAX,samples)
# SUBSETTING phyloseq obejct
carbom <- subset_samples(carbom, (date %in% c("t0","t1","t2","t3","t4","t5", "t6","t7","t8","t9","t10")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)

gt_diversity_samples <- plot_richness(carbom, measures=c("Chao1","Shannon", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), 
                                      color="date", x="date") +
  guides(colour = guide_legend(nrow = 1))+
  theme(legend.position="top")

######### plotting above results in a different way: 
# focus on three measures;
# whisker plots instead 

my_comparisons <- data.frame(my_comparisons=c("t0_t2", "t2_t4", "t4_t6", "t6_t8", "t8_t10"))


Chao1_data <- gt_diversity_samples$data %>%
  filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="t10") %>%
  filter(variable=="Chao1")
Chao1_data_test <- as.data.frame(compare_means(value ~ date,
                                               p.adjust.method = "bonferroni", 
                                               method='t.test', data = Chao1_data) %>%
                                   dplyr::mutate(my_comparisons=paste0(group1,"_",group2))) %>%
  inner_join(., my_comparisons) %>%
  dplyr::mutate(y.position=c(300,320,340,360,380))
Chao1_plot <- Chao1_data %>%
  ggplot(., aes(x=date, y=value, color=date)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=1)+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "time point",
       y = "Chao1") +
  theme_minimal() +
  stat_pvalue_manual(Chao1_data_test, label = "p.adj")


Shannon_data <- gt_diversity_samples$data %>%
  dplyr::filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="t10") %>%
  dplyr::filter(variable=="Shannon")
Shannon_data_test <- as.data.frame(compare_means(value ~ date, 
                                                 p.adjust.method = "bonferroni", 
                                                 method='t.test', data = Shannon_data) %>%
                                     dplyr::mutate(my_comparisons=paste0(group1,"_",group2))) %>%
  inner_join(., my_comparisons) %>%
  dplyr::mutate(y.position=c(5.2,5.4,5.6,5.8,6.0))
Shannon_plot <- Shannon_data %>%
  ggplot(., aes(x=date, y=value, color=date)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=1)+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "time point",
       y = "Shannon") +
  theme_minimal() +
  stat_pvalue_manual(Shannon_data_test, label = "p.adj")


Simpson_data <- gt_diversity_samples$data %>%
  dplyr::filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="t10") %>%
  dplyr::filter(variable=="Simpson")
Simpson_data_test <- as.data.frame(compare_means(value ~ date, 
                                                 p.adjust.method = "bonferroni", 
                                                 method='t.test', data = Simpson_data) %>%
                                     mutate(my_comparisons=paste0(group1,"_",group2))) %>%
  inner_join(., my_comparisons) %>%
  dplyr::mutate(y.position=c(1.02,1.04,1.06,1.08,1.1))
Simpson_plot <- Simpson_data %>%
  ggplot(., aes(x=date, y=value, color=date)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=1)+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "time point",
       y = "Simpson") +
  theme_minimal() +
  stat_pvalue_manual(Simpson_data_test, label = "p.adj")

tosave <- ggarrange(Chao1_plot, 
                    Shannon_plot,
                    Simpson_plot,
                    ncol = 3, 
                    nrow=1,
                    labels=c("A","B","C"), 
                    common.legend = TRUE)

ggsave(filename = paste0(out_dir,"gt_phylo_diversity_boxplot.pdf"), plot = tosave)


######################################################################
######################################################################


# ordination - CO-HOUSING : 


# Control 
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("Control")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)

# Remove taxa not seen more than 10 times in at least a 50% of the samples 
carbom <- filter_taxa(carbom, function(x) sum(x > 10) > (0.5*length(x)), TRUE)


carbom.ord <- ordinate(carbom, "PCoA", "bray")
ctrlt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="Control") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("Control")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
# Remove taxa not seen more than 10 times in at least a 50% of the samples 
carbom <- filter_taxa(carbom, function(x) sum(x > 10) > (0.5*length(x)), TRUE)

carbom.ord <- ordinate(carbom, "PCoA", "bray")
ctrlt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="Control") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# DScour 
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("D-Scour")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
# Remove taxa not seen more than 10 times in at least a 50% of the samples 
carbom <- filter_taxa(carbom, function(x) sum(x > 10) > (0.5*length(x)), TRUE)

carbom.ord <- ordinate(carbom, "PCoA", "bray")
dscourt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                            title="D-Scour") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("D-Scour")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
# Remove taxa not seen more than 10 times in at least a 50% of the samples 
carbom <- filter_taxa(carbom, function(x) sum(x > 10) > (0.5*length(x)), TRUE)

carbom.ord <- ordinate(carbom, "PCoA", "bray")
dscourt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                            title="D-Scour") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# ColiGuard 
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("ColiGuard")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
# Remove taxa not seen more than 10 times in at least a 50% of the samples 
carbom <- filter_taxa(carbom, function(x) sum(x > 10) > (0.5*length(x)), TRUE)

carbom.ord <- ordinate(carbom, "PCoA", "bray")
coligt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                           title="ColiGuard") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("ColiGuard")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
# Remove taxa not seen more than 10 times in at least a 50% of the samples 
carbom <- filter_taxa(carbom, function(x) sum(x > 10) > (0.5*length(x)), TRUE)

carbom.ord <- ordinate(carbom, "PCoA", "bray")
coligt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                           title="ColiGuard") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# Neomycin
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("Neomycin")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
# Remove taxa not seen more than 10 times in at least a 50% of the samples 
carbom <- filter_taxa(carbom, function(x) sum(x > 10) > (0.5*length(x)), TRUE)

carbom.ord <- ordinate(carbom, "PCoA", "bray")
neot0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                         title="Neomycin") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("Neomycin")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
# Remove taxa not seen more than 10 times in at least a 50% of the samples 
carbom <- filter_taxa(carbom, function(x) sum(x > 10) > (0.5*length(x)), TRUE)

carbom.ord <- ordinate(carbom, "PCoA", "bray")
neot8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                         title="Neomycin") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# NeoD 
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("NeoD")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
# Remove taxa not seen more than 10 times in at least a 50% of the samples 
carbom <- filter_taxa(carbom, function(x) sum(x > 10) > (0.5*length(x)), TRUE)

carbom.ord <- ordinate(carbom, "PCoA", "bray")
neoDt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="NeoD") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("NeoD")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
# Remove taxa not seen more than 10 times in at least a 50% of the samples 
carbom <- filter_taxa(carbom, function(x) sum(x > 10) > (0.5*length(x)), TRUE)

carbom.ord <- ordinate(carbom, "PCoA", "bray")
neoDt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="NeoD") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
# NeoC 
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t0")))
carbom <- subset_samples(carbom, (cohort %in% c("NeoC")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
# Remove taxa not seen more than 10 times in at least a 50% of the samples 
carbom <- filter_taxa(carbom, function(x) sum(x > 10) > (0.5*length(x)), TRUE)

carbom.ord <- ordinate(carbom, "PCoA", "bray")
neoCt0 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="NeoC") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################
carbom <- phyloseq(gOTU,TAX,samples)
carbom <- subset_samples(carbom, (date %in% c("t8")))
carbom <- subset_samples(carbom, (cohort %in% c("NeoC")))
# RAREFY
carbom = rarefy_even_depth(carbom,
                           replace=TRUE, 
                           rngseed = 42)
# Remove taxa not seen more than 10 times in at least a 50% of the samples 
carbom <- filter_taxa(carbom, function(x) sum(x > 10) > (0.5*length(x)), TRUE)

carbom.ord <- ordinate(carbom, "PCoA", "bray")

neoCt8 <- plot_ordination(carbom, carbom.ord, type="samples", color="pen", shape="pen",
                          title="NeoC") + 
  geom_point(size=2) +
  facet_wrap(~date, scales = c("free"))
######################


ctrl <- ggarrange(ctrlt0,ctrlt8)
neo <- ggarrange(neot0,neot8)
dscour <- ggarrange(dscourt0,dscourt8)
colig <- ggarrange(coligt0,coligt8)
neoD <- ggarrange(neoDt0,neoDt8)
neoC <- ggarrange(neoCt0,neoCt8)

pdf(paste0(out_dir,"gt_phylo_ordination_cohousing.pdf"))
ggarrange(
  ctrl, neo, nrow=2, labels=c("A","B")
)
ggarrange(
  dscour, colig, nrow=2, labels=c("A","B")
)
ggarrange(
  neoD, neoC, nrow=2, labels=c("A","B")
)
dev.off()


# function to test correlation with Dunn.test

pen_function <- function(myd){
  
  x_axis1 <- dunnTest(myd$Axis.1~pen,data=myd,method = "bonferroni")
  x_axis1 <- as.data.frame(x_axis1$res)
  x_axis1$axis = "Axis.1"
  
  x_axis2 <- dunnTest(myd$Axis.2~pen,data=myd,method = "bonferroni")
  x_axis2 <- as.data.frame(x_axis2$res)
  x_axis2$axis = "Axis.2"
  
  x <- rbind(x_axis1,x_axis2)
  
}

out1 <- pen_function(ctrlt0$data)
out1$type = "ctrlt0"
out2 <- pen_function(ctrlt8$data)
out2$type = "ctrlt8"
out3 <- pen_function(dscourt0$data)
out3$type = "dscourt0"
out4 <- pen_function(dscourt8$data)
out4$type = "dscourt8"
out5 <- pen_function(coligt0$data)
out5$type = "coligt0"
out6 <- pen_function(coligt8$data)
out6$type = "coligt8"
out7 <- pen_function(neot0$data)
out7$type = "neot0"
out8 <- pen_function(neot8$data)
out8$type = "neot8"
out9 <- pen_function(neoDt0$data)
out9$type = "neoDt0"
out10 <- pen_function(neoDt8$data)
out10$type = "neoDt8"
out11 <- pen_function(neoCt0$data)
out11$type = "neoCt0"
out12 <- pen_function(neoCt8$data)
out12$type = "neoCt8"


pen_stats <- rbind(out1,out2,out3,out4,out5,out6,
                   out7,out8,out9,out10,out11,out12)
pen_stats$test = "Dunn.test"
pen_stats$correction = "Bonferroni"

# save 
fwrite(x=pen_stats, file = paste0(out_dir_git,"gtdb_cohousing.csv"), sep = ",")


######################################################################
######################################################################


# ORDINATION  - effect of weight 



my.theme <- theme(axis.title.x = element_text(size=8),
                  axis.title.y = element_text(size=8),
                  axis.text.x = element_text(size=7),
                  axis.text.y = element_text(size=7),
                  legend.text = element_blank(),
                  legend.position = "top",
                  title=element_text(size=7))


# ORDINATION  - effect of weight 


weight_corr_plot <- function(ccc) {
  
  # RAREFY
  ccc = rarefy_even_depth(ccc,
                          replace=TRUE, 
                          rngseed = 42)
  carbom.ord <- ordinate(ccc, "PCoA", "bray")
  mydf <- as.data.frame(carbom.ord$vectors)
  mydf$sample <- rownames(mydf)
  mydf <- mydf %>%
    pivot_longer(., cols=-sample, names_to = "component", values_to = "value")
  mydf <- cSplit(indt = mydf, "sample", sep = "_", drop = NA)
  colnames(mydf) <- c("sample","component","value","date","pig")
  mydf$pig <- as.character(mydf$pig)
  mydf1 <- left_join(mydf,weights)
  mydf2 <- na.omit(mydf1) 
  
  mid<-median(mydf2$weight)
  
  # Spearman! 
  z <- mydf2 %>% 
    nest(-component) %>% 
    dplyr::mutate(cor=map(data,~cor.test(.x$value, .x$weight, method = "sp")),
                  tidied = map(cor, tidy)) %>% 
    unnest(tidied, .drop = T) %>%
    dplyr::arrange(p.value) %>%
    dplyr::select(component,estimate,p.value,method,alternative)
  
  # plot it! 
  the_weight_plot <- plot_ordination(ccc, carbom.ord, type="samples", color="weight", 
                                     axes = c(as.numeric(gsub("Axis.","",z$component[1])),
                                              as.numeric(gsub("Axis.","",z$component[2])))) + 
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
  
  return(the_weight_plot)
}

weight_corr_output <- function(ccc) {
  
  # RAREFY
  ccc = rarefy_even_depth(ccc,
                          replace=TRUE, 
                          rngseed = 42)
  carbom.ord <- ordinate(ccc, "PCoA", "bray")
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
    dplyr::mutate(cor=map(data,~cor.test(.x$value, .x$weight, method = "sp")),
                  tidied = map(cor, tidy)) %>% 
    unnest(tidied, .drop = T) %>%
    dplyr::arrange(p.value) %>%
    dplyr::select(component,estimate,p.value,method,alternative)
  
  return(z)
}



all_timepoints <- ggarrange(weight_corr_plot(subset_samples(phyloseq(gOTU,TAX,samples), date %in% c("t0"))),
                            weight_corr_plot(subset_samples(phyloseq(gOTU,TAX,samples), date %in% c("t2"))),
                            weight_corr_plot(subset_samples(phyloseq(gOTU,TAX,samples), date %in% c("t4"))),
                            weight_corr_plot(subset_samples(phyloseq(gOTU,TAX,samples), date %in% c("t6"))),
                            weight_corr_plot(subset_samples(phyloseq(gOTU,TAX,samples), date %in% c("t8"))),
                            weight_corr_plot(subset_samples(phyloseq(gOTU,TAX,samples), date %in% c("t10"))),
                            ncol=2,nrow=3,common.legend = TRUE)


pdf(paste0(out_dir,"gt_phylo_ordination_weight.pdf"))
all_timepoints
dev.off()

# save data of corr with weight
w_t0 <- as.data.frame(weight_corr_output(subset_samples(phyloseq(gOTU,TAX,samples), date %in% c("t0"))))
w_t0$date="t0"
w_t2 <- as.data.frame(weight_corr_output(subset_samples(phyloseq(gOTU,TAX,samples), date %in% c("t2"))))
w_t2$date="t2"
w_t4 <- as.data.frame(weight_corr_output(subset_samples(phyloseq(gOTU,TAX,samples), date %in% c("t4"))))
w_t4$date="t4"
w_t6 <- as.data.frame(weight_corr_output(subset_samples(phyloseq(gOTU,TAX,samples), date %in% c("t6"))))
w_t6$date="t6"
w_t8 <- as.data.frame(weight_corr_output(subset_samples(phyloseq(gOTU,TAX,samples), date %in% c("t8"))))
w_t8$date="t8"
w_t10 <- as.data.frame(weight_corr_output(subset_samples(phyloseq(gOTU,TAX,samples), date %in% c("t10"))))
w_t10$date="t10"

weight_stats <- rbind(w_t0,w_t2,w_t4,w_t6,w_t8,w_t10)

# save 
fwrite(x=weight_stats, file = paste0(out_dir_git,"gtdb_weight.csv"), sep = ",")



######################################################################
######################################################################




# BREED

breed_corr_plot <- function(ccc) {
  
  # RAREFY
  ccc = rarefy_even_depth(ccc,
                          replace=TRUE, 
                          rngseed = 42)
  
  # Remove taxa not seen more than 2 times in at least a 20% of the samples 
  ccc_abund <- filter_taxa(ccc, function(x) sum(x > 2) > (0.2*length(x)), TRUE)
  
  ccc.ord <- ordinate(ccc_abund, "PCoA", "bray")
  
  the_breed_plot <- plot_ordination(ccc_abund, ccc.ord, type="samples", color="breed") + 
    geom_point(size=1)  +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x=element_text(size = 6),
          axis.text.y=element_text(size = 6),
          axis.title.x=element_text(size = 8),
          axis.title.y=element_text(size = 8),
          strip.text = element_text(size=8)) +
    ggtitle(paste0(sample_data(ccc_abund)$date)[1])+
    facet_wrap(~birth_day) +
    scale_color_brewer(palette = "Dark2")+
    guides(col = guide_legend(nrow = 4, title.position = "top"))
  return(the_breed_plot)
}


breed_corr_stats <- function(ccc) {
  
  # RAREFY
  ccc = rarefy_even_depth(ccc,
                          replace=TRUE, 
                          rngseed = 42)
  
  # Remove taxa not seen more than 2 times in at least a 20% of the samples 
  ccc_abund <- filter_taxa(ccc, function(x) sum(x > 2) > (0.2*length(x)), TRUE)
  
  ccc.ord <- ordinate(ccc_abund, "PCoA", "bray")
  
  myd <- plot_ordination(ccc_abund, ccc.ord, type="samples", color="breed") 
  
  x_axis1 <- dunnTest(myd$data$Axis.1~myd$data$breed,data=myd$data,method = "bonferroni")
  x_axis1 <- as.data.frame(x_axis1$res)
  x_axis1$axis = "Axis.1"
  
  x_axis2 <- dunnTest(myd$data$Axis.2~myd$data$breed,data=myd$data,method = "bonferroni")
  x_axis2 <- as.data.frame(x_axis2$res)
  x_axis2$axis = "Axis.2"
  
  x <- rbind(x_axis1,x_axis2) 
  x$group <- as.character(myd$data$date[1])
  
  x <- as.data.frame(x)
  
  return(x)
}

########

age_breed_plot_t0 <- breed_corr_plot(subset_samples(phyloseq(gOTU,TAX,samples), (date %in% c("t0"))))
age_breed_plot_t2 <- breed_corr_plot(subset_samples(phyloseq(gOTU,TAX,samples), (date %in% c("t2"))))
age_breed_plot_t4 <- breed_corr_plot(subset_samples(phyloseq(gOTU,TAX,samples), (date %in% c("t4"))))

all_breed_and_age_plot <- ggarrange(age_breed_plot_t0, 
                                    age_breed_plot_t2, 
                                    age_breed_plot_t4, 
                                    common.legend=TRUE,
                                    nrow=3)

########
breed_subset <- subset_samples(carbom, (breed %in% c("DxL","DxLW") & birth_day %in% c("2017-01-08")))
# plots
breed_plot_t0 <- breed_corr_plot(subset_samples(breed_subset, (date %in% c("t0"))))
breed_plot_t2 <- breed_corr_plot(subset_samples(breed_subset, (date %in% c("t2"))))
breed_plot_t4 <- breed_corr_plot(subset_samples(breed_subset, (date %in% c("t4"))))

breed_plot_all <- ggarrange(breed_plot_t0, 
                            breed_plot_t2,
                            breed_plot_t4,
                            common.legend=TRUE,
                            nrow=3)
########

# stats
breed_stats <- rbind(breed_corr_stats(subset_samples(breed_subset, (date %in% c("t0")))),
                     breed_corr_stats(subset_samples(breed_subset, (date %in% c("t2")))),
                     breed_corr_stats(subset_samples(breed_subset, (date %in% c("t4")))))

breed_stats$subgroup <- "bday_20170108"
breed_stats$test = "Dunn.test"
breed_stats$correction = "Bonferroni"

fwrite(x=breed_stats, file = paste0(out_dir_git,"gtdb_breed.csv"), sep = ",")


######################################################################



# AGE

age_corr_plot <- function(ccc) {
  
  # RAREFY
  ccc = rarefy_even_depth(ccc,
                          replace=TRUE, 
                          rngseed = 42)
  
  # Remove taxa not seen more than 2 times in at least a 20% of the samples 
  ccc_abund <- filter_taxa(ccc, function(x) sum(x > 2) > (0.2*length(x)), TRUE)
  
  ccc.ord <- ordinate(ccc_abund, "PCoA", "bray")
  
  the_age_plot <- plot_ordination(ccc_abund, ccc.ord, type="samples", color="birth_day") + 
    geom_point(size=1)  +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x=element_text(size = 6),
          axis.text.y=element_text(size = 6),
          axis.title.x=element_text(size = 8),
          axis.title.y=element_text(size = 8),
          strip.text = element_text(size=8)) +
    ggtitle(paste0(sample_data(ccc_abund)$date)[1])+
    facet_wrap(~breed) +
    scale_color_jco()+
    guides(col = guide_legend(nrow = 4, title.position = "top"))
  
  return(the_age_plot)
}

age_corr_stats <- function(ccc) {
  
  # RAREFY
  ccc = rarefy_even_depth(ccc,
                          replace=TRUE, 
                          rngseed = 42)
  
  # Remove taxa not seen more than 2 times in at least a 20% of the samples 
  ccc_abund <- filter_taxa(ccc, function(x) sum(x > 2) > (0.2*length(x)), TRUE)
  
  ccc.ord <- ordinate(ccc_abund, "PCoA", "bray")
  
  myd <- plot_ordination(ccc_abund, ccc.ord, type="samples", color="birth_day") 
  
  x_axis1 <- dunnTest(myd$data$Axis.1~myd$data$birth_day,data=myd$data,method = "bonferroni")
  x_axis1 <- as.data.frame(x_axis1$res)
  x_axis1$axis = "Axis.1"
  
  x_axis2 <- dunnTest(myd$data$Axis.2~myd$data$birth_day,data=myd$data,method = "bonferroni")
  x_axis2 <- as.data.frame(x_axis2$res)
  x_axis2$axis = "Axis.2"
  
  x <- rbind(x_axis1,x_axis2) 
  x$group <- as.character(myd$data$date[1])
  
  x <- as.data.frame(x)
  
  return(x)
}


age_subset <- subset_samples(carbom, (birth_day %in% c("2017-01-08","2017-01-11") & breed %in% c("DxL")))
# plots
age_plot_t0 <- age_corr_plot(subset_samples(age_subset, (date %in% c("t0"))))
age_plot_t2 <- age_corr_plot(subset_samples(age_subset, (date %in% c("t2"))))
age_plot_t4 <- age_corr_plot(subset_samples(age_subset, (date %in% c("t4"))))

age_plot_all <- ggarrange(age_plot_t0,
                          age_plot_t2,
                          age_plot_t4,
                          common.legend=TRUE,
                          nrow=3)

########

# stats
age_stats <- rbind(age_corr_stats(subset_samples(age_subset, (date %in% c("t0")))),
                   age_corr_stats(subset_samples(age_subset, (date %in% c("t2")))),
                   age_corr_stats(subset_samples(age_subset, (date %in% c("t4")))))

age_stats$subgroup <- "breed_DxL"
age_stats$test = "Dunn.test"
age_stats$correction = "Bonferroni"

fwrite(x=age_stats, file = paste0(out_dir_git,"gtdb_age.csv"), sep = ",")


########

# arrange all plots in one pdf

all_plotties <- plot_grid(
  all_breed_and_age_plot, 
  age_plot_all,
  breed_plot_all,
  ncol=3,nrow=1,
  rel_widths = c(5,2.5,2.5), labels = c("A","B","C")
)

pdf(paste0(out_dir,"gt_phylo_ordination_breed&age.pdf"))
all_plotties
dev.off()

#############################################

age_stats_sign <- age_stats %>%
  dplyr::select(P.adj,group_analyzed,group,Comparison) %>%
  dplyr::filter(P.adj<0.05) %>%
  dplyr::arrange(P.adj)
  
breed_stats_sign <- breed_stats %>%
  dplyr::select(P.adj,group_analyzed,group,Comparison) %>%
  dplyr::filter(P.adj<0.05) %>%
  dplyr::arrange(P.adj)

sink(paste0(out_dir,"breed_age_significant.txt"), type = "output", append = FALSE)
start_message <- " ########################## AGE - significant ########################## "
start_message
age_stats_sign
start_message <- " ########################## BREED - significant ########################## "
start_message
breed_stats_sign
sink()

############################################################################################################
############################################################################################################
############################################################################################################

# Presence of PROBIOTIC species in dataset - plotting the taxa that correspond most closely to what the prob formulations contain


carbom_ProbioticCheck <- phyloseq(gOTU,TAX,samples)
carbom_ProbioticCheck <- subset_samples(carbom_ProbioticCheck, (date %in% c("t0","t1","t2","t3","t4","t5", "t6","t7","t8","t9")))
carbom_ProbioticCheck <- subset_taxa(carbom_ProbioticCheck, (species %in% c("Lactobacillus_B salivarius",
                                                                            "Lactobacillus_F plantarum",
                                                                            "Enterococcus_B faecium_B",
                                                                            "Bifidobacterium bifidum",
                                                                            "Streptococcus thermophilus",
                                                                            "Lactobacillus_F plantarum", 
                                                                            "Lactobacillus helveticus",
                                                                            "Lactobacillus delbreuckii",
                                                                            "Lactobacillus_C rhamnosus")))
# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(carbom_ProbioticCheck))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_ProbioticCheck = transform_sample_counts(carbom_ProbioticCheck, standf)
sample_variables(carbom_ProbioticCheck)

# HEATMAP

sampleOrder = sort(sample_names(carbom_ProbioticCheck))
pdf(paste0(out_dir,"gt_phylo_heatmap_ProbioticCheck.pdf"))
# HEATMAP time - cohorts
plot_heatmap(carbom_ProbioticCheck, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "species", taxa.order = "species", sample.order = sampleOrder,
             trans=NULL, low="blue", high="red", na.value="blue") +
  facet_grid(~ cohort, switch = "x", scales = "free_x", space = "free_x")+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "Probiotic Species - Presence in cohorts")
dev.off()



######################################################################
######################################################################


# PROBIOTIC species check without phyloseq (so the lib size normalization is done before plotting)


# lib size normalisation

#################################
# STEP 1.

# normalization for library size 
df1 <- df0 %>%
  dplyr::group_by(pig,date) %>% # a pig sample is: pig+date
  dplyr::mutate(norm_value=value/sum(value))    
head(df1)
sort(unique(df1$species))


df2 <- df1 %>%
  dplyr::filter(!date=="tM") %>%
  dplyr::filter(!date=="t10") %>%
  group_by(cohort,pig,date,species) %>%
  dplyr::summarise(disp = mean(norm_value), sd = sd(norm_value), n=n()) %>%
  dplyr::filter(species=="Lactobacillus_B salivarius"|
           species=="Lactobacillus_F plantarum"|
           species=="Bifidobacterium bifidum"|
           species=="Enterococcus_B faecium_B"|
           species=="Lactobacillus delbrueckii"|
           species=="Lactobacillus helveticus"|
           species=="Lactobacillus_C rhamnosus"|
           species=="Lactobacillus_F plantarum"|
           species=="Streptococcus thermophilus")

# reorder cohorts 
df2$cohort  = factor(df2$cohort, levels=c("Control", 
                                          "D-Scour",
                                          "ColiGuard", 
                                          "Neomycin",
                                          "NeoD",
                                          "NeoC"))

multiple_DFs <- split( df2 , f = df2$species )



pdf(paste0(out_dir,"gt_ProbioticCheck.pdf"))
for (single_DF in multiple_DFs) {
  
  
  DF <- as.data.frame(single_DF)
  
  species <- DF$species[1]
  DF$species <- NULL
  
  print(ggplot(DF, aes(x=date,y=disp, color=cohort)) +
          #geom_bar(stat="identity",position=position_dodge())+
          geom_point(aes(x=date,y=disp, group=pig),size=2)+
          geom_line(aes(x=date,y=disp, group=pig))+
          theme_minimal()+
          stat_n_text(size = 2) +
          facet_wrap(~cohort, scales = c("fixed")) +
          ggtitle(paste0(species)))
  
  
}
dev.off()



############################################################################################################
############################################################################################################


