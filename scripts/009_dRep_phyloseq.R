

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
#library(ComplexHeatmap)
library(cluster)
library(circlize)
library(readxl)
library(data.table)
library(FSA)
library(openxlsx)



source_dir = "/Users/12705859/metapigs_dry/source_data/" # git 
middle_dir = "/Users/12705859/metapigs_dry/middle_dir/" # git 
out_dir = "/Users/12705859/Desktop/metapigs_dry/dRep/"  # local
out_dir_git = "/Users/12705859/metapigs_dry/out/" # git 

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

# ###### run this part if you want to only retain piglets that 
# # where present throughput trial (not euthanised) 
# 
# to_keep <- no_reps_all %>%
#   filter(date=="t8") %>%
#   dplyr::select(pig) %>%
#   distinct()
# 
# to_keep <- as.character(to_keep$pig)
# 
# no_reps_all <- subset(no_reps_all, (pig %in% to_keep))
# NROW(unique(no_reps_all$pig))
# 
# ############################################################

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
unique(mdat$Cohort)

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

# merge bday info : 
bday <- breed_bday %>%
  dplyr::select(pig,birth_day)

weights <- left_join(weights, bday)
NROW(weights)
unique(weights$birth_day)

weights_rest <- weights %>% 
  dplyr::filter(!birth_day=="2017-01-06") %>%
  group_by(birth_day,date) %>%
  dplyr::mutate(weight_category=cut(weight, breaks=c(summary(weight)[1], summary(weight)[2], summary(weight)[5], summary(weight)[6]), 
                             labels=c("under","normal","over"))) 

weights_rest <- as.data.frame(weights_rest)

# quickly visualize weight category distribution
ggplot(weights_rest,aes(x=date,fill=weight_category)) +
  geom_bar() +
  facet_wrap(~birth_day)

weights6 <- weights %>% 
  dplyr::filter(birth_day=="2017-01-06")
weights6$weight_category <- NA

weights <- rbind(weights_rest,weights6) %>%
  dplyr::select(pig,date,weight_category)
head(weights)


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
  dplyr::select(secondary_cluster,species,genus,family,order,class,phylum,domain) %>%
  group_by(secondary_cluster) %>%
  slice(1) %>%
  dplyr::filter(!secondary_cluster=="no_cluster")

taxa_mat_df <- as.data.frame(taxa_mat)

NROW(taxa_mat_df)
NROW(unique(taxa_mat_df$secondary_cluster))

# to matrix 
taxa_mat <- taxa_mat_df
rownames(taxa_mat) <- taxa_mat[,1]
taxa_mat[,1] <- NULL
taxa_mat <- as.matrix(taxa_mat)

# ready 

NROW(unique(rownames(taxa_mat)))
head(taxa_mat_df)


######################################################################

# species  

# columns to be kept 
keep <- c("cohort","pig","bin","date","value","secondary_cluster")
df1 <- df0[ , (names(df0) %in% keep)]

# NA to zeros 
df1$value[is.na(df1$value)] <- 0

# as dates with NA were giving problems, change to class character and swap date NAs with "no-t"
df1$date <- as.character(df1$date)
df1$date[is.na(df1$date)] <- "no-t"

NROW(df1)

# sum up all the counts from the same sample (pig and date) that belong to the same OTU
df2 <- df1 %>%
  group_by(pig,date,secondary_cluster) %>%
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
  dplyr::select(sample,pig,date,cohort,pen,birth_day,breed,weight_category) %>%
  group_by(sample) %>%
  dplyr::slice(1)

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


# Prepare function for rarefaction: 

myrarefy_fun <- function(your_phyloseq_obj) {
  
  # removal of samples with low count (any sample with a count lower than 10k)
  r <- which(colSums(otu_table(your_phyloseq_obj))<10000) 
  to_remove <- rownames(as.data.frame(r))
  
  carbom_noFailSamples <- prune_samples(!(sample_names(your_phyloseq_obj) %in% to_remove), your_phyloseq_obj)
  
  # RAREFY
  carbom_rarefied = rarefy_even_depth(carbom_noFailSamples, 
                                      sample.size = min(sample_sums(carbom_noFailSamples)), 
                                      rngseed = 42)
  return(carbom_rarefied)
  
}

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


# # Remove taxa if not seen in at least 30% of the samples 
# carbom_abund <- filter_taxa(carbom, function(x) sum(x > 1) > (0.3*length(x)), TRUE)

carbom_abund <- carbom 

carbom_abund.ord <- ordinate(carbom_abund, "NMDS", "bray")

dRep_ordination_plot <- plot_ordination(carbom_abund, carbom_abund.ord, type="samples", color="date") + 
  geom_point(size=1) +
  theme_bw()+
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=7))+
  guides(colour = guide_legend(nrow = 1))+
  theme(legend.position="top")

pdf(paste0(out_dir,"dRep_phylo_ordination.pdf"))
dRep_ordination_plot #+
  #facet_wrap(~cohort)
dev.off()


########################################################################################

# NETWORK ANALYSIS 

ig = make_network(carbom_abund, type = "samples", distance = "bray", max.dist = 0.65) #0.55
dRep_network_plot <- plot_network(ig, carbom_abund, color = "date", shape = "cohort", line_weight = 0.3, 
                                label = NULL, point_size = 1)+
  theme(legend.position = "bottom")+
  guides(shape = guide_legend(nrow = 1))+
  guides(size = "legend", colour = "none")

pdf(paste0(out_dir,"dRep_phylo_network.pdf"))
dRep_network_plot
dev.off()

########################################################################################


# BAR PLOT

# BAR GRAPH - by time point
pdf(paste0(out_dir,"dRep_phylo_barplot_time.pdf"))
plot_bar(carbom_abund, fill = "phylum") +
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack") +
  facet_grid(~date,scales="free_x") +
  theme(axis.text.x = element_blank())
plot_bar(carbom_abund, fill = "class") +
  geom_bar(aes(color=class, fill=class), stat="identity", position="stack") +
  facet_grid(~date,scales="free_x") +
  theme(axis.text.x = element_blank())
dev.off()


######################

# HEATMAP

# NORMALIZATION BY RAREFACTION
carbom <- phyloseq(gOTU,TAX,samples)
# SUBSETTING phyloseq obejct
carbom <- subset_samples(carbom, (date %in% c("t0","t2","t4","t6","t8","t10")))


# cut out samples with extremely low counts and RAREFY: 
carbom_rarefied <- myrarefy_fun(carbom)


c <- carbom_rarefied

c

c1.2 <- filter_taxa(c, function(x) sum(x > 100) > (0.2*length(x)), TRUE)
c1.2


# out of these, take the ones with the highest inter-samples variance 
c2 = filter_taxa(c1.2, function(x) var(x) > 40000000, TRUE)
c2


carbom_abund <- c2

random_tree = rtree(ntaxa(carbom_abund), rooted=TRUE, tip.label=taxa_names(carbom_abund))
physeq1 = merge_phyloseq(carbom_abund,random_tree)


# HEATMAP time - genus, family, order, etc ...
pdf(paste0(out_dir,"dRep_phylo_heatmap.pdf"))
plot_heatmap(physeq1, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "species", taxa.order="species",
             sample.order = "date", trans=identity_trans(),
             low="blue", high="red", na.value="white") +
  ggtitle(label = "Microbe Species Diversity")
plot_heatmap(physeq1, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "genus", taxa.order="species",
             sample.order = "date", trans=identity_trans(),
             low="blue", high="red", na.value="white") +
  ggtitle(label = "Microbe Genus Diversity") 
plot_heatmap(physeq1, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "family", taxa.order="species",
             sample.order = "date", trans=identity_trans(),
             low="blue", high="red", na.value="white") +
  ggtitle(label = "Microbe Family Diversity") 
plot_heatmap(physeq1, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "order", taxa.order="species",
             sample.order = "date", trans=identity_trans(),
             low="blue", high="red", na.value="white") +
  ggtitle(label = "Microbe Order Diversity") 
plot_heatmap(physeq1, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "class", taxa.order="species",
             sample.order = "date", trans=identity_trans(),
             low="blue", high="red", na.value="white") +
  ggtitle(label = "Microbe Class Diversity") 
plot_heatmap(physeq1, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "phylum", taxa.order="species",
             sample.order = "date", trans=identity_trans(),
             low="blue", high="red", na.value="white") +
  ggtitle(label = "Microbe Phylum Diversity") 
dev.off()


# HEATMAP time - cohorts
pdf(paste0(out_dir,"dRep_phylo_heatmap_cohorts.pdf"))
plot_heatmap(physeq1, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "species", taxa.order="species",
             sample.order = "date", trans=identity_trans(),
             low="blue", high="red", na.value="white") +
  facet_grid(~ cohort, switch = "x", scales = "free_x", space = "free_x")+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "Microbe Species Diversity")
dev.off()


##############################################

# this is a small break in between the phyloseq analysis: 

# DISPLAYING MOST ABUNDANT SPECIES TIME-TREND

# this is a zoom in on the species shown in phylo_heatmap because most abundant 
# (taking gOTUs that represent at least 3% of the sample and present in at least 40 samples)

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
                                      "NeoC"))


# split df by species
multiple_DFs <- split( z , f = z$species )

pdf(paste0(out_dir,"dRep_zoomIN_on_phylo_heatmap.pdf"))
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


##############################################


# DIVERSITY 

# NORMALIZATION BY RAREFACTION
carbom <- phyloseq(gOTU,TAX,samples)
# SUBSETTING phyloseq obejct
carbom <- subset_samples(carbom, (date %in% c("t0","t1","t2","t3","t4","t5", "t6","t7","t8","t9","t10")))

# cut out samples with extremely low counts and RAREFY: 
carbom_rarefied <- myrarefy_fun(carbom)

c <- carbom_rarefied

dRep_diversity_samples <- plot_richness(c, measures=c("Chao1","Shannon", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), 
                                      color="date", x="date") +
  guides(colour = guide_legend(nrow = 1))+
  theme(legend.position="top")

######### plotting above results in a different way: 
# focus on three measures;
# whisker plots instead 

comparison <- data.frame(comparison=c("t2-t0", "t4-t2", "t6-t4", "t8-t6", "t10-t8"))

Chao1_data <- dRep_diversity_samples$data %>%
  filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="t10") %>%
  filter(variable=="Chao1")
Shannon_data <- dRep_diversity_samples$data %>%
  dplyr::filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="t10") %>%
  dplyr::filter(variable=="Shannon")
Simpson_data <- dRep_diversity_samples$data %>%
  dplyr::filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="t10") %>%
  dplyr::filter(variable=="Simpson")



estimates_CI_ttest_fun <- function(xxx) {
  
  # estimates and CIs: 
  apc <- pairwiseCI(value ~ date, data=xxx,
                    method="Param.diff")
  s <- summary(apc)
  s_estimate <- round(s$estimate,2)
  s_CI <- as.data.frame(s$conf.int)
  s_CI$comparison <- rownames(s_CI)
  s <- cbind(s_CI,s_estimate)
  s$se <- round(s$estimate-s$lower,2)
  
  # t-test with Bonferroni adjust: 
  apcTest <- pairwiseTest(value ~ date, data=xxx,
                          method="t.test")
  t <- summary(apcTest, p.adjust.method = "bonferroni")
  
  st <- inner_join(s,t)
  st$variable <- unique(xxx$variable)
  colnames(st) <- c("lower","upper","comparison","estimate","se","p.val.adj","p.val.raw",
                    "group1","group2","variable")
  return(st)
}

Chao1_data_stats <- estimates_CI_ttest_fun(Chao1_data)
Shannon_data_stats <- estimates_CI_ttest_fun(Shannon_data)
Simpson_data_stats <- estimates_CI_ttest_fun(Simpson_data)


# plotting: 

Chao1_data_stats <- Chao1_data_stats %>%
  inner_join(., comparison, by="comparison") %>%
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
  stat_pvalue_manual(Chao1_data_stats, label = "p.val.adj")


Shannon_data_stats <- Shannon_data_stats %>%
  inner_join(., comparison, by="comparison") %>%
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
  stat_pvalue_manual(Shannon_data_stats, label = "p.val.adj")


Simpson_data_stats <- Simpson_data_stats %>%
  inner_join(., comparison, by="comparison") %>%
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
  stat_pvalue_manual(Simpson_data_stats, label = "p.val.adj")


tosave <- ggarrange(Chao1_plot, 
                    Shannon_plot,
                    Simpson_plot,
                    ncol = 3, 
                    nrow=1,
                    labels=c("A","B","C"), 
                    common.legend = TRUE)

ggsave(filename = paste0(out_dir,"dRep_phylo_diversity_boxplot.pdf"), plot = tosave)


# save stats : 

dRep_diversity <- rbind(Shannon_data_stats,
                        Simpson_data_stats,
                        Chao1_data_stats)

# save 
fwrite(x=dRep_diversity, file = paste0(out_dir_git,"dRep_diversity.csv"), sep = ",")

############################################################################################################
############################################################################################################

