
# upload libraries
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

######################################################################

middle_dir = "/Users/12705859/metapigs_dry/middle_dir/" # git 
out_dir = "/Users/12705859/Desktop/metapigs_dry/checkm/"  # local


# input files: 
# checkm/checkm_all_nearly
# no_reps_all.csv (BINS COUNTS)

# OUTPUTS:

# cm_phylo_barplot.pdf
# cm_phylo_heatmap.pdf
# cm_phylo_ordination.pdf
# cm_phylo_diversity.pdf
# cm_phylo_network.pdf


######################################################################

# counts data 

# upload bins with counts (sample-dereplicated- output of 7.R, now 003_depths_to_counts.R)

no_reps_all <- read.csv(paste0(middle_dir,"no_reps_all.csv"), 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)

# remove .fa extension to match bins in checkm df 
no_reps_all$bin <- gsub(".fa","", no_reps_all$bin)


######################################################################


# checkM output of NEARLY COMPLETE bins 

# upload file
# careful cause we don't have the pigID to distinguish which bins to which sample
checkm_all_nearly <- read_delim(paste0(middle_dir,"checkm_all_nearly"), 
                                "\t", escape_double = FALSE, col_types = cols(pigid = col_character()), 
                                trim_ws = TRUE)


# some formatting 
checkm_all_nearly$Completeness <- as.numeric(checkm_all_nearly$Completeness)
checkm_all_nearly$Contamination <- as.numeric(checkm_all_nearly$Contamination)
# remove rows containing NAs as in the original file these rows where headers
checkm_all_nearly <- na.omit(checkm_all_nearly)


# filter >90 <5 bins 
checkm_all_nearly <- dplyr::filter(checkm_all_nearly, !grepl("Completeness",Completeness))
newdata <- subset(checkm_all_nearly, Completeness >= 90 & Contamination <= 5)

# rename cols to matching colnames between dataframes to merge  
colnames(newdata)[colnames(newdata)=="pigid"] <- "pig"
colnames(newdata)[colnames(newdata)=="Bin Id"] <- "bin"
colnames(newdata)[colnames(newdata)=="Taxonomy (contained)"] <- "taxa"

newdata <- newdata %>%
  dplyr::select(pig,bin,taxa)

######################################################################


# create "gOTU" table: a separate "gOTU" identifier for each different taxa

new <- newdata

# split taxa column into several (kingdom, phylum, etc ...) 
new$all_taxa <- new$taxa
new <- cSplit(new, "taxa", sep=";")
head(new)
colnames(new) <- c("pig","bin","all_taxa","kingdom","phylum","class","order","family","genus","species")
head(new)
NROW(new)

new <- new %>%
  dplyr::group_by(all_taxa) %>%
  dplyr::mutate(gOTU = paste0("gOTU_",group_indices())) 

taxa_mat <- new[4:11] 

NROW(taxa_mat)
NROW(unique(taxa_mat$gOTU))
taxa_mat <- taxa_mat %>%
  dplyr::distinct()
NROW(taxa_mat)
NROW(unique(taxa_mat$gOTU))

taxa_mat_df <- as.data.frame(taxa_mat)
# to matrix 
taxa_mat <- taxa_mat_df
rownames(taxa_mat) <- taxa_mat[,8]
taxa_mat[,8] <- NULL
taxa_mat <- as.matrix(taxa_mat)
# ready 

NROW(unique(rownames(taxa_mat)))
head(taxa_mat_df)

############################################################################################################

df <- right_join(no_reps_all, new, by=c("pig","bin"))
head(df)
NROW(df)

NROW(unique(paste0(df$pig,df$date)))
NROW(unique(paste0(df$gOTU)))

# columns to be kept 
keep <- c("cohort","pig","bin","date","value","gOTU")
df1 <- df[ , (names(df) %in% keep)]

# NA to zeros 
df1$value[is.na(df1$value)] <- 0

# as dates with NA was giving problems, change to class character and swap date NAs with "no-t"
df1$date <- as.character(df1$date)
df1$date[is.na(df1$date)] <- "no-t"

NROW(df1)
# sum up all the counts from the same sample (pig and date) that belong to the same gOTU
df2 <- df1 %>%
  group_by(pig,date,gOTU) %>%
  dplyr::summarise(all_bins_value = sum(value))

NROW(df2)
NROW(unique(paste0(df2$pig,df2$date)))


# assign a unique sample name 
df2$sample <- paste0(df2$pig,"_",df2$date)
# remove now pig and date (redundant)
df2 <- df2[3:5]

# long to wide 
df3 <- df2 %>%
  pivot_wider(names_from = sample, values_from = all_bins_value, values_fill = list(all_bins_value = 0)) 

# to matrix 
gOTU_mat <- as.data.frame(df3)
rownames(gOTU_mat) <- gOTU_mat[,1]
gOTU_mat[,1] <- NULL
gOTU_mat <- as.matrix(gOTU_mat)
mode(gOTU_mat) <- "integer"

# ready to transform to gOTU table

NROW(unique(rownames(gOTU_mat)))
NROW(unique(colnames(gOTU_mat)))

############################################################################################################


# create a sample table

df_dummy <- df1
df_dummy$sample_name <- paste0(df_dummy$pig,"_",df_dummy$date,"_")
df_dummy <- df_dummy %>%
  dplyr::select(sample_name, everything()) 
df_dummy <- df_dummy[,1:2]
head(df_dummy)

df_dummy$sample_name2 <- df_dummy$sample_name
df_dummy <- cSplit(df_dummy, "sample_name2", sep="_")

colnames(df_dummy) <- c("sample_name","cohort","pig","date")
NROW(df_dummy)

df_dummy$sample_name <- df_dummy$sample_name %<>%
  gsub('_$', '', .)  # removes the last _

df_dummy <- unique(df_dummy)
sample_df <- as.data.frame(df_dummy)

rownames(sample_df) <- sample_df[,1]

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
                                    "t10"))

# reorder cohorts 
sample_df$cohort  = factor(sample_df$cohort, levels=c("Control", 
                                        "D-Scour",
                                        "ColiGuard", 
                                        "Neomycin",
                                        "NeoD",
                                        "NeoC"))

############################################################################################################

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

cm_ordination_plot <- plot_ordination(carbom_abund, carbom_abund.ord, type="samples", color="date") + 
  geom_point(size=1) +
  theme_bw()+
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=7))+
  guides(colour = guide_legend(nrow = 1))+
  theme(legend.position="top") 

pdf(paste0(out_dir,"cm_phylo_ordination.pdf"))
cm_ordination_plot #+
  #facet_wrap(~cohort)
dev.off()

######################

# NETWORK ANALYSIS 

ig = make_network(carbom_abund, type = "samples", distance = "bray", max.dist = 0.35)
cm_network_plot <- plot_network(ig, carbom_abund, color = "date", shape = "cohort", line_weight = 0.3, 
                                  label = NULL, point_size = 1)+
  theme(legend.position = "bottom")+
  guides(shape = guide_legend(nrow = 1))+
  guides(size = "legend", colour = "none")

pdf(paste0(out_dir,"cm_phylo_network.pdf"))
cm_network_plot
dev.off()

######################


# BAR PLOT

# BAR GRAPH - by time point
pdf(paste0(out_dir,"cm_phylo_barplot_time.pdf"))
plot_bar(carbom_abund, fill = "phylum") +
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack") +
  facet_grid(~date,scales="free_x") +
  theme(axis.text.x = element_blank())
plot_bar(carbom_abund, fill = "class") +
  geom_bar(aes(color=class, fill=class), stat="identity", position="stack") +
  facet_grid(~date,scales="free_x") +
  theme(axis.text.x = element_blank())
plot_bar(carbom_abund, fill = "order") +
  geom_bar(aes(color=order, fill=order), stat="identity", position="stack") +
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
pdf(paste0(out_dir,"cm_phylo_heatmap.pdf"))
plot_heatmap(physeq1, 
             taxa.label = "species", 
             sample.order = "date") +
  facet_grid(~ date, switch = "x", scales = "free_x", space = "free_x")+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = paste0("Diversity of GTDB-predicted species in the piglet population"))
plot_heatmap(physeq1, 
             taxa.label = "genus", 
             sample.order = "date") +
  facet_grid(~ date, switch = "x", scales = "free_x", space = "free_x")+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = paste0("Diversity of GTDB-predicted species in the piglet population"))
dev.off()

######################
