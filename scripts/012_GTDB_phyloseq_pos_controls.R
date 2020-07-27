
library(readr)
library(splitstackshape)
library(tidyr)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(viridis)


source_dir = "/Users/12705859/metapigs_dry/source_data/" # git 
middle_dir = "/Users/12705859/metapigs_dry/middle_dir/" # git 
out_dir_git = "/Users/12705859/metapigs_dry/out/" # git 
out_dir = "/Users/12705859/Desktop/metapigs_dry/gtdbtk/"  # local


######################################################################


# input files: 
# gtdbtk_bins_completeTaxa
# no_reps_pos_controls.csv (BINS COUNTS)


# OUTPUTS:
# gt_phylo_PosControls_barplot.pdf
# gt_phylo_PosControls_heatmap.pdf



######################################################################

# counts data 

no_reps_pos_controls <- read.csv(paste0(middle_dir,"no_reps_pos_controls.csv"), 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)


# remove .fa extension to match bins in checkm df 
no_reps_pos_controls$bin <- gsub(".fa","", no_reps_pos_controls$bin)
head(no_reps_pos_controls)
NROW(no_reps_pos_controls)

colnames(no_reps_pos_controls)[colnames(no_reps_pos_controls) == 'isolation_source'] <- 'pig'

######################################################################

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv(paste0(middle_dir,"gtdb_bins_completeTaxa"),
                        col_types = cols(node = col_character(),
                                         pig = col_character()))


head(gtdbtk_bins)

######################################################################


# merge info 

NROW(gtdbtk_bins)
NROW(no_reps_pos_controls)
head(gtdbtk_bins)
head(no_reps_pos_controls)
df <- merge(no_reps_pos_controls, gtdbtk_bins, by=c("pig","bin"))

# rename node as gOTU and place "gOTU_" in front of node number: a separate genomic OTU identifier for each different genome

colnames(df)[colnames(df) == 'node'] <- 'gOTU'
df$gOTU <- paste0("gOTU_",df$gOTU)

NROW(unique(df$gOTU))
NROW(df)

# change Protexin with D-Scour
df$pig <- gsub("Protexin","D-Scour", df$pig)


######################################################################


# TAXA


taxa_mat <- df %>%
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

# gOTU  

# columns to be kept 
keep <- c("cohort","pig","bin","variable","value","gOTU")
df1 <- df[ , (names(df) %in% keep)]

# NA to zeros 
df1$value[is.na(df1$value)] <- 0

# as dates with NA was giving problems, change to class character and swap date NAs with "no-t"
df1$date <- as.character(df1$variable)

NROW(df1)
# sum up all the counts from the same sample (pig and date) that belong to the same OTU
df2 <- df1 %>%
  group_by(pig,variable,gOTU) %>%
  dplyr::summarise(all_bins_value = sum(value))

NROW(df2)
NROW(unique(paste0(df2$pig,df2$variable)))


# assign a unique sample name 
df2$sample <- paste0(df2$pig,"_",df2$variable)
# remove now pig and date (redundant)
df2$pig <- NULL
df2$variable <- NULL

# long to wide 
df3 <- df2 %>%
  pivot_wider(names_from = sample, values_from = all_bins_value, values_fill = list(all_bins_value = 0)) 

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

sample_df <- df

sample_df$sample <- paste0(sample_df$pig,"_",sample_df$variable)
NROW(unique(sample_df$sample))

sample_df <- sample_df %>%
  dplyr::select(sample,pig,variable,cohort) %>%
  group_by(sample) %>%
  dplyr::slice(1)

sample_df$gOTU <- NULL
sample_df <- as.data.frame(sample_df)

NROW(sample_df)
head(sample_df)

rownames(sample_df) <- sample_df[,1]
# ready


######################################################################


# create phyloseq object

OTU = otu_table(gOTU_mat, taxa_are_rows = TRUE)
TAX = tax_table(taxa_mat)
samples = sample_data(sample_df)

carbom <- phyloseq(OTU,TAX,samples)


############################################################################################################

# NORMALIZATION 

# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)

############################################################################################################

# PLOT

######################

# HEATMAP
sampleOrder = sort(sample_names(carbom))

pos_controls_heatmap <- plot_heatmap(carbom, method = "MDS", distance="unifrac",weighted=TRUE, 
             taxa.label = "species", taxa.order = "species", sample.order = sampleOrder,
             trans = NULL, low="blue", high="red", na.value="blue") +
  theme(axis.text.x = element_text(size=9),
        axis.title.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.text = element_text(size=10, vjust=0),
        legend.title = element_text(size=11, vjust=0),
        legend.position = "right")+
  facet_grid(~pig,scales="free")+
  labs(x="Replicate",
       y="Species", 
       fill = "Relative\nabundance")

# HEATMAP time - genus, family, order, etc ...
pdf(paste0(out_dir,"gt_phylo_PosControls_heatmap.pdf"))
pos_controls_heatmap
dev.off()

######################


# BAR PLOT

pos_controls_barplot <- plot_bar(carbom, fill = "species") +
  facet_wrap(~pig,scales="free") +
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=7, angle=90, vjust=1),
        legend.text = element_text(size=9),
        legend.position="right")+
  labs(x="Replicate",
       y="Relative abundance")+
  guides(fill=guide_legend(ncol=1,bycol=TRUE))

# BAR GRAPH - by time point
pdf(paste0(out_dir,"gt_phylo_PosControls_barplot.pdf"))
pos_controls_barplot
dev.off()


######################################################################
######################################################################


# function to get barplots : one for each positive control 
# (played around with phyloseq subset data but couldn t get levels to drop!)

myfun <- function(single_pos_control_df,palette_type) {
  
  df <- single_pos_control_df
  mypalette <- as.character(palette_type)
  
  # TAXA
  taxa_mat <- df %>%
    dplyr::select(gOTU,species,genus,family,order,class,phylum,domain) %>%
    group_by(gOTU) %>%
    dplyr::slice(1)
  
  taxa_mat_df <- as.data.frame(taxa_mat)
  taxa_mat <- taxa_mat_df
  rownames(taxa_mat) <- taxa_mat[,1]
  taxa_mat[,1] <- NULL
  taxa_mat <- as.matrix(taxa_mat)
  
  # gOTU  
  keep <- c("cohort","pig","bin","variable","value","gOTU")
  df1 <- df[ , (names(df) %in% keep)]
  df1$value[is.na(df1$value)] <- 0
  df1$date <- as.character(df1$variable)
  df2 <- df1 %>%
    group_by(pig,variable,gOTU) %>%
    dplyr::summarise(all_bins_value = sum(value))
  df2$sample <- paste0(df2$pig,"_",df2$variable)
  df2$pig <- NULL
  df2$variable <- NULL
  df3 <- df2 %>%
    pivot_wider(names_from = sample, values_from = all_bins_value, values_fill = list(all_bins_value = 0)) 
  gOTU_mat <- as.data.frame(df3)
  rownames(gOTU_mat) <- gOTU_mat[,1]
  gOTU_mat[,1] <- NULL
  gOTU_mat <- as.matrix(gOTU_mat)
  # SAMPLES 
  sample_df <- df
  sample_df$sample <- paste0(sample_df$pig,"_",sample_df$variable)
  sample_df <- sample_df %>%
    dplyr::select(sample,pig,variable,cohort) %>%
    group_by(sample) %>%
    dplyr::slice(1)
  sample_df$gOTU <- NULL
  sample_df <- as.data.frame(sample_df)
  rownames(sample_df) <- sample_df[,1]
  
  # create phyloseq object
  
  OTU = otu_table(gOTU_mat, taxa_are_rows = TRUE)
  TAX = tax_table(taxa_mat)
  samples = sample_data(sample_df)
  
  carbom <- phyloseq(OTU,TAX,samples)
  
  ##########################################################################
  
  # NORMALIZATION 
  
  # Normalize number of reads in each sample using median sequencing depth.
  total = median(sample_sums(carbom))
  standf = function(x, t=total) round(t * (x / sum(x)))
  carbom = transform_sample_counts(carbom, standf)
  sample_variables(carbom)
  
  ##########################################################################
  
  # PLOT
  
  ######################
  
  # BAR PLOT
  
  barplot <- plot_bar(carbom, fill = "species") +
    facet_wrap(~pig,scales="free") +
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=7, angle=90, vjust=1),
          legend.text = element_text(size=7),
          legend.title=element_blank(),
          legend.position="top")+
    labs(x="Replicate",
         y="Relative abundance")+
    scale_fill_brewer(palette = as.character(mypalette))+
    guides(fill=guide_legend(ncol=1,bycol=TRUE))
  
  return(barplot)
  
}

df_sub1 <- df %>% dplyr::filter(pig=="MockCommunity")
df_sub2 <- df %>% dplyr::filter(pig=="D-Scour")
df_sub3 <- df %>% dplyr::filter(pig=="ColiGuard")


plot1 <- myfun(df_sub1,as.character("Paired"))  
plot2 <- myfun(df_sub2,as.character("Set1"))  
plot3 <- myfun(df_sub3,as.character("Set2"))  

p123 <- ggarrange(plot1,
          plot2,
          plot3,
          ncol=3)

pdf(paste0(out_dir,"gt_phylo_PosControls_barplotS.pdf"))
p123
dev.off()

######################################################################



