
library(readr)
library(splitstackshape)
library(treemap)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(robCompositions)
library(ggbiplot)


source_dir = "/Users/12705859/metapigs_dry/source_data/" # git 
middle_dir = "/Users/12705859/metapigs_dry/middle_dir/" # git 
out_dir = "/Users/12705859/Desktop/metapigs_dry/gtdbtk/"  # local


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

# taxa overall - based on gtdbtk assignment


phylum_counts <- setDT(gtdbtk_bins)[, .(Freq = .N), by = .(phylum)]

# most abundant 
phylum_counts_most_ab <- phylum_counts %>%
  dplyr::filter(!Freq<2) %>% # two bins in the whole dataset matching a phylum are not worth taking along
  dplyr::mutate(perc=round(Freq/sum(Freq)*100,2)) %>%
  # if percentage lower than 1.4 group as "others"
  dplyr::mutate(phylum = ifelse(perc < 0.7, "others", phylum)) %>%
  group_by(phylum) %>%
  dplyr::summarise(perc=sum(perc))
  
# least abundant
phylum_counts_least_ab <- phylum_counts %>%
  dplyr::filter(!Freq<2) %>% # two bins in the whole dataset matching a phylum are not worth taking along
  dplyr::mutate(perc=round(Freq/sum(Freq)*100,2)) %>%
  dplyr::filter(perc<0.7)

phylum_counts_most_ab$label <- paste(paste(phylum_counts_most_ab$phylum,
                                           phylum_counts_most_ab$perc,sep = "\n"),"%")
phylum_counts_least_ab$label <- paste(paste(phylum_counts_least_ab$phylum,
                                            phylum_counts_least_ab$perc,sep = "\n"),"%")


pdf(paste0(out_dir,"gt_treemap_phyla.pdf"))
# most abundant 
treemap(phylum_counts_most_ab, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "perc",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Phyla distribution from all MAGs (GTDB) - most common", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
# least abundant 
treemap(phylum_counts_least_ab, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "perc",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Phyla distribution from all MAGs (GTDB) - least common", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
dev.off()




######################################################################



# Species overall - based on gtdbtk assignment



species_counts <- setDT(gtdbtk_bins)[, .(Freq = .N), by = .(family,species)]

# all 
species_counts <- species_counts %>%
  #filter(!Freq<10) %>% # two bins in the whole dataset matching a phylum are not worth taking along
  dplyr::mutate(perc=round(Freq/sum(Freq)*100,2)) %>%
  dplyr::arrange(desc(perc))
species_counts$label <- paste(paste(species_counts$species,
                                    species_counts$perc),"%")

# most abundant 
species_counts_most_ab <- species_counts[1:50]


pdf(paste0(out_dir,"gt_treemap_species.pdf"))
# most abundant 
treemap(species_counts_most_ab, #Your data frame object
        index=c("family","label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Species distribution from all MAGs (GTDB) - 50 most common", #Customize your title
        fontsize.title = 15, #Change the font size of the title
        overlap.labels=0,
        fontsize.labels=c(13,10),                # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
        fontcolor.labels=c("white","black"),    # Color of labels
        fontface.labels=c(2,1),                  # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
        bg.labels=c("transparent"),              # Background color of labels
        align.labels=list(
          c("left", "top"), 
          c("right", "bottom")
        )
        #fontsize.labels = 8
)
dev.off()


######################################################################



# Genus overall - based on gtdbtk assignment



genus_counts <- setDT(gtdbtk_bins)[, .(Freq = .N), by = .(order,genus)]

# all 
genus_counts <- genus_counts %>%
  #filter(!Freq<10) %>% # two bins in the whole dataset matching a phylum are not worth taking along
  dplyr::mutate(perc=round(Freq/sum(Freq)*100,2)) %>%
  dplyr::arrange(desc(perc))
genus_counts$label <- paste(paste(genus_counts$genus,
                                  genus_counts$perc),"%")

# most abundant 
genus_counts_most_ab <- genus_counts[1:50]


pdf(paste0(out_dir,"gt_treemap_genus.pdf"))
# most abundant 
treemap(genus_counts_most_ab, #Your data frame object
        index=c("order","label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Genus distribution from all MAGs (GTDB) - 50 most common", #Customize your title
        fontsize.title = 15, #Change the font size of the title
        overlap.labels=0,
        fontsize.labels=c(13,10),                # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
        fontcolor.labels=c("white","black"),    # Color of labels
        fontface.labels=c(2,1),                  # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
        bg.labels=c("transparent"),              # Background color of labels
        align.labels=list(
          c("left", "top"), 
          c("right", "bottom")
        )
        #fontsize.labels = 8
)
dev.off()



######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################



# Principal component analysis: clustering of bins based on phyla with time (labels per cohort)




# merge info 

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

##############################################

# how many species per sample do we have ? 

head(df0)

z <- df0 %>%
  group_by(pig,date,species) %>%
  distinct() %>%
  group_by(pig,date) %>%
  tally()

sink(file = paste0(out_dir,"gt_species_per_sample.txt"), 
     append = FALSE, type = c("output"))
paste0("Mean - number of gOTUs (GTDB species) per sample")
mean(z$n)  
paste0("Median - number of gOTUs (GTDB species) per sample")
median(z$n)
sink()

##############################################

df <- df0 %>%
  filter(!cohort=="Mothers")

#################################

# CREATE COUNTS TABLE 

df1 <- df
head(df1)
NROW(df1)

# PROCEED to all: 

# for each sample (pig,date), sum up the counts that fall within one species (same species assigned to distinct bins)
df2 <- df1 %>%
  group_by(pig,gOTU,date) %>%
  dplyr::summarize(sum_value = sum(value)) 
head(df2)
sum(df2$sum_value)

# normalize by library size 
df3 <- df2 %>% 
  group_by(pig,date) %>% 
  dplyr::mutate(norm_value = sum_value/sum(sum_value)) %>% 
  dplyr::select(-sum_value)
head(df3)

# if your total sum is equal to the total number of samples, 
# it means that the sum within each sample (pig,date) is 1, and that's correct  
NROW(unique(paste0(df3$pig,df3$date)))==sum(df3$norm_value)


df3 <- as.data.frame(df3)
df3$sample = paste0(df3$date,"_",df3$pig)
head(df3)

# pivot wider
df3 <- df3 %>%
  dplyr::select(sample,gOTU,norm_value) %>%
  pivot_wider(names_from = gOTU, values_from = norm_value, values_fill = list(norm_value = 0))

feat <- as.data.frame(df3)
which(is.na(feat[,1]))

# rownames(feat) <- feat[,1]
# feat[,1] <- NULL
# 
# head(feat)
# dim(feat)
# 
# # is the sum of each columns 1? 
# colSums(feat)
# # yes 

# ready! 




# get a quick cohorts to pig table
cohorts <- df %>% dplyr::select(cohort,pig,date) %>% distinct()
cohorts$sample <- paste0(cohorts$date,"_",cohorts$pig)
cohorts <- as.data.frame(cohorts)


df5 <- inner_join(cohorts,feat) 
df5$sample <- paste0(df5$date,"_",df5$cohort)

df5$pig <- NULL
df5$date <- NULL
df5$cohort <- NULL



df6 <- df5 %>%
  group_by(sample) %>%
  dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)


df6 <- as.data.frame(df6)
rowSums(df6[,-1])



rownames(df6) <- df6$sample
df6$sample <- NULL
m <- as.matrix(df6)

df6.pca <- prcomp(m, center = FALSE,scale. = FALSE)
summary(df6.pca)

# to get samples info showing on PCA plot
this_mat_samples <- data.frame(sample=rownames(m)) 
this_mat_samples <- cSplit(indt = this_mat_samples, "sample", sep = "_", drop = NA)

# reorder dates 
this_mat_samples$sample_1  = factor(this_mat_samples$sample_1, levels=c("t0",
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

gt_PC12 <- ggbiplot(df6.pca,
                    labels=this_mat_samples$sample_2,
                    groups=this_mat_samples$sample_1,
                    ellipse=TRUE,
                    var.axes = FALSE,
                    labels.size = 2,
                    choices = (1:2)) +
  theme_bw() +
  xlim(c(-2,1)) +
  scale_colour_discrete(name="timepoint")+
  guides(color = guide_legend(ncol = 1))

gt_PC34 <- ggbiplot(df6.pca,
                    labels=this_mat_samples$sample_2,
                    groups=this_mat_samples$sample_1,
                    ellipse=TRUE,
                    var.axes = FALSE,
                    labels.size = 2,
                    choices = (3:4)) +
  theme_bw() +
  scale_colour_discrete(name="timepoint")+
  guides(color = guide_legend(ncol = 1))


gt_PCA <- ggarrange(gt_PC12,gt_PC34,
                    ncol=2,legend = "right",
                    common.legend=TRUE)

pdf(paste0(out_dir,"gt_PCA.pdf"), width=7,height=4)
gt_PCA
dev.off()
