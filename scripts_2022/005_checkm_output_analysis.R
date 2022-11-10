

# upload libraries
library(tidyverse)
library(gplots)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(splitstackshape)
library(treemap)
library(data.table)
library(readxl)
library(pheatmap)
library(compositions) # <- this is the good one for clr transform
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
library(EnvStats)
library(treemapify)
library(scales)
library(compositions)

middle_dir = "/Users/dgaio/cloudstor/Gaio/github/metapigs_dry/middle_dir/" # git 
out_dir = "/Users/dgaio/cloudstor/Desktop/metapigs_dry/checkm/"  # local

# input files (in middle_dir)
# all_checkm_output.tsv 
# checkm_all_nearly
# all_bin_depths.csv (new bin depths - old file name: no_reps_all.csv)

# OUTPUTS:

# about bins in general :
#  - numbers_2022.txt

# from all_checkm_output :
#  - cm_Compl_vs_contam_2022.pdf
#  - cm_numbers_2022.txt
#  - cm_contigs_distribution_2022.pdf
#  - cm_scaffolds_predicted_genes_2022.pdf

# from BINS COUNTS + checkm_all_nearly: 
# based on bins frequency: 
## - cm_treemap_phylum_2022.pdf
## - cm_treemap_class_2022.pdf
## - cm_treemap_order_2022.pdf
## - numbers_2022.txt
# based on (counts) log10 relative abundance: 
## - cm_parallel_coordinates_phyla_2022.pdf
# based on (counts) relative abundance: 
## - cm_rel_ab_phyla_2022.pdf
## - cm_rel_ab_phyla_cohorts_2022.pdf
## - cm_balloonplot_phyla_2022.pdf
## - cm_balloonplot_phyla_cohorts_2022.pdf
# based on (counts) cenLR relative abundance: 
## - cm_PCA_2022.pdf



######################################################################

# template to collect and store bins info 

sink(paste0(out_dir,"numbers_2022.txt"))
start_message <- " ########################## BINS NUMBERS ########################## "
start_message
sink()

# template to collect and store checkM info 

sink(paste0(out_dir,"cm_numbers_2022.txt"))
start_message <- " ########################## CHECKM ANALYSIS ########################## "
start_message
sink()

######################################################################

# checkM output of ALL bins 

# upload file
# careful cause we don't have the pigID to distinguish which bins to which sample
all_checkm_output <- read_delim(paste0(middle_dir,"all_checkm_output.tsv"),
                                "\t", escape_double = FALSE, trim_ws = TRUE)
all_checkm_output <- dplyr::filter(all_checkm_output, !grepl("Completeness",Completeness))
all_checkm_output$Completeness <- as.numeric(all_checkm_output$Completeness)
all_checkm_output$Contamination <- as.numeric(all_checkm_output$Contamination)
NROW(all_checkm_output)


######################################################################


# checkM output of NEARLY COMPLETE bins 

# upload file
# careful cause we don't have the pigID to distinguish which bins to which sample
checkm_all_nearly <- read_delim(paste0(middle_dir,"checkm_all_nearly"), 
                                "\t", escape_double = FALSE, col_types = cols(pigid = col_character()), 
                                trim_ws = TRUE)


######################################################################


# upload bins with counts (sample-dereplicated- output of 7.R)

all_bins_depths <- read.csv(paste0(middle_dir,"all_bin_depths.csv"), 
                        na.strings=c("","NA"),
                        check.names = FALSE, 
                        header = TRUE)

# remove .fa extension to match bins in checkm df 
all_bins_depths$bin <- gsub(".fa","", all_bins_depths$bin)
head(all_bins_depths)


######################################################################

# merge cohort info 


#input files
cohorts <- read_excel("/Users/dgaio/cloudstor/Gaio/github/metapigs_dry/source_data/cohorts.xlsx")
colnames(cohorts)[colnames(cohorts)=="Animal ID"] <- "pig"
colnames(cohorts)[colnames(cohorts)=="Cohort Name"] <- "cohort"
cohorts <- as.data.frame(cohorts)

head(all_bins_depths)
all_bins_depths <- merge(cohorts,all_bins_depths, by="pig")


######################################################################
######################################################################

# Nearly complete bins (>=90 <=5) taxonomic categorization: proportions explained at each taxonomic level: 

# -- proportion of >=90 and <=5 bins assigned to phylum/order/etc... : 

# how many of >90 <5 bins are classified at the each level? 
nn <- checkm_all_nearly %>% 
  dplyr::rename(taxa = `Taxonomy (contained)`)
NROW(nn)

archaea <- dplyr::filter(nn, grepl('k__Archaea', taxa))
bacteria <- dplyr::filter(nn, grepl('k__Bacteria', taxa))
kingdom <- dplyr::filter(nn, grepl('k__', taxa))
phylum <- dplyr::filter(nn, grepl('p__', taxa))
class <- dplyr::filter(nn, grepl('c__', taxa))
order <- dplyr::filter(nn, grepl('o__', taxa))
family <- dplyr::filter(nn, grepl('f__', taxa))
genus <- dplyr::filter(nn, grepl('g__', taxa))
species <- dplyr::filter(nn, grepl('s__', taxa))

######################################################################



# MERGE checkM info of nearly complete bins to bins counts : 

# rename cols of checkm nearly complete bins to match colnames of no_reps_all for dataframes to merge  
colnames(nn)[colnames(nn)=="pigid"] <- "pig"
colnames(nn)[colnames(nn)=="Bin Id"] <- "bin"
colnames(nn)[colnames(nn)=="Taxonomy (contained)"] <- "taxa"

nn <- nn %>%
  dplyr::select(pig,bin,taxa)

# merge 
df <- right_join(all_bins_depths, nn, by=c("pig","bin"))
head(df)
NROW(df)


# split taxa column into several (kingdom, phylum, etc ...) 
df <- cSplit(df, "taxa", sep=";")

# reorder dates 
df$date  = factor(df$date, levels=c("tM", 
                                    "t0",
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
                                    "tNONE"))

# reorder cohorts 
unique(df$date)
df$cohort  = factor(df$cohort, levels=c("Control", 
                                        "D-Scour",
                                        "ColiGuard", 
                                        "Neomycin",
                                        "NeoD",
                                        "NeoC"))


######################################################################


# treemap PHYLA (piglet samples, mothers exlcuded)

# keep only rows that contain phyla info, discard others; also exclude mothers
for_treemap <- df %>%
  dplyr::select(cohort,pig,bin,taxa_2) %>%
  dplyr::filter(!cohort=="Mothers")

# remove "p__" before phylum 
for_treemap[,4] <- lapply(
  for_treemap[,4], 
  gsub, 
  pattern = "p__", 
  replacement = "", 
  fixed = TRUE)

counts <- setDT(for_treemap)[, .(Freq = .N), by = .(taxa_2)]

counts <- counts %>%
  dplyr::mutate(perc=round(Freq/sum(Freq)*100,2))

counts$label <- paste(paste(counts$taxa_2,counts$perc,sep = "\n"),"%")

pdf(paste0(out_dir,"cm_treemap_phyla_2022.pdf"))
treemap(counts, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Phyla distribution - from nearly complete MAGs - CheckM", #Customize your title
        fontsize.title = 15 #Change the font size of the title
)
dev.off()

# treemap CLASS

# keep only rows that contain phyla info, discard others; also exclude mothers
for_treemap <- df %>%
  dplyr::select(cohort,pig,bin,taxa_3) %>%
  dplyr::filter(!cohort=="Mothers")

# remove "p__" before phylum 
for_treemap[,4] <- lapply(
  for_treemap[,4], 
  gsub, 
  pattern = "c__", 
  replacement = "", 
  fixed = TRUE)

counts <- setDT(for_treemap)[, .(Freq = .N), by = .(taxa_3)]

counts <- counts %>%
  dplyr::mutate(perc=round(Freq/sum(Freq)*100,2))

counts$label <- paste(paste(counts$taxa_3,counts$perc,sep = "\n"),"%")

pdf(paste0(out_dir,"cm_treemap_class_2022.pdf"))
treemap(counts, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Class distribution from nearly complete MAGs - CheckM", #Customize your title
        fontsize.title = 15 #Change the font size of the title
)
dev.off()

# treemap ORDER

# keep only rows that contain phyla info, discard others; also exclude mothers
for_treemap <- df %>%
  dplyr::select(cohort,pig,bin,taxa_4) %>%
  dplyr::filter(!cohort=="Mothers")

# remove "p__" before phylum 
for_treemap[,4] <- lapply(
  for_treemap[,4], 
  gsub, 
  pattern = "o__", 
  replacement = "", 
  fixed = TRUE)

counts <- setDT(for_treemap)[, .(Freq = .N), by = .(taxa_4)]

counts <- counts %>%
  dplyr::mutate(perc=round(Freq/sum(Freq)*100,2))

counts$label <- paste(paste(counts$taxa_4,counts$perc,sep = "\n"),"%")

pdf(paste0(out_dir,"cm_treemap_order_2022.pdf"))
treemap(counts, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Order distribution from nearly complete MAGs - CheckM", #Customize your title
        fontsize.title = 15 #Change the font size of the title
)
dev.off()




######################################################################################################

# parallel coordinates & relative abundance - PHYLA

# keep only rows that contain phyla info, discard others; also exclude mothers
for_parallel_coo <- df %>%
  dplyr::select(cohort,pig,bin,date,taxa_2,value) %>%
  dplyr::filter(!cohort=="Mothers")

# remove "p__" before phylum 
for_parallel_coo[,5] <- lapply(
  for_parallel_coo[,5], 
  gsub, 
  pattern = "p__", 
  replacement = "", 
  fixed = TRUE)

# general time change 
summs_for_parallel_coo <- for_parallel_coo %>% group_by(taxa_2,date) %>% 
  dplyr::summarise(min = min(value)
                   ,max = max(value)
                   ,mean = mean(value)
                   ,n = n()
                   ,sd = sd(value)
                   ,q25 = quantile(value, .25)
                   ,q75 = quantile(value, .75)) 



pdf(paste0(out_dir,"cm_rel_ab_phyla.pdf"))
ggplot(summs_for_parallel_coo, aes(fill=taxa_2, y=mean, x=date)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.title.y = element_text(),
        legend.title = element_text()) +
  labs(x = "collection date",
       y = "relative abundance (ratio)",
       fill = "Phylum")  
dev.off()


summs_for_parallel_coo_cohorts <- for_parallel_coo %>% group_by(taxa_2,date,cohort) %>% 
  dplyr::summarise(min = min(value)
                   ,max = max(value)
                   ,mean = mean(value)
                   ,n = n()
                   ,sd = sd(value)
                   ,q25 = quantile(value, .25)
                   ,q75 = quantile(value, .75)) 

pdf(paste0(out_dir,"cm_rel_ab_phyla_cohorts.pdf"))
ggplot(summs_for_parallel_coo_cohorts, aes(fill=taxa_2, y=mean, x=date)) + 
  geom_bar(position="fill", stat="identity")+
  facet_wrap(~cohort) +
  theme(axis.title.y = element_text(),
        legend.title = element_text()) +
  labs(x = "collection date",
       y = "relative abundance (ratio)",
       fill = "Phylum") 
dev.off()

# higher proportion of Tenericutes in Neomycin from t5 to t10 compared to other cohorts?
# Euryarchaeota "disapper" from t5 to t10. 
# Actinobacteria slowly increase in proportion with time
# igher bacteroidetes proportion between t3 and t4 in Neo than other cohorts? 
# Synergistetes prettu much disappear after t1 in Control, Dscour and ColiGuard,
# whereas they disappera later (after t2) in Neo and NeoD. 
# In NeoC in t4 still high proportion before gone. 
# Proteobacterai from t4 on higher proportion in Dscour, low in any other cohort. 


######################################################################################################

# PHYLA composition change through time: PARALLEL COORDINATES 

# keep only rows that contain phyla info, discard others; also exclude mothers
df1 <- df %>%
  dplyr::select(cohort,pig,bin,date,taxa_2,value) %>%
  dplyr::filter(!cohort=="Mothers")

# remove "p__" before phylum 
df1[,5] <- lapply(
  df1[,5], 
  gsub, 
  pattern = "p__", 
  replacement = "", 
  fixed = TRUE)

df1 <- na.omit(df1)

df2 <- df1

# lib size normalization
df2 <- df2 %>% 
  group_by(pig,date) %>% 
  dplyr::mutate(norm_value = (value/sum(value))) %>% 
  dplyr::select(pig, date, taxa_2, value, norm_value)

# sum all the norm values that fall within same pig,date,phylum:
df2 <- df2 %>%  
  group_by(pig,date,taxa_2) %>%
  dplyr::summarise(indiv_sum = sum(norm_value))

# take the mean of each phylum by date:
df2 <- df2 %>%
  group_by(taxa_2,date) %>%
  dplyr::summarize(mean = mean(indiv_sum, na.rm=TRUE)) %>%
  dplyr::mutate(perc=mean*100) %>%
  dplyr::select(taxa_2,date,perc)


df2 <- as.data.frame(df2)

pdf(paste0(out_dir,"cm_parallel_coordinates_phyla.pdf"))
ggplot(df2, aes(x=date, y=log10(perc), group=taxa_2, color=taxa_2)) + 
  geom_line() + geom_point(size=0.8)+
  theme_bw()+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "collection date",
       y = "relative abundance (log10)",
       color = "Phylum")  +
  theme(legend.title = element_text()) 
dev.off()


######################################################################################################

# PHYLA composition change through time: BALLOONPLOTS - all


df3 <- df2
# pivot wider
df3 <- df3 %>%
  dplyr::select(taxa_2,date,perc) %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))

df4 <- df3
df4 <- as.data.frame(df4)
rownames(df4) <- df4[,1]
df4[,1] <- NULL
m <- as.table(as.matrix(df4))

pdf(paste0(out_dir,"cm_balloonplot_phyla.pdf"))
balloonplot(t(m), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
dev.off()


######################################################################################################

# PHYLA composition change through time: HEATMAP - all

phyla_pheatmap <- pheatmap(t(m), display_numbers = T,
                           cluster_rows = F, cluster_cols = F, fontsize_number = 8,
                           fontsize_row = 8,show_rownames = TRUE)

pdf(paste0(out_dir,"cm_heatmap_phyla.pdf"))
phyla_pheatmap
dev.off()




######################################################################################################

# PHYLA composition change through time: BALLOONPLOTS - by cohort

df2 <- df1

# splitting into multiple dataframes (by cohort)
multi_df <- split( df2 , f = df2$cohort )

# construct an empty dataframe to build on 
final <- data.frame(
  pig = character(),
  bin = character(),
  date = character(),
  taxa_2 = character(),
  value = character(),
  stringsAsFactors = FALSE
)


for (single_df in multi_df) {
  
  single_df <- as.data.frame(single_df)
  coho <- as.character(single_df$cohort[1])
  
  # lib size normalization
  df2 <- single_df %>% 
    group_by(pig,date) %>% 
    dplyr::mutate(norm_value = (value/sum(value))) %>% 
    dplyr::select(pig, date, taxa_2, value, norm_value)
  
  # sum all the norm values that fall within same pig,date,phylum:
  df2 <- df2 %>%  
    group_by(pig,date,taxa_2) %>%
    dplyr::summarise(indiv_sum = sum(norm_value))
  
  # take the mean of each phylum by date:
  df2 <- df2 %>%
    group_by(taxa_2,date) %>%
    dplyr::summarize(mean = mean(indiv_sum, na.rm=TRUE)) %>%
    dplyr::mutate(perc=mean*100) %>%
    dplyr::select(taxa_2,date,perc)
  
  df2 <- as.data.frame(df2)
  df2$cohort <- coho
  
  final <- rbind(final,df2)
  
}

pdf(paste0(out_dir,"cm_balloonplot_phyla_cohorts.pdf"))
layout(matrix(c(1:3), 3, 1) )
par(mar = c(0, 4.1,2.5, 2.1),oma = c(2, 0, 2, 0))
# Control
fin <- final %>%
  dplyr::filter(cohort=="Control") %>%
  dplyr::select(taxa_2,date,perc) %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))
fin <- as.data.frame(fin)
rownames(fin) <- fin[,1]
fin[,1] <- NULL
m <- as.table(as.matrix(fin))
balloonplot(t(m), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
# DScour
fin <- final %>%
  dplyr::filter(cohort=="D-Scour") %>%
  dplyr::select(taxa_2,date,perc) %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))
fin <- as.data.frame(fin)
rownames(fin) <- fin[,1]
fin[,1] <- NULL
m <- as.table(as.matrix(fin))
balloonplot(t(m), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
# ColiGuard
fin <- final %>%
  dplyr::filter(cohort=="ColiGuard") %>%
  dplyr::select(taxa_2,date,perc) %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))
fin <- as.data.frame(fin)
rownames(fin) <- fin[,1]
fin[,1] <- NULL
m <- as.table(as.matrix(fin))
balloonplot(t(m), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
layout(matrix(c(1:3), 3, 1) )
par(mar = c(0, 4.1,2.5, 2.1),oma = c(2, 0, 2, 0))
# Neomycin
fin <- final %>%
  dplyr::filter(cohort=="Neomycin") %>%
  dplyr::select(taxa_2,date,perc) %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))
fin <- as.data.frame(fin)
rownames(fin) <- fin[,1]
fin[,1] <- NULL
m <- as.table(as.matrix(fin))
balloonplot(t(m), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
# NeoD
fin <- final %>%
  dplyr::filter(cohort=="NeoD") %>%
  dplyr::select(taxa_2,date,perc) %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))
fin <- as.data.frame(fin)
rownames(fin) <- fin[,1]
fin[,1] <- NULL
m <- as.table(as.matrix(fin))
balloonplot(t(m), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
# NeoC
fin <- final %>%
  dplyr::filter(cohort=="NeoC") %>%
  dplyr::select(taxa_2,date,perc) %>%
  pivot_wider(names_from = taxa_2, values_from = perc, values_fill = list(perc = 0))
fin <- as.data.frame(fin)
rownames(fin) <- fin[,1]
fin[,1] <- NULL
m <- as.table(as.matrix(fin))
balloonplot(t(m), main =NULL, xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=40, colmar=2,
            text.size=0.7)
dev.off()



######################################################################################################
######################################################################################################

# Principal component analysis: clustering of bins based on full taxa with time 

df1 <- df 
NROW(df1)
head(df1)

# define full_taxa 
df1$full_taxa <- paste0(df1$taxa_1,"..",
                        df1$taxa_2,"..",
                        df1$taxa_3,"..",
                        df1$taxa_4,"..",
                        df1$taxa_5,"..",
                        df1$taxa_6,"..",
                        df1$taxa_7,"..")

df1 <- df1 %>%
  dplyr::select(pig,bin,date,value,full_taxa,cohort)

df1 <- as.data.frame(na.omit(df1))
NCOL(df1)
NROW(df1)
head(df1)

unique(df1$full_taxa)


#################################
# STEP 1.


# for each sample (pig,date), sum up the counts that fall within one species (same species assigned to distinct bins)
df2 <- df1 %>%
  group_by(pig,full_taxa,date) %>%
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

##############################################################################################
# CLR TRANSFORM AND pivot wider 

# pivot wider
df3 <- df3 %>%
  dplyr::select(sample,full_taxa,norm_value) %>%
  dplyr::mutate(norm_value=as.numeric(clr(norm_value))) %>%
  pivot_wider(names_from = full_taxa, values_from = norm_value, values_fill = list(norm_value = 0))
##############################################################################################

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


cm_PC12 <- ggbiplot(df6.pca,
                    labels=this_mat_samples$sample_2,
                    groups=this_mat_samples$sample_1,
                    ellipse=TRUE,
                    var.axes = FALSE,
                    labels.size = 2,
                    choices = (1:2)) +
  theme_bw() +
  #xlim(c(-2,1)) +
  scale_colour_discrete(name="timepoint")+
  guides(color = guide_legend(ncol = 1))

cm_PC34 <- ggbiplot(df6.pca,
                    labels=this_mat_samples$sample_2,
                    groups=this_mat_samples$sample_1,
                    ellipse=TRUE,
                    var.axes = FALSE,
                    labels.size = 2,
                    choices = (3:4)) +
  theme_bw() +
  scale_colour_discrete(name="timepoint")+
  guides(color = guide_legend(ncol = 1))


cm_PCA <- ggarrange(cm_PC12,cm_PC34,
                    ncol=2,legend = "right",
                    common.legend=TRUE)

pdf(paste0(out_dir,"cm_PCA.pdf"), width=7,height=4)
cm_PCA
dev.off()












# normalization for library size 
df2 <- df1 %>%
  dplyr::group_by(pig,date) %>%
  dplyr::mutate(norm_value=value/sum(value)) 
NROW(df2)
head(df2)

# # test:
# test <- df2 %>%
#   filter(pig=="14159") %>%
#   filter(date=="t0") %>%
#   dplyr::mutate(norm_value=value/sum(value))
# head(test)
# sum(test$norm_value)

#################################
# STEP 2.

# sum all the norm values that fall within same pig,date,taxa_2
df3 <- df2 %>%
  dplyr::group_by(pig,date,full_taxa) %>%
  dplyr::summarise(indiv_sum = sum(norm_value))
head(df3)

# # test:
# test2 <- test %>%
#   group_by(full_taxa) %>%
#   dplyr::summarise(indiv_sum = sum(norm_value))
# head(test2)
# sum(test2$indiv_sum)

#################################
# STEP 3.



##############################################################################################
# CLR TRANSFORM AND pivot wider 

# pivot wider
# df4 <- df3 %>%
#   dplyr::select(pig,date,full_taxa,indiv_sum) %>%
#   dplyr::mutate(indiv_sum=as.numeric(clr(indiv_sum))) %>%
#   pivot_wider(names_from = full_taxa, values_from = indiv_sum, values_fill = list(indiv_sum = 0))
##############################################################################################


# long to wide format
df4 <- df3 %>%
  pivot_wider(names_from = full_taxa, values_from = indiv_sum, values_fill = list(indiv_sum = 0)) 
head(df4)

# # test:
# test3 <- test2 %>%
#   pivot_wider(names_from = full_taxa, values_from = indiv_sum, values_fill = list(indiv_sum = 0))
# head(test3)
# sum(test3[1,])


#################################



# get a quick cohorts to pig table
cohorts <- df %>% dplyr::select(cohort,pig,date) %>% distinct()
cohorts$sample <- paste0(cohorts$date,"_",cohorts$pig)
cohorts <- as.data.frame(cohorts)


df5 <- inner_join(cohorts,df4) 
df5$sample <- paste0(df5$date,"_",df5$cohort)

df5$pig <- NULL
df5$date <- NULL
df5$cohort <- NULL



df6 <- df5 %>%
  group_by(sample) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)


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

cm_PC12 <- ggbiplot(df6.pca,
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
cm_PC34 <- ggbiplot(df6.pca,
                    labels=this_mat_samples$sample_2,
                    groups=this_mat_samples$sample_1,
                    ellipse=TRUE,
                    var.axes = FALSE,
                    labels.size = 2,
                    choices = (3:4)) +
  theme_bw() +
  scale_colour_discrete(name="timepoint")+
  guides(color = guide_legend(ncol = 1))


cm_PCA <- ggarrange(cm_PC12,cm_PC34,
                    ncol=2,legend = "right",
                    common.legend=TRUE)

pdf(paste0(out_dir,"cm_PCA.pdf"), width=7,height=4)
cm_PCA
dev.off()

