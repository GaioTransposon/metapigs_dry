

library(readr)
library(dplyr)
library(tidyr)
library(data.table)
library(splitstackshape)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(treemapify)
library(ggrepel)
library(factoextra)
library(stringr)


setwd("~/Desktop/metapigs_dry/dbcan")
basedir = "~/Desktop/metapigs_dry/"

########################

# load dbcan_diamond data (output of mapping our-bins-AAs against CAZy db)
hmmer <- read_table2("hmmer.out", col_names = FALSE)

########################

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv(paste0(basedir,"gtdbtk/gtdbtk_bins_completeTaxa"),
                        col_types = cols(node = col_character(),
                                         pig = col_character()))

########################

# counts data 

no_reps_all <- read.csv(paste0(basedir,"no_reps_all.csv"), 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)


# remove .fa extension to match bins in checkm df 
no_reps_all$bin <- gsub(".fa","", no_reps_all$bin)
head(no_reps_all)
NROW(no_reps_all)

########################

# load checkM assignments of the bins
checkm_all_nearly <- read_delim(paste0(basedir,"checkm/checkm_all_nearly"), 
                                "\t", escape_double = FALSE, col_types = cols(pigid = col_character()), 
                                trim_ws = TRUE)

##########################################################
##########################################################


# PART 1

head(hmmer)
colnames(hmmer)

# hmmer data cleaning 

# remove .hmm from enzyme ID
hmmer <- cSplit(hmmer,"X1",".")

# separate by _
hmmer <- cSplit(hmmer,"X1_1","_")

# copy enzyme ID into new column 
hmmer$enzymeID <- hmmer$X1_1_1

# extra enzyme ID - suffixes together 
hmmer <- hmmer  %>%
  mutate_at(c('X1_1_2', 'X1_1_3','X1_1_4','X1_1_5'), ~str_replace_na(., "")) %>%
  mutate(combo_var = paste0(X1_1_2,".", X1_1_3,".",X1_1_4,".",X1_1_5))

unique(hmmer$enzymeID)

# subject/bin/contig
hmmer <- cSplit(hmmer,"X3","_")

# columns selection
hmmer <- hmmer %>%
  dplyr::select(enzymeID,combo_var,
                X3_1,X3_2, X3_3, # pig, bin, contig number 
                X2,X4, # length enzyme, length subject protein
                X5,X6,X7,X8,X9,X10) # evalue, start query, end query, start sample, end sample, coverage  

# get enzyme class
hmmer$getclass <- hmmer$enzymeID
hmmer <- hmmer %>%
  separate(getclass, 
           into = c("text", "num"), 
           sep = "(?<=[A-Za-z])(?=[0-9])"
  )

hmmer$num <- NULL

# cols renaming
colnames(hmmer) <- c("enzymeID","enzymeID_suffixes",
                     "pig","bin","contig",
                     "HMM_length","query_lenght",
                     "evalue",
                     "start HMM","end HMM",
                     "start sample","end sample",
                     "coverage", "enzymeNAME")
                     
hmmer$bin <- gsub(".fa","", hmmer$bin)

head(hmmer)

########################################################################

# QUALITY of hmmer data : 
# how well our predicted proteins - matched against the CAZy database 
sink(file="dbcan_HMMER_numbers.txt", append=FALSE)
paste0("total number of proteins matching the CAZy db with dbcan-hmmer")
NROW(hmmer)
paste0("total number of unique enzymeIDs with dbcan-hmmer")
NROW(unique(hmmer$enzymeID))
paste0("coverage")
summary(hmmer$coverage)
paste0("e-values")
summary(hmmer$evalue)
paste0("e-values")
hmmer %>%
  group_by(enzymeNAME) %>%
  tally() %>%
  mutate(perc=n/sum(n))
sink()


# re-rder enzymes to follow the same aesthetics as Stewart et al 2019
unique(hmmer$enzymeNAME)
levels(hmmer$enzymeNAME)

hmmer$enzymeNAME <- as.character(hmmer$enzymeNAME)
hmmer$enzymeNAME  = factor(
  hmmer$enzymeNAME, levels=c("AA",
                               "CE",
                               "SLH",
                               "GH",
                               "GT",
                               "CBM",
                               "cohesin",
                               "PL"))

# display.brewer.pal(8, "Spectral")

# PLOT 
hmmer_perc_identity_plot <- hmmer %>%
  group_by(enzymeNAME) %>%
  ggplot(.,aes(enzymeNAME,coverage*100,fill=enzymeNAME))+
  ylab("Percentage identity")+
  geom_boxplot(coef=1e30, lwd=0.1)+ # high value cutoff to deact outliers; pdf becomes otherwise heavy to move around
  ylim(0,100)+
  scale_fill_brewer(palette="Spectral")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=9),
        axis.text.y=element_text(size=8),
        axis.title.y = element_text(size=9),
        axis.title.x=element_blank(),
        legend.position="none",legend.text = element_text(size=8),
        legend.title = element_blank(),title=element_text(size=9))+
  ggtitle("Percentage identity against CAZy")

enzymes_proportion <- hmmer %>% 
  group_by(enzymeNAME) %>% 
  tally() %>%
  mutate(perc = paste0(round(n/sum(n)*100,2),"%"))

hmmer_perc_identity_plot <- hmmer_perc_identity_plot + 
  geom_text(data = enzymes_proportion,
            aes(enzymeNAME, 0, label = perc), vjust="inward",size=2)

pdf("dbcan_HMMER_perc_identity.pdf")
hmmer_perc_identity_plot
dev.off()

########################################################################

# # clean checkM all nearly data:
# 
# # filter >90 <5 bins 
# checkm_all_nearly <- dplyr::filter(checkm_all_nearly, !grepl("Completeness",Completeness))
# checkm_all_nearly <- subset(checkm_all_nearly, Completeness >= 90 & Contamination <= 5)
# # some formatting 
# checkm_all_nearly$Completeness <- as.numeric(checkm_all_nearly$Completeness)
# checkm_all_nearly$Contamination <- as.numeric(checkm_all_nearly$Contamination)
# # rename cols to matching colnames between dataframes to merge  
# colnames(checkm_all_nearly)[colnames(checkm_all_nearly)=="pigid"] <- "pig"
# colnames(checkm_all_nearly)[colnames(checkm_all_nearly)=="Bin Id"] <- "bin"
# colnames(checkm_all_nearly)[colnames(checkm_all_nearly)=="Taxonomy (contained)"] <- "taxa"
# 
# checkm_clean <- checkm_all_nearly %>%
#   select(pig,bin,taxa)
# checkm_clean <- cSplit(checkm_clean, "taxa", sep=";")
# colnames(checkm_clean) <- c("pig","bin","kingdom","phylum","class","order","family","genus","species")
# head(checkm_clean)
# 
# checkm_clean$phylum <- as.character(checkm_clean$phylum)
# checkm_clean$phylum[is.na(checkm_clean$phylum)] <- "Unknown"
# 
# checkm_clean$phylum <- gsub("p__","", checkm_clean$phylum)
# 
# ########################################################################
# 
# # merge diamond data with checkM- nearly complete bins - taxonomic assignments 
# 
# NROW(checkm_clean)
# NROW(hmmer)
# 
# cm_hmmer <- inner_join(checkm_clean,hmmer)
# NROW(cm_hmmer)
# 
# 
# cm_hmmer_taxa <- cm_hmmer %>%
#   group_by(phylum,enzymeNAME) %>%
#   tally() %>%
#   mutate(Proportion = n) 
# 
# cm_hmmer_prop <- cm_hmmer %>%
#   group_by(phylum) %>%
#   tally() %>%
#   mutate(proportion=n/sum(n)*100)
# 
# a <- ggplot(cm_hmmer_prop, aes(y=proportion, x=phylum)) + 
#   geom_bar(position="dodge", lwd=0.1,stat="identity",colour="black",fill="snow3")+ #lightskyblue happier
#   theme(axis.text.x=element_text(angle=90))+
#   theme_pubr()+
#   ylab("Proportion of CAZymes (%)")+
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank())
# b <- ggplot(cm_hmmer_taxa, aes(fill=enzymeNAME, y=Proportion, x=phylum)) + 
#   geom_bar(position="fill", lwd=0.1,stat="identity",colour="black")+
#   theme(axis.text.x=element_text(angle=90))+
#   scale_fill_brewer(palette="Spectral")+
#   theme_pubr()+
#   theme(axis.text.x=element_text(angle=45,hjust=1),
#         axis.title.x = element_blank(),
#         legend.position="right",
#         legend.title = element_blank())
# 
# # RColorBrewer::display.brewer.all(n=8,select = "Spectral")
# 
# CAZ_HMMER_CMphyla_plot <- plot_grid(
#   plot_grid(
#     a + theme(legend.position = "none"),
#     b + theme(legend.position = "none"),
#     ncol = 1,
#     align = "v"),
#   plot_grid(
#     get_legend(a),
#     get_legend(b),
#     ncol =1),
#   rel_widths = c(8,2)
# )
# 
# pdf("dbcan_HMMER_CAZ_CMphyla.pdf")
# CAZ_HMMER_CMphyla_plot
# dev.off()
# 
# 
# # perc identity plot + phyla distr plot in one
# pdf("dbcan_HMMER_CAZ_PercID_Phyla.pdf", width=7,height=5)
# plot_grid(hmmer_perc_identity_plot,
#           CAZ_HMMER_GTphyla_plot,
#           nrow=1,ncol=2,
#           rel_width=c(4,6),
#           labels=c("A","B"))
# dev.off()

########################################################################

# merge diamond data with gtdbtk- bins taxonomic assignments 


NROW(gtdbtk_bins)
NROW(hmmer)

gt_hmmer <- inner_join(gtdbtk_bins,hmmer)
NROW(gt_hmmer)

gt_hmmer$phylum[is.na(gt_hmmer$phylum)] <- "Unknown"

gt_hmmer_taxa <- gt_hmmer %>%
  group_by(phylum,enzymeNAME) %>%
  tally() %>%
  mutate(Proportion = n) 

gt_hmmer_prop <- gt_hmmer %>%
  group_by(phylum) %>%
  tally() %>%
  mutate(proportion=n/sum(n)*100)

a_gt <- ggplot(gt_hmmer_prop, aes(y=proportion, x=phylum)) + 
  geom_bar(position="dodge", stat="identity",colour="black",fill="snow3")+ #lightskyblue happier
  theme(axis.text.x=element_text(angle=90))+
  theme_pubr()+
  ylab("Proportion of CAZymes (%)")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=8),
        axis.text.x = element_blank())
b_gt <- ggplot(gt_hmmer_taxa, aes(fill=enzymeNAME, y=Proportion, x=phylum)) + 
  geom_bar(position="fill", stat="identity",colour="black")+
  theme(axis.text.x=element_text(angle=90))+
  scale_fill_brewer(palette="Spectral")+
  theme_pubr()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=6),
        axis.text.y=element_text(size=6),
        axis.title.y = element_text(size=8),
        axis.title.x=element_blank(),
        legend.position="right",legend.text = element_text(size=8),
        legend.title = element_blank())


CAZ_HMMER_GTphyla_plot <- plot_grid(
  plot_grid(
    a_gt + theme(legend.position = "none"),
    b_gt + theme(legend.position = "none"),
    ncol = 1,
    rel_heights=c(4,6),
    align = "v"),
  plot_grid(
    get_legend(a_gt),
    get_legend(b_gt),
    ncol =1,
    rel_heights = c(3,7)),
  rel_widths = c(8,2)
)

pdf("dbcan_HMMER_CAZ_GTphyla.pdf")
CAZ_HMMER_GTphyla_plot
dev.off()


##########################################################
##########################################################


# PART 2


gt_hmmer <- as.data.frame(gt_hmmer)

# normalize by lib size: 
no_reps_all_norm <- no_reps_all %>%
  group_by(cohort,pig,date) %>%
  mutate(norm_value=value/sum(value))

hmmer_count <- hmmer %>%
  group_by(pig,bin,enzymeID,enzymeNAME)%>%
  dplyr::summarise(enz_count=n())


df <- left_join(no_reps_all_norm,hmmer_count) 


df1 <- df %>%
  mutate(enz_countBYmap=norm_value*enz_count)

df1 <- df1[!is.na(df1$enzymeID),]

df2 <- df1 %>%
  group_by(cohort,pig,date,enzymeID,enzymeNAME) %>%
  dplyr::summarise(norm_count=sum(enz_countBYmap))


# reorder dates 
df2$date  = factor(df2$date, levels=c("t0",
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
                                      "tM"))



df2 <- as.data.frame(df2)

df_piggies <- df2 %>%
  filter(!cohort=="Mothers")

# split df by enzymeID
multiple_DFs <- split( df_piggies , f = df_piggies$enzymeID ,drop = TRUE)

# total number of unique CAZ IDs
NROW(multiple_DFs)

# empty df
all_pvalues <- data.frame(group1=character(),
                          group2=character(),
                          p_value=double(),
                          enzID=character(),
                          test=character(),
                          p_adj_method=character())


for (single_DF in multiple_DFs) {
  
  single <- as.data.frame(single_DF)
  
  
  if (NROW(unique(single$date))>10) {
    
    enzID <- single$enzymeID[1]
    
    out <- pairwise.wilcox.test(single$norm_count, single$date,
                                p.adjust.method = "bonferroni")
    
    out2 <- as.data.frame(out$p.value)
    out2$group1 <- rownames(out2)
    out2 <- out2 %>% pivot_longer(cols = -group1,names_to = "group2", values_to="p_value")
    out2 <- as.data.frame(out2[!is.na(out2$p_value),])
    out2$enzID = paste0(enzID)
    out2$test <- "pairwise.wilcox.test"
    out2$p_adj_method = "Bonferroni"
    
    all_pvalues = rbind(all_pvalues,out2)
    
  }
  
  else (print("not enough observations"))
}


write.csv(all_pvalues,
          "dbcan_HMMER_pvalues.csv", 
          row.names = FALSE)

significant <- all_pvalues %>%
  filter(p_value<0.05) %>%
  arrange(p_value)
tail(significant)
mylist <- unique(significant$enzID)

sink(file="dbcan_HMMER_numbers.txt", append=TRUE)
paste0("number of unique enzymes that showed a significant change between timepoints")
NROW(unique(significant$enzID))
sink()

NROW(unique(hmmer$enzymeID)) # tot enzymes
NROW(unique(significant$enzID)) # tot enzymes with significant diffs when comparing any time point
##########################################################
##########################################################


# PART 3

# now I will use these IDs to plot interesting stuff
df_part <- subset(df2, (enzymeID %in% mylist))

# renaming this column 
colnames(df_part)[6] <- "tot"

# getting a tally of number of piglets carrying each specific enzyme
piglets_with_enzyme <- df_part %>%
  filter(!date=="tM") %>%
  ungroup() %>%
  select(enzymeID,pig) %>%
  distinct() %>%
  group_by(enzymeID) %>%
  tally() 

# getting a tally of number of moms carrying each specific enzyme
moms_with_enzyme <- df_part %>%
  filter(date=="tM") %>%
  ungroup() %>%
  select(enzymeID,pig) %>%
  distinct() %>%
  group_by(enzymeID) %>%
  tally() %>%
  mutate(n_moms=n) %>%
  dplyr::select(enzymeID,n_moms)

df_part <- as.data.frame(inner_join(df_part,piglets_with_enzyme))
df_part <- as.data.frame(inner_join(df_part,moms_with_enzyme))

df_part$enzymeNAME  = factor(
  df_part$enzymeNAME, levels=c("AA",
                               "CE",
                               "SLH",
                               "GH",
                               "GT",
                               "CBM",
                               "cohesin",
                               "PL"))

# set defined colors (releveling )
scale_fill_gaio8 <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(c("#D53E4F","#F46D43","#FDAE61","#FEE08B","#E6F598","#ABDDA4","#66C2A5","#3288BD"), 
                      c("AA","CE","SLH","GH","GT","CBM","cohesin","PL")), 
    ...
  )
}

# levels(df_part$enzymeNAME)
# brewer.pal(8, name="Spectral")

# put df rows in order of list (list is ordered by p-value (descending))
require(gdata)
df_part$enzymeID <- reorder.factor(df_part$enzymeID, new.order=mylist)
df_part <- df_part %>%
  arrange(enzymeID)
# now they are plotted in order of significance! 
# that means that the first pages will be most interesting
# 

## commented this chunk out as redudant with the ones below, 
## however this one is nicer (larger) 
# pdf("dbcan_HMMER_CAZ_time_boxplots.pdf")
# for (i in seq(1, length(mylist), 12)) {    # can also use: length(unique(df_part$enzymeID))
#   
#   print(ggplot(df_part[df_part$enzymeID %in% mylist[i:(i+11)], ], 
#                aes(date, log(tot),fill=enzymeNAME)) + 
#           geom_boxplot(outlier.size = 1) +
#           facet_wrap(~ enzymeID, scales = "free_y") +
#           theme_bw()+
#           scale_fill_gaio8()+
#           theme(legend.position="none",
#                 axis.text.y=element_text(size = 4),
#                 axis.ticks.length.y = unit(.05, "cm"),
#                 axis.text.x=element_text(size=6))+
#           geom_text(aes(x="t1", y=max(log(tot)), label=paste0("n=(",n,")")),
#                     size=2,colour="black", inherit.aes=TRUE, parse=FALSE,check_overlap = TRUE)+
#           geom_text(aes(x="t10", y=max(log(tot)), label=paste0("n=(",n_moms,")")),
#                     size=2,colour="black", inherit.aes=TRUE, parse=FALSE,check_overlap = TRUE))
# }
# dev.off()
# 

##########################################################


# HEATMAP of enzymes trend over time 

df_part$sample <- paste0(df_part$date,"_",df_part$pig)

# reorder dates 
df_part$date  = factor(df_part$date, levels=c("t0",
                                              "t1", 
                                              "t2",
                                              "t3",
                                              "t4",
                                              "t5",
                                              "t6",
                                              "t7",
                                              "t8",
                                              "t9",
                                              "tM"))

df_part <- df_part[!is.na(df_part$date),]
head(df_part)

# subsetting whole

# AA & CE
df_part_AA_CE <- subset(df_part, (enzymeNAME %in% c("AA","CE")))
list_AA_CE <- unique(df_part_AA_CE$enzymeID)

# PL 
df_part_PL <- subset(df_part, (enzymeNAME %in% "PL"))
list_PL <- unique(df_part_PL$enzymeID)

# GT
df_part_GT <- subset(df_part, (enzymeNAME %in% "GT"))
NROW(unique(df_part_GT$enzymeID))
list_GT_1 <- unique(df_part_GT$enzymeID)[1:15]
list_GT_2 <- unique(df_part_GT$enzymeID)[16:30]
list_GT_3 <- unique(df_part_GT$enzymeID)[31:45]
df_part1_GT <- subset(df_part_GT, (enzymeID %in% list_GT_1))
df_part2_GT <- subset(df_part_GT, (enzymeID %in% list_GT_2))
df_part3_GT <- subset(df_part_GT, (enzymeID %in% list_GT_3))

# GH
df_part_GH <- subset(df_part, (enzymeNAME %in% "GH"))
NROW(unique(df_part_GH$enzymeID))
list_GH_1 <- unique(df_part_GH$enzymeID)[1:25]
list_GH_2 <- unique(df_part_GH$enzymeID)[26:50]
list_GH_3 <- unique(df_part_GH$enzymeID)[51:75]
list_GH_4 <- unique(df_part_GH$enzymeID)[76:106]
df_part1_GH <- subset(df_part_GH, (enzymeID %in% list_GH_1))
df_part2_GH <- subset(df_part_GH, (enzymeID %in% list_GH_2))
df_part3_GH <- subset(df_part_GH, (enzymeID %in% list_GH_3))
df_part4_GH <- subset(df_part_GH, (enzymeID %in% list_GH_4))

# CBM
df_part_CBM <- subset(df_part, (enzymeNAME %in% "CBM"))
NROW(unique(df_part_CBM$enzymeID))
list_CBM_1 <- unique(df_part_CBM$enzymeID)[1:20]
list_CBM_2 <- unique(df_part_CBM$enzymeID)[21:40]
df_part1_CBM <- subset(df_part_CBM, (enzymeID %in% list_CBM_1))
df_part2_CBM <- subset(df_part_CBM, (enzymeID %in% list_CBM_2))

# SLH
df_part_SLH <- subset(df_part, (enzymeNAME %in% "SLH"))
NROW(unique(df_part_SLH$enzymeID))
list_SLH <- unique(df_part_SLH$enzymeID)
df_part_SLH <- subset(df_part_SLH, (enzymeID %in% list_SLH))

# cohesin
df_part_cohesin <- subset(df_part, (enzymeNAME %in% "cohesin"))
NROW(unique(df_part_cohesin$enzymeID))
list_cohesin <- unique(df_part_cohesin$enzymeID)
df_part_cohesin <- subset(df_part_cohesin, (enzymeID %in% list_cohesin))


######
# get species count for each enzymeID and make subsets of the data 

t <- gt_hmmer %>%
  dplyr::select(pig,enzymeID,enzymeNAME,species,family) %>%   # dplyr::select(pig,contig,enzymeID,enzymeNAME,species,family) %>%
  distinct() %>%
  group_by(pig,enzymeID,species) %>%
  add_tally() %>%
  group_by(enzymeNAME,enzymeID,species) %>%
  dplyr::summarise(n_sum_species=sum(n)) %>%
  mutate(perc_species=round(n_sum_species/sum(n_sum_species)*100,2)) %>% 
  drop_na()
s <- as.data.frame(t)

t_genus <- gt_hmmer %>%
  dplyr::select(pig,enzymeID,enzymeNAME,genus,family) %>%   # dplyr::select(pig,contig,enzymeID,enzymeNAME,species,family) %>%
  distinct() %>%
  group_by(pig,enzymeID,genus) %>%
  add_tally() %>%
  group_by(enzymeNAME,enzymeID,genus) %>%
  dplyr::summarise(n_sum_genus=sum(n)) %>%
  mutate(perc_genus=round(n_sum_genus/sum(n_sum_genus)*100,2)) %>% 
  drop_na()
s_genus <- as.data.frame(t_genus)


# subsets (species)
s_part1_GH <- subset(s, (enzymeID %in% list_GH_1))
s_part2_GH <- subset(s, (enzymeID %in% list_GH_2))
s_part3_GH <- subset(s, (enzymeID %in% list_GH_3))
s_part4_GH <- subset(s, (enzymeID %in% list_GH_4))
s_part1_GT <- subset(s, (enzymeID %in% list_GT_1))
s_part2_GT <- subset(s, (enzymeID %in% list_GT_2))
s_part3_GT <- subset(s, (enzymeID %in% list_GT_3))
s_part1_CBM <- subset(s, (enzymeID %in% list_CBM_1))
s_part2_CBM <- subset(s, (enzymeID %in% list_CBM_2))
s_part_AA_CE <- subset(s, (enzymeID %in% list_AA_CE))
s_part_PL <- subset(s, (enzymeID %in% list_PL))
s_part_SLH <- subset(s, (enzymeID %in% list_SLH))
s_part_cohesin <- subset(s, (enzymeID %in% list_cohesin))
######

# subsets (genus)
s_genus_part1_GH <- subset(s_genus, (enzymeID %in% list_GH_1))
s_genus_part2_GH <- subset(s_genus, (enzymeID %in% list_GH_2))
s_genus_part3_GH <- subset(s_genus, (enzymeID %in% list_GH_3))
s_genus_part4_GH <- subset(s_genus, (enzymeID %in% list_GH_4))
s_genus_part1_GT <- subset(s_genus, (enzymeID %in% list_GT_1))
s_genus_part2_GT <- subset(s_genus, (enzymeID %in% list_GT_2))
s_genus_part3_GT <- subset(s_genus, (enzymeID %in% list_GT_3))
s_genus_part1_CBM <- subset(s_genus, (enzymeID %in% list_CBM_1))
s_genus_part2_CBM <- subset(s_genus, (enzymeID %in% list_CBM_2))
s_genus_part_AA_CE <- subset(s_genus, (enzymeID %in% list_AA_CE))
s_genus_part_PL <- subset(s_genus, (enzymeID %in% list_PL))
s_genus_part_SLH <- subset(s_genus, (enzymeID %in% list_SLH))
s_genus_part_cohesin <- subset(s_genus, (enzymeID %in% list_cohesin))
######

##########################################################
##########################################################


# PART 4


# # function for making enzyme heatmaps (whole large heatmap)
# make_enzyme_heatmap <- function(x) {
#   p <- x %>%
#     filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="t10"|date=="tM") %>%
#     group_by(pig, enzymeID, date) %>%
#     dplyr::summarise(tot = mean(tot, na.rm = TRUE)) %>%
#     group_by(enzymeID) %>%
#     mutate(tot = tot/max(tot)) %>%
#     ungroup() %>%
#     ggplot(aes(x = factor(pig), y = reorder(enzymeID, tot, FUN = mean), fill = tot)) +
#     geom_tile() +
#     #scale_fill_gaio8()+
#     scale_fill_distiller(type = "div", palette = "Spectral") +
#     facet_grid(~date, scales = "free") +
#     labs(x = "date", y = "enzymeID", fill = "normalized abundance")+
#     theme(axis.text.x=element_blank(),
#           axis.title.x=element_blank(),
#           legend.position="top")
#   return(p)
# }

# # function for making species heatmaps 
# make_species_heatmap <- function(species_df,CAZ_heatmap) {
#   
#   CAZ_order_vector <- CAZ_heatmap$data %>% group_by(enzymeID) %>% dplyr::summarise(mean=mean(tot)) %>% 
#     group_by(enzymeID) %>%
#     arrange(desc(mean))
#   
#   g <- species_df %>% 
#     group_by(enzymeID,species) %>% 
#     dplyr::summarise(n_sum_species = mean(n_sum_species, na.rm = TRUE)) %>% 
#     group_by(enzymeID) %>%
#     mutate(n_sum_species = n_sum_species/sum(n_sum_species)) %>%
#     top_n(n = 4, wt = n_sum_species) %>%
#     arrange(desc(n)) %>%
#     slice(1:4) %>%
#     pivot_wider(names_from = enzymeID,values_from=n_sum_species) %>% 
#     pivot_longer(cols=-species,names_to="enzymeID",values_to = "n_sum_species",values_drop_na = FALSE) %>%
#     ungroup() 
#   
#   #require(gdata)
#   g$enzymeID <- reorder.factor(g$enzymeID, new.order=CAZ_order_vector$enzymeID)
#   
#   p <- ggplot(g, aes(x = reorder(enzymeID, n_sum_species, FUN = mean), y = species, fill = n_sum_species)) +
#     geom_tile(size = 0.5, color = "black") +
#     #scale_fill_gaio8()+
#     scale_fill_distiller(palette = "Spectral",na.value = "black") +
#     labs(x = "date", y = "enzymeID", fill = "normalized abundance")+
#     theme_bw()+
#     theme(legend.position="top",
#           axis.text.x=element_text(angle=90,size=6),
#           axis.text.y=element_text(size=5),
#           axis.title.x=element_blank(),
#           axis.title.y=element_blank())
#   #scale_y_discrete(labels=function(x){sub("\\s", "\n", x)})
#   
#   return(p)
# }

# function for making top n species plots (per enzymeID) 
make_species_CAZ_plots <- function(x) {
  
  p <- x %>% 
    group_by(enzymeID) %>%
    top_n(n = 3, wt = perc_species) %>%
    arrange(desc(perc_species)) %>%
    slice(1:3) %>%
    drop.levels() %>% 
    ggplot(aes(x = reorder(enzymeID, perc_species, FUN = mean), y = reorder(species, perc_species, FUN = mean), fill = perc_species)) +
    geom_tile() +
    #scale_fill_gaio8()+
    scale_fill_distiller(type = "div", palette = "Spectral") +
    geom_text(aes(y=species,label=round(perc_species)), size=2)+
    facet_wrap(~enzymeID, scales = "free",ncol = 3)+
    theme(axis.title.x=element_blank(),axis.text.x = element_blank(),
          axis.text.y=element_text(size=5),
          axis.title.y=element_blank(),
          strip.text.x = element_text(size = 7, colour = "black"),
          legend.key.height = unit(.5, "cm"),
          legend.key.width = unit(.5, "cm"),
          legend.key.size = unit(.5, "cm"),
          legend.text.align = 1,
          legend.title = element_text(size=7),
          legend.text = element_text(size=6),
          legend.position="bottom",
          legend.justification="bottom",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-6,-6,-6,-6), # get closer farther from zero to get it closer to edge 
          axis.ticks.x = element_blank())+
    labs(fill="percentage of species carrying enzyme")+
    scale_y_discrete(labels=function(x){sub("\\s", "\n", x)})
  return(p)
}


# function for making top n genus plots (per enzymeID) 
make_genus_CAZ_plots <- function(x) {
  
  p <- x %>% 
    group_by(enzymeID) %>%
    top_n(n = 3, wt = perc_genus) %>%
    arrange(desc(perc_genus)) %>%
    slice(1:3) %>%
    drop.levels() %>% 
    ggplot(aes(x = reorder(enzymeID, perc_genus, FUN = mean), y = reorder(genus, perc_genus, FUN = mean), fill = perc_genus)) +
    geom_tile() +
    #scale_fill_gaio8()+
    scale_fill_distiller(type = "div", palette = "Spectral") +
    geom_text(aes(y=genus,label=round(perc_genus)), size=2)+
    facet_wrap(~enzymeID, scales = "free",ncol = 3)+
    theme(axis.title.x=element_blank(),axis.text.x = element_blank(),
          axis.text.y=element_text(size=5),
          axis.title.y=element_blank(),
          strip.text.x = element_text(size = 7, colour = "black"),
          legend.key.height = unit(.5, "cm"),
          legend.key.width = unit(.5, "cm"),
          legend.key.size = unit(.5, "cm"),
          legend.text.align = 1,
          legend.title = element_text(size=7),
          legend.text = element_text(size=6),
          legend.position="bottom",
          legend.justification="bottom",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-6,-6,-6,-6), # get closer farther from zero to get it closer to edge 
          axis.ticks.x = element_blank())+
    labs(fill="percentage of genera carrying enzyme")+
    scale_y_discrete(labels=function(x){sub("\\s", "\n", x)})
  return(p)
}


# function for making CAZ time trend boxplots : 
# these (different from previous) are adapted to be plotted togetehr with species distribution 
make_enzyme_boxplots <- function(x) {
  
  p <- x %>% 
    filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="t10"|date=="tM") %>%
    group_by(pig, enzymeID, enzymeNAME,date,n,n_moms) %>% 
    dplyr::summarise(tot = mean(tot, na.rm = TRUE)) %>% 
    mutate(tot=log(tot)) %>%
    arrange(desc(enzymeID)) %>%
    drop.levels() %>%
    ggplot(., aes(x = date, y = tot, fill = enzymeNAME)) +
    geom_boxplot(lwd=0.1, outlier.size = 0.2)+
    scale_fill_gaio8()+
    #scale_fill_distiller(type = "div", palette = "Spectral") + # this was for the heatmap 
    facet_wrap(~enzymeID, scales = "free_y", ncol = 3) +
    labs(x = "date", y = "enzymeID", fill = "avg. abundance (log)")+
    theme_bw()+
    theme(legend.position="none",
          axis.text.y=element_text(size = 4),
          axis.ticks.length.y = unit(.05, "cm"),
          axis.text.x=element_text(size=6,angle=90))+
    geom_text(aes(x="ss_piglets",y=Inf, label=paste0("n=",n)),
              size=2.3,colour="black", hjust = 1.5, angle=90, inherit.aes=TRUE, parse=FALSE,check_overlap = TRUE)+
    geom_text(aes(x="ss_sows",y=Inf, label=paste0("n=",n_moms)),
              size=2.3,colour="black", hjust = 1.5, angle=90, inherit.aes=TRUE, parse=FALSE,check_overlap = TRUE)
  
  return(p)
  
}

a_AA_CE <- make_enzyme_boxplots(df_part_AA_CE)
a_PL <- make_enzyme_boxplots(df_part_PL)
a_GT_1 <- make_enzyme_boxplots(df_part1_GT)
a_GT_2 <- make_enzyme_boxplots(df_part2_GT)
a_GT_3 <- make_enzyme_boxplots(df_part3_GT)
a_GH_1 <- make_enzyme_boxplots(df_part1_GH)
a_GH_2 <- make_enzyme_boxplots(df_part2_GH)
a_GH_3 <- make_enzyme_boxplots(df_part3_GH)
a_GH_4 <- make_enzyme_boxplots(df_part4_GH)
a_CBM_1 <- make_enzyme_boxplots(df_part1_CBM)
a_CBM_2 <- make_enzyme_boxplots(df_part2_CBM)
a_SLH <- make_enzyme_boxplots(df_part_SLH)
a_cohesin <- make_enzyme_boxplots(df_part_cohesin)
a_SLH_cohesin <- make_enzyme_boxplots(rbind(df_part_SLH,df_part_cohesin))

b_AA_CE <- make_species_CAZ_plots(s_part_AA_CE)
b_PL <- make_species_CAZ_plots(s_part_PL)
b_GT_1 <- make_species_CAZ_plots(s_part1_GT)
b_GT_2 <- make_species_CAZ_plots(s_part2_GT)
b_GT_3 <- make_species_CAZ_plots(s_part3_GT)
b_GH_1 <- make_species_CAZ_plots(s_part1_GH)
b_GH_2 <- make_species_CAZ_plots(s_part2_GH)
b_GH_3 <- make_species_CAZ_plots(s_part3_GH)
b_GH_4 <- make_species_CAZ_plots(s_part4_GH)
b_CBM_1 <- make_species_CAZ_plots(s_part1_CBM)
b_CBM_2 <- make_species_CAZ_plots(s_part2_CBM)
b_SLH <- make_species_CAZ_plots(s_part_SLH)
b_cohesin <- make_species_CAZ_plots(s_part_cohesin)
b_SLH_cohesin <- make_species_CAZ_plots(rbind(s_part_SLH,s_part_cohesin))

b2_AA_CE <- make_genus_CAZ_plots(s_genus_part_AA_CE)
b2_PL <- make_genus_CAZ_plots(s_genus_part_PL)
b2_GT_1 <- make_genus_CAZ_plots(s_genus_part1_GT)
b2_GT_2 <- make_genus_CAZ_plots(s_genus_part2_GT)
b2_GT_3 <- make_genus_CAZ_plots(s_genus_part3_GT)
b2_GH_1 <- make_genus_CAZ_plots(s_genus_part1_GH)
b2_GH_2 <- make_genus_CAZ_plots(s_genus_part2_GH)
b2_GH_3 <- make_genus_CAZ_plots(s_genus_part3_GH)
b2_GH_4 <- make_genus_CAZ_plots(s_genus_part4_GH)
b2_CBM_1 <- make_genus_CAZ_plots(s_genus_part1_CBM)
b2_CBM_2 <- make_genus_CAZ_plots(s_genus_part2_CBM)
b2_SLH <- make_genus_CAZ_plots(s_genus_part_SLH)
b2_cohesin <- make_genus_CAZ_plots(s_genus_part_cohesin)
b2_SLH_cohesin <- make_genus_CAZ_plots(rbind(s_genus_part_SLH,
                                         s_genus_part_cohesin))

# h_a_AA_CE <- make_enzyme_heatmap(df_part_AA_CE)
# h_a_PL <- make_enzyme_heatmap(df_part_PL)
# h_a_GT_1 <- make_enzyme_heatmap(df_part1_GT)
# h_a_GT_2 <- make_enzyme_heatmap(df_part2_GT)
# h_a_GT_3 <- make_enzyme_heatmap(df_part3_GT)
# h_a_GH_1 <- make_enzyme_heatmap(df_part1_GH)
# h_a_GH_2 <- make_enzyme_heatmap(df_part2_GH)
# h_a_GH_3 <- make_enzyme_heatmap(df_part3_GH)
# h_a_GH_4 <- make_enzyme_heatmap(df_part4_GH)
# h_a_CBM_1 <- make_enzyme_heatmap(df_part1_CBM)
# h_a_CBM_2 <- make_enzyme_heatmap(df_part2_CBM)
# h_a_SLH <- make_enzyme_heatmap(df_part_SLH)
# h_a_cohesin <- make_enzyme_heatmap(df_part_cohesin)



####
# combining the plots: CAZ time trend boxplots with species per enzymeID info: 
AA_CE <- plot_grid(a_AA_CE,b_AA_CE,ncol=2)

PL <- plot_grid(a_PL,b_PL,ncol=2)

CBM1 <- plot_grid(a_CBM_1,b_CBM_1,ncol=2)
CBM2 <- plot_grid(a_CBM_2,b_CBM_2,ncol=2)

GT1 <- plot_grid(a_GT_1,b_GT_1,ncol=2)
GT2 <- plot_grid(a_GT_2,b_GT_2,ncol=2)
GT3 <- plot_grid(a_GT_3,b_GT_3,ncol=2)

GH1 <- plot_grid(a_GH_1,b_GH_1,ncol=2)
GH2 <- plot_grid(a_GH_2,b_GH_2,ncol=2)
GH3 <- plot_grid(a_GH_3,b_GH_3,ncol=2)
GH4 <- plot_grid(a_GH_4,b_GH_4,ncol=2)

SLH <- plot_grid(a_SLH,b_SLH,ncol=2)

cohesin <- plot_grid(a_cohesin,b_cohesin,ncol=2)

SLH_cohesin <- plot_grid(a_SLH_cohesin,
                         b_SLH_cohesin,
                         ncol=2)
####

####
# combining the plots: CAZ time trend boxplots with genus per enzymeID info: 
g_AA_CE <- plot_grid(a_AA_CE,b2_AA_CE,ncol=2)

g_PL <- plot_grid(a_PL,b2_PL,ncol=2)

g_CBM1 <- plot_grid(a_CBM_1,b2_CBM_1,ncol=2)
g_CBM2 <- plot_grid(a_CBM_2,b2_CBM_2,ncol=2)

g_GT1 <- plot_grid(a_GT_1,b2_GT_1,ncol=2)
g_GT2 <- plot_grid(a_GT_2,b2_GT_2,ncol=2)
g_GT3 <- plot_grid(a_GT_3,b2_GT_3,ncol=2)

g_GH1 <- plot_grid(a_GH_1,b2_GH_1,ncol=2)
g_GH2 <- plot_grid(a_GH_2,b2_GH_2,ncol=2)
g_GH3 <- plot_grid(a_GH_3,b2_GH_3,ncol=2)
g_GH4 <- plot_grid(a_GH_4,b2_GH_4,ncol=2)

g_SLH <- plot_grid(a_SLH,b2_SLH,ncol=2)

g_cohesin <- plot_grid(a_cohesin,b2_cohesin,ncol=2)

g_SLH_cohesin <- plot_grid(a_SLH_cohesin,
                         b2_SLH_cohesin,
                         ncol=2)
####


####
# Two heatmaps per page: one for the CAZ time trend, one for the species corresponding to the bin where the CAZ was found 
# sh_AA_CE <- make_species_heatmap(s_part_AA_CE,h_a_AA_CE)
# H1 <- plot_grid(h_a_AA_CE,sh_AA_CE,ncol=2)
# 
# sh_PL <- make_species_heatmap(s_part_PL,h_a_PL)
# H2 <- plot_grid(h_a_PL,sh_PL,ncol=2)
# 
# sh_GT1 <- make_species_heatmap(s_part1_GT,h_a_GT_1)
# H3 <- plot_grid(h_a_GT_1,sh_GT1,ncol=2)
# 
# sh_GT2 <- make_species_heatmap(s_part2_GT,h_a_GT_2)
# H4 <- plot_grid(h_a_GT_2,sh_GT2,ncol=2)
# 
# sh_GT3 <- make_species_heatmap(s_part3_GT,h_a_GT_3)
# H5 <- plot_grid(h_a_GT_3,sh_GT3,ncol=2)
# 
# sh_GH1 <- make_species_heatmap(s_part1_GH,h_a_GH_1)
# H6 <- plot_grid(h_a_GH_1,sh_GH1,ncol=2)
# 
# sh_GH2 <- make_species_heatmap(s_part2_GH,h_a_GH_2)
# H7 <- plot_grid(h_a_GH_2,sh_GH2,ncol=2)
# 
# sh_GH3 <- make_species_heatmap(s_part3_GH,h_a_GH_3)
# H8 <- plot_grid(h_a_GH_3,sh_GH3,ncol=2)
# 
# sh_CBM1 <- make_species_heatmap(s_part1_CBM,h_a_CBM_1)
# H9 <- plot_grid(h_a_CBM_1,sh_CBM1,ncol=2)
# 
# sh_CBM2 <- make_species_heatmap(s_part2_CBM,h_a_CBM_2)
# H10 <- plot_grid(h_a_CBM_2,sh_CBM2,ncol=2)
# 
# sh_SLH <- make_species_heatmap(s_part_SLH,h_a_SLH)
# H11 <- plot_grid(h_a_SLH,sh_SLH,ncol=2)
# 
# sh_cohesin <- make_species_heatmap(s_part_cohesin,h_a_cohesin)
# H12 <- plot_grid(h_a_cohesin,sh_cohesin,ncol=2)
####




# pdf("dbcan_HMMER_CAZ_time_heatmaps.pdf")
# h_a_AA_CE 
# h_a_PL 
# h_a_GT_1 
# h_a_GT_2 
# h_a_GT_3 
# h_a_GH_1 
# h_a_GH_2 
# h_a_GH_3 
# h_a_CBM_1 
# h_a_CBM_2 
# h_a_SLH
# h_a_cohesin
# dev.off()
# 
# pdf("dbcan_CAZ_ALL_heatmaps.pdf")
# H1
# H2
# H3
# H4
# H5
# H6
# H7
# H8
# H9
# H10
# H11
# H12
# dev.off()

pdf("dbcan_HMMER_CAZ_time_species.pdf")
AA_CE
PL
CBM1
CBM2
GT1
GT2
GT3
GH1
GH2
GH3
GH4
#SLH
#cohesin
SLH_cohesin
dev.off()

pdf("dbcan_HMMER_CAZ_time_genus.pdf")
g_AA_CE
g_PL
g_CBM1
g_CBM2
g_GT1
g_GT2
g_GT3
g_GH1
g_GH2
g_GH3
g_GH4
#g_SLH
#g_cohesin
g_SLH_cohesin
dev.off()






# PLOT showing all enzymes that significantly changed over time, 
# with popularity (y) and "singularity" (x) meaning how uniquely an enzyme is expressed by one (right) or multiple (left) taxa

CAZy_species_representation <- gt_hmmer %>% 
  group_by(enzymeID,species) %>%
  tally() %>%
  group_by(enzymeID) %>%
  dplyr::mutate(perc_spec=round(n/sum(n)*100,2)) %>%
  drop_na() %>%
  group_by(enzymeID) %>%
  top_n(n=1,wt=perc_spec) %>%
  arrange(enzymeID) %>%
  dplyr::mutate(top_species = perc_spec) %>%
  dplyr::select(enzymeID, species, top_species)

head(CAZy_species_representation)
tail(CAZy_species_representation)

CAZy_prevalence_in_population <- gt_hmmer %>%
  filter(!pig=="Protexin"&!pig=="ColiGuard"&!pig=="MockCommunity"&!pig=="NegativeControl") %>%
  dplyr::select(pig,enzymeID, enzymeNAME) %>%
  distinct() %>%
  group_by(enzymeID, enzymeNAME) %>%
  tally() %>%
  dplyr::mutate(popularity = n) %>%
  dplyr::select(enzymeID, enzymeNAME, popularity)

head(CAZy_prevalence_in_population)
tail(CAZy_prevalence_in_population)

df <- full_join(CAZy_species_representation, CAZy_prevalence_in_population) %>% distinct()

pdf("dbcan_HMMER_all_CAZ_plot.pdf")
ggplot(df) +
  geom_point(aes(top_species, popularity, fill=enzymeNAME), 
             colour="black",pch=21, size=4)+
  labs(x="percentage of top species",
       y="prevalence in the pig population",
       fill="enzyme class")+
  scale_fill_gaio8()+
  theme_bw(base_size = 10) +
  xlim(0,130) +
  ggrepel::geom_label_repel(data = subset(df, top_species > 55 & popularity > 5),
                            aes(
                              x = top_species,
                              y = popularity,
                              fill = factor(enzymeNAME),
                              label = paste0(enzymeID,"\n",species)
                            ),
                            box.padding   = 0.3,label.padding = 0.1,
                            point.padding = 0.5,alpha=0.7,
                            force         = 10,
                            segment.size  = 0.2,
                            segment.color = "grey50", nudge_x = 30,
                            direction     = "y", nudge_y = 10,
                            size=2.5,
                            color="black"
  ) 
dev.off()



# EDIT THIS BELOW WITH NEW DF ("BOTH")


# CE

sink(file = "dbcan_HMMER_numbers.text",append = TRUE)

paste0(" ##############################  CE  ############################## ")
paste0("number of all enzymeIDs within this class")
NROW(subset(df, enzymeNAME == "CE")) 
paste0("number of all enzymeIDs within this class, that are present in multiple species")
NROW(subset(df, enzymeNAME == "CE" & top_species<25)) 
CE_df <- df %>%
  filter(enzymeNAME=="CE") %>%
  filter(top_species<25) 
summary(CE_df)

paste0(" ##############################  PL  ############################## ")
paste0("number of all enzymeIDs within this class")
NROW(subset(df, enzymeNAME == "PL")) 
paste0("number of all enzymeIDs within this class, that are present in multiple species")
NROW(subset(df, enzymeNAME == "PL" & top_species<25)) 
PL_df <- df %>%
  filter(enzymeNAME=="PL") %>%
  filter(top_species<25) 
summary(PL_df)

paste0(" ##############################  GT  ############################## ")
paste0("number of all enzymeIDs within this class")
NROW(subset(df, enzymeNAME == "GT")) 
paste0("number of all enzymeIDs within this class, that are present in multiple species")
NROW(subset(df, enzymeNAME == "GT" & top_species<25 & popularity >50)) 
34*100/68
GT_df <- df %>%
  filter(enzymeNAME=="GT") %>%
  filter(top_species<25) 
summary(GT_df)

paste0(" ##############################  GH  ############################## ")
paste0("number of all enzymeIDs within this class")
NROW(subset(df, enzymeNAME == "GH")) 
paste0("number of all enzymeIDs within this class, that are present in multiple species")
NROW(subset(df, enzymeNAME == "GH" & top_species<25 & popularity >50)) 
94*100/128
GH_df <- df %>%
  filter(enzymeNAME=="GH") %>%
  filter(top_species<25) 
summary(GH_df)

paste0(" ##############################  CBM  ############################## ")
paste0("number of all enzymeIDs within this class")
NROW(subset(df, enzymeNAME == "CBM")) 
CBM_df <- df %>%
  filter(enzymeNAME=="CBM") 
summary(CBM_df)

paste0(" ##############################  AA  ############################## ")
paste0("number of all enzymeIDs within this class")
NROW(subset(df, enzymeNAME == "AA")) 
AA_df <- df %>%
  filter(enzymeNAME=="AA") 
AA_df
summary(AA_df)

paste0(" ##############################  cohesin  ############################## ")
paste0("number of all enzymeIDs within this class")
NROW(subset(df, enzymeNAME == "cohesin")) 
cohesin_df <- df %>%
  filter(enzymeNAME=="cohesin") 
cohesin_df



paste0(" ##############################  SLH  ############################## ")
paste0("number of all enzymeIDs within this class")
NROW(subset(df, enzymeNAME == "SLH")) 
SLH_df <- df %>%
  filter(enzymeNAME=="SLH") 
SLH_df


paste0(" ##############################  enzymes present majorly in one species ############################## ")
rightmost <- df %>%
  filter(top_species>70) %>%
  filter(popularity>15) %>%
  arrange(desc(top_species))
NROW(rightmost)
rightmost$species
rightmost$enzymeID
sink()



# Prevotella vs Bacteroidetes enzyme specificity

test <- subset(df_part, (enzymeID %in% mylist))


# 1 
# plot all prevotella and all bacteroidetes genus, all enzymes, no labels 


# 2 
# plot only Bacteroidetes A and Prevotella, 50 most abundant enzymes; no labels 

# 3 
# plot only Bacteroidetes A and Prevotella, 60th-100th most abundant enzymes; labels 




gt_hmmer_sub <- gt_hmmer %>% 
  filter(genus=="Prevotella"|genus=="Bacteroides_A"|genus=="Bacteroides"|genus=="Bacteroides_B") 

# reorder  
gt_hmmer_sub$genus  = factor(gt_hmmer_sub$genus, levels=c("Bacteroides_A",
                                                              "Bacteroides_B", 
                                                              "Bacteroides",
                                                              "Prevotella"))


# all Bact. and Prev
test2 <- gt_hmmer_sub %>% 
  filter(genus=="Prevotella"|genus=="Bacteroides_A"|genus=="Bacteroides"|genus=="Bacteroides_B") %>%
  dplyr::select(genus,enzymeID) %>%
  distinct()
test4 <- inner_join(test,test2, by="enzymeID")
NROW(test4)
test4$sample=paste0(".",test4$genus,"__",test4$sample)
df3 <- test4 %>% dplyr::select(sample,enzymeID,tot) %>% 
  pivot_wider(id_cols = sample, names_from = enzymeID, 
              values_from=tot, values_fill = list(tot = 0))
x <- as.data.frame(df3)
rownames(x) <- x$sample
x$sample <- NULL
# order left to right in descending order 
x <- x[,names(sort(colSums(x), decreasing = TRUE))]
#############
# PCA
mtcars.pca2 <- prcomp(x, center = TRUE,scale. = TRUE)   # [,1:100]
genus <- as.character(qdapRegex::ex_between(rownames(x),".", "__"))
dates <- as.character(qdapRegex::ex_between(rownames(x), "__", "_"))
p1 <- fviz_pca_ind(mtcars.pca2, 
                   geom.ind="point",
                   #fill.ind = dates, 
                   #col.ind = rainbow(n = 11),
                   pointsize = 2, 
                   habillage = genus, 
                   pointshape=21,
                   #geom.ind = "point", # show points only (nbut not "text") 
                   col.ind = genus, # color by groups
                   palette = c("#8000FFFF", "#FF0000FF","#00FFFFFF","#80FF00FF"),
                   addEllipses = FALSE, # Concentration ellipses
                   title="")+
  theme(legend.position="top",
        legend.title = element_blank())+
  guides(color = guide_legend(nrow = 1))

myleg <- get_legend(p1)

p1 <- p1 + theme(legend.position="none")


# only Bact.A and Prev
test2 <- gt_hmmer_sub %>% 
  filter(genus=="Prevotella"|genus=="Bacteroides_A") %>%
  dplyr::select(genus,enzymeID) %>%
  distinct()
test4 <- inner_join(test,test2, by="enzymeID")
NROW(test4)
test4$sample=paste0(".",test4$genus,"__",test4$sample)
df3 <- test4 %>% dplyr::select(sample,enzymeID,tot) %>% 
  pivot_wider(id_cols = sample, names_from = enzymeID, 
              values_from=tot, values_fill = list(tot = 0))
x <- as.data.frame(df3)
rownames(x) <- x$sample
x$sample <- NULL
# order left to right in descending order 
x <- x[,names(sort(colSums(x), decreasing = TRUE))]

# first 60
mtcars.pca2 <- prcomp(x[,1:60], center = TRUE,scale. = TRUE)   
genus <- as.character(qdapRegex::ex_between(rownames(x[,1:60]),".", "__"))
dates <- as.character(qdapRegex::ex_between(rownames(x[,1:60]), "__", "_"))
p2 <- fviz_pca_ind(mtcars.pca2, 
                   geom.ind="point",
                   #fill.ind = dates, 
                   #col.ind = rainbow(n = 11),
                   pointsize = 2, 
                   habillage = genus, 
                   pointshape=21,
                   #geom.ind = "point", # show points only (nbut not "text") 
                   col.ind = genus, # color by groups
                   palette = c("#FF0000FF","#80FF00FF"),
                   addEllipses = FALSE, # Concentration ellipses
                   title="")+
  theme(legend.position="none")+
  guides(color = guide_legend(nrow = 1))


# 60-100
mtcars.pca2 <- prcomp(x[,60:100], center = TRUE,scale. = TRUE)   
genus <- as.character(qdapRegex::ex_between(rownames(x[,60:100]),".", "__"))
dates <- as.character(qdapRegex::ex_between(rownames(x[,60:100]), "__", "_"))
p3 <- fviz_pca_biplot(mtcars.pca2,
                      geom.ind="point",
                      pointsize = 2, label = c("var"), #fill.ind = dates,
                      habillage = genus, 
                      pointshape=21, 
                      col.ind = genus, # select.var = list(contrib = 20), # selecting 20 most contributing vars
                      alpha.var ="contrib", labelsize=4,
                      palette = c("#FF0000FF","#80FF00FF"),
                      repel = TRUE,
                      title="") +
  ggtitle("")+
  theme(legend.position="none") # panel.border = element_rect(colour = "black", fill=NA, size=1)




p1 <- ggarrange(p1,labels = "A")
p23 <- ggarrange(p2,p3, 
                 widths = c(1,2), labels=c("B","C"))
tosave <- ggarrange(p1,p23, widths = c(1,3))
tosave <- ggarrange(myleg,tosave,heights=c(1,9),nrow=2)

pdf("dbcan_HMMER_specificity_genera.pdf")
tosave
dev.off()


##################################

# ENZYMES CO-OCCURRENCE in the same GENUS - piglets 

z <- gt_hmmer %>% dplyr::select(enzymeID,pig,genus)
piggiesIDs <- no_reps_all %>% filter(!cohort=="Mothers") %>% dplyr::select(pig) %>% distinct()

# collect mothers enzyme info only 
z_piggies <- left_join(piggiesIDs,z)

z1 <- z_piggies %>% 
  group_by(genus,enzymeID) %>% 
  tally() 

z2 <- z_piggies %>% group_by(enzymeID) %>% tally() 

# set limits 
firstQu <- as.numeric(summary(z2$n)[2])
median <- as.numeric(summary(z2$n)[3])
thirdQu <- as.numeric(summary(z2$n)[5])

# PLOT 
pdf("dbcan_HMMER_occurr_piglets_same_genus.pdf")
# below 1st Qu.
z3 <- subset(z2,n < firstQu)
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
        main="Co-occurrence - Frequency: below 1st Qu.") 
# 1st Qu. - median
z3 <- subset(z2, n >= firstQu & n < median)
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
        main="Co-occurrence - Frequency: 1st Qu. - median") 
# median  - 3rd Qu. 
z3 <- subset(z2,n >= median & n < thirdQu)
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
        main="Co-occurrence - Frequency: median - 3rd Qu.") 
# above 3rd Qu. 
z3 <- subset(z2,n >= thirdQu)
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
        main="Co-occurrence - Frequency: above 3rd Qu.") 
# without any subsetting
z3 <- z2
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none",cexRow = 0.06, cexCol=0.06,
        main="Co-occurrence - all") 
dev.off()


##################################

# ENZYMES CO-OCCURRENCE in the same GENUS - mothers 

z <- gt_hmmer %>% dplyr::select(enzymeID,pig,genus)

mothersIDs <- no_reps_all %>% filter(cohort=="Mothers") %>% dplyr::select(pig) %>% distinct()
# collect mothers enzyme info only 
z_moms <- left_join(mothersIDs,z)

z1 <- z_moms %>% 
  group_by(genus,enzymeID) %>% 
  tally() 

z2 <- z_moms %>% group_by(enzymeID) %>% tally() 

# set limits 
firstQu <- as.numeric(summary(z2$n)[2])
median <- as.numeric(summary(z2$n)[3])
thirdQu <- as.numeric(summary(z2$n)[5])

# PLOT 
pdf("dbcan_HMMER_occurr_mothers_same_genus.pdf")
# below 1st Qu.
z3 <- subset(z2,n < firstQu)
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
        main="Co-occurrence - Frequency: below 1st Qu.") 
# 1st Qu. - median
z3 <- subset(z2, n >= firstQu & n < median)
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
        main="Co-occurrence - Frequency: 1st Qu. - median") 
# median  - 3rd Qu. 
z3 <- subset(z2,n >= median & n < thirdQu)
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
        main="Co-occurrence - Frequency: median - 3rd Qu.") 
# above 3rd Qu. 
z3 <- subset(z2,n >= thirdQu)
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
        main="Co-occurrence - Frequency: above 3rd Qu.") 
# without any subsetting
z3 <- z2
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none",cexRow = 0.06, cexCol=0.06,
        main="Co-occurrence - all")
dev.off()

##################################
##################################

# ENZYMES CO-OCCURRENCE in the same GENUS - mothers AND piglets 

z <- gt_hmmer %>% dplyr::select(enzymeID,pig,genus)

IDs <- no_reps_all %>% dplyr::select(pig) %>% distinct()
# collect mothers enzyme info only 
z <- left_join(IDs,z)

z1 <- z %>% 
  group_by(genus,enzymeID) %>% 
  tally() 

z2 <- z %>% group_by(enzymeID) %>% tally() 

# PLOT 
pdf("dbcan_HMMER_occurr_moms&piglets_same_genus.pdf")
# without any subsetting
z3 <- z2
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none",cexRow = 0.06, cexCol=0.06,
        main="Co-occurrence of CAZy enzymes in the same genus")
dev.off()


##################################
##################################

# ENZYMES CO-OCCURRENCE in the same SPECIES - mothers AND piglets 

z <- gt_hmmer %>% dplyr::select(enzymeID,pig,species)

IDs <- no_reps_all %>% dplyr::select(pig) %>% distinct()
# collect mothers enzyme info only 
z <- left_join(IDs,z)

z1 <- z %>% 
  group_by(species,enzymeID) %>% 
  tally() 

z2 <- z %>% group_by(enzymeID) %>% tally() 

# PLOT 
pdf("dbcan_HMMER_occurr_moms&piglets_same_species.pdf")
# without any subsetting
z3 <- z2
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none",cexRow = 0.06, cexCol=0.06,
        main="Co-occurrence of CAZy enzymes in the same species")
dev.off()

##################################
##################################


# ENZYMES CO-OCCURRENCE in the same subject 

# create an adjancy matrix: rows are species, columns are enzyme IDs
z <- gt_hmmer %>% dplyr::select(enzymeID,pig)
z

z1 <- z %>% 
  group_by(pig,enzymeID) %>% 
  tally() 

z2 <- z %>% group_by(enzymeID) %>% tally() 

# set limits 
firstQu <- as.numeric(summary(z2$n)[2])
median <- as.numeric(summary(z2$n)[3])
thirdQu <- as.numeric(summary(z2$n)[5])

# PLOT 
pdf("dbcan_HMMER_occurr_same_host.pdf")
# below 1st Qu.
z3 <- subset(z2,n < firstQu)
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
        main="Co-occurrence - Frequency: below 1st Qu.") 
# 1st Qu. - median
z3 <- subset(z2, n >= firstQu & n < median)
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
        main="Co-occurrence - Frequency: 1st Qu. - median") 
# median  - 3rd Qu. 
z3 <- subset(z2,n >= median & n < thirdQu)
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
        main="Co-occurrence - Frequency: median - 3rd Qu.") 
# above 3rd Qu. 
z3 <- subset(z2,n >= thirdQu)
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
        main="Co-occurrence - Frequency: above 3rd Qu.") 
# without any subsetting
z3 <- z2
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none",cexRow = 0.06, cexCol=0.06,
        main="Co-occurrence of CAZy enzymes in the same host") 
dev.off()
