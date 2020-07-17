

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


setwd("~/Desktop/metapigs_dry/dbcan")
basedir = "~/Desktop/metapigs_dry/"

########################

# load dbcan_diamond data (output of mapping our-bins-AAs against CAZy db)
diamond <- read_table2("diamond.out", col_names = FALSE)

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



# DIAMOND data cleaning 

# get the enzyme ID clean
diamond <- cSplit(diamond,"X2","|")
diamond <- cSplit(diamond,"X2_2","_")
diamond <- cSplit(diamond,"X2_2_1"," ")

# subject/bin/contig
diamond <- cSplit(diamond,"X1","_")

diamond <- diamond %>%
  dplyr::select(X1_1,X1_2,X1_3,
                X2_1,X2_2_1_1,
                X3,X4,X5,X6,X7,X8,X9,X10,X11,X12)

# save column with full enzyme ID 
diamond$enzyme <- diamond$X2_2_1_1

# split the enzyme ID from the enzyme number 
diamond <- diamond %>%
  separate(X2_2_1_1, 
           into = c("text", "num"), 
           sep = "(?<=[A-Za-z])(?=[0-9])")

# cols renaming
colnames(diamond) <- c("pig","bin","contig","CAZ_accession","enzymeNAME","enzymeDIGIT",
                       "perc_identity","align_length","mismatches","gap_openings",
                       "start_query","end_query","start_sample","end_sample",
                       "evalue","bitscore","enzymeID")
diamond$bin <- gsub(".fa","", diamond$bin)

head(diamond)

########################################################################

# QUALITY of diamond data : 
# how well our predicted proteins - matched against the CAZy database 

# our evalues 
summary(diamond$evalue)

# our percentages of identity 
summary(diamond$perc_identity)

# re-rder enzymes to follow the same aesthetics as Stewart et al 2019
diamond$enzymeNAME  = factor(
  diamond$enzymeNAME, levels=c("AA",
                               "CE",
                               "GH",
                               "GT",
                               "CBM",
                               "PL"))

diamond$enzymeNAME <- as.character(diamond$enzymeNAME)

# PLOT 
diamond_perc_identity_plot <- diamond %>%
  group_by(enzymeNAME) %>%
  ggplot(.,aes(enzymeNAME,perc_identity,fill=enzymeNAME))+
  ylab("Percentage identity")+
  geom_boxplot()+
  ylim(0,100)+
  scale_fill_brewer(palette="Spectral")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  ggtitle("Percentage identity against CAZy")

enzymes_proportion <- diamond %>% 
  group_by(enzymeNAME) %>% 
  tally() %>%
  mutate(perc = paste0(round(n/sum(n)*100,2),"%"))

#RColorBrewer::display.brewer.all(n=6)

diamond_perc_identity_plot <- diamond_perc_identity_plot + 
  geom_text(data = enzymes_proportion,
            aes(enzymeNAME, 0, label = perc), vjust="inward",size=3)

pdf("dbcan_perc_identity.pdf")
diamond_perc_identity_plot
dev.off()

########################################################################

# clean checkM all nearly data:

# filter >90 <5 bins 
checkm_all_nearly <- dplyr::filter(checkm_all_nearly, !grepl("Completeness",Completeness))
checkm_all_nearly <- subset(checkm_all_nearly, Completeness >= 90 & Contamination <= 5)
# some formatting 
checkm_all_nearly$Completeness <- as.numeric(checkm_all_nearly$Completeness)
checkm_all_nearly$Contamination <- as.numeric(checkm_all_nearly$Contamination)
# rename cols to matching colnames between dataframes to merge  
colnames(checkm_all_nearly)[colnames(checkm_all_nearly)=="pigid"] <- "pig"
colnames(checkm_all_nearly)[colnames(checkm_all_nearly)=="Bin Id"] <- "bin"
colnames(checkm_all_nearly)[colnames(checkm_all_nearly)=="Taxonomy (contained)"] <- "taxa"

checkm_clean <- checkm_all_nearly %>%
  select(pig,bin,taxa)
checkm_clean <- cSplit(checkm_clean, "taxa", sep=";")
colnames(checkm_clean) <- c("pig","bin","kingdom","phylum","class","order","family","genus","species")
head(checkm_clean)

checkm_clean$phylum <- as.character(checkm_clean$phylum)
checkm_clean$phylum[is.na(checkm_clean$phylum)] <- "Unknown"

checkm_clean$phylum <- gsub("p__","", checkm_clean$phylum)

########################################################################

# merge diamond data with checkM- nearly complete bins - taxonomic assignments 


NROW(checkm_clean)
NROW(diamond)

cm_diamond <- inner_join(checkm_clean,diamond)
NROW(cm_diamond)


cm_diamond_taxa <- cm_diamond %>%
  group_by(phylum,enzymeNAME) %>%
  tally() %>%
  mutate(Proportion = n) 

cm_diamond_prop <- cm_diamond %>%
  group_by(phylum) %>%
  tally() %>%
  mutate(proportion=n/sum(n)*100)

a <- ggplot(cm_diamond_prop, aes(y=proportion, x=phylum)) + 
  geom_bar(position="dodge", stat="identity",colour="black",fill="snow3")+ #lightskyblue happier
  theme(axis.text.x=element_text(angle=90))+
  theme_pubr()+
  ylab("Proportion of CAZymes (%)")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
b <- ggplot(cm_diamond_taxa, aes(fill=enzymeNAME, y=Proportion, x=phylum)) + 
  geom_bar(position="fill", stat="identity",colour="black")+
  theme(axis.text.x=element_text(angle=90))+
  scale_fill_brewer(palette="Spectral")+
  theme_pubr()+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x = element_blank(),
        legend.position="right",
        legend.title = element_blank())

# RColorBrewer::display.brewer.all(n=6)

CAZ_CMphyla_plot <- plot_grid(
  plot_grid(
    a + theme(legend.position = "none"),
    b + theme(legend.position = "none"),
    ncol = 1,
    align = "v"),
  plot_grid(
    get_legend(a),
    get_legend(b),
    ncol =1),
  rel_widths = c(8,2)
)

pdf("dbcan_CAZ_CMphyla.pdf")
CAZ_CMphyla_plot
dev.off()


########################################################################

# merge diamond data with gtdbtk- bins taxonomic assignments 


NROW(gtdbtk_bins)
NROW(diamond)

gt_diamond <- inner_join(gtdbtk_bins,diamond)
NROW(gt_diamond)

gt_diamond$phylum[is.na(gt_diamond$phylum)] <- "Unknown"

gt_diamond_taxa <- gt_diamond %>%
  group_by(phylum,enzymeNAME) %>%
  tally() %>%
  mutate(Proportion = n) 

gt_diamond_prop <- gt_diamond %>%
  group_by(phylum) %>%
  tally() %>%
  mutate(proportion=n/sum(n)*100)

a_gt <- ggplot(gt_diamond_prop, aes(y=proportion, x=phylum)) + 
  geom_bar(position="dodge", stat="identity",colour="black",fill="snow3")+ #lightskyblue happier
  theme(axis.text.x=element_text(angle=90))+
  theme_pubr()+
  ylab("Proportion of CAZymes (%)")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
b_gt <- ggplot(gt_diamond_taxa, aes(fill=enzymeNAME, y=Proportion, x=phylum)) + 
  geom_bar(position="fill", stat="identity",colour="black")+
  theme(axis.text.x=element_text(angle=90))+
  scale_fill_brewer(palette="Spectral")+
  theme_pubr()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.title.x = element_blank(),
        legend.position="right",
        legend.title = element_blank())


CAZ_GTphyla_plot <- plot_grid(
  plot_grid(
    a_gt + theme(legend.position = "none"),
    b_gt + theme(legend.position = "none"),
    ncol = 1,
    align = "v"),
  plot_grid(
    get_legend(a_gt),
    get_legend(b_gt),
    ncol =1,
    rel_heights = c(3,7)),
  rel_widths = c(8,2)
)

pdf("dbcan_CAZ_GTphyla.pdf")
CAZ_GTphyla_plot
dev.off()


##########################################################
##########################################################


# PART 2


gt_diamond <- as.data.frame(gt_diamond)

# normalize by lib size: 
no_reps_all_norm <- no_reps_all %>%
  group_by(cohort,pig,date) %>%
  mutate(norm_value=value/sum(value))

diamond_count <- diamond %>%
  group_by(pig,bin,enzymeID,enzymeNAME)%>%
  dplyr::summarise(enz_count=n())


df <- left_join(no_reps_all_norm,diamond_count) 


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


NROW(all_pvalues)

write.csv(all_pvalues,
          "dbcan_pvalues", 
          row.names = FALSE)

significant <- all_pvalues %>%
  filter(p_value<0.05) %>%
  arrange(p_value)
tail(significant)
mylist <- unique(significant$enzID)
head(mylist)
head(significant)




##########################################################
##########################################################


# PART 3

# now I will use these IDs to immediately plot interesting stuff
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

# set defined colors (releveling )
scale_fill_gaio <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(c('#C72A3F', '#EC7746', '#FADA78', '#E0F686','#8ACE81','#2872AF'), 
                      c("AA","CE","GH","GT","CBM","PL")), 
    ...
  )
}

# put df rows in order of list (list is ordered by p-value (descending))
require(gdata)
df_part$enzymeID <- reorder.factor(df_part$enzymeID, new.order=mylist)
df_part <- df_part %>%
  arrange(enzymeID)
# now they are plotted in order of significance! 
# that means that the first pages will be most interesting


pdf("dbcan_CAZ_time_boxplots.pdf")
for (i in seq(1, length(mylist), 12)) {    # can also use: length(unique(df_part$enzymeID))
  
  print(ggplot(df_part[df_part$enzymeID %in% mylist[i:(i+11)], ], 
               aes(date, log(tot),fill=enzymeNAME)) + 
          geom_boxplot(outlier.size = 1) +
          facet_wrap(~ enzymeID, scales = "free_y") +
          theme_bw()+
          scale_fill_gaio()+
          theme(legend.position="none",
                axis.text.y=element_text(size = 4),
                axis.ticks.length.y = unit(.05, "cm"),
                axis.text.x=element_text(size=6))+
          geom_text(aes(x="t1", y=max(log(tot)), label=paste0("n=(",n,")")),
                    size=2,colour="black", inherit.aes=TRUE, parse=FALSE,check_overlap = TRUE)+
          geom_text(aes(x="t10", y=max(log(tot)), label=paste0("n=(",n_moms,")")),
                    size=2,colour="black", inherit.aes=TRUE, parse=FALSE,check_overlap = TRUE))
}
dev.off()


# use list of enzymes that turned out to be significant
# to retrieve taxonomic info of those enzymes
gt_diamond_sub <- subset(gt_diamond, (enzymeID %in% mylist)) 
# plotting treemaps representing bins (species level) that carry a specific enzyme


# as this one takes long, uncomment when you want to plot again 

# pdf("dbcan_CAZ_species_treemaps_new.pdf")
# for (i in seq(1, length(mylist), by = 1)) {    # can also use: length(unique(df_part$enzymeID))
# 
#   this <- gt_diamond_sub[gt_diamond_sub$enzymeID %in% mylist[i], ] %>%
#     dplyr::select(species,family) %>%
#     distinct() %>%
#     group_by(family) %>%
#     add_tally() %>%
#     mutate(perc=round(n/sum(n)*100,2)) %>%
#     drop_na()
# 
#   ID <- paste0(mylist[i])
#   this <- as.data.frame(this)
# 
#   print(ggplot(this, aes(area = perc, fill = family, label = species,
#                          subgroup = family, subgroup2=species)) +
#           geom_treemap() +
#           geom_treemap_subgroup_border() +
#           geom_treemap_subgroup2_border(colour="black",size=1) +
#           # geom_treemap_subgroup_text(place = "top", grow = T, alpha = .9, colour =
#           #                              "White", fontface = "italic", min.size = 0) +
#           geom_treemap_subgroup2_text(place = "topleft", colour ="Black", fontface = "italic",
#                                       min.size = 1,reflow=T,grow=T)+
#           theme(legend.position="none")+
#           ggtitle(ID))
# }
# dev.off()


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
list_GT_1 <- unique(df_part_GT$enzymeID)[1:20]
list_GT_2 <- unique(df_part_GT$enzymeID)[21:39]
df_part1_GT <- subset(df_part_GT, (enzymeID %in% list_GT_1))
df_part2_GT <- subset(df_part_GT, (enzymeID %in% list_GT_2))

# GH
df_part_GH <- subset(df_part, (enzymeNAME %in% "GH"))
list_GH_1 <- unique(df_part_GH$enzymeID)[1:25]
list_GH_2 <- unique(df_part_GH$enzymeID)[26:50]
list_GH_3 <- unique(df_part_GH$enzymeID)[51:75]
list_GH_4 <- unique(df_part_GH$enzymeID)[76:102]
df_part1_GH <- subset(df_part_GH, (enzymeID %in% list_GH_1))
df_part2_GH <- subset(df_part_GH, (enzymeID %in% list_GH_2))
df_part3_GH <- subset(df_part_GH, (enzymeID %in% list_GH_3))
df_part4_GH <- subset(df_part_GH, (enzymeID %in% list_GH_4))

# CBM
df_part_CBM <- subset(df_part, (enzymeNAME %in% "CBM"))
list_CBM_1 <- unique(df_part_CBM$enzymeID)[1:22]
list_CBM_2 <- unique(df_part_CBM$enzymeID)[23:44]
df_part1_CBM <- subset(df_part_CBM, (enzymeID %in% list_CBM_1))
df_part2_CBM <- subset(df_part_CBM, (enzymeID %in% list_CBM_2))



######
# get species count for each enzymeID and make subsets of the data 

t <- gt_diamond %>%
  dplyr::select(pig,enzymeID,enzymeNAME,species,family) %>%   # dplyr::select(pig,contig,enzymeID,enzymeNAME,species,family) %>%
  distinct() %>%
  group_by(pig,enzymeID,species) %>%
  add_tally() %>%
  group_by(enzymeNAME,enzymeID,species) %>%
  dplyr::summarise(n_sum_species=sum(n)) %>%
  mutate(perc_species=round(n_sum_species/sum(n_sum_species)*100,2)) %>% 
  drop_na()
s <- as.data.frame(t)

# ###
# # it makes sense: little test: 
# test4 <- gt_diamond %>% 
#   filter(enzymeID=="AA1") %>%
#   group_by(pig,enzymeID,species) %>% 
#   tally() %>%
#   group_by(enzymeID,species) %>%
#   dplyr::summarise(allpiggies=sum(n)) %>%
#   mutate(allpiggiesperc=allpiggies/sum(allpiggies))
# head(test4)
# ###

# subsets
s_part1_GH <- subset(s, (enzymeID %in% list_GH_1))
s_part2_GH <- subset(s, (enzymeID %in% list_GH_2))
s_part3_GH <- subset(s, (enzymeID %in% list_GH_3))
s_part4_GH <- subset(s, (enzymeID %in% list_GH_4))
s_part1_GT <- subset(s, (enzymeID %in% list_GT_1))
s_part2_GT <- subset(s, (enzymeID %in% list_GT_2))
s_part1_CBM <- subset(s, (enzymeID %in% list_CBM_1))
s_part2_CBM <- subset(s, (enzymeID %in% list_CBM_2))
s_part_AA_CE <- subset(s, (enzymeID %in% list_AA_CE))
s_part_PL <- subset(s, (enzymeID %in% list_PL))
######

##########################################################
##########################################################


# PART 4


# function for making enzyme heatmaps
make_enzyme_heatmap <- function(x) {
  p <- x %>% 
    filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="t10"|date=="tM") %>%
    group_by(pig, enzymeID, date) %>% 
    dplyr::summarise(tot = mean(tot, na.rm = TRUE)) %>% 
    group_by(enzymeID) %>%
    mutate(tot = tot/max(tot)) %>%
    ungroup() %>% 
    ggplot(aes(x = factor(pig), y = reorder(enzymeID, tot, FUN = mean), fill = tot)) +
    geom_tile() +
    scale_fill_distiller(type = "div", palette = "Spectral") +
    facet_grid(~date, scales = "free") +
    labs(x = "date", y = "enzymeID", fill = "normalized abundance")+
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          legend.position="top")
  return(p)
}

# function for making species heatmaps 
make_species_heatmap <- function(species_df,CAZ_heatmap) {
  
  CAZ_order_vector <- CAZ_heatmap$data %>% group_by(enzymeID) %>% dplyr::summarise(mean=mean(tot)) %>% 
    group_by(enzymeID) %>%
    arrange(desc(mean))
  
  g <- species_df %>% 
    group_by(enzymeID,species) %>% 
    dplyr::summarise(n_sum_species = mean(n_sum_species, na.rm = TRUE)) %>% 
    group_by(enzymeID) %>%
    mutate(n_sum_species = n_sum_species/sum(n_sum_species)) %>%
    top_n(n = 4, wt = n_sum_species) %>%
    slice(1:4) %>%
    pivot_wider(names_from = enzymeID,values_from=n_sum_species) %>% 
    pivot_longer(cols=-species,names_to="enzymeID",values_to = "n_sum_species",values_drop_na = FALSE) %>%
    ungroup() 
  
  #require(gdata)
  g$enzymeID <- reorder.factor(g$enzymeID, new.order=CAZ_order_vector$enzymeID)
  
  p <- ggplot(g, aes(x = reorder(enzymeID, n_sum_species, FUN = mean), y = species, fill = n_sum_species)) +
    geom_tile(size = 0.5, color = "black") +
    scale_fill_distiller(palette = "Spectral",na.value = "black") +
    labs(x = "date", y = "enzymeID", fill = "normalized abundance")+
    theme_bw()+
    theme(legend.position="top",
          axis.text.x=element_text(angle=90,size=6),
          axis.text.y=element_text(size=5),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  #scale_y_discrete(labels=function(x){sub("\\s", "\n", x)})
  
  return(p)
}

# function for making top n species plots (per enzymeID) 
make_species_CAZ_plots <- function(x) {
  
  p <- x %>% 
    group_by(enzymeID) %>%
    top_n(n = 3, wt = perc_species) %>%
    slice(1:3) %>%
    drop.levels() %>% 
    ggplot(aes(x = reorder(enzymeID, perc_species, FUN = mean), y = reorder(species, perc_species, FUN = mean), fill = perc_species)) +
    geom_tile() +
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
    scale_fill_gaio()+
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
a_GH_1 <- make_enzyme_boxplots(df_part1_GH)
a_GH_2 <- make_enzyme_boxplots(df_part2_GH)
a_GH_3 <- make_enzyme_boxplots(df_part3_GH)
a_GH_4 <- make_enzyme_boxplots(df_part4_GH)
a_CBM_1 <- make_enzyme_boxplots(df_part1_CBM)
a_CBM_2 <- make_enzyme_boxplots(df_part2_CBM)

b_AA_CE <- make_species_CAZ_plots(s_part_AA_CE)
b_PL <- make_species_CAZ_plots(s_part_PL)
b_GT_1 <- make_species_CAZ_plots(s_part1_GT)
b_GT_2 <- make_species_CAZ_plots(s_part2_GT)
b_GH_1 <- make_species_CAZ_plots(s_part1_GH)
b_GH_2 <- make_species_CAZ_plots(s_part2_GH)
b_GH_3 <- make_species_CAZ_plots(s_part3_GH)
b_GH_4 <- make_species_CAZ_plots(s_part4_GH)
b_CBM_1 <- make_species_CAZ_plots(s_part1_CBM)
b_CBM_2 <- make_species_CAZ_plots(s_part2_CBM)

h_a_AA_CE <- make_enzyme_heatmap(df_part_AA_CE)
h_a_PL <- make_enzyme_heatmap(df_part_PL)
h_a_GT_1 <- make_enzyme_heatmap(df_part1_GT)
h_a_GT_2 <- make_enzyme_heatmap(df_part2_GT)
h_a_GH_1 <- make_enzyme_heatmap(df_part1_GH)
h_a_GH_2 <- make_enzyme_heatmap(df_part2_GH)
h_a_GH_3 <- make_enzyme_heatmap(df_part3_GH)
h_a_GH_4 <- make_enzyme_heatmap(df_part4_GH)
h_a_CBM_1 <- make_enzyme_heatmap(df_part1_CBM)
h_a_CBM_2 <- make_enzyme_heatmap(df_part2_CBM)







####
# combining the plots: CAZ time trend boxplots with species per enzymeID info: 
AA_CE <- plot_grid(a_AA_CE,b_AA_CE,ncol=2)

PL <- plot_grid(a_PL,b_PL,ncol=2)

CBM1 <- plot_grid(a_CBM_1,b_CBM_1,ncol=2)
CBM2 <- plot_grid(a_CBM_2,b_CBM_2,ncol=2)

GT1 <- plot_grid(a_GT_1,b_GT_1,ncol=2)
GT2 <- plot_grid(a_GT_2,b_GT_2,ncol=2)

GH1 <- plot_grid(a_GH_1,b_GH_1,ncol=2)
GH2 <- plot_grid(a_GH_2,b_GH_2,ncol=2)
GH3 <- plot_grid(a_GH_3,b_GH_3,ncol=2)
GH4 <- plot_grid(a_GH_4,b_GH_4,ncol=2)
####



####
# Two heatmaps per page: one for the CAZ time trend, one for the species corresponding to the bin where the CAZ was found 
sh_AA_CE <- make_species_heatmap(s_part_AA_CE,h_a_AA_CE)
H1 <- plot_grid(h_a_AA_CE,sh_AA_CE,ncol=2)

sh_PL <- make_species_heatmap(s_part_PL,h_a_PL)
H2 <- plot_grid(h_a_PL,sh_PL,ncol=2)

sh_GT1 <- make_species_heatmap(s_part1_GT,h_a_GT_1)
H3 <- plot_grid(h_a_GT_1,sh_GT1,ncol=2)

sh_GT2 <- make_species_heatmap(s_part2_GT,h_a_GT_2)
H4 <- plot_grid(h_a_GT_2,sh_GT2,ncol=2)

sh_GH1 <- make_species_heatmap(s_part1_GH,h_a_GH_1)
H5 <- plot_grid(h_a_GH_1,sh_GH1,ncol=2)

sh_GH2 <- make_species_heatmap(s_part2_GH,h_a_GH_2)
H6 <- plot_grid(h_a_GH_2,sh_GH2,ncol=2)

sh_GH3 <- make_species_heatmap(s_part3_GH,h_a_GH_3)
H7 <- plot_grid(h_a_GH_3,sh_GH3,ncol=2)

sh_CBM1 <- make_species_heatmap(s_part1_CBM,h_a_CBM_1)
H8 <- plot_grid(h_a_CBM_1,sh_CBM1,ncol=2)

sh_CBM2 <- make_species_heatmap(s_part2_CBM,h_a_CBM_2)
H9 <- plot_grid(h_a_CBM_2,sh_CBM2,ncol=2)
####




pdf("dbcan_CAZ_time_heatmaps.pdf")
h_a_AA_CE 
h_a_PL 
h_a_GT_1 
h_a_GT_2 
h_a_GH_1 
h_a_GH_2 
h_a_GH_3 
h_a_CBM_1 
h_a_CBM_2 
dev.off()

pdf("dbcan_CAZ_ALL_heatmaps.pdf")
H1
H2
H3
H4
H5
H6
H7
H8
H9
dev.off()

pdf("dbcan_CAZ_time_species.pdf")
AA_CE
PL
CBM1
CBM2
GT1
GT2
GH1
GH2
GH3
GH4
dev.off()



# PLOT showing all enzymes that significantly changed over time, 
# with popularity (y) and "singularity" (x) meaning how uniquely an enzyme is expressed by one (right) or multiple (left) taxa


df <- gt_diamond %>%
  dplyr::select(pig,enzymeID,enzymeNAME,species,family) %>%
  distinct() %>%
  group_by(pig,enzymeID,species) %>%
  add_tally() %>%
  group_by(enzymeNAME,enzymeID,species) %>%
  dplyr::summarise(n_sum_species=sum(n)) %>%
  mutate(perc_species=round(n_sum_species/sum(n_sum_species)*100,2)) %>% 
  drop_na()

df2 <- as.data.frame(df)
df2 <- df2 %>%
  group_by(enzymeID) %>%
  top_n(n = 1, wt = perc_species) 

# subsetting to contain only sigbificantly changing enzymes over time
df2 <- subset(df2, (enzymeID %in% df_part$enzymeID))

# to get the species showing in the labels 
df2$species <- gsub(" ","\n", df2$species)

df2$enzymeNAME  = factor(
  df2$enzymeNAME, levels=c("AA",
                            "CE",
                            "GH",
                            "GT",
                            "CBM",
                            "PL"))

pdf("dbcan_all_CAZ_plot.pdf")
ggplot(df2) +
  geom_point(aes(perc_species, n_sum_species, fill = factor(enzymeNAME))) +
  #ylim(0,10000)+
  #xlim(-25,100)+
  geom_point(aes(perc_species, n_sum_species, fill=enzymeNAME), 
             colour="black",pch=21, size=4)+
  labs(x="percentage of top species per enzymeID",
       y="number of subjects carrying the enzyme",
       fill="enzyme class")+
  scale_fill_manual(values = setNames(c('#C72A3F', '#EC7746', '#FADA78', '#E0F686','#8ACE81','#2872AF'), 
                                      c("AA","CE","GH","GT","CBM","PL")))+
  theme_bw(base_size = 10)+
  ggrepel::geom_label_repel(data = subset(df2, perc_species >= 50),
                            aes(
                              x = perc_species,
                              y = n_sum_species,
                              fill = factor(enzymeNAME),
                              label = paste0(enzymeID,"\n",species)
                            ),
                            #nudge_y       = 100 - subset(test, perc_species >= 55)$perc_species, #rep(c(4000,45000),7),
                            box.padding   = 0.3,label.padding = 0.1,
                            point.padding = 0.5,
                            #force         = 70,
                            segment.size  = 0.2,
                            segment.color = "grey50",
                            #direction     = "x",
                            size=2.5,
                            color="black"
  ) 
dev.off()



# CE

sink(file = "enzymes_popularity and species-uniqueness",append = FALSE)

paste0(" ##############################  CE  ############################## ")
paste0("number of all enzymeIDs within this class")
NROW(subset(df2, enzymeNAME == "CE")) 
paste0("number of all enzymeIDs within this class, that are present in multiple species")
NROW(subset(df2, enzymeNAME == "CE" & perc_species<25)) 
CE_df2 <- df2 %>%
  filter(enzymeNAME=="CE") %>%
  filter(perc_species<25) 
summary(CE_df2)

paste0(" ##############################  PL  ############################## ")
paste0("number of all enzymeIDs within this class")
NROW(subset(df2, enzymeNAME == "PL")) 
paste0("number of all enzymeIDs within this class, that are present in multiple species")
NROW(subset(df2, enzymeNAME == "PL" & perc_species<25)) 
PL_df2 <- df2 %>%
  filter(enzymeNAME=="PL") %>%
  filter(perc_species<25) 
summary(PL_df2)

paste0(" ##############################  GT  ############################## ")
paste0("number of all enzymeIDs within this class")
NROW(subset(df2, enzymeNAME == "GT")) 
paste0("number of all enzymeIDs within this class, that are present in multiple species")
NROW(subset(df2, enzymeNAME == "GT" & perc_species<25)) 
GT_df2 <- df2 %>%
  filter(enzymeNAME=="GT") %>%
  filter(perc_species<25) 
summary(GT_df2)

paste0(" ##############################  GH  ############################## ")
paste0("number of all enzymeIDs within this class")
NROW(subset(df2, enzymeNAME == "GH")) 
paste0("number of all enzymeIDs within this class, that are present in multiple species")
NROW(subset(df2, enzymeNAME == "GH" & perc_species<25)) 
GH_df2 <- df2 %>%
  filter(enzymeNAME=="GH") %>%
  filter(perc_species<25) 
summary(GH_df2)

paste0(" ##############################  CBM  ############################## ")
paste0("number of all enzymeIDs within this class")
NROW(subset(df2, enzymeNAME == "CBM")) 
CBM_df2 <- df2 %>%
  filter(enzymeNAME=="CBM") 
summary(CBM_df2)

paste0(" ##############################  AA  ############################## ")
paste0("number of all enzymeIDs within this class")
NROW(subset(df2, enzymeNAME == "AA")) 
AA_df2 <- df2 %>%
  filter(enzymeNAME=="AA") 
summary(AA_df2)


paste0(" ##############################  enzymes present majorly in one species ############################## ")
rightmost <- df2 %>%
  filter(perc_species>70) %>%
  arrange(n_sum_species)
NROW(rightmost)
rightmost

sink()





# Prevotella vs Bacteroidetes enzyme specificity

test <- subset(df_part, (enzymeID %in% mylist))


# 1 
# plot all prevotella and all bacteroidetes genus, all enzymes, no labels 


# 2 
# plot only Bacteroidetes A and Prevotella, 50 most abundant enzymes; no labels 

# 3 
# plot only Bacteroidetes A and Prevotella, 60th-100th most abundant enzymes; labels 




gt_diamond_sub <- gt_diamond_sub %>% 
  filter(genus=="Prevotella"|genus=="Bacteroides_A"|genus=="Bacteroides"|genus=="Bacteroides_B") 

# reorder  
gt_diamond_sub$genus  = factor(gt_diamond_sub$genus, levels=c("Bacteroides_A",
                                                              "Bacteroides_B", 
                                                              "Bacteroides",
                                                              "Prevotella"))


# all Bact. and Prev
test2 <- gt_diamond_sub %>% 
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
test2 <- gt_diamond_sub %>% 
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




p1 <- ggarrange(p1)
p23 <- ggarrange(p2,p3, 
                 widths = c(1,2))
tosave <- ggarrange(p1,p23, widths = c(1,3))
tosave <- ggarrange(myleg,tosave,heights=c(1,9),nrow=2)

pdf("dbcan_specificity_genera.pdf")
tosave
dev.off()
