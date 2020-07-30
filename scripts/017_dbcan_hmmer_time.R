
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
library(zip)


source_dir = "/Users/12705859/metapigs_dry/source_data/" # git 
middle_dir = "/Users/12705859/metapigs_dry/middle_dir/" # git 
out_dir_git = "/Users/12705859/metapigs_dry/out/" # git 
out_dir = "/Users/12705859/Desktop/metapigs_dry/dbcan/"  # local

######################################################################

# load CleanCounts - hmmer data 

hmmerCleanCounts <- read.csv(paste0(middle_dir, "hmmer.out_CleanCounts.csv.gz"))
hmmerCleanCounts$enzymeNAME <- as.character(hmmerCleanCounts$enzymeNAME)
hmmerCleanCounts$enzymeNAME  = factor(hmmerCleanCounts$enzymeNAME, levels=c("AA","CE","SLH","GH","GT","CBM","cohesin","PL"))

##########################################################

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv(paste0(middle_dir,"gtdb_bins_completeTaxa"),
                        col_types = cols(node = col_character(),
                                         pig = col_character()))

##########################################################

# load counts data 

no_reps_all <- read.csv(paste0(middle_dir,"no_reps_all.csv"), 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)

# remove .fa extension to match bins in checkm df 
no_reps_all$bin <- gsub(".fa","", no_reps_all$bin)
# normalize by lib size: 
no_reps_all_norm <- no_reps_all %>%
  group_by(cohort,pig,date) %>%
  dplyr::mutate(norm_value=value/sum(value))

##########################################################
##########################################################
##########################################################


# CAZy enzymes (HMMER) join with GTDB

NROW(gtdbtk_bins)
NROW(hmmerCleanCounts)
# merge gtdbtk- bins taxonomic assignments with hmmer data
gt_hmmer <- inner_join(gtdbtk_bins,hmmerCleanCounts)
NROW(gt_hmmer)

gt_hmmer$phylum[is.na(gt_hmmer$phylum)] <- "Unknown"

gt_hmmer <- as.data.frame(gt_hmmer)

##########################################################

# create an enzyme to genus & species profile 

gt_hmmer_genus <- gt_hmmer %>% 
  group_by(enzymeID,genus) %>%
  dplyr::summarise(sum=sum(enz_count)) %>%
  group_by(enzymeID) %>%
  dplyr::mutate(percentage_genus=(sum/sum(sum))*100)

gt_hmmer_species <- gt_hmmer %>% 
  group_by(enzymeID,species) %>%
  dplyr::summarise(sum=sum(enz_count)) %>%
  group_by(enzymeID) %>%
  dplyr::mutate(percentage_species=(sum/sum(sum))*100)

##########################################################


# join bins abundance with enzyme data 

df <- left_join(no_reps_all_norm,hmmerCleanCounts) 

df1 <- df %>%
  dplyr::mutate(enz_countBYmap=norm_value*enz_count)

df1 <- df1[!is.na(df1$enzymeID),]

# within the same sample, sum the noralized values of enzymes that fall under the same enzyme ID 
df2 <- df1 %>%
  group_by(cohort,pig,date,enzymeID,enzymeNAME) %>%
  dplyr::summarise(norm_count=sum(enz_countBYmap))

DF <- as.data.frame(df2)
head(DF)
NROW(DF)

##########################################################
##########################################################


# PART 2 : X-Y plot: popularity (y) and "singularity" (x) 


##########################################################

# Set aesthetics 


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

##########################################################


# PLOT showing all enzymes 
# with popularity (y) and "singularity" (x) 
# meaning how uniquely an enzyme is expressed by one (right) or multiple (left) taxa
# and how popular among the pig population (top=popular; bottom=unpopular)

CAZy_species_representation <- gt_hmmer_species %>%
  drop_na() %>%
  group_by(enzymeID) %>%
  top_n(n=1,wt=percentage_species) %>%
  dplyr::arrange(enzymeID) %>%
  dplyr::mutate(top_species = percentage_species) %>%
  dplyr::select(enzymeID, species, top_species) %>%
  group_by(enzymeID) %>%
  arrange(desc(top_species)) %>%
  slice(1)

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

df_XY <- full_join(CAZy_species_representation, CAZy_prevalence_in_population) %>% distinct()

pdf(paste0(out_dir,"dbcan_HMMER_all_plot.pdf"))
ggplot(df_XY) +
  geom_point(aes(top_species, popularity, fill=enzymeNAME), 
             colour="black",pch=21, size=4)+
  labs(x="percentage of top species bearing the enzyme",
       y="prevalence in the pig population",
       fill="enzyme class")+
  scale_fill_gaio8()+
  theme_bw(base_size = 10) +
  xlim(0,130) +
  ggrepel::geom_label_repel(data = subset(df_XY, top_species > 55 & popularity > 5),
                            aes(
                              x = top_species,
                              y = popularity,
                              fill = factor(enzymeNAME),
                              label = paste0(enzymeID,"\n",species)
                            ), show.legend = FALSE,
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

###
# other way to represent species specificity and prevalence in pig population
popul_plot <- ggplot(df_XY,aes(x=enzymeNAME,y=popularity))+
  geom_boxplot()+
  geom_jitter(alpha=0.5, width = 0.1)+
  theme(axis.text.x = element_text(angle=90))

top_sp_plot <- ggplot(df_XY,aes(x=enzymeNAME,y=top_species))+
  geom_boxplot()+
  geom_jitter(alpha=0.5, width = 0.1)+
  theme(axis.text.x = element_text(angle=90))

pdf(paste0(out_dir,"dbcan_HMMER_popul&top_species.pdf"))
ggarrange(popul_plot,top_sp_plot)
dev.off()
###

##########################################################

# Some numbers for manuscript:

sink(file = paste0(out_dir,"dbcan_hmmerClean_numbers.txt"),append = FALSE)

z <- df_XY %>%
  group_by(enzymeNAME) %>%
  dplyr::summarise(median_species=round(median(top_species),1),
                   mean_species=round(mean(top_species),1),
                   median_popul=round(median(popularity),1),
                   mean_popul=round(mean(popularity),1),
                   number_of_enzIDs=n())
z

paste0(" ############################################################################################################ ")
paste0(" ##############################  AA enzymes species ######################################################### ")
df_XY %>% dplyr::filter(enzymeNAME=="AA")

paste0(" ############################################################################################################ ")
paste0(" ##############################  SLH,cohesin  ############################################################### ")
df_XY %>% dplyr::filter(enzymeNAME=="SLH"|enzymeNAME=="cohesin")

paste0(" ############################################################################################################ ")
paste0(" ##############################  Found majorly represented by a single species ############################## ")
df_XY %>%
  dplyr::filter(top_species>70) %>%
  dplyr::filter(popularity>15) %>%
  dplyr::arrange(desc(top_species))

sink()


##########################################################
##########################################################


# PART 3


##########################################################


# Test which enzymes significantly differ in abundance between time points 


DF <- as.data.frame(DF)

df_piggies <- DF %>%
  dplyr::filter(!cohort=="Mothers")

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
          paste0(out_dir_git,"dbcan_HMMER_pvalues.csv"), 
          row.names = FALSE)

significant <- all_pvalues %>%
  dplyr::filter(p_value<0.05) %>%
  dplyr::arrange(p_value)
tail(significant)
mylist <- unique(significant$enzID)

sink(file=paste0(out_dir,"dbcan_hmmerClean_numbers.txt"), append=TRUE)
paste0("number of total, unique CAZy enzymes in our dataset")
NROW(unique(hmmerCleanCounts$enzymeID)) # tot enzymes
paste0("number of unique enzymes that showed a significant change between timepoints")
NROW(unique(significant$enzID)) # tot enzymes with significant diffs when comparing any time point
sink()



##########################################################

###########
# create column to specify if subject ID is mother or piglet 

# create list of mother IDs
moms_IDs <- no_reps_all %>%
  dplyr::filter(date=="tM") %>%
  dplyr::select(pig) %>%
  distinct() %>%
  mutate(is_mom="yes")

# join list of mothers IDs to gt_hmmer
gt_hmmer <- left_join(gt_hmmer,moms_IDs)
gt_hmmer <- gt_hmmer %>% 
  dplyr::mutate(is_mom = replace_na(is_mom, "no")) 
###########

piglets_with_enzyme <- gt_hmmer %>%
  filter(!pig=="Protexin"&!pig=="ColiGuard"&!pig=="MockCommunity"&!pig=="NegativeControl") %>%
  dplyr::filter(!is_mom=="yes") %>%
  dplyr::select(pig,enzymeID, enzymeNAME) %>%
  distinct() %>%
  group_by(enzymeID, enzymeNAME) %>%
  tally() %>%
  dplyr::mutate(n_piggies=n) %>%
  dplyr::select(enzymeID, n_piggies) 

moms_with_enzyme <- gt_hmmer %>%
  filter(!pig=="Protexin"&!pig=="ColiGuard"&!pig=="MockCommunity"&!pig=="NegativeControl") %>%
  dplyr::filter(is_mom=="yes") %>%
  dplyr::select(pig,enzymeID, enzymeNAME) %>%
  distinct() %>%
  group_by(enzymeID, enzymeNAME) %>%
  tally() %>%
  dplyr::mutate(n_moms=n) %>%
  dplyr::select(enzymeID, n_moms) 
  

piglets_with_enzyme <- as.data.frame(piglets_with_enzyme)
moms_with_enzyme <- as.data.frame(moms_with_enzyme)

###########

# now I will use these IDs to plot interesting stuff
df_part <- subset(DF, (enzymeID %in% mylist))

# renaming this column 
colnames(df_part)[6] <- "tot"

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
                                              "t10",
                                              "tM"))

df_part <- df_part[!is.na(df_part$date),]
# put df rows in order of list (list is ordered by p-value (descending))
require(gdata)
df_part$enzymeID <- reorder.factor(df_part$enzymeID, new.order=mylist)
# now they will be plotted in order of significance!
# that means that the first pages will be most interesting

##########################################################

# PLOTS of enzyme trends alone 

## can be commented out as I won't be using it in the manuscript, 
## but it's nicer as it has time points (boxplots only)

# pdf(paste0(out_dir,"dbcan_HMMER_time_boxplots.pdf"))
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
#           geom_text(aes(x="t1", y=max(log(tot)), label=paste0("n=(",n_piggies,")")),
#                     size=2,colour="black", inherit.aes=TRUE, parse=FALSE,check_overlap = TRUE)+
#           geom_text(aes(x="t10", y=max(log(tot)), label=paste0("n=(",n_moms,")")),
#                     size=2,colour="black", inherit.aes=TRUE, parse=FALSE,check_overlap = TRUE))
# }
# dev.off()


##########################################################

#  filter out smaller collection date - time points 
df_part <- df_part %>%
  filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="tM")
#  re-level 
df_part$date  = factor(df_part$date, levels=c("t0","t2","t4","t6","t8","tM"))

##########################################################

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
list_GH_4 <- unique(df_part_GH$enzymeID)[76:105]
df_part1_GH <- subset(df_part_GH, (enzymeID %in% list_GH_1))
df_part2_GH <- subset(df_part_GH, (enzymeID %in% list_GH_2))
df_part3_GH <- subset(df_part_GH, (enzymeID %in% list_GH_3))
df_part4_GH <- subset(df_part_GH, (enzymeID %in% list_GH_4))

# CBM
df_part_CBM <- subset(df_part, (enzymeNAME %in% "CBM"))
NROW(unique(df_part_CBM$enzymeID))
list_CBM_1 <- unique(df_part_CBM$enzymeID)[1:14]
list_CBM_2 <- unique(df_part_CBM$enzymeID)[15:29]
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
# formula is the same as the one used before, except the top_n(1), 
# as we need the full profile of proportions and not just the top 1 


species <- gt_hmmer_species %>%
  drop_na() %>%
  group_by(enzymeID) %>%
  top_n(n=3,wt=percentage_species) %>%
  dplyr::arrange(enzymeID) %>%
  dplyr::mutate(top_species = percentage_species) %>%
  dplyr::select(enzymeID, species, top_species) %>%
  group_by(enzymeID) %>%
  arrange(desc(top_species)) %>%
  slice(1:3)
species <- as.data.frame(species)

# # enzyme order must follow the order of significance
# require(gdata)
# species$enzymeID <- reorder.factor(species$enzymeID, new.order=mylist)
# species <- species %>%
#   dplyr::arrange(enzymeID)
# View(species)

genus <- gt_hmmer_genus %>%
  drop_na() %>%
  group_by(enzymeID) %>%
  top_n(n=3,wt=percentage_genus) %>%
  dplyr::arrange(enzymeID) %>%
  dplyr::mutate(top_genus = percentage_genus) %>%
  dplyr::select(enzymeID, genus, top_genus) %>%
  group_by(enzymeID) %>%
  arrange(desc(top_genus)) %>%
  slice(1:3)
genus <- as.data.frame(genus)
# # enzyme order must follow the order of significance
# require(gdata)
# genus$enzymeID <- reorder.factor(genus$enzymeID, new.order=mylist)
# genus <- genus %>%
#   dplyr::arrange(enzymeID)


# subsets (species)
s_part1_GH <- subset(species, (enzymeID %in% list_GH_1))
s_part2_GH <- subset(species, (enzymeID %in% list_GH_2))
s_part3_GH <- subset(species, (enzymeID %in% list_GH_3))
s_part4_GH <- subset(species, (enzymeID %in% list_GH_4))
s_part1_GT <- subset(species, (enzymeID %in% list_GT_1))
s_part2_GT <- subset(species, (enzymeID %in% list_GT_2))
s_part3_GT <- subset(species, (enzymeID %in% list_GT_3))
s_part1_CBM <- subset(species, (enzymeID %in% list_CBM_1))
s_part2_CBM <- subset(species, (enzymeID %in% list_CBM_2))
s_part_AA_CE <- subset(species, (enzymeID %in% list_AA_CE))
s_part_PL <- subset(species, (enzymeID %in% list_PL))
s_part_SLH <- subset(species, (enzymeID %in% list_SLH))
s_part_cohesin <- subset(species, (enzymeID %in% list_cohesin))
######

# subsets (genus)
s_genus_part1_GH <- subset(genus, (enzymeID %in% list_GH_1))
s_genus_part2_GH <- subset(genus, (enzymeID %in% list_GH_2))
s_genus_part3_GH <- subset(genus, (enzymeID %in% list_GH_3))
s_genus_part4_GH <- subset(genus, (enzymeID %in% list_GH_4))
s_genus_part1_GT <- subset(genus, (enzymeID %in% list_GT_1))
s_genus_part2_GT <- subset(genus, (enzymeID %in% list_GT_2))
s_genus_part3_GT <- subset(genus, (enzymeID %in% list_GT_3))
s_genus_part1_CBM <- subset(genus, (enzymeID %in% list_CBM_1))
s_genus_part2_CBM <- subset(genus, (enzymeID %in% list_CBM_2))
s_genus_part_AA_CE <- subset(genus, (enzymeID %in% list_AA_CE))
s_genus_part_PL <- subset(genus, (enzymeID %in% list_PL))
s_genus_part_SLH <- subset(genus, (enzymeID %in% list_SLH))
s_genus_part_cohesin <- subset(genus, (enzymeID %in% list_cohesin))
######

##########################################################
##########################################################


# PART 4


##########################################################

# # function for making enzyme heatmaps (whole large heatmap)
# make_enzyme_heatmap <- function(x) {
#   p <- x %>%
#     dplyr::filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="t10"|date=="tM") %>%
#     group_by(pig, enzymeID, date) %>%
#     dplyr::summarise(tot = mean(tot, na.rm = TRUE)) %>%
#     group_by(enzymeID) %>%
#     dplyr::mutate(tot = tot/max(tot)) %>%
#     dplyr::ungroup() %>%
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

# # # function for making species heatmaps
# make_species_heatmap <- function(species_df,CAZ_heatmap) {
# 
#   CAZ_order_vector <- CAZ_heatmap$data %>% 
#     group_by(enzymeID) %>% 
#     dplyr::summarise(mean=mean(tot)) %>%
#     group_by(enzymeID) %>%
#     dplyr::arrange(desc(mean))
# 
#   g <- species_df %>%
#     group_by(enzymeID,species) %>%
#     dplyr::summarise(n_sum_species = mean(n_sum_species, na.rm = TRUE)) %>%
#     group_by(enzymeID) %>%
#     dplyr::mutate(n_sum_species = n_sum_species/sum(n_sum_species)) %>%
#     top_n(n = 4, wt = n_sum_species) %>%
#     dplyr::arrange(desc(n)) %>%
#     dplyr::slice(1:4) %>%
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
  
  # enzyme order must follow the order of significance
  require(gdata)
  x$enzymeID <- reorder.factor(x$enzymeID, new.order=mylist)
  
  p <- x %>% 
    ggplot(aes(x = reorder(enzymeID, top_species, FUN = mean), y = reorder(species, top_species, FUN = mean), fill = top_species)) +
    geom_tile() +
    #scale_fill_gaio8()+
    scale_fill_distiller(type = "div", palette = "Spectral") +
    geom_text(aes(y=species,label=round(top_species)), size=2)+
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
# WARNING (does not keep the same order of plots as make_enzyme_boxplots function below)
make_genus_CAZ_plots <- function(x) {
  
  # enzyme order must follow the order of significance
  require(gdata)
  x$enzymeID <- reorder.factor(x$enzymeID, new.order=mylist)
  
  p <- x %>% 
    ggplot(aes(x = reorder(enzymeID, top_genus, FUN = mean), y = reorder(genus, top_genus, FUN = mean), fill = top_genus)) +
    geom_tile() +
    #scale_fill_gaio8()+
    scale_fill_distiller(type = "div", palette = "Spectral") +
    geom_text(aes(y=genus,label=round(top_genus)), size=2)+
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
  
  x1 <- x %>% 
    dplyr::filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="t10"|date=="tM") %>%
    group_by(pig, enzymeID, enzymeNAME,date,n_piggies,n_moms) %>% 
    dplyr::summarise(tot = mean(tot, na.rm = TRUE)) %>% 
    dplyr::mutate(tot=log(tot)) %>%
    dplyr::arrange(desc(enzymeID)) %>%
    drop.levels() 
  
  # enzyme order must follow the order of significance
  require(gdata)
  x1$enzymeID <- reorder.factor(x1$enzymeID, new.order=mylist)
  
  p <- x1 %>% 
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
    geom_text(aes(x="ss_piglets",y=Inf, label=paste0("n=",n_piggies)),
              size=2.3,colour="black", hjust = 1.5, angle=90, inherit.aes=TRUE, parse=FALSE,check_overlap = TRUE)+
    geom_text(aes(x="ss_sows",y=Inf, label=paste0("n=",n_moms)),
              size=2.3,colour="black", hjust = 1.5, angle=90, inherit.aes=TRUE, parse=FALSE,check_overlap = TRUE)
  
  return(p)
  
}

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

# pdf("dbcan_HMMER_time_heatmaps.pdf")
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
# pdf("dbcan_ALL_heatmaps.pdf")
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

pdf(paste0(out_dir,"dbcan_HMMER_time_species.pdf"))
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

pdf(paste0(out_dir,"dbcan_HMMER_time_genus.pdf"))
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



##########################################################

