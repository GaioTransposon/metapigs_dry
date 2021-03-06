
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


##########################################################
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



##########################################################


unique(df_part$date)

significant <- all_pvalues %>%
  dplyr::filter(    enzID=="GH54"|
                      #enzID=="GT53"|
                      enzID=="CBM71"|
                      enzID=="GT94"|
                      #enzID=="GH70"|
                      enzID=="GH110"|
                      #enzID=="GH1"|
                    enzID=="CE13"|
                    enzID=="GH25"|
                    enzID=="GH68"|
                  enzID=="GT85"|
                  enzID=="AA3") %>% #CBM25
  dplyr::filter(p_value<0.05) #%>%
  #dplyr::arrange(p_value)
tail(significant)
mylist <- unique(significant$enzID)




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

##########################################################

# subsetting whole

df_part <- df_part
list <- unique(df_part$enzymeID)


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
s_part <- subset(species, (enzymeID %in% list))

######


##########################################################
##########################################################


# PART 4


##########################################################

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
          legend.position="top",
          #legend.justification="right",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-6,-6,-6,-6), # get closer farther from zero to get it closer to edge 
          axis.ticks.x = element_blank(),
          strip.text = element_text(angle = 90))+
    labs(fill="percentage of species \n carrying the enzyme")+
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
    dplyr::filter(!date=="t10") %>%
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
    labs(x = "sample size and timepoint", y = "enzyme abundance (log)", fill = "avg. abundance (log)")+
    theme_bw()+
    theme(legend.position="none",
          axis.text.y=element_text(size = 5),
          axis.ticks.length.y = unit(.05, "cm"),
          axis.text.x=element_text(size=5,angle=90,vjust=0.5))+
    geom_text(aes(x="ss_piglets",y=Inf, label=paste0("n=",n_piggies)),
              size=1.5,colour="black", hjust = 1.5, angle=90, inherit.aes=TRUE, parse=FALSE,check_overlap = TRUE)+
    geom_text(aes(x="ss_sows",y=Inf, label=paste0("n=",n_moms)),
              size=1.5,colour="black", hjust = 1.5, angle=90, inherit.aes=TRUE, parse=FALSE,check_overlap = TRUE)
  
  return(p)
  
}

b <- make_species_CAZ_plots(s_part) 

a <- make_enzyme_boxplots(df_part)


####
# combining the plots: CAZ time trend boxplots with species per enzymeID info: 
A <- plot_grid(a,b,ncol=2, labels = c("A","B"))
####


pdf(paste0(out_dir,"dbcan_HMMER_time_species_sub.pdf"), width=6.5, height=4.5)
A
dev.off()

