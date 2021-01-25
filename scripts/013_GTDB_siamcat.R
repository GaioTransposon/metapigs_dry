

library(readr)
library(splitstackshape)
library(tidyr)
library(dplyr)
library(ggplot2)
library(SIAMCAT)
library(matrixStats)
library(data.table)
library(pheatmap)
library(readxl)
library(ggpubr)
library(forcats)
library(data.table)
library(gmodels)




source_dir = "/Users/12705859/metapigs_dry/source_data/" # git 
middle_dir = "/Users/12705859/metapigs_dry/middle_dir/" # git 
out_dir_git = "/Users/12705859/metapigs_dry/out/" # git 
out_dir = "/Users/12705859/Desktop/metapigs_dry/gtdbtk/"  # local

######################################################################

# input files: 
# gtdbtk_bins_completeTaxa
# no_reps_all.csv (BINS COUNTS)
# dictionary

# OUTPUTS:

# lots!!


######################################################################

# a look at SIAMCAT dummy dataset:

#these are the counts and metadata from siamcat as they should look 
head(feat.crc.zeller)
head(meta.crc.zeller)
class(feat.crc.zeller)
class(meta.crc.zeller)

# the sum of each column (1 column: 1 sample) always equals 1. 
# this implies SIAMCAT uses the simplex (ratios)
colSums(feat.crc.zeller)

######################################################################

# counts data 

no_reps_all <- read.csv(paste0(middle_dir,"no_reps_all.csv"), 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)


# remove .fa extension 
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
unique(mdat$Cohort)


###
# metadata mdat2 is used below for the cohousing effect

mdat2 <- mdat %>%
  dplyr::filter(!Cohort=="Sows")  %>%
  dplyr::filter(!`*collection_date`=="2017-01-31"|
                  `*collection_date`=="2017-02-01"|
                  `*collection_date`=="2017-02-03") %>%
  dplyr::select(isolation_source,PigPen)

colnames(mdat2) <- c("pig","pen")
mdat2 <- as.data.frame(mdat2)

mdat2 <- mdat2 %>%
  group_by(pig) %>%
  dplyr::distinct()
NROW(mdat2)

mdat2$pen <- gsub("nan",NA,mdat2$pen)
mdat2 <- na.omit(mdat2)

# we need to keep only record of pigs that were not relocated. 
mdat2 <- setDT(mdat2)[,if(.N ==1) .SD,by=pig]
###

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
head(weights)

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

weights_rest<- as.data.frame(weights_rest)

# quickly visualize weight category distribution
# ggplot(weights_rest,aes(x=date,fill=weight_category)) +
#   geom_bar() +
#   facet_wrap(~birth_day)

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
# df0$gOTU <- paste0("gOTU_",df0$gOTU)

df0$gOTU <- as.character(df0$gOTU)


NROW(unique(df0$gOTU))
NROW(df0)

df0$gOTU <- paste0(df0$species,"__",df0$gOTU)


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




######################################################################

# CREATE COUNTS TABLE (like feat.crc.zeller)

df1 <- df0
head(df1)
NROW(df1)

### minitest: 
# filter a sample (pig,date)
test <- df1 %>% dplyr::filter(pig=="14159") %>% dplyr::filter(date=="t2")
sum(test$value)
NROW(test)
# for each sample (pig,date), sum up together the counts that fall within one species (same species assigned to distinct bins)
test2 <- test %>%
  group_by(pig,gOTU,date) %>%
  dplyr::summarize(sum_value = sum(value)) 
sum(test2$sum_value)
NROW(test2)
# normalize by library size 
test3 <- test2 %>% 
  group_by(pig,date) %>% 
  dplyr::mutate(norm_value = sum_value/sum(sum_value)) %>% 
  dplyr::select(-sum_value)
sum(test3$norm_value)
NROW(test3)
###


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
  pivot_wider(names_from = sample, values_from = norm_value, values_fill = list(norm_value = 0))

feat <- as.data.frame(df3)
which(is.na(feat[,1]))

rownames(feat) <- feat[,1]
feat[,1] <- NULL

head(feat)
dim(feat)

# is the sum of each columns 1? 
colSums(feat)
# yes 

# ready! 

######################################################################

# CREATE METADATA TABLE (like meta.crc.zeller)

theseAREtheSamples <- as.data.frame(colnames(feat))
colnames(theseAREtheSamples) <- "sample"

df1$sample <- paste0(df1$date,"_",df1$pig)

df1 <- df1 %>%
  dplyr::select(sample,cohort,pig,date,breed,birth_day) %>%
  distinct()

# add anther grouping: date+cohort:
df1$group <- paste0(df1$date,"_",df1$cohort)

# add another grouping: date+breed+birth_day:
df1$group2 <- paste0(df1$date,"_",df1$breed,"_",df1$birth_day)

head(theseAREtheSamples)

meta <- left_join(theseAREtheSamples,df1)

rownames(meta) <- meta[,1]
meta[,1] <- NULL

# ready! 

######################################################################
######################################################################

# SIAMCAT starts! 

class(feat)
class(feat.crc.zeller)
class(meta)
class(meta.crc.zeller)

head(feat)
head(feat.crc.zeller)
head(meta)
head(meta.crc.zeller)

colnames(feat)==rownames(meta)


######################################################################
######################################################################

# create function to compare time points 

# prepare empty df to be filled
empty_df <- data.frame(
  fc = numeric(),
  p.val = numeric(),
  auc = character(),
  auc.ci.l = character(),
  auc.ci.h = character(),
  pr.shift = character(),
  pr.n = character(),
  pr.p = character(),
  bcol = character(),
  p.adj = numeric(),
  comparison = character(),
  species = character(),
  stringsAsFactors = FALSE
)

# save (now empty; to fill inside the loop)
fwrite(x=empty_df, file = paste0(out_dir_git,"gt_siamcat_time.csv"), sep = ",", append = FALSE)

comparetimepoints_full <- function(t1,t2) { # where c is the cohort of interest, t0 is the start time point, t1 is the next time point)
  
  label.normalized <- create.label(meta=meta,
                                   label='date', 
                                   case= paste0(t2),
                                   control= paste0(t1))
  
  siamcat <- siamcat(feat=feat,label=label.normalized,meta=meta)
  siamcat <- filter.features(siamcat,filter.method = 'abundance',cutoff = 0.0001)
  
  # check for significant associations
  siamcat <- check.associations(
    siamcat,
    fn.plot = paste0(out_dir,"gt_siamcatA_",t1,t2,".pdf"),
    sort.by = 'fc',
    alpha = 0.05, 
    mult.corr = "fdr",
    detect.lim = 10 ^-6,
    plot.type = "quantile.box",
    panels = c("fc", "auroc"))
  
  
  # save the data (significantly associated hits)
  
  mydata <- associations(siamcat,verbose=1)
  mydata$comparison <- paste0(t1,"_",t2)
  mydata$species <- rownames(mydata)
  rownames(mydata) <- NULL
  
  fwrite(x=mydata, file=paste0(out_dir_git,"gt_siamcat_time.csv"), sep = ",",
         append = TRUE)
  
  mydata_sub50 <- mydata %>%
    dplyr::arrange(p.adj) %>%
    slice(1:50) %>%
    arrange(desc(fc)) %>%
    dplyr::select(fc,p.adj)
  
  fwrite(x=mydata_sub50, file=paste0(out_dir_git,"gt_siamcat_time_sub50.csv"), sep = ",",
         append = TRUE)
  
  # Model building
  
  siamcat <- normalize.features(
    siamcat,
    norm.method = "log.unit",
    norm.param = list(
      log.n0 = 1e-06, 
      n.p = 2,
      norm.margin = 1
    )
  )
  
  siamcat <-  create.data.split(
    siamcat,
    num.folds = 5,
    num.resample = 2
  )
  
  siamcat <- train.model(
    siamcat,
    method = "lasso"
  )
  
  siamcat <- make.predictions(siamcat)
  pred_matrix <- pred_matrix(siamcat)
  siamcat <-  evaluate.predictions(siamcat)
  #model.evaluation.plot(siamcat)
  
  model.interpretation.plot(
    siamcat,
    fn.plot = paste0(out_dir,"gt_siamcatH_",t1,t2,".pdf"),
    consens.thres = 0.5, 
    limits = c(-3, 3),
    heatmap.type = 'zscore')
  
}


# run the timepoint comparisons you want: (outputs are two plots each, automatically saved)
# 1 week interval: 
comparetimepoints_full("t0","t2")
comparetimepoints_full("t2","t4")
comparetimepoints_full("t4","t6")
comparetimepoints_full("t6","t8")
comparetimepoints_full("t8","t10")
# 2 weeks interval:
comparetimepoints_full("t0","t4")
comparetimepoints_full("t4","t8")
# 4 weeks interval:
comparetimepoints_full("t0","t8")


############################################################################################################################################
############################################################################################################################################

# retrieve the associations we just obtained 
TimeAssociations <- read_csv(paste0(out_dir_git,"gt_siamcat_time.csv"))

# extract the significant hits
significant_with_time <- TimeAssociations %>%
  dplyr::filter(comparison=="t0_t2"|comparison=="t2_t4"|comparison=="t4_t6"|comparison=="t6_t8"|comparison=="t8_t10") %>%
  dplyr::filter(p.adj<=0.05) 

# save
fwrite(x=significant_with_time, file=paste0(out_dir,"gt_siamcat_time_sign.csv"), sep = ",",
       append = FALSE)

###############################################################################

# what species show signif associations more than once ? 
often_found_associating_with_time <- significant_with_time %>%
  group_by(species) %>%
  tally() %>%
  dplyr::arrange(desc(n))

fwrite(x=often_found_associating_with_time, file=paste0(out_dir,"gt_siamcat_time_sign_species.csv"), sep = ",",
       append = FALSE)

###############################################################################

# how many associations found per each time intervals comparison
sink(file = paste0(out_dir,"gt_siamcat_time_sign.txt"), 
     append = FALSE, type = c("output"))
TimeAssociations %>%
  dplyr::filter(comparison=="t0_t2"|comparison=="t2_t4"|comparison=="t4_t6"|comparison=="t6_t8"|comparison=="t8_t10") %>%
  dplyr::filter(p.adj<=0.05) %>%
  group_by(comparison) %>%
  tally() %>%
  dplyr::mutate(perc=n/sum(n)*100)
sink()

############################################################################################################################################
############################################################################################################################################


###### majorly shifting species with time - different plotting (1 plot for all sign shifts at any time interval)

gt_siamcat_time_sign <- read_csv(paste0(out_dir,"gt_siamcat_time_sign.csv"))

# some filtering and parse
gt_siamcat_time_sign$comparison <- as.factor(gt_siamcat_time_sign$comparison)

z <- gt_siamcat_time_sign %>% 
  dplyr::filter(p.adj < 0.05)

# split to get the node number 
z <- cSplit(z, "species", "__", drop = FALSE)

# # function to plot groups
# plot_genera_groups <- function(df) {
#   print(ggplot(df,aes(y=fct_reorder(species, species_2),x=fc,color=comparison))+
#           geom_point(alpha=0.8) +
#           geom_vline(xintercept = 0, linetype="dashed", 
#                      color = "black", size=0.5)+
#           scale_color_discrete(drop=FALSE) +
#           labs(color="interval")+
#           theme(axis.title.x=element_blank(),
#                 axis.text.y=element_text(size=6), 
#                 axis.title.y=element_blank()))
# }

# function to plot all sign species
plot_all <- function(df) {
  print(ggplot(df,aes(y=fct_reorder(species, species_2),x=fc,color=comparison))+
          geom_point(alpha=0.8, size=1) +
          geom_vline(xintercept = 0, linetype="dashed", 
                     color = "black", size=0.5)+
          scale_color_discrete(drop=FALSE) +
          labs(color="interval")+
          theme(axis.title.x=element_blank(),
                axis.text.y=element_text(size=3), 
                axis.title.y=element_blank()))
}


# CAG_plots <- plot_genera_groups(subset(z, grepl("^CAG", z$species))) # CAG : co-abundance genomes (Lesker et al, 2020)
# UBA_plots <- plot_genera_groups(subset(z, grepl("^UBA", z$species))) # UBA : uncultured bacteria and archaea 
# Prevo_plots <- plot_genera_groups(subset(z, grepl("^Prev", z$species)))
# Agath_plots <- plot_genera_groups(subset(z, grepl("^Agath", z$species)))
all_species_plots <- plot_all(z)

# all <- ggarrange(CAG_plots,
#                  UBA_plots,
#                  Prevo_plots,labels=c("A","B","C","D"),align="hv",
#                  Agat_plots, common.legend=TRUE)

pdf(paste0(out_dir,"gt_siamcat_time_sign_major.pdf"))
# all
all_species_plots
dev.off()

#######

# positive versus negative fold shifts plot at each time point: 

setDT(z)    
z <- z[, ":="(positive = sum(fc > 0), negative = sum(fc < 0)), by = comparison]


fold_changes_plot <- z %>%
  dplyr::select(comparison,positive,negative) %>%
  dplyr::distinct() %>%
  pivot_longer(., cols=c("negative","positive")) %>%
  dplyr::mutate(name = factor(name, levels = c("positive", "negative"))) %>%
  group_by(comparison) %>%
  dplyr::mutate(prop=paste0(round(value/sum(value)*100),"%")) %>%
  ggplot(., aes(x=comparison,y=value,fill=name))+
  geom_bar(stat="identity", position="dodge")+
  theme_minimal()+
  labs(fill="fold change",
       y="significant abundance shifts (counts)",
       x="time interval") +
  theme(legend.position = c(0.8, 0.9))+
  geom_text(aes(label=prop), position=position_dodge(width=0.9), vjust=-0.25, size=3)


pdf(paste0(out_dir,"gt_siamcat_time_pos_vs_neg.pdf"))
fold_changes_plot
dev.off()


############################################################################################################################################
############################################################################################################################################


# comparing time point (top 10 hits; for figures in main manuscript)

comparetimepoints_mini <- function(t1,t2) { 
  
  label.normalized <- create.label(meta=meta,
                                   label='date', 
                                   case= paste0(t2),
                                   control= paste0(t1))
  
  siamcat <- siamcat(feat=feat,label=label.normalized,meta=meta)
  siamcat <- filter.features(siamcat,filter.method = 'abundance',cutoff = 0.0001)
  
  # check for significant associations
  siamcat <- check.associations(
    siamcat,
    fn.plot = paste0(out_dir,"mini_gt_siamcatA_",t1,t2,".pdf"),
    sort.by = 'fc', 
    alpha = 0.05, max.show = 10,
    mult.corr = "fdr",
    detect.lim = 10 ^-6,
    plot.type = "quantile.box",
    panels = c("fc", "auroc"))
  
  
}

# run the timepoint comparisons you want: (outputs are two plots each, automatically saved)
# 1 week interval: 
comparetimepoints_mini("t0","t2")
comparetimepoints_mini("t2","t4")
comparetimepoints_mini("t4","t6")
comparetimepoints_mini("t8","t10")
# 2 weeks interval:
comparetimepoints_mini("t0","t4")
comparetimepoints_mini("t4","t8")
# 4 weeks interval:
comparetimepoints_mini("t0","t8")

############################################################################################################################################
############################################################################################################################################


# comparing piglets from the same breed, two birthday groups


# t0

label.normalized <- create.label(meta=meta,
                                 label='group2', 
                                 case= "t0_Duroc x Landrace_2017-01-11",
                                 control= "t0_Duroc x Landrace_2017-01-08")

siamcat <- siamcat(feat=feat,label=label.normalized,meta=meta)
siamcat <- filter.features(siamcat, filter.method = 'abundance',cutoff = 0.001)

# check for significant associations
siamcat <- check.associations(
  siamcat,
  sort.by = 'fc',
  fn.plot = paste0(out_dir,"gt_siamcatA_age_","t0_Duroc x Landrace.pdf"),
  alpha = 0.07, #0.6
  mult.corr = "fdr",
  detect.lim = 10 ^-30,
  prompt = FALSE,
  panels = c("fc", "auroc"))

# save the data (significantly associated hits with bday - t0)
mydata <- associations(siamcat,verbose=1)
mydata$comparison <- paste0("breed_DxL_bday08vs11_t0")
mydata$species <- rownames(mydata)
rownames(mydata) <- NULL

# t2

label.normalized <- create.label(meta=meta,
                                 label='group2', 
                                 case= "t2_Duroc x Landrace_2017-01-11",
                                 control= "t2_Duroc x Landrace_2017-01-08")

siamcat <- siamcat(feat=feat,label=label.normalized,meta=meta)
siamcat <- filter.features(siamcat,filter.method = 'abundance',cutoff = 0.001)

# check for significant associations
siamcat <- check.associations(
  siamcat,
  sort.by = 'fc',
  fn.plot = paste0(out_dir,"gt_siamcatA_age_","t2_Duroc x Landrace.pdf"),
  alpha = 0.16, #0.13
  mult.corr = "fdr",
  detect.lim = 10 ^-30,
  prompt = FALSE,
  panels = c("fc", "auroc"))

# save the data (significantly associated hits with bday - t2)
mydata2 <- associations(siamcat,verbose=1)
mydata2$comparison <- paste0("breed_DxL_bday08vs11_t2")
mydata2$species <- rownames(mydata2)
rownames(mydata2) <- NULL


both <- rbind(mydata,mydata2)
fwrite(x=both, file=paste0(out_dir_git,"gt_siamcat_age.csv"), sep = ",")


############################################################################################################################################
############################################################################################################################################

# comparing piglets from the same birth day (2017-01-08), two breeds (DxL vs DxLW)

# t0

label.normalized <- create.label(meta=meta,
                                 label='group2', 
                                 case= "t0_Duroc x Landrace_2017-01-08",
                                 control= "t0_Duroc x Large white_2017-01-08")

siamcat <- siamcat(feat=feat,label=label.normalized,meta=meta)
siamcat <- filter.features(siamcat, filter.method = 'abundance',cutoff = 0.001)

# check for significant associations
siamcat <- check.associations(
  siamcat,
  sort.by = 'fc',
  fn.plot = paste0(out_dir,"gt_siamcatA_breed_","t0_bday08_DxL_vs_DxLW.pdf"),
  alpha = 0.2,
  mult.corr = "fdr",
  detect.lim = 10 ^-30,
  prompt = FALSE,
  panels = c("fc", "auroc"))

# save the data (significantly associated hits with bday - t0)
mydata <- associations(siamcat,verbose=1)
mydata$comparison <- paste0("t0_bday08_DxL_vs_DxLW")
mydata$species <- rownames(mydata)
rownames(mydata) <- NULL

# t2

label.normalized <- create.label(meta=meta,
                                 label='group2', 
                                 case= "t2_Duroc x Landrace_2017-01-08",
                                 control= "t2_Duroc x Large white_2017-01-08")

siamcat <- siamcat(feat=feat,label=label.normalized,meta=meta)
siamcat <- filter.features(siamcat, filter.method = 'abundance',cutoff = 0.001)

# check for significant associations
siamcat <- check.associations(
  siamcat,
  sort.by = 'fc',
  fn.plot = paste0(out_dir,"gt_siamcatA_breed_","t2_bday08_DxL_vs_DxLW.pdf"),
  alpha = 0.2,
  mult.corr = "fdr",
  detect.lim = 10 ^-30,
  prompt = FALSE,
  panels = c("fc", "auroc"))

# save the data (significantly associated hits with bday - t0)
mydata2 <- associations(siamcat,verbose=1)
mydata2$comparison <- paste0("t2_bday08_DxL_vs_DxLW")
mydata2$species <- rownames(mydata2)
rownames(mydata2) <- NULL

both <- rbind(mydata,mydata2)
fwrite(x=both, file=paste0(out_dir_git,"gt_siamcat_breed.csv"), sep = ",")


############################################################################################################################################
############################################################################################################################################

# as Amy suggested to add confidence intervals for major shifts mentioned in manuscript ... 


df4 <- cSplit(df3, "gOTU", "__", drop = FALSE)

keep <- c("Blautia_A wexlerae",
          "CAG-83 sp003487665",
          "Lactobacillus amylovorus", 
          "CAG-45 sp002299665",
          "CAG-110 sp002437585",
          "CAG-110 sp002437585 ",
          "Methanobrevibacter_A smithii",
          "Phil1 sp001940855",
          "Prevotella copri",
          "Prevotella sp000434515",
          "Lactobacillus amylovorus",
          "Corynebacterium xerosis",
          "Blautia_A obeum",
          "Blautia_A sp000285855",
          "Clostridium sp000435835",
          "Clostridium_P ventriculi",
          "UBA7748 sp900314535",
          "Clostridium sp000435835",
          "Prevotella copri",
          "Bifidobacterium boum",
          "Clostridium sp000435835",
          "Corynebacterium xerosis",
          "Prevotella copri_A",
          "Prevotella copri",
          "Prevotella sp000434515")

keep <- data.frame(keep)
colnames(keep) <- "gOTU_1"

df_CI <-inner_join(keep,df4)

df_CI <- df_CI %>%
  group_by(gOTU,date) %>%
  summarise(lowCI = ci(norm_value)[2],
            hiCI = ci(norm_value)[3])

fwrite(x=df_CI, file=paste0(out_dir_git,"gt_siamcat_time_major_species_shifts.csv"), sep = ",")


