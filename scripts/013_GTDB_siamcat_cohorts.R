
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

source_dir = "/Users/12705859/Desktop/metapigs_dry/source_data/" # git 
middle_dir = "/Users/12705859/Desktop/metapigs_dry/middle_dir/" # git 
out_dir_git = "/Users/12705859/Desktop/metapigs_dry/out/" # git 
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



# comparing cohorts 


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

# save it
fwrite(x=empty_df, file=paste0(out_dir,"gt_siamcat_cohortsbetween.csv"), sep = ",",
       append = FALSE)

######################################################################

# function to compare cohorts: 
comparebetween <- function(t,c2,c1) { 
  
  label.normalized <- create.label(meta=meta,
                                   label='group', 
                                   case= paste0(t,"_",c2),
                                   control= paste0(t,"_",c1))
  
  siamcat <- siamcat(feat=feat,label=label.normalized,meta=meta)
  siamcat <- filter.features(siamcat,filter.method = 'abundance',cutoff = 0.0001)
  
  # check for significant associations
  siamcat <- check.associations(
    siamcat,
    sort.by = 'fc',
    alpha = 0.05,
    fn.plot = paste0(out_dir,"gt_siamcatA_",t,"_",c1,"_vs_",c2,".pdf"),
    mult.corr = "fdr",
    detect.lim = 10 ^-6)
  
  # save the data 
  mydata_cohorts_between <- associations(siamcat,verbose=1)
  mydata_cohorts_between$comparison <- paste0(t,"_",c1,"_vs_",c2)
  mydata_cohorts_between$species <- rownames(mydata_cohorts_between)
  rownames(mydata_cohorts_between) <- NULL
  
  fwrite(x=mydata_cohorts_between, file=paste0(out_dir,"gt_siamcat_cohortsbetween.csv"), sep = ",",
         append = TRUE)
  
}


# run the association tests and plot the output
comparebetween("t0","DScour","Control")
comparebetween("t2","DScour","Control")
comparebetween("t4","DScour","Control")
comparebetween("t6","DScour","Control")
comparebetween("t8","DScour","Control")
comparebetween("t10","DScour","Control")
comparebetween("t0","ColiGuard","Control")
comparebetween("t2","ColiGuard","Control")
comparebetween("t4","ColiGuard","Control")
comparebetween("t6","ColiGuard","Control")
comparebetween("t8","ColiGuard","Control")
comparebetween("t10","ColiGuard","Control")
comparebetween("t0","Neomycin","Control")
comparebetween("t2","Neomycin","Control")
comparebetween("t4","Neomycin","Control")
comparebetween("t6","Neomycin","Control")
comparebetween("t8","Neomycin","Control")
comparebetween("t10","Neomycin","Control")
comparebetween("t0","NeoD","Neomycin")
comparebetween("t2","NeoD","Neomycin")
comparebetween("t4","NeoD","Neomycin")
comparebetween("t6","NeoD","Neomycin")
comparebetween("t8","NeoD","Neomycin")
comparebetween("t10","NeoD","Neomycin")
comparebetween("t0","NeoC","Neomycin")
comparebetween("t2","NeoC","Neomycin")
comparebetween("t4","NeoC","Neomycin")
comparebetween("t6","NeoC","Neomycin")
comparebetween("t8","NeoC","Neomycin")
comparebetween("t10","NeoC","Neomycin")


######################################################################

# retrieve the significant hits from the data we just obtained 
CohortsAssociations <- read_csv(paste0(out_dir,"gt_siamcat_cohortsbetween.csv"))

significant_between_cohorts <- CohortsAssociations %>%
  dplyr::filter(p.adj<=0.05) 

fwrite(x=significant_between_cohorts, file=paste0(out_dir_git,"gt_siamcat_cohortsbetween_sign.csv"), sep = ",",
       append = FALSE)

# how many associations found per each time intervals comparison
sink(file = paste0(out_dir,"gt_siamcat_cohortsbetween_sign.txt"), 
     append = FALSE, type = c("output"))
significant_between_cohorts %>%
  dplyr::filter(p.adj<=0.05) %>%
  group_by(comparison) %>%
  tally() %>%
  dplyr::mutate(perc=n/sum(n)*100)
sink()

######################################################################
######################################################################
# 
# # Function to plot the significant changes per cohort for each time interval : 
# 
# comparewithin_plots <- function(c,t1,t2) { # where c is the cohort of interest, t0 is the start time point, t1 is the next time point)
#   
#   label.normalized <- create.label(meta=meta,
#                                    label='group', 
#                                    case= paste0(t2,"_",c),
#                                    control= paste0(t1,"_",c))
#   
#   siamcat <- siamcat(feat=feat,label=label.normalized,meta=meta)
#   siamcat <- filter.features(siamcat,filter.method = 'abundance',cutoff = 0.0001)
#   
#   
#   # check for significant associations
#   siamcat <- check.associations(
#     siamcat,
#     fn.plot = paste0("gt_siamcatA_",c,"_",t1,"_",t2,".pdf"),
#     sort.by = 'fc',
#     alpha = 0.05, 
#     mult.corr = "fdr",
#     detect.lim = 10 ^-6,
#     plot.type = "quantile.box",
#     panels = c("fc", "prevalence", "auroc"))
#   
#   
# }
# 
# 
# # Plot : 
# 
# comparewithin_plots("Control","t0","t2")
# comparewithin_plots("Control","t2","t4")
# comparewithin_plots("Control","t4","t6")
# comparewithin_plots("Control","t6","t8")
# comparewithin_plots("Control","t8","t10")
# comparewithin_plots("Control","t4","t8")
# comparewithin_plots("Control","t4","t10")
# comparewithin_plots("Control","t2","t8")
# #
# comparewithin_plots("DScour","t0","t2")
# comparewithin_plots("DScour","t2","t4")
# comparewithin_plots("DScour","t4","t6")
# comparewithin_plots("DScour","t6","t8")
# comparewithin_plots("DScour","t8","t10")
# comparewithin_plots("DScour","t4","t8")
# comparewithin_plots("DScour","t4","t10")
# comparewithin_plots("DScour","t2","t8")
# #
# comparewithin_plots("ColiGuard","t0","t2")
# comparewithin_plots("ColiGuard","t2","t4")
# comparewithin_plots("ColiGuard","t4","t6")
# comparewithin_plots("ColiGuard","t6","t8")
# comparewithin_plots("ColiGuard","t8","t10")
# comparewithin_plots("ColiGuard","t4","t8")
# comparewithin_plots("ColiGuard","t4","t10")
# comparewithin_plots("ColiGuard","t2","t8")
# #
# comparewithin_plots("Neomycin","t0","t2")
# comparewithin_plots("Neomycin","t2","t4")
# comparewithin_plots("Neomycin","t4","t6")
# comparewithin_plots("Neomycin","t6","t8")
# comparewithin_plots("Neomycin","t8","t10")
# comparewithin_plots("Neomycin","t4","t8")
# comparewithin_plots("Neomycin","t4","t10")
# comparewithin_plots("Neomycin","t2","t8")
# #
# comparewithin_plots("NeoD","t0","t2")
# comparewithin_plots("NeoD","t2","t4")
# comparewithin_plots("NeoD","t4","t6")
# comparewithin_plots("NeoD","t6","t8")
# comparewithin_plots("NeoD","t8","t10")
# comparewithin_plots("NeoD","t4","t8")
# comparewithin_plots("NeoD","t4","t10")
# comparewithin_plots("NeoD","t2","t8")
# #
# comparewithin_plots("NeoC","t0","t2")
# comparewithin_plots("NeoC","t2","t4")
# comparewithin_plots("NeoC","t4","t6")
# comparewithin_plots("NeoC","t6","t8")
# comparewithin_plots("NeoC","t8","t10")
# comparewithin_plots("NeoC","t4","t8")
# comparewithin_plots("NeoC","t4","t10")
# comparewithin_plots("NeoC","t2","t8")


######################################################################
# 
# # Function to obtain the (data only) significant changes per cohort for each time interval : 
# 
# dfz <- data.frame(
#   fc = numeric(),
#   p.adj = numeric(),
#   gOTU = character(),
#   comparison = character(),
#   stringsAsFactors = FALSE
# )
# 
# comparewithin <- function(c,t1,t2) { # where c is the cohort of interest, t0 is the start time point, t1 is the next time point)
#   
#   label.normalized <- create.label(meta=meta,
#                                    label='group', 
#                                    case= paste0(t2,"_",c),
#                                    control= paste0(t1,"_",c))
#   
#   siamcat <- siamcat(feat=feat,label=label.normalized,meta=meta)
#   siamcat <- filter.features(siamcat,filter.method = 'abundance',cutoff = 0.0001)
#   
#   # check for significant associations
#   siamcat <- check.associations(
#     siamcat,
#     sort.by = 'fc',
#     alpha = 0.05,
#     mult.corr = "fdr",
#     detect.lim = 10 ^-6,
#     prompt = FALSE,
#     panels = c("fc", "prevalence", "auroc"))
#   
#   yy <- SIAMCAT::associations(siamcat)
#   yy$gOTU <- rownames(yy)
#   yy <- yy %>% filter(p.adj<0.05)
#   
#   
#   if (NROW(yy)>0) {
#     yy$comparison=paste0(c,"_",t1,"_",t2)
#     yy$feat <- rownames(yy)
#     yy <- yy %>%
#       dplyr::select(comparison,gOTU,fc,p.adj)
#     dfz <- rbind(dfz,yy)
#   } else {
#     print ("No significant associations found. No plot will be produced.")
#   }
#   
#   return(dfz)
# }
# 
# 
# 
# last <- rbind(comparewithin("Control","t0","t2"),
#               comparewithin("Control","t2","t4"),
#               comparewithin("Control","t4","t6"),
#               comparewithin("Control","t6","t8"),
#               comparewithin("Control","t8","t10"),
#               comparewithin("Control","t4","t8"),
#               comparewithin("Control","t4","t10"),
#               comparewithin("Control","t2","t8"),
#               #
#               comparewithin("DScour","t0","t2"),
#               comparewithin("DScour","t2","t4"),
#               comparewithin("DScour","t4","t6"),
#               comparewithin("DScour","t6","t8"),
#               comparewithin("DScour","t8","t10"),
#               comparewithin("DScour","t4","t8"),
#               comparewithin("DScour","t4","t10"),
#               comparewithin("DScour","t2","t8"),
#               #
#               comparewithin("ColiGuard","t0","t2"),
#               comparewithin("ColiGuard","t2","t4"),
#               comparewithin("ColiGuard","t4","t6"),
#               comparewithin("ColiGuard","t6","t8"),
#               comparewithin("ColiGuard","t8","t10"),
#               comparewithin("ColiGuard","t4","t8"),
#               comparewithin("ColiGuard","t4","t10"),
#               comparewithin("ColiGuard","t2","t8"),
#               #
#               comparewithin("Neomycin","t0","t2"),
#               comparewithin("Neomycin","t2","t4"),
#               comparewithin("Neomycin","t4","t6"),
#               comparewithin("Neomycin","t6","t8"),
#               comparewithin("Neomycin","t8","t10"),
#               comparewithin("Neomycin","t4","t8"),
#               comparewithin("Neomycin","t4","t10"),
#               comparewithin("Neomycin","t2","t8"),
#               #
#               comparewithin("NeoD","t0","t2"),
#               comparewithin("NeoD","t2","t4"),
#               comparewithin("NeoD","t4","t6"),
#               comparewithin("NeoD","t6","t8"),
#               comparewithin("NeoD","t8","t10"),
#               comparewithin("NeoD","t4","t8"),
#               comparewithin("NeoD","t4","t10"),
#               comparewithin("NeoD","t2","t8"),
#               #
#               comparewithin("NeoC","t0","t2"),
#               comparewithin("NeoC","t2","t4"),
#               comparewithin("NeoC","t4","t6"),
#               comparewithin("NeoC","t6","t8"),
#               comparewithin("NeoC","t8","t10"),
#               comparewithin("NeoC","t4","t8"),
#               comparewithin("NeoC","t4","t10"),
#               comparewithin("NeoC","t2","t8"))
# 
# last_to_use <- last 
# 
# NROW(last_to_use)
# last_to_use$p.adj <- NULL
# 
# toplot <- last_to_use %>%
#   dplyr::select(gOTU,comparison,fc)
# 
# 
# 
# # save hits as text file 
# sink(file = "gt_siamcatA_cohortswithin.txt", 
#      append = TRUE, type = c("output"))
# last_to_use
# sink()
# 
# 
# 
# toplot <- cSplit(toplot, "comparison", "_")
# toplot$interval <- paste0(toplot$comparison_2,"_",toplot$comparison_3)
# 
# toplot <- toplot %>%
#   dplyr::select(gOTU,comparison_1,interval,fc)
# 
# colnames(toplot) <- c("gOTU","cohort","interval","fc")
# 
# 
# # I am adding this here because the pdf below often doesn t get printed
# closeAllConnections()
# 
# 
# # split df by time interval 
# multiple_DFs <- split( toplot , f = toplot$interval )
# 
# 
# pdf("gt_siamcatA_cohortswithin.pdf")
# for (single_DF in multiple_DFs) {
#   
#   
#   DF <- as.data.frame(single_DF)
#   
#   interval <- DF$interval[1]
#   
#   DF_wide <- DF %>%
#     pivot_wider(names_from = cohort, values_from = fc, values_fill = list(fc = 0))
#   
#   if (NCOL(DF_wide) > 1 & NROW(DF_wide) >1 ) {
#     
#     DF_wide <- as.data.frame(DF_wide)
#     
#     rownames(DF_wide) <- DF_wide[,1]
#     DF_wide[,1] <- NULL
#     DF_wide$interval <- NULL
#     
#     m <- as.matrix(DF_wide)
#     
#     pheatmap(m, fontsize_row = 5, main = interval, display_numbers = TRUE, cluster_cols = FALSE)
#     
#   }
#   
#   
# }
# dev.off()
# 



