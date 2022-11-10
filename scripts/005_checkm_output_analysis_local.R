

# upload libraries
library(tidyverse)
library(gplots)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(splitstackshape)
library(treemap)
library(dplyr)
library(data.table)
library(readxl)
library(pheatmap)
#library(robCompositions)
library(compositions) # <- this is the good one for clr transform
library(ggbiplot)
library(EnvStats)
library(treemapify)
library(scales)
library(compositions)

middle_dir = "/Users/danielagaio/Gaio/github/metapigs_dry/middle_dir/" # git 
out_dir = "/Users/danielagaio/Desktop/metapigs_dry/checkm/"  # local

# input files (in middle_dir)
# all_checkm_output.tsv 
# checkm/checkm_all_nearly
# no_reps_all.csv (BINS COUNTS)

# OUTPUTS:

# about bins in general :
#  - numbers.txt

# from all_checkm_output :
#  - cm_Compl_vs_contam.pdf
#  - cm_numbers.txt
#  - cm_contigs_distribution.pdf
#  - cm_scaffolds_predicted_genes.pdf

# from BINS COUNTS + checkm_all_nearly: 
# based on bins frequency: 
## - cm_treemap_phylum.pdf
## - cm_treemap_class.pdf
## - cm_treemap_order.pdf
## - numbers.txt
# based on (counts) log10 relative abundance: 
## - cm_parallel_coordinates_phyla.pdf
# based on (counts) relative abundance: 
## - cm_rel_ab_phyla.pdf
## - cm_rel_ab_phyla_cohorts.pdf
## - cm_balloonplot_phyla.pdf
## - cm_balloonplot_phyla_cohorts.pdf
# based on (counts) cenLR relative abundance: 
## - cm_PCA.pdf



######################################################################

# template to collect and store bins info 

sink(paste0(out_dir,"numbers.txt"))
start_message <- " ########################## BINS NUMBERS ########################## "
start_message
sink()

# template to collect and store checkM info 

sink(paste0(out_dir,"cm_numbers.txt"))
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

no_reps_all <- read.csv(paste0(middle_dir,"no_reps_all.csv"), 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)

# remove .fa extension to match bins in checkm df 
no_reps_all$bin <- gsub(".fa","", no_reps_all$bin)
head(no_reps_all)


######################################################################


# NUMBERS 


##########################

# Bins distribution:
# how many bins each subject has:
az <- no_reps_all %>%
  dplyr::select(bin, pig, cohort) %>%
  distinct() 

az <- az %>%
  dplyr::mutate(group = ifelse(cohort == "Mothers", "mothers (n=42)","piglets (n=126)"))

az <- az %>%
  group_by(group,pig) %>%
  dplyr::summarise(`number of metagenomes obtained`= n()) 
tail(az)

pdf(paste0(out_dir,"#bins_subject.pdf"))
ggplot(data=az, mapping=aes(x=group, y=`number of metagenomes obtained`)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1) +
  theme_bw()
dev.off()


##########################

# Bins per cohort: 

# calculate number of MAGs obtaiend per subject
az2 <- no_reps_all %>%
  dplyr::select(date,bin,pig) %>%
  group_by(date,pig) %>%
  dplyr::summarise(`number of metagenomes obtained`= n()) 
tail(az2)
az2$date <- NULL

az2 <- az2 %>%
  group_by(pig,`number of metagenomes obtained`) %>%
  dplyr::summarise(mean_bins= mean(`number of metagenomes obtained`),
                   median_bins= median(`number of metagenomes obtained`),
                   tot_bins= sum(`number of metagenomes obtained`))

# TIMEPOINTS PER SUBJ
az3 <- no_reps_all %>%
  dplyr::select(date,pig,cohort) %>%
  distinct() %>%
  group_by(pig,cohort) %>%
  dplyr::summarise(sampling_frequency= n())
tail(az3)


final_pig <- left_join(az2,az3)

# reorder cohorts 
final_pig$cohort  = factor(final_pig$cohort, levels=c("Control", 
                                                      "D-Scour",
                                                      "ColiGuard", 
                                                      "Neomycin",
                                                      "NeoD",
                                                      "NeoC",
                                                      "Mothers"))


# mean - Bins per cohort: 
final_cohort <- final_pig %>%
  group_by(cohort) %>%
  dplyr::summarise(mean_bins= mean(`number of metagenomes obtained`),
                   median_bins= median(`number of metagenomes obtained`),
                   tot_bins= sum(`number of metagenomes obtained`),
                   mean_sampling_frequency=mean(sampling_frequency),
                   median_sampling_frequency=median(sampling_frequency))

sink(file = paste0(out_dir,"numbers.txt"), 
     append = TRUE, type = c("output"))
paste0("Total bins and sampling frequency per cohort ")
as.data.frame(final_cohort)
sink()

# means as labels in plot
tot <- final_cohort %>%
  dplyr::select(cohort,tot_bins)
labs <- round(tot$tot_bins)

bins_per_cohort_plot <- final_pig %>%
  ggplot(data=., mapping=aes(x=cohort, y=`number of metagenomes obtained`, color=sampling_frequency)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, size=0.5,aes(colour = sampling_frequency)) +
  stat_n_text(size = 3) +
  theme_bw() +
  theme(legend.position="top",
        axis.text.x=element_text(angle=90))+
  scale_color_gradient(low = "blue", 
                       high = "green",
                       name = "sampling frequency")+
  guides(color=guide_colorbar(direction = "horizontal",
                              title.position = "top"))+
  annotate("text", label = labs[1], x = 1, y = 835, size = 3, colour = "black") +
  annotate("text", label = labs[2], x = 2, y = 835, size = 3, colour = "black") +
  annotate("text", label = labs[3], x = 3, y = 835, size = 3, colour = "black") +
  annotate("text", label = labs[4], x = 4, y = 835, size = 3, colour = "black") +
  annotate("text", label = labs[5], x = 5, y = 835, size = 3, colour = "black") +
  annotate("text", label = labs[6], x = 6, y = 835, size = 3, colour = "black") +
  annotate("text", label = labs[7], x = 7, y = 835, size = 3, colour = "black")


pdf(paste0(out_dir,"#bins_cohort.pdf"))
bins_per_cohort_plot 
dev.off()




##########################


# Get the numbers: 
az_piglets <- no_reps_all %>%
  filter(!cohort=="Mothers") %>%
  dplyr::select(bin,pig)
bz_piglets <- unique(az_piglets)
cz_piglets <- data.frame(table(bz_piglets$pig))

az_mothers <- no_reps_all %>%
  filter(cohort=="Mothers") %>%
  dplyr::select(bin,pig)
bz_mothers <- unique(az_mothers)
cz_mothers <- data.frame(table(bz_mothers$pig))


sink(file = paste0(out_dir,"numbers.txt"), 
     append = TRUE, type = c("output"))
paste0("total of bins from Piglets: ", sum(cz_piglets$Freq), " Percentage: ", round(sum(cz_piglets$Freq)/51170*100,2) )
paste0("total of bins from Mothers: ", sum(cz_mothers$Freq), " Percentage: ", round(sum(cz_mothers$Freq)/51170*100,2) )
paste0("total of bins from pos and neg controls: ", 51175-NROW(bz_piglets)-NROW(bz_mothers) , " Percentage: ", round((51175-NROW(bz_piglets)-NROW(bz_mothers))/51170*100,2) )
sink()


##########################


# number of samples per time point: 

az <- no_reps_all %>%
  filter(!cohort=="Mothers") %>%
  dplyr::select(pig, date) %>%
  distinct() 

az <- az %>%
  group_by(date) %>%
  dplyr::summarise(`number of samples`= n()) 
tail(az)

# reorder dates 
az$date  = factor(az$date, levels=c("t0",
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
az <- az %>%
  dplyr::mutate(sampling = ifelse(date == "t0" | date == "t2" | 
                                    date == "t4" | date == "t6" | 
                                    date == "t8" | date == "t10" , "all subjects", "subset"))


sample_per_timepoint_plot <- ggplot(data=az, mapping=aes(x=date, y=`number of samples`, fill=sampling)) + 
  geom_point(aes(group=sampling, shape=sampling),size=3) +
  geom_line(aes(group = sampling), linetype = 2) +
  theme_bw() + 
  theme(legend.position="top")+
  labs(x = "timepoint",
       y = "number of samples",
       subtitle=NULL)

pdf(paste0(out_dir,"#samples_per_timepoint.pdf"))
sample_per_timepoint_plot
dev.off()


##########################
# Correlation between timepoints available and bins obtained per subject: 
# How many timepoints scale to how many bins? 

# BINS PER SUBJ
az1 <- no_reps_all %>%
  dplyr::select(bin, pig) %>%
  distinct() %>%
  group_by(pig) %>%
  dplyr::summarise(`number of metagenomes per subject`= n()) 
tail(az1)
NROW(az1)
head(az1)

sink(file = paste0(out_dir,"numbers.txt"), 
     append = TRUE, type = c("output"))
paste0("Bins per subject: ")
summary(az1$`number of metagenomes per subject`)
sink()

# TIMEPOINTS PER SUBJ
az2 <- no_reps_all %>%
  dplyr::select(date, pig) %>%
  distinct() %>%
  group_by(pig) %>%
  dplyr::summarise(`number of timepoints per subject`= n()) 
tail(az2)
NROW(az2)

sink(file = paste0(out_dir,"numbers.txt"), 
     append = TRUE, type = c("output"))
paste0("Timepoints per subject: ")
summary(az2$`number of timepoints per subject`)
sink()


az3 <- cbind(az1,az2)
az3$pig <- NULL

az3$`number of timepoints per subject` <- as.numeric(az3$`number of timepoints per subject`)

pdf(paste0(out_dir,"#bins_vs_#timepoints.pdf"), height = 4, width = 7)
ggplot(data=az3, mapping=aes(x=`number of timepoints per subject`, y=`number of metagenomes per subject`)) + 
  geom_boxplot(outlier.shape=NA, inherit.aes = TRUE, aes(group=`number of timepoints per subject`))+
  geom_point(color='blue', size=0.6) +
  geom_smooth(method = "loess", se=TRUE, aes(group=1))+
  labs(x = "timepoints available per subject",
       y = "MAGs obtained per subject",
       subtitle=NULL) +
  theme_bw()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14))
dev.off()

###
# correlation between time points and number of bins obtained per subject 
# using spearman as variables are not normally distributed 
res <- cor.test(az3$`number of metagenomes per subject`, az3$`number of timepoints per subject`, 
                method = c("spearman")) #, "kendall", "spearman"
res$estimate^2
###

az4 <- az3 %>%
  group_by(`number of timepoints per subject`) %>%
  dplyr::summarise(`mean_bins`= mean(`number of metagenomes per subject`),
                   `median_bins`= median(`number of metagenomes per subject`)) 

sink(file = paste0(out_dir,"numbers.txt"), 
     append = TRUE, type = c("output"))
paste0("Correlation timepoint sample/ mean bins obtained: ")
as.data.frame(az4)
sink()


# TIMEPOINTS PER PIGGIES SUBJECT
az2 <- no_reps_all %>%
  dplyr::filter(!cohort=="Mothers") %>%
  dplyr::select(date, pig) %>%
  distinct() %>%
  group_by(pig) %>%
  dplyr::summarise(`number of timepoints per subject`= n()) 
tail(az2)
NROW(az2)

sink(file = paste0(out_dir,"numbers.txt"), 
     append = TRUE, type = c("output"))
paste0("Timepoints per subject (piggies): ")
summary(az2$`number of timepoints per subject`)
sink()

######################################################################

# bins distribution - cohorts & piglets

test <- no_reps_all %>%
  dplyr::filter(!cohort=="Mothers") %>%
  dplyr::select(pig,bin,cohort) %>%
  distinct() %>%
  group_by(cohort,pig) %>%
  tally()  %>%
  dplyr::mutate(perc=round(n/sum(n)*100,2)) 

# reorder cohorts 
test$cohort  = factor(test$cohort, levels=c("Control", 
                                            "D-Scour",
                                            "ColiGuard", 
                                            "Neomycin",
                                            "NeoD",
                                            "NeoC"))

treemap_ggplot <- ggplot(test, aes(area = perc, fill = cohort, label = pig,
                                   subgroup = cohort, subgroup2=pig)) +
  geom_treemap() +
  geom_treemap_subgroup_border() +
  geom_treemap_subgroup2_border(colour="black",size=1) +
  geom_treemap_subgroup_text(place = "centre", grow = T, alpha = .9, colour =
                               "White", fontface = "italic", min.size = 0) +
  #geom_treemap_text(colour = "black", place = "topleft", reflow = T)+
  theme(legend.position="none")+
  ggtitle("       ")

pdf(paste0(out_dir,"bins_distribution_over_cohorts_piglets.pdf"))
treemap_bins <- treemap(test, 
                        index=c("cohort","pig"), 
                        vSize="perc", 
                        type="index",                            # How you color the treemap. type help(treemap) for more info
                        palette = "Pastel1",                        # Select your color palette from the RColorBrewer presets or make your own.
                        fontsize.title=14,                       # Size of the title
                        title=paste0("Bins distribution over cohorts and piglets (",NROW(unique(paste0(no_reps_all$pig,no_reps_all$bin))),")")
)
dev.off()

# library(RColorBrewer)
# brewer.pal(n=6, "Pastel1")
# display.brewer.pal(n=6, "Pastel1")

######################################################################

bins_to_sampling_3_figures <- ggarrange(sample_per_timepoint_plot,
                                        bins_per_cohort_plot,
                                        treemap_ggplot,
                                        nrow=1,
                                        ncol=3,
                                        labels=c("A","B","C"))

#pdf("sample_effort_to_MAGs_plots.pdf",width = 9,height = 5)
#bins_to_sampling_3_figures
#dev.off()

# preparing data for Completeness vs Contamination dot plot

# nearly complete genomes 
df_90100 <- subset(all_checkm_output, Completeness>= 90 & Completeness <= 100)
df_90100_5 <- subset(df_90100, Contamination <=5)
df_90100_5$type = paste0(paste0("Nearly complete: ",
                                round(NROW(df_90100_5)/NROW(all_checkm_output)*100,2),
                                "%",
                                " (n=",
                                NROW(df_90100_5),
                                ")"))

# Medium quality
df_8090 <- subset(all_checkm_output, Completeness>= 80 & Completeness < 90)
df_8090_10 <- subset(df_8090, Contamination <=10)
# add completeness between 90 and 100, where contamination is between 5 and 10 
df_90100 <- subset(all_checkm_output, Completeness>= 90 & Completeness <= 100)
df_90100_510 <- subset(df_90100, Contamination >5 & Contamination <= 10)
medium_quality <- rbind(df_8090_10, df_90100_510)

medium_quality$type = paste0(paste0("Medium quality: ",
                                    round(NROW(medium_quality)/NROW(all_checkm_output)*100,2),
                                    "%",
                                    " (n=",
                                    NROW(medium_quality),
                                    ")"))


best <- rbind(df_90100_5,medium_quality)

# Assessment of genome quality
cm_Compl_vs_Contam <- ggplot(best, aes(x=Completeness, y=Contamination)) + 
  geom_point(aes(color=factor(type)),size=0.1,shape=21)+
  ylab("Contamination (%)")+
  xlab("Completeness (%)")+
  xlim(75,100)+
  theme_minimal()+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(size=11),
        axis.text.y=element_text(size=11),
        legend.title = element_blank(),
        legend.position = "none")

# similar to https://www.nature.com/articles/s41564-017-0012-7#Sec2


# library(scales)
# scales::hue_pal()(3)


# perfect bins
df_99100 <- subset(all_checkm_output, Completeness>= 99 & Completeness <= 100)
df_99100_01 <- subset(df_99100, Contamination <=0.1)


# perfect bins
df_70100 <- subset(all_checkm_output, Completeness>= 70 & Completeness <= 100)
df_70100_10 <- subset(df_70100, Contamination <=10)


sink(file = paste0(out_dir,"cm_numbers.txt"), 
     append = TRUE, type = c("output"))
paste0("Total number of bins is:   ",
       NROW(all_checkm_output))
paste0("Of which: >=70 <=10):   ",
       NROW(df_70100_10)," ", round(NROW(df_70100_10)/NROW(all_checkm_output),4)*100,"%")
paste0("Of which: nearly complete (>=90 <=5):   ",
       NROW(df_90100_5)," ", round(NROW(df_90100_5)/NROW(all_checkm_output),4)*100,"%")
paste0("Of which: =>99% completeness and <=0.1% contamination: ",
       NROW(df_99100_01)," ", round(NROW(df_99100_01)/NROW(all_checkm_output),4)*100,"%")
sink()


######################################################################


# Contig number and N50 distribution 

# as numeric 
df_90100_5$`N50 (contigs)` <- as.numeric(df_90100_5$`N50 (contigs)`)
df_90100_5$`# contigs` <- as.numeric(df_90100_5$`# contigs`)


sink(file = paste0(out_dir,"cm_numbers.txt"), 
     append = TRUE, type = c("output"))
paste0("Nearly complete genomes, contigs stats :   ")
summary(df_90100_5$`# contigs`)
paste0("Nearly complete genomes, contigs N50 stats :   ")
summary(df_90100_5$`N50 (contigs)`)
sink()


pdf(paste0(out_dir,"cm_contigs_distribution.pdf"), height=4, width=5)
# Distribution of contig number across bins 
ggplot(df_90100_5, 
       aes(x = `N50 (contigs)`)) + 
  geom_histogram(bins = 100, colour="black",fill="white") + 
  scale_x_log10() +
  labs(x = "N50 (log10)",
       y = "Frequency")+
  theme_bw()+
  theme(axis.title.x=element_text(size=14),
        axis.text.x=element_text(size=12),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=12), title = element_text(size=8))
dev.off()


# create another category: partial quality: 
# Completeness => 60 <80
# Contamination <= 10
df_6080 <- subset(all_checkm_output, Completeness>= 60 & Completeness < 80)
df_6080_10 <- subset(df_6080, Contamination <=10)
df_6080_10$type <- paste0("Partial quality: ",
                          round(NROW(df_6080_10)/NROW(all_checkm_output)*100,2),
                          "%",
                          " (n=",
                          NROW(df_6080_10),
                          ")")

toplot <- rbind(df_90100_5,medium_quality,df_6080_10)

toplot$type <- as.factor(toplot$type)
types <- as.character(levels(toplot$type))

# reorder type 
toplot$type  = factor(toplot$type, levels=c(
  types[3],
  types[1],
  types[2]))

toplot$`# contigs` <- as.numeric(toplot$`# contigs`)
toplot$`# predicted genes` <- as.numeric(toplot$`# predicted genes`)
toplot$`N50 (contigs)` <- as.numeric(toplot$`N50 (contigs)`)

plot_contigs <- ggplot(data= toplot, aes(x=`# contigs`, fill=type)) +
  geom_histogram(alpha=0.6, position = 'stack', binwidth=15) +
  scale_fill_manual(values=c("#6C6C6C", #grey
                             "#F8766D", #red
                             "#45B4B8" #blue
  )) + 
  xlab("number of contigs") +
  ylab("Frequency")+
  theme_minimal()+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(size=11),
        axis.text.y=element_text(size=11))+
  labs(fill="") +
  lims(x=c(0,750))

plot_N50 <- ggplot(data= toplot, aes(x=`N50 (contigs)`, fill=type)) +
  geom_histogram(bins = 50, alpha=0.6, position = 'stack') + 
  scale_x_log10() +
  scale_fill_manual(values=c("#6C6C6C", #grey
                             "#F8766D", #red
                             "#45B4B8" #blue
  )) + 
  xlab("N50 (log10)") +
  ylab("Frequency")+
  theme_minimal()+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(size=11),
        axis.text.y=element_text(size=11))+
  labs(fill="") 
  


plot_pred_genes <- ggplot(data= toplot, aes(x=`# predicted genes`, fill=type)) +
  geom_histogram(alpha=0.6, position = 'stack', binwidth=15) +
  scale_fill_manual(values=c("#6C6C6C", #grey
                             "#F8766D", #red
                             "#45B4B8" #blue
  )) + 
  scale_x_continuous(breaks=seq(1e3,6e3, by = 2e3), label=scientific) +
  xlab("predicted genes") +
  ylab("Frequency")+
  theme_minimal()+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(size=11),
        axis.text.y=element_text(size=11))+
  #scale_x_continuous(breaks=pretty(`# predicted genes`, length(`# predicted genes`)/2))+
  labs(fill="")

plots <- ggarrange(plot_contigs,plot_pred_genes,cm_Compl_vs_Contam,nrow=1,
                   labels=c("A","B","C"),
                   common.legend=TRUE)

plots_new <- ggarrange(plot_contigs,plot_N50,cm_Compl_vs_Contam,nrow=1,
                   labels=c("A","B","C"),
                   common.legend=TRUE)

pdf(paste0(out_dir,"cm_Compl_vs_Contam_contigs_genes.pdf"), width=7, height=4)
plots
dev.off()


pdf(paste0(out_dir,"cm_Compl_vs_Contam_contigs_genes_new.pdf"), width=7, height=4)
plots_new
dev.off()

summary(as.numeric(toplot$`Genome size (bp)`))
summary(as.numeric(toplot$`N50 (contigs)`))
summary(as.numeric(toplot$`# marker sets`))

summary(toplot$Completeness)
summary(toplot$Contamination)


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

sink(file = paste0(out_dir,"cm_numbers.txt"), 
     append = TRUE, type = c("output"))
paste0("######################################################################")
paste0("Of the nearly complete bins: ")
paste0("Assigned to the Kingdom of Bacteria:   ",
       NROW(bacteria),"  ",round(NROW(bacteria)/NROW(nn),4)*100,"%")
paste0("Assigned to the Kingdom of Archaea:   ",
       NROW(archaea),"  ",round(NROW(archaea)/NROW(nn),4)*100,"%")
paste0("Resolved at Phylum level:   ",
       NROW(phylum),"  ",round(NROW(phylum)/NROW(nn),4)*100,"%")
paste0("Resolved at Class level:   ",
       NROW(class),"  ",round(NROW(class)/NROW(nn),4)*100,"%")
paste0("Resolved at Order level:   ",
       NROW(order),"  ",round(NROW(order)/NROW(nn),4)*100,"%")
paste0("Resolved at Family level:   ",
       NROW(family),"  ",round(NROW(family)/NROW(nn),4)*100,"%")
paste0("Resolved at Genus level:   ",
       NROW(genus),"  ",round(NROW(genus)/NROW(nn),4)*100,"%")
paste0("Resolved at Species level:   ",
       NROW(species),"  ",round(NROW(species)/NROW(nn),4)*100,"%")
sink()


######################################################################

# -- Proportions of ALL TAXA in the >90 <5 bins: 


k <- kingdom %>%
  dplyr::mutate(taxa = sub(".*k__ *(.*?) *;p__.*", "\\1", taxa)) %>%
  dplyr::mutate(taxa = gsub(".*k__", "", taxa)) %>%
  group_by(taxa) %>%
  tally() %>%
  dplyr::mutate(perc=round(n/sum(n)*100,2)) %>%
  dplyr::arrange(desc(perc))

p <- phylum %>%
  dplyr::mutate(taxa = sub(".*p__ *(.*?) *;c__.*", "\\1", taxa)) %>%
  dplyr::mutate(taxa = gsub(".*p__", "", taxa)) %>%
  group_by(taxa) %>%
  tally() %>%
  dplyr::mutate(perc=round(n/sum(n)*100,2)) %>%
  dplyr::arrange(desc(perc))

c <- class %>%
  dplyr::mutate(taxa = sub(".*c__ *(.*?) *;o__.*", "\\1", taxa)) %>%
  dplyr::mutate(taxa = gsub(".*c__", "", taxa)) %>%
  group_by(taxa) %>%
  tally() %>%
  dplyr::mutate(perc=round(n/sum(n)*100,2)) %>%
  dplyr::arrange(desc(perc))

o <- order %>%
  dplyr::mutate(taxa = sub(".*o__ *(.*?) *;f__.*", "\\1", taxa)) %>%
  dplyr::mutate(taxa = gsub(".*o__", "", taxa)) %>%
  group_by(taxa) %>%
  tally() %>%
  dplyr::mutate(perc=round(n/sum(n)*100,2)) %>%
  dplyr::arrange(desc(perc))

f <- family %>%
  dplyr::mutate(taxa = sub(".*f__ *(.*?) *;g__.*", "\\1", taxa)) %>%
  dplyr::mutate(taxa = gsub(".*f__", "", taxa)) %>%
  group_by(taxa) %>%
  tally() %>%
  dplyr::mutate(perc=round(n/sum(n)*100,2)) %>%
  dplyr::arrange(desc(perc))

g <- genus %>%
  dplyr::mutate(taxa = sub(".*g__ *(.*?) *;s__.*", "\\1", taxa)) %>%
  dplyr::mutate(taxa = gsub(".*g__", "", taxa)) %>%
  group_by(taxa) %>%
  tally() %>%
  dplyr::mutate(perc=round(n/sum(n)*100,2)) %>%
  dplyr::arrange(desc(perc))

s <- species %>%
  dplyr::mutate(taxa = sub(".*s__ *(.*?) *", "\\1", taxa)) %>%
  dplyr::mutate(taxa = gsub(".*s__", "", taxa)) %>%
  group_by(taxa) %>%
  tally() %>%
  dplyr::mutate(perc=round(n/sum(n)*100,2)) %>%
  dplyr::arrange(desc(perc))




p <- as.data.frame(p)
c <- as.data.frame(c)
o <- as.data.frame(o)
f <- as.data.frame(f)
g <- as.data.frame(g)
s <- as.data.frame(s)

sink(file = paste0(out_dir,"cm_numbers.txt"), 
     append = TRUE, type = c("output"))
paste0("################################ PIGGIES  #############################")
paste0(" Kingdom distribution : ")
as.data.frame(k)
paste0("######################################################################")
paste0(" Phyla distribution : ")
as.data.frame(p)
paste0("######################################################################")
paste0(" Class distribution : ")
as.data.frame(c)
paste0("######################################################################")
paste0(" Order distribution : ")
as.data.frame(o)
paste0("######################################################################")
paste0(" Family distribution : ")
as.data.frame(f)
paste0("######################################################################")
paste0(" Genus distribution : ")
as.data.frame(g)
paste0("######################################################################")
paste0(" Species distribution : ")
as.data.frame(s)
sink()


######################################################################



# MERGE checkM info of nearly complete bins to bins counts : 

# rename cols of checkm nearly complete bins to match colnames of no_reps_all for dataframes to merge  
colnames(nn)[colnames(nn)=="pigid"] <- "pig"
colnames(nn)[colnames(nn)=="Bin Id"] <- "bin"
colnames(nn)[colnames(nn)=="Taxonomy (contained)"] <- "taxa"

nn <- nn %>%
  dplyr::select(pig,bin,taxa)

# merge 
df <- right_join(no_reps_all, nn, by=c("pig","bin"))
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
                                    "t10"))

# reorder cohorts 
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

pdf(paste0(out_dir,"cm_treemap_phyla.pdf"))
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

pdf(paste0(out_dir,"cm_treemap_class.pdf"))
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

pdf(paste0(out_dir,"cm_treemap_order.pdf"))
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

