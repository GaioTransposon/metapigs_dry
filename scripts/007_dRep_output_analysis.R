
# dRep_output_analysis.R
# analysis of dRep output 

library(readxl)
library(data.table)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(splitstackshape)
library(pheatmap)
library(ggpubr)
library(robCompositions)
library(ggbiplot)

source_dir = "/Users/12705859/metapigs_dry/source_data/" # git 
middle_dir = "/Users/12705859/metapigs_dry/middle_dir/" # git 
out_dir = "/Users/12705859/Desktop/metapigs_dry/dRep/"  # local


# upload dRep output 
Cdb <- read_csv(paste0(middle_dir,"Cdb.csv"))

# upload cohorts info
cohorts <- read_xlsx(paste0(source_dir,"cohorts.xlsx"))

C1 <- separate(data = Cdb, col = genome, into = c("pig", "bin"), sep = "_")
C1 <- C1[,c("pig","bin","primary_cluster","secondary_cluster")]
C1$primary_cluster <- as.character(C1$primary_cluster)
head(C1)

# upload bins with counts (from output of 7.R)
no_reps_all <- read.csv(paste0(middle_dir,"no_reps_all.csv"), 
                      na.strings=c("","NA"),
                      check.names = FALSE,
                      header = TRUE)

# remove .fa extension to match bins in checkm df 
no_reps_all$bin <- gsub(".fa","", no_reps_all$bin)
head(no_reps_all)
NROW(no_reps_all)

no_reps_all$primary_cluster <- paste0(no_reps_all$secondary_cluster)
no_reps_all <- cSplit(no_reps_all,"primary_cluster","_")
no_reps_all$primary_cluster_2 <- NULL
colnames(no_reps_all)[colnames(no_reps_all)=="primary_cluster_1"] <- "primary_cluster"

######################################################################

# load gtdbtk assignments of the bins

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv(paste0(middle_dir,"gtdb_bins_completeTaxa"),
                        col_types = cols(node = col_character(),
                                         pig = col_character()))


######################################################################


# create text file to contain dRep text output

sink(file = paste0(out_dir,"dRep_numbers.txt"), 
     append = FALSE, type = c("output"))
sink()


###########################


sink(file = paste0(out_dir,"dRep_numbers.txt"), 
     append = TRUE, type = c("output"))
paste0("dRep-clustered bins: ", 
       round(NROW(C1)/NROW(gtdbtk_bins)*100,2),
       "%",
       " (n=",NROW(C1),")" )
paste0("of which primary clusters: ", length(unique(C1$primary_cluster)) )
paste0("of which secondary clusters ", length(unique(C1$secondary_cluster)) )
sink()



########################################################################################################


# Extent of agreement between dRep and GTDBTK classification: 


df <- merge(no_reps_all, gtdbtk_bins, by=c("pig","bin"))



# Primary clusters: 


b <- df %>%
  dplyr::filter(!primary_cluster=="no") %>%
  dplyr::select(node,domain,phylum,class,order,family,genus,species,primary_cluster) 


b_species <- b %>%
  dplyr::group_by(primary_cluster,species) %>%
  dplyr::summarise(freq= n())  %>%
  dplyr::mutate(num2= freq/sum(freq))  %>%
  dplyr::group_by(primary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(primary_cluster,freq,num2) %>%
  dplyr::rename(., agree = num2) 

b_genus <- b %>%
  dplyr::group_by(primary_cluster,genus) %>%
  dplyr::summarise(freq= n())  %>%
  dplyr::mutate(num2= freq/sum(freq))  %>%
  dplyr::group_by(primary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(primary_cluster,freq,num2) %>%
  dplyr::rename(., agree = num2) 

b_family <- b %>%
  dplyr::group_by(primary_cluster,family) %>%
  dplyr::summarise(freq= n())  %>%
  dplyr::mutate(num2= freq/sum(freq))  %>%
  dplyr::group_by(primary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(primary_cluster,freq,num2) %>%
  dplyr::rename(., agree = num2) 

b_order <- b %>%
  dplyr::group_by(primary_cluster,order) %>%
  dplyr::summarise(freq= n())  %>%
  dplyr::mutate(num2= freq/sum(freq))  %>%
  dplyr::group_by(primary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(primary_cluster,freq,num2) %>%
  dplyr::rename(., agree = num2) 

b_class <- b %>%
  dplyr::group_by(primary_cluster,class) %>%
  dplyr::summarise(freq= n())  %>%
  dplyr::mutate(num2= freq/sum(freq))  %>%
  dplyr::group_by(primary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(primary_cluster,freq,num2) %>%
  dplyr::rename(., agree = num2) 

b_phylum <- b %>%
  dplyr::group_by(primary_cluster,phylum) %>%
  dplyr::summarise(freq= n())  %>%
  dplyr::mutate(num2= freq/sum(freq))  %>%
  dplyr::group_by(primary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(primary_cluster,freq,num2) %>%
  dplyr::rename(., agree = num2) 



b_phylum$GTDB_taxa_level="phylum"
b_class$GTDB_taxa_level="class"
b_order$GTDB_taxa_level="order"
b_family$GTDB_taxa_level="family"
b_genus$GTDB_taxa_level="genus"
b_species$GTDB_taxa_level="species"


prim_clu_agree <- rbind(b_phylum,
                        b_class,
                        b_order,
                        b_family,
                        b_genus,
                        b_species)

# reorder taxa levels 
prim_clu_agree$GTDB_taxa_level  = factor(prim_clu_agree$GTDB_taxa_level, levels=c("phylum",
                                                                                  "class",
                                                                                  "order",
                                                                                  "family",
                                                                                  "genus",
                                                                                  "species"))

means <- aggregate(agree ~  GTDB_taxa_level, prim_clu_agree, mean)

primary_clusters_agreement_plot <- ggplot(prim_clu_agree, aes(agree,GTDB_taxa_level))+
  geom_boxplot()+
  geom_point(aes(size = freq^3,color=freq)) +
  guides(size=FALSE)+
  xlim(0,1.2)+
  geom_text(data = means, aes(label = paste0(round(agree*100,2),"%"), x = 1.1))+
  theme_bw() +
  ggtitle("GTDB- MAGs assignments agreement with primary clusters (95% ANI)")


sink(file = paste0(out_dir,"dRep_numbers.txt"), 
     append = TRUE, type = c("output"))
paste0("Extent of agreement between dRep classification and gtdbtk assignment of bins")
paste0("Primary clusters: ")
tapply(prim_clu_agree$agree, prim_clu_agree$GTDB_taxa_level, summary)
sink()


######################################################################


# Secondary clusters: 


b <- df %>%
  dplyr::filter(!secondary_cluster=="no_cluster") %>%
  dplyr::select(node,domain,phylum,class,order,family,genus,species,secondary_cluster) 


b_species <- b %>%
  dplyr::group_by(secondary_cluster,species) %>%
  dplyr::summarise(freq= n())  %>%
  dplyr::mutate(num2= freq/sum(freq))  %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(secondary_cluster,freq,num2) %>%
  dplyr::rename(., agree = num2) 

b_genus <- b %>%
  dplyr::group_by(secondary_cluster,genus) %>%
  dplyr::summarise(freq= n())  %>%
  dplyr::mutate(num2= freq/sum(freq))  %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(secondary_cluster,freq,num2) %>%
  dplyr::rename(., agree = num2) 

b_family <- b %>%
  dplyr::group_by(secondary_cluster,family) %>%
  dplyr::summarise(freq= n())  %>%
  dplyr::mutate(num2= freq/sum(freq))  %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(secondary_cluster,freq,num2) %>%
  dplyr::rename(., agree = num2) 

b_order <- b %>%
  dplyr::group_by(secondary_cluster,order) %>%
  dplyr::summarise(freq= n())  %>%
  dplyr::mutate(num2= freq/sum(freq))  %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(secondary_cluster,freq,num2) %>%
  dplyr::rename(., agree = num2) 

b_class <- b %>%
  dplyr::group_by(secondary_cluster,class) %>%
  dplyr::summarise(freq= n())  %>%
  dplyr::mutate(num2= freq/sum(freq))  %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(secondary_cluster,freq,num2) %>%
  dplyr::rename(., agree = num2) 

b_phylum <- b %>%
  dplyr::group_by(secondary_cluster,phylum) %>%
  dplyr::summarise(freq= n())  %>%
  dplyr::mutate(num2= freq/sum(freq))  %>%
  dplyr::group_by(secondary_cluster) %>%
  dplyr::top_n(1, num2) %>% 
  dplyr::select(secondary_cluster,freq,num2) %>%
  dplyr::rename(., agree = num2) 



b_phylum$GTDB_taxa_level="phylum"
b_class$GTDB_taxa_level="class"
b_order$GTDB_taxa_level="order"
b_family$GTDB_taxa_level="family"
b_genus$GTDB_taxa_level="genus"
b_species$GTDB_taxa_level="species"


sec_clu_agree <- rbind(b_phylum,
                       b_class,
                       b_order,
                       b_family,
                       b_genus,
                       b_species)

# reorder taxa levels 
sec_clu_agree$GTDB_taxa_level  = factor(sec_clu_agree$GTDB_taxa_level, levels=c("phylum",
                                                                                "class",
                                                                                "order",
                                                                                "family",
                                                                                "genus",
                                                                                "species"))

means <- aggregate(agree ~  GTDB_taxa_level, sec_clu_agree, mean)

secondary_clusters_agreement_plot <- ggplot(sec_clu_agree, aes(agree,GTDB_taxa_level))+
  geom_boxplot()+
  geom_point(aes(size = freq^3,color=freq)) +
  guides(size=FALSE)+
  xlim(0,1.2)+
  geom_text(data = means, aes(label = paste0(round(agree*100,2),"%"), x = 1.1))+
  theme_bw() +
  ggtitle("GTDB- MAGs assignments agreement with secondary clusters (99% ANI)")


sink(file = paste0(out_dir,"dRep_numbers.txt"), 
     append = TRUE, type = c("output"))
paste0("Extent of agreement between dRep classification and gtdbtk assignment of bins")
paste0("Secondary clusters: ")
tapply(sec_clu_agree$agree, sec_clu_agree$GTDB_taxa_level, summary)
sink()


agreement_plots <- ggarrange(primary_clusters_agreement_plot,
          secondary_clusters_agreement_plot,
          nrow=2,
          common.legend = TRUE)

pdf(paste0(out_dir,"dRep_GTDB_extent_of_agreement.pdf"))
agreement_plots
dev.off()


######################################################################
######################################################################

########################################################################################################



# Display amount of shared vs unique clusters 


# primary clusters, piglets
primary_piglets <- no_reps_all %>%
  dplyr::filter(!cohort=="Mothers") %>%
  dplyr::select(pig,primary_cluster) %>% 
  dplyr::distinct() %>%
  group_by(primary_cluster) %>%
  dplyr::mutate(type = ifelse(n() > 1, "common","unique")) %>%
  ggplot(., aes(pig)) +
  geom_bar(aes(fill=type), width = 0.5) + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  labs(title="Shared versus common primary clusters (95% ANI)", 
       subtitle="distribution among piglets",
       x = "piglets",
       y = "clustered bins") +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_text(),
        title=element_text(size=7))

# secondary clusters, piglets
secondary_piglets <- no_reps_all %>%
  dplyr::filter(!cohort=="Mothers") %>%
  dplyr::select(pig,secondary_cluster) %>% 
  dplyr::distinct() %>%
  group_by(secondary_cluster) %>%
  dplyr::mutate(type = ifelse(n() > 1, "common","unique")) %>%
  ggplot(., aes(pig)) +
  geom_bar(aes(fill=type), width = 0.5) + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  labs(title="Shared versus common secondary clusters (99% ANI)", 
       subtitle="distribution among piglets",
       x = "piglets",
       y = "clustered bins") +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_text(),
        title=element_text(size=7))




# primary clusters, mothers
primary_mothers <- no_reps_all %>%
  dplyr::filter(cohort=="Mothers") %>%
  dplyr::select(pig,primary_cluster) %>% 
  dplyr::distinct() %>%
  group_by(primary_cluster) %>%
  dplyr::mutate(type = ifelse(n() > 1, "common","unique")) %>%
  ggplot(., aes(pig)) +
  geom_bar(aes(fill=type), width = 0.5) + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  labs(title="Shared versus common primary clusters (95% ANI)", 
       subtitle="distribution among mothers",
       x = "mothers",
       y = "clustered bins") +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_text(),
        title=element_text(size=7))

# secondary clusters, mothers
secondary_mothers <- no_reps_all %>%
  dplyr::filter(cohort=="Mothers") %>%
  dplyr::select(pig,secondary_cluster) %>% 
  dplyr::distinct() %>%
  group_by(secondary_cluster) %>%
  dplyr::mutate(type = ifelse(n() > 1, "common","unique")) %>%
  ggplot(., aes(pig)) +
  geom_bar(aes(fill=type), width = 0.5) + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  labs(title="Shared versus common secondary clusters (99% ANI)", 
       subtitle="distribution among mothers",
       x = "mothers",
       y = "clustered bins") +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_text(),
        title=element_text(size=7))



tosave <- ggarrange(primary_piglets,
          secondary_piglets,
          primary_mothers,
          secondary_mothers,
          ncol=2,
          nrow=2,
          common.legend = TRUE,
          widths = c(1,1),
          heights = c(1,1))


pdf(paste0(out_dir,"dRep_common_vs_unique_clusters.pdf"))
tosave
dev.off()


sink(file = paste0(out_dir,"dRep_numbers.txt"), 
     append = TRUE, type = c("output"))
paste0("Percentage of shared primary clusters among piglets:")
NROW(which(primary_piglets$data$type=="common"))/NROW(primary_piglets$data)*100
paste0("Percentage of shared secondary clusters among piglets:")
NROW(which(secondary_piglets$data$type=="common"))/NROW(secondary_piglets$data)*100
paste0("#######")
paste0("Percentage of shared primary clusters among mothers:")
NROW(which(primary_mothers$data$type=="common"))/NROW(primary_mothers$data)*100
paste0("Percentage of shared secondary clusters among mothers:")
NROW(which(secondary_mothers$data$type=="common"))/NROW(secondary_mothers$data)*100
sink()



######################################################################################################
######################################################################################################

# Principal component analysis: clustering of bins based on phyla with time (labels per cohort)




# I hereby select taxa_2 only (corresponds to phylum) and remove any row where no phylum was resolved
df1 <- no_reps_all %>%
  dplyr::select(pig,bin,date,value,secondary_cluster,cohort)

df1 <- as.data.frame(na.omit(df1))

unique(df1$secondary_cluster)


#################################
#################################

# CREATE COUNTS TABLE 

df1 <- df1 %>%
  dplyr::filter(!secondary_cluster=="no_cluster")
head(df1)
NROW(df1)

# PROCEED to all: 

# for each sample (pig,date), sum up the counts that fall within one species (same species assigned to distinct bins)
df2 <- df1 %>%
  group_by(pig,secondary_cluster,date) %>%
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
  dplyr::select(sample,secondary_cluster,norm_value) %>%
  pivot_wider(names_from = secondary_cluster, values_from = norm_value, values_fill = list(norm_value = 0))

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
dRep_PC12 <- ggbiplot(df6.pca,
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
dRep_PC34 <- ggbiplot(df6.pca,
                      labels=this_mat_samples$sample_2,
                      groups=this_mat_samples$sample_1,
                      ellipse=TRUE,
                      var.axes = FALSE,
                      labels.size = 2,
                      choices = (3:4)) +
  theme_bw() +
  scale_colour_discrete(name="timepoint")+
  guides(color = guide_legend(ncol = 1))


dRep_PCA <- ggarrange(dRep_PC12,dRep_PC34,
                      ncol=2,legend = "right",
                      common.legend=TRUE)

pdf(paste0(out_dir,"dRep_PCA.pdf"), width=7,height=4)
dRep_PCA
dev.off()

