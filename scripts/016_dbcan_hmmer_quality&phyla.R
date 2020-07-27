
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
library(VennDiagram)


source_dir = "/Users/12705859/metapigs_dry/source_data/" # git 
middle_dir = "/Users/12705859/metapigs_dry/middle_dir/" # git 
out_dir_git = "/Users/12705859/metapigs_dry/out/" # git 
out_dir = "/Users/12705859/Desktop/metapigs_dry/dbcan/"  # local

######################################################################

# load Clean & CleanCounts - hmmer data 

hmmerClean <- read.csv(paste0(middle_dir, "hmmer.out_Clean.csv.gz"))
# re-order enzymes to follow the same aesthetics as Stewart et al 2019
hmmerClean$enzymeNAME <- as.character(hmmerClean$enzymeNAME)
hmmerClean$enzymeNAME  = factor(hmmerClean$enzymeNAME, levels=c("AA","CE","SLH","GH","GT","CBM","cohesin","PL"))

hmmerCleanCounts <- read.csv(paste0(middle_dir, "hmmer.out_CleanCounts.csv.gz"))
hmmerCleanCounts$enzymeNAME <- as.character(hmmerCleanCounts$enzymeNAME)
hmmerCleanCounts$enzymeNAME  = factor(hmmerCleanCounts$enzymeNAME, levels=c("AA","CE","SLH","GH","GT","CBM","cohesin","PL"))

##########################################################

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv(paste0(middle_dir,"gtdb_bins_completeTaxa"),
                        col_types = cols(node = col_character(),
                                         pig = col_character()))

##########################################################
##########################################################
##########################################################


# PART 1 : QUALITY of mapping


##########################################################


# PLOT quality

enzymes_proportion <- hmmerClean %>% 
  group_by(enzymeNAME) %>% 
  tally() %>%
  dplyr::mutate(perc = paste0(round(n/sum(n)*100,2),"%"))

hmmerClean_perc_identity_plot <- hmmerClean %>%
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
  ggtitle("Percentage identity against CAZy") + 
  geom_text(data = enzymes_proportion,
            aes(enzymeNAME, 0, label = perc), vjust="inward",size=2)

pdf(paste0(out_dir,"dbcan_HMMER_identity.pdf"))
hmmerClean_perc_identity_plot
dev.off()


# how well our predicted proteins - matched against the CAZy database 
sink(file=paste0(out_dir,"dbcan_hmmerClean_quality.txt"), append=FALSE)
paste0("total number of proteins matching the CAZy db with dbcan-hmmer")
NROW(hmmerClean)
paste0("total number of unique enzymeIDs with dbcan-hmmer")
NROW(unique(hmmerClean$enzymeID))
paste0("coverage")
summary(hmmerClean$coverage)
paste0("e-values")
summary(hmmerClean$evalue)
paste0("Enzyme classes - Proportions")
hmmerClean %>%
  group_by(enzymeNAME) %>%
  tally() %>%
  dplyr::mutate(perc=n/sum(n))
sink()

##########################################################


# PART 2 : join GTDB data with hmmer CleanCounts data


##########################################################

NROW(gtdbtk_bins)
NROW(hmmerClean)
# merge gtdbtk- bins taxonomic assignments with hmmer data
gt_hmmer <- inner_join(gtdbtk_bins,hmmerCleanCounts)
NROW(gt_hmmer)
head(gt_hmmer)

gt_hmmer$phylum[is.na(gt_hmmer$phylum)] <- "Unknown"

gt_hmmer <- as.data.frame(gt_hmmer)

##########################################################


# PART 3 : Distribution over Phyla 


##########################################################

head(gt_hmmer)

prop_CAZ_per_phylum <- gt_hmmer %>%
  group_by(phylum) %>%
  dplyr::summarise(sum=sum(enz_count)) %>%
  dplyr::mutate(value=(sum/sum(sum))*100)

prop_CAZclass_per_phylum <- gt_hmmer %>%
  group_by(enzymeNAME,phylum) %>%
  dplyr::summarise(sum=sum(enz_count)) %>%
  group_by(phylum) %>%
  dplyr::mutate(value=(sum/sum(sum))*100)

a_gt <- ggplot(prop_CAZ_per_phylum, aes(y=value, x=phylum)) + 
  geom_bar(position="dodge", stat="identity",colour="black",fill="snow3")+ #lightskyblue happier
  theme(axis.text.x=element_text(angle=90))+
  theme_pubr()+
  ylab("Proportion of CAZymes (%)")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=8),
        axis.text.x = element_blank())
b_gt <- ggplot(prop_CAZclass_per_phylum, aes(fill=enzymeNAME, y=value, x=phylum)) + 
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

pdf(paste0(out_dir,"dbcan_HMMER_GTphyla.pdf"))
CAZ_HMMER_GTphyla_plot
dev.off()


##########################################################

# merge figures : quality + phyla distribution 

pdf(paste0(out_dir,"dbcan_HMMER_identity&phyla.pdf"))
ggarrange(hmmerClean_perc_identity_plot,
          CAZ_HMMER_GTphyla_plot,
          ncol=2)
dev.off()

##########################################################
##########################################################


# PART 4


##########################################################

# create an enzymeID to genus & enzymeID to species profile 

gt_hmmer_genus <- gt_hmmer %>% 
  group_by(enzymeID,genus) %>%
  dplyr::summarise(sum=sum(enz_count),
                   sd=sd(enz_count)) %>%
  group_by(enzymeID) %>%
  dplyr::mutate(percentage_genus=(sum/sum(sum))*100)

gt_hmmer_species <- gt_hmmer %>% 
  group_by(enzymeID,species) %>%
  dplyr::summarise(sum=sum(enz_count),
                   sd=sd(enz_count)) %>%
  group_by(enzymeID) %>%
  dplyr::mutate(percentage_species=(sum/sum(sum))*100)

# save as tsv
fwrite(gt_hmmer_genus, file=paste0(out_dir_git,"gt_hmmer_genus.tsv"), sep = "\t")
fwrite(gt_hmmer_species, file=paste0(out_dir_git,"gt_hmmer_species.tsv"), sep = "\t")


##########################################################
##########################################################


# PART 5


##########################################################


# Venn diagram to show, for each enzyme class, the top three genera
# and the extent of share of unique enzymeIDs


# function to generate Venn plot
top_genera_per_enzymeCLASS <- function(df_enzymeNAME_selection) {
  
  df1 <- df_enzymeNAME_selection %>% 
    group_by(enzymeNAME,genus) %>%
    dplyr::summarise(sum=sum(enz_count)) %>%
    dplyr::arrange(desc(sum))
  
  top_genera <- as.list(df1[1:3,2])
  
  set1 <- df %>% 
    dplyr::filter(genus==top_genera$genus[1]) %>% 
    dplyr::select(enzymeID) %>%
    distinct()
  set1 <- as.character(set1$enzymeID)
  
  set2 <- df %>% 
    dplyr::filter(genus==top_genera$genus[2]) %>% 
    dplyr::select(enzymeID) %>%
    distinct()
  set2 <- as.character(set2$enzymeID)
  
  set3 <- df %>% 
    dplyr::filter(genus==top_genera$genus[3]) %>% 
    dplyr::select(enzymeID) %>%
    distinct()
  set3 <- as.character(set3$enzymeID)
  
  # Prepare a palette of 3 colors with R colorbrewer:
  myCol <- brewer.pal(3, "Pastel2")
  
  # Chart
  v1 <- venn.diagram(
    x = list(set1, set2, set3),
    category.names = c(as.character(top_genera$genus[1]), 
                       as.character(top_genera$genus[2]), 
                       as.character(top_genera$genus[3])),
    
    main=paste0(as.character(df_enzymeNAME_selection$enzymeNAME[1]),
                " (n=",NROW(unique(df_enzymeNAME_selection$enzymeID)),")"),
    
    #filename = paste0(out_dir,'zzzz_diagramm.tiff'),
    output=FALSE,
    filename = NULL,
    
    # Output features
    imagetype="tiff" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontgenus = "sans",
    
    # Set names
    cat.cex = 0.6,
    cat.fontface = "italic",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontgenus = "sans",
    rotation = 1
  )
  
  return(v1)
}


flog.threshold(ERROR)  # to omit the log message from venn.diagram


pdf(paste0(out_dir,"dbcan_hmmer_Venn.pdf"))
grid.arrange(gTree(children=top_genera_per_enzymeCLASS(as.data.frame(gt_hmmer %>% dplyr::filter(enzymeNAME=="GH")))),
             gTree(children=top_genera_per_enzymeCLASS(as.data.frame(gt_hmmer %>% dplyr::filter(enzymeNAME=="CE")))),
             gTree(children=top_genera_per_enzymeCLASS(as.data.frame(gt_hmmer %>% dplyr::filter(enzymeNAME=="AA")))),
             gTree(children=top_genera_per_enzymeCLASS(as.data.frame(gt_hmmer %>% dplyr::filter(enzymeNAME=="CBM")))),
             gTree(children=top_genera_per_enzymeCLASS(as.data.frame(gt_hmmer %>% dplyr::filter(enzymeNAME=="SLH")))),
             gTree(children=top_genera_per_enzymeCLASS(as.data.frame(gt_hmmer %>% dplyr::filter(enzymeNAME=="GT")))),
             gTree(children=top_genera_per_enzymeCLASS(as.data.frame(gt_hmmer %>% dplyr::filter(enzymeNAME=="cohesin")))),
             ncol=2,nrow=4)
dev.off()

##########################################################
##########################################################
