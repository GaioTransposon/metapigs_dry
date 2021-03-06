
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

# NROW(unique(gt_hmmer$species))
##########################################################


# PART 3 : Distribution over Phyla 


##########################################################

head(gt_hmmer)

#######
# merge gtdbtk- bins taxonomic assignments with hmmer data (here I need the number of proteins)
gt_hmmer1 <- inner_join(gtdbtk_bins,hmmerClean)
head(gt_hmmer1)

gt_hmmer2 <- gt_hmmer1 %>%
  dplyr::select(predicted_protein_ID,phylum) %>%
  distinct() %>%
  add_tally() %>%
  group_by(phylum) %>%
  dplyr::summarise(sum=sum(n)) %>%
  dplyr::arrange(desc(sum))

this_order <- as.list(gt_hmmer2$phylum)
#######


prop_CAZ_per_phylum <- gt_hmmer %>%
  group_by(phylum) %>%
  dplyr::summarise(sum=sum(enz_count)) %>%
  dplyr::mutate(value=(sum/sum(sum))*100)

prop_CAZclass_per_phylum <- gt_hmmer %>%
  group_by(enzymeNAME,phylum) %>%
  dplyr::summarise(sum=sum(enz_count)) %>%
  group_by(phylum) %>%
  dplyr::mutate(value=(sum/sum(sum))*100)

aa_gt <- ggplot(gt_hmmer2, aes(y=sum, x=factor(phylum, levels = this_order))) + 
  geom_bar(position="dodge", stat="identity",colour="black",fill="snow1")+ #lightskyblue happier
  theme_pubr()+
  ylab("Total predicted proteins")+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.y = element_text(size=11),
        axis.title.x=element_blank())

a_gt <- ggplot(prop_CAZ_per_phylum, aes(y=value, x=factor(phylum, levels = this_order))) + 
  geom_bar(position="dodge", stat="identity",colour="black",fill="snow3")+ #lightskyblue happier
  theme_pubr()+
  ylab("Proportion of CAZymes (%)")+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.y = element_text(size=11),
        axis.title.x=element_blank())


b_gt <- ggplot(prop_CAZclass_per_phylum, aes(fill=enzymeNAME, y=value, x=factor(phylum, levels = this_order))) + 
  geom_bar(position="fill", stat="identity",colour="black")+
  scale_fill_brewer(palette="Spectral")+
  theme_pubr()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=6),
        axis.text.y=element_text(size=8),
        axis.title.y = element_text(size=10),
        axis.title.x=element_blank(),
        legend.position="right",legend.text = element_text(size=8),
        legend.title = element_blank())

CAZ_HMMER_GTphyla_plot <- plot_grid(
  plot_grid(
    aa_gt + theme(legend.position = "none"),
    a_gt + theme(legend.position = "none"),
    b_gt + theme(legend.position = "none"),
    ncol = 1,
    rel_heights=c(3,3,4),
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
    dplyr::arrange(desc(sum)) %>% 
    dplyr::mutate(perc=round((sum/sum(sum))*100,1))
  
  top_genera <- as.list(df1[1:3,2])
  top_perc <- as.list(df1[1:3,4])
  
  set1 <- df_enzymeNAME_selection %>% 
    dplyr::filter(genus==top_genera$genus[1]) %>% 
    dplyr::select(enzymeID) %>%
    distinct()
  set1 <- as.character(set1$enzymeID)
  
  set2 <- df_enzymeNAME_selection %>% 
    dplyr::filter(genus==top_genera$genus[2]) %>% 
    dplyr::select(enzymeID) %>%
    distinct()
  set2 <- as.character(set2$enzymeID)
  
  set3 <- df_enzymeNAME_selection %>% 
    dplyr::filter(genus==top_genera$genus[3]) %>% 
    dplyr::select(enzymeID) %>%
    distinct()
  set3 <- as.character(set3$enzymeID)
  
  # Prepare a palette of 3 colors with R colorbrewer:
  myCol <- brewer.pal(3, "Pastel2")
  
  # Chart
  v1 <- venn.diagram(
    x = list(set1, set2, set3),
    category.names = c(as.character(paste0(top_genera$genus[1]," ",round(top_perc$perc[1],1),"%")), 
                       as.character(paste0(top_genera$genus[2]," ",round(top_perc$perc[2],1),"%")), 
                       as.character(paste0(top_genera$genus[3]," ",round(top_perc$perc[3],1),"%"))),
    
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
             gTree(children=top_genera_per_enzymeCLASS(as.data.frame(gt_hmmer %>% dplyr::filter(enzymeNAME=="PL")))),
             gTree(children=top_genera_per_enzymeCLASS(as.data.frame(gt_hmmer %>% dplyr::filter(enzymeNAME=="GT")))),
             gTree(children=top_genera_per_enzymeCLASS(as.data.frame(gt_hmmer %>% dplyr::filter(enzymeNAME=="SLH")))),
             gTree(children=top_genera_per_enzymeCLASS(as.data.frame(gt_hmmer %>% dplyr::filter(enzymeNAME=="cohesin")))),
             ncol=2,nrow=4)
dev.off()
##########################################################
##########################################################

# Gene Frequency and distribution across species 

gt_hmmer2 <- inner_join(gtdbtk_bins,hmmerClean)


unique(gt_hmmer2$phylum)


from = c("Firmicutes","Firmicutes_A",
         "Bacteroidota","Actinobacteriota",
         "Patescibacteria","Spirochaetota",
         "Firmicutes_C","Planctomycetota",  
         "Myxococcota","Proteobacteria", 
         "Cyanobacteria","Euryarchaeota",
         "Campylobacterota","Verrucomicrobiota_A",
         "Fibrobacterota","Desulfobacterota_A",
         "Synergistota","Thermoplasmatota",
         "Firmicutes_B","Verrucomicrobiota",  
         "Elusimicrobiota","Fusobacteriota",    
         "Deferribacterota","Eremiobacterota",    
         "Thermotogota",NA)

to = c("Firmicutes","Firmicutes",
       "Bacteroidota","Actinobacteriota",
       "Other Phyla","Other Phyla",
       "Firmicutes","Other Phyla",  
       "Other Phyla","Proteobacteria", 
       "Other Phyla","Other Phyla",
       "Other Phyla","Other Phyla",
       "Other Phyla","Other Phyla",
       "Other Phyla","Other Phyla",
       "Firmicutes","Other Phyla",  
       "Other Phyla","Other Phyla",    
       "Other Phyla","Other Phyla",    
       "Other Phyla","Other Phyla")

# replace 
gt_hmmer2$phylum <- plyr::mapvalues(as.character(gt_hmmer2$phylum), from, to)
unique(gt_hmmer2$phylum)


z <- gt_hmmer2 %>% 
  dplyr::select(pig,bin,predicted_protein_ID,enzymeID,enzymeNAME,species, phylum) %>% 
  distinct() %>%
  dplyr::select(enzymeID,enzymeNAME,species,phylum) 


# prepare empty file to capture genes per genome info 
sink(file=paste0(out_dir,"dbcan_hmmerClean_genesXgenome.txt"), append=FALSE)
paste0("summary of : genes per genome")
sink()

# prepare empty to catch csv data 
fwrite(x = save, 
       file=paste0(out_dir_git,"dbcan_hmmer_genesXgenomes.csv"), 
       sep = ",", append = FALSE)


# function to plot frequency of genes per genome and number of diff enzymes within each genome
gene_frequency_fun <- function(df_enzymeNAME_selected) {
  
  x <- df_enzymeNAME_selected %>% 
    dplyr::mutate(n=1) %>%
    group_by(enzymeID,species,phylum) %>%
    dplyr::summarise(sum_count=sum(n)) %>% 
    group_by(species,phylum) %>%
    dplyr::summarise(x=mean(sum_count)) 
  
  
  y <- df_enzymeNAME_selected %>% 
    dplyr::select(enzymeID,species) %>%
    distinct() %>% 
    dplyr::mutate(n=1) %>%
    group_by(species) %>%
    dplyr::summarise(y=sum(n)) 
  
  xy <- inner_join(x,y) 
  
  xy$enzyme_class <- unique(df_enzymeNAME_selected$enzymeNAME)
  
  sink(file=paste0(out_dir,"dbcan_hmmerClean_genesXgenome.txt"), append=TRUE)
  print(paste0("enzyme class : ",unique(xy$enzyme_class)))
  print(summary(xy$x))
  sink()
  
  save <- xy %>% 
    dplyr::mutate(genes_per_genome=x) %>%
    dplyr::mutate(distinct_enzymes=y) %>%
    dplyr::select(species,phylum,enzyme_class,genes_per_genome,distinct_enzymes)
  
  fwrite(x = save, 
         file=paste0(out_dir_git,"dbcan_hmmer_genesXgenomes.csv"), 
         sep = ",", append = TRUE)
  
  if (NROW(unique(y$y)) == 1) {
    
    # plot for enzyme class with only one enzymeID
    print(ggplot(xy,aes(x=x,y=y,color=phylum)) +
            geom_point()+
            labs(x = paste0("Number of genes encoding ",
                            unique(df_enzymeNAME_selected$enzymeNAME),
                            " per genome"),
                 y = paste0("Number of different ",
                            unique(df_enzymeNAME_selected$enzymeNAME),
                            " class enzymes within a genome"))+
            geom_text_repel(
              data          = subset(xy, x > 5 ), # max(xy$x)/2 ), 
              aes(x,y,label=species), 
              size=2,
              segment.size  = 0.2,
              segment.color = "grey50",
              direction     = "both"
            ) +
            theme_bw())
    
    
  } else {
    
    print(ggplot(xy,aes(x=x,y=y,color=phylum)) +
            geom_point()+
            labs(x = paste0("Number of genes encoding ",
                            unique(df_enzymeNAME_selected$enzymeNAME),
                            " per genome"),
                 y = paste0("Number of different ",
                            unique(df_enzymeNAME_selected$enzymeNAME),
                            " class enzymes within a genome"))+
            geom_text_repel(
              data          = subset(xy, x > max(xy$x)/2 | y > max(xy$y)/2+max(xy$y)/4), 
              aes(x,y,label=species), 
              size=2,
              segment.size  = 0.2,
              segment.color = "grey50",
              direction     = "both"
            ) +
            theme_bw())
    
  }
  
  
  
}


p1 <- gene_frequency_fun(as.data.frame(z %>% dplyr::filter(enzymeNAME=="GT")))
p2 <- gene_frequency_fun(as.data.frame(z %>% dplyr::filter(enzymeNAME=="GH")))
p3 <- gene_frequency_fun(as.data.frame(z %>% dplyr::filter(enzymeNAME=="CBM")))
p4 <- gene_frequency_fun(as.data.frame(z %>% dplyr::filter(enzymeNAME=="CE")))
p5 <- gene_frequency_fun(as.data.frame(z %>% dplyr::filter(enzymeNAME=="PL")))
p6 <- gene_frequency_fun(as.data.frame(z %>% dplyr::filter(enzymeNAME=="AA")))
p7 <- gene_frequency_fun(as.data.frame(z %>% dplyr::filter(enzymeNAME=="SLH")))
p8 <- gene_frequency_fun(as.data.frame(z %>% dplyr::filter(enzymeNAME=="cohesin")))


pdf(paste0(out_dir,"dbcan_hmmer_geneFrequency.pdf"))
p1
p2
p3
p4
p5
p6
p7
p8
dev.off()



