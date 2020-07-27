

##########################################################


# Prevotella vs Bacteroidetes enzyme specificity

test <- subset(df_part, (enzymeID %in% mylist))


# 1 
# plot all prevotella and all bacteroidetes genus, all enzymes, no labels 


# 2 
# plot only Bacteroidetes A and Prevotella, 50 most abundant enzymes; no labels 

# 3 
# plot only Bacteroidetes A and Prevotella, 60th-100th most abundant enzymes; labels 


View(gt_hmmer)

gt_hmmer_sub <- gt_hmmer %>% 
  dplyr::filter(genus=="Prevotella"|genus=="Bacteroides_A"|genus=="Bacteroides"|genus=="Bacteroides_B") 

# reorder  
gt_hmmer_sub$genus  = factor(gt_hmmer_sub$genus, levels=c("Bacteroides_A",
                                                          "Bacteroides_B", 
                                                          "Bacteroides",
                                                          "Prevotella"))


# all Bact. and Prev
test2 <- gt_hmmer_sub %>% 
  dplyr::filter(genus=="Prevotella"|genus=="Bacteroides_A"|genus=="Bacteroides"|genus=="Bacteroides_B") %>%
  dplyr::select(genus,enzymeID) %>%
  dplyr::distinct()
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
  dplyr::filter(genus=="Prevotella"|genus=="Bacteroides_A") %>%
  dplyr::select(genus,enzymeID) %>%
  dplyr::distinct()
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

pdf(paste0(out_dir,"dbcan_HMMER_specificity_genera.pdf"))
tosave
dev.off()


##################################

# # ENZYMES CO-OCCURRENCE in the same GENUS - piglets 
# 
# z <- gt_hmmer %>% dplyr::select(enzymeID,pig,genus)
# piggiesIDs <- no_reps_all %>% 
#   dplyr::filter(!cohort=="Mothers") %>% 
#   dplyr::select(pig) %>% 
#   dplyr::distinct()
# 
# # collect mothers enzyme info only 
# z_piggies <- left_join(piggiesIDs,z)
# 
# z1 <- z_piggies %>% 
#   group_by(genus,enzymeID) %>% 
#   tally() 
# 
# z2 <- z_piggies %>% group_by(enzymeID) %>% tally() 
# 
# # set limits 
# firstQu <- as.numeric(summary(z2$n)[2])
# median <- as.numeric(summary(z2$n)[3])
# thirdQu <- as.numeric(summary(z2$n)[5])
# 
# # PLOT 
# pdf(paste0(out_dir,"dbcan_HMMER_occurr_piglets_same_genus.pdf"))
# # below 1st Qu.
# z3 <- subset(z2,n < firstQu)
# thelist <- as.character(z3$enzymeID)
# z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
# V <- crossprod(table(z1_sub[1:2]))
# heatmap(V, Rowv = NULL,Colv=NULL,
#         na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
#         main="Co-occurrence - Frequency: below 1st Qu.") 
# # 1st Qu. - median
# z3 <- subset(z2, n >= firstQu & n < median)
# thelist <- as.character(z3$enzymeID)
# z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
# V <- crossprod(table(z1_sub[1:2]))
# heatmap(V, Rowv = NULL,Colv=NULL,
#         na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
#         main="Co-occurrence - Frequency: 1st Qu. - median") 
# # median  - 3rd Qu. 
# z3 <- subset(z2,n >= median & n < thirdQu)
# thelist <- as.character(z3$enzymeID)
# z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
# V <- crossprod(table(z1_sub[1:2]))
# heatmap(V, Rowv = NULL,Colv=NULL,
#         na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
#         main="Co-occurrence - Frequency: median - 3rd Qu.") 
# # above 3rd Qu. 
# z3 <- subset(z2,n >= thirdQu)
# thelist <- as.character(z3$enzymeID)
# z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
# V <- crossprod(table(z1_sub[1:2]))
# heatmap(V, Rowv = NULL,Colv=NULL,
#         na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
#         main="Co-occurrence - Frequency: above 3rd Qu.") 
# # without any subsetting
# z3 <- z2
# thelist <- as.character(z3$enzymeID)
# z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
# V <- crossprod(table(z1_sub[1:2]))
# heatmap(V, Rowv = NULL,Colv=NULL,
#         na.rm = TRUE, scale="none",cexRow = 0.06, cexCol=0.06,
#         main="Co-occurrence - all") 
# dev.off()
# 
# 
# ##################################
# 
# # ENZYMES CO-OCCURRENCE in the same GENUS - mothers 
# 
# z <- gt_hmmer %>% 
#   dplyr::select(enzymeID,pig,genus)
# 
# mothersIDs <- no_reps_all %>% 
#   dplyr::filter(cohort=="Mothers") %>% 
#   dplyr::select(pig) %>% 
#   dplyr::distinct()
# 
# # collect mothers enzyme info only 
# z_moms <- left_join(mothersIDs,z)
# 
# z1 <- z_moms %>% 
#   group_by(genus,enzymeID) %>% 
#   tally() 
# 
# z2 <- z_moms %>% group_by(enzymeID) %>% tally() 
# 
# # set limits 
# firstQu <- as.numeric(summary(z2$n)[2])
# median <- as.numeric(summary(z2$n)[3])
# thirdQu <- as.numeric(summary(z2$n)[5])
# 
# # PLOT 
# pdf(paste0(out_dir,"dbcan_HMMER_occurr_mothers_same_genus.pdf"))
# # below 1st Qu.
# z3 <- subset(z2,n < firstQu)
# thelist <- as.character(z3$enzymeID)
# z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
# V <- crossprod(table(z1_sub[1:2]))
# heatmap(V, Rowv = NULL,Colv=NULL,
#         na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
#         main="Co-occurrence - Frequency: below 1st Qu.") 
# # 1st Qu. - median
# z3 <- subset(z2, n >= firstQu & n < median)
# thelist <- as.character(z3$enzymeID)
# z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
# V <- crossprod(table(z1_sub[1:2]))
# heatmap(V, Rowv = NULL,Colv=NULL,
#         na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
#         main="Co-occurrence - Frequency: 1st Qu. - median") 
# # median  - 3rd Qu. 
# z3 <- subset(z2,n >= median & n < thirdQu)
# thelist <- as.character(z3$enzymeID)
# z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
# V <- crossprod(table(z1_sub[1:2]))
# heatmap(V, Rowv = NULL,Colv=NULL,
#         na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
#         main="Co-occurrence - Frequency: median - 3rd Qu.") 
# # above 3rd Qu. 
# z3 <- subset(z2,n >= thirdQu)
# thelist <- as.character(z3$enzymeID)
# z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
# V <- crossprod(table(z1_sub[1:2]))
# heatmap(V, Rowv = NULL,Colv=NULL,
#         na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
#         main="Co-occurrence - Frequency: above 3rd Qu.") 
# # without any subsetting
# z3 <- z2
# thelist <- as.character(z3$enzymeID)
# z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
# V <- crossprod(table(z1_sub[1:2]))
# heatmap(V, Rowv = NULL,Colv=NULL,
#         na.rm = TRUE, scale="none",cexRow = 0.06, cexCol=0.06,
#         main="Co-occurrence - all")
# dev.off()
# 
# ##################################
# ##################################
# 
# # ENZYMES CO-OCCURRENCE in the same GENUS - mothers AND piglets 
# 
# z <- gt_hmmer %>% dplyr::select(enzymeID,pig,genus)
# 
# IDs <- no_reps_all %>% 
#   dplyr::select(pig) %>% 
#   distinct()
# 
# # collect mothers enzyme info only 
# z <- left_join(IDs,z)
# 
# z1 <- z %>% 
#   group_by(genus,enzymeID) %>% 
#   tally() 
# 
# z2 <- z %>% 
#   group_by(enzymeID) %>% 
#   tally() 
# 
# # PLOT 
# pdf(paste0(out_dir,"dbcan_HMMER_occurr_moms&piglets_same_genus.pdf"))
# # without any subsetting
# z3 <- z2
# thelist <- as.character(z3$enzymeID)
# z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
# V <- crossprod(table(z1_sub[1:2]))
# heatmap(V, Rowv = NULL,Colv=NULL,
#         na.rm = TRUE, scale="none",cexRow = 0.06, cexCol=0.06,
#         main="Co-occurrence of CAZy enzymes in the same genus")
# dev.off()
# 
# 
# ##################################
# ##################################
# 
# # ENZYMES CO-OCCURRENCE in the same SPECIES - mothers AND piglets 
# 
# z <- gt_hmmer %>% 
#   dplyr::select(enzymeID,pig,species)
# 
# IDs <- no_reps_all %>% 
#   dplyr::select(pig) %>% 
#   dplyr::distinct()
# 
# # collect mothers enzyme info only 
# z <- left_join(IDs,z)
# 
# z1 <- z %>% 
#   group_by(species,enzymeID) %>% 
#   tally() 
# 
# z2 <- z %>% 
#   group_by(enzymeID) %>% 
#   tally() 
# 
# # PLOT 
# pdf(paste0(out_dir,"dbcan_HMMER_occurr_moms&piglets_same_species.pdf"))
# # without any subsetting
# z3 <- z2
# thelist <- as.character(z3$enzymeID)
# z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
# V <- crossprod(table(z1_sub[1:2]))
# heatmap(V, Rowv = NULL,Colv=NULL,
#         na.rm = TRUE, scale="none",cexRow = 0.06, cexCol=0.06,
#         main="Co-occurrence of CAZy enzymes in the same species")
# dev.off()
# 
# ##################################
# ##################################
# 
# 
# # ENZYMES CO-OCCURRENCE in the same subject 
# 
# # create an adjancy matrix: rows are species, columns are enzyme IDs
# z <- gt_hmmer %>% 
#   dplyr::select(enzymeID,pig)
# 
# 
# z1 <- z %>% 
#   group_by(pig,enzymeID) %>% 
#   tally() 
# 
# z2 <- z %>% 
#   group_by(enzymeID) %>% 
#   tally() 
# 
# # set limits 
# firstQu <- as.numeric(summary(z2$n)[2])
# median <- as.numeric(summary(z2$n)[3])
# thirdQu <- as.numeric(summary(z2$n)[5])
# 
# # PLOT 
# pdf(paste0(out_dir,"dbcan_HMMER_occurr_same_host.pdf"))
# # below 1st Qu.
# z3 <- subset(z2,n < firstQu)
# thelist <- as.character(z3$enzymeID)
# z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
# V <- crossprod(table(z1_sub[1:2]))
# heatmap(V, Rowv = NULL,Colv=NULL,
#         na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
#         main="Co-occurrence - Frequency: below 1st Qu.") 
# # 1st Qu. - median
# z3 <- subset(z2, n >= firstQu & n < median)
# thelist <- as.character(z3$enzymeID)
# z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
# V <- crossprod(table(z1_sub[1:2]))
# heatmap(V, Rowv = NULL,Colv=NULL,
#         na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
#         main="Co-occurrence - Frequency: 1st Qu. - median") 
# # median  - 3rd Qu. 
# z3 <- subset(z2,n >= median & n < thirdQu)
# thelist <- as.character(z3$enzymeID)
# z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
# V <- crossprod(table(z1_sub[1:2]))
# heatmap(V, Rowv = NULL,Colv=NULL,
#         na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
#         main="Co-occurrence - Frequency: median - 3rd Qu.") 
# # above 3rd Qu. 
# z3 <- subset(z2,n >= thirdQu)
# thelist <- as.character(z3$enzymeID)
# z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
# V <- crossprod(table(z1_sub[1:2]))
# heatmap(V, Rowv = NULL,Colv=NULL,
#         na.rm = TRUE, scale="none",cexRow = 0.3, cexCol=0.3,
#         main="Co-occurrence - Frequency: above 3rd Qu.") 
# # without any subsetting
# z3 <- z2
# thelist <- as.character(z3$enzymeID)
# z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
# V <- crossprod(table(z1_sub[1:2]))
# heatmap(V, Rowv = NULL,Colv=NULL,
#         na.rm = TRUE, scale="none",cexRow = 0.06, cexCol=0.06,
#         main="Co-occurrence of CAZy enzymes in the same host") 
# dev.off()

