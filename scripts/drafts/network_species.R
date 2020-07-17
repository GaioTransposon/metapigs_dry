
library(readr)
library(dplyr)
library(tidyr)
library(data.table)
library(splitstackshape)
library(ggplot2)
library(RColorBrewer)
library(forcats)
library(igraph)
library(visNetwork)


# instructions from https://kateto.net/network-visualization


setwd("~/Desktop/metapigs_dry/network_analysis")
basedir = "~/Desktop/metapigs_dry/"

########################

# load secondary_clusters to pig/bin file: 

# upload input file
df <- read.csv(paste0(basedir,"merged_all_clustered_wa_bins_with_cohorts.csv"),
               na.strings=c("","NA"),
               check.names = FALSE,
               header = TRUE)
# bins that don't have secondary cluster assigned are filled with "no_cluster"
df <- df %<>% mutate(secondary_cluster = fct_explicit_na(secondary_cluster, na_level = "no_cluster")) %>% dplyr::select(pig,bin,secondary_cluster)
# remove .fa extension to match bins in checkm df 
df$bin <- gsub(".fa","", df$bin)
head(df)

########################

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv(paste0(basedir,"gtdbtk/gtdbtk_bins_completeTaxa"),
                        col_types = cols(node = col_character(),
                                         pig = col_character()))

########################

# load edges  


edges_moms <- read_table2("fastspar/edges_mothers_99_significant_grouped.csv",
                          col_types=cols(Source=col_double(),
                                         Target=col_double(),
                                         NegPos=col_character(),
                                         Weight=col_double()))

########################

# load nodes

nodes_moms <- read_table2("fastspar/nodes_mothers_99_significant_grouped.csv",
                          col_types=cols(Id=col_double(),
                                         Label=col_character(),
                                         group=col_character()))


########################


# get dataframe with species and upper taxa info of clusters 
z0 <- inner_join(df,gtdbtk_bins) %>%
  filter(!secondary_cluster=="no_cluster") # don't need bins that were not assigned any cluster by dRep 

NROW(unique(z0$secondary_cluster))
head(z0)

z <- z0 %>%
  dplyr::select(secondary_cluster,phylum,class,order,family,genus,species) %>% # selecting cols of interest
  group_by(secondary_cluster,phylum,class,order,family,genus,species) %>%
  tally() %>% # get number of occurrences of species per secondary cluster
  group_by(secondary_cluster) %>%
  top_n(1) %>% # get the top species per secondary_cluster
  slice(1) %>%
  mutate(Label=secondary_cluster)

# do we now have distinct secondary clusters with one and only one species assigned (the most common) ? 
NROW(z)==NROW(unique(z0$secondary_cluster))

NROW(z)
head(z)


########################

# join taxa info to nodes (nodes IDs are my secondary clusters: 

# run after dbcan_hmmer.R

z <- gt_hmmer %>% dplyr::select(enzymeID,pig,genus)
piggiesIDs <- no_reps_all %>% filter(!cohort=="Mothers") %>% dplyr::select(pig) %>% distinct()

# collect mothers enzyme info only 
z_piggies <- left_join(piggiesIDs,z)

z1 <- z_piggies %>% 
  group_by(genus,enzymeID) %>% 
  tally() 

z2 <- z_piggies %>% group_by(enzymeID) %>% tally() 

# set limits 
firstQu <- as.numeric(summary(z2$n)[2])
median <- as.numeric(summary(z2$n)[3])
thirdQu <- as.numeric(summary(z2$n)[5])


# 1st Qu. - median
z3 <- subset(z2, n >= firstQu & n < median)
thelist <- as.character(z3$enzymeID)
z1_sub <- subset(z1, (enzymeID %in% thelist)) %>% drop.levels()
V <- crossprod(table(z1_sub[1:2]))
heatmap(V, Rowv = NULL,Colv=NULL,
        na.rm = TRUE, scale="none", cexRow = 0.3, cexCol=0.3,
        main="Co-occurrence - Frequency: 1st Qu. - median") 



# graphing
g <- graph.adjacency(V, weighted=TRUE, mode ='undirected')
g <- simplify(g)
#plot(g)

# set labels and degrees of vertices
V(g)$label <- V(g)$name
V(g)$degree <- degree(g)
#plot(g)

# Set edge width based on weight:
E(g)$width <- E(g)$weight  # log(2+E(g)$weight)
hist(E(g)$width) # visualize distribution
#plot(g)

# keeping only the most important ties and discarding the rest
#plot(g,
#     vertex.label=NA)
cut.off <- mean(V(g)$degree)
g <- delete_edges(g, V(g)[degree<cut.off])
plot(g,
     vertex.label=NA)

# Generate colors based on enzyme class:
enzymeNAME <- (str_extract(V(g)$label, "[aA-zZ]+")) # extract characters only
colrs <- brewer.pal(NROW(unique(enzymeNAME)), "Spectral")
V(g)$color <- colrs[as.factor(enzymeNAME)]
plot(g)



# l <- layout_with_fr(g, dim=3)
# plot(g, layout=l,vertex.label=NA)
# 
# l <- layout_with_kk(g)
# plot(g, layout=l,vertex.label=NA)
# 
# l <- layout_with_graphopt(g)
# plot(g, layout=l,vertex.label=NA)


# The charge parameter below changes node repulsion:
l2 <- layout_with_graphopt(g,
                           charge=0.00001)   # charge=0.00000001
plot(g,
     layout=l2)












V_long <- data.frame(from=rownames(V)[row(V)], to=colnames(V)[col(V)],
                     weight=c(V))
mylinks <- V_long 

mynodes <- z %>% 
  group_by(enzymeID) %>% 
  tally() %>%
  mutate(enzymeFreq=n) %>%
  mutate(getclass = enzymeID) %>%
  separate(getclass, into = c("text", "num"),
           sep = "(?<=[A-Za-z])(?=[0-9])") %>% 
  mutate(enzymeNAME=text) %>% 
  dplyr::select(enzymeID,enzymeNAME,enzymeFreq)


vis.nodes <- mynodes
vis.links <- mylinks


vis.nodes$shape  <- "dot"  
vis.nodes$shadow <- TRUE # Nodes will drop shadow
vis.nodes$title  <- vis.nodes$enzymeFreq # Text on click
vis.nodes$label  <- vis.nodes$enzymeID # Node label
vis.nodes$size   <- log(vis.nodes$enzymeFreq) # Node size
vis.nodes$borderWidth <- 2 # Node border width
colrs <- brewer.pal(NROW(unique(vis.nodes$enzymeNAME)), "Spectral")
vis.nodes$color   <- colrs[as.factor(vis.nodes$enzymeNAME)]

vis.links$width <- log(10+(vis.links$weight)) # line width
hist(log(vis.links$weight+10))
vis.links$color <- "black"    # line color 
vis.links$arrows <- "middle" # arrows: 'from', 'to', or 'middle'
vis.links$smooth <- FALSE    # should the edges be curved?
vis.links$shadow <- TRUE    # edge shadow


# # For more information, you can also check out:
#   
# ?visOptions # available options 
# ?visLayout  # available layouts
# ?visGroups  # using node groups
# ?visLegend  # adding a legend

# links still not showing 
visNetwork(vis.nodes, vis.links, height = "700px") %>% # width = "100%"
  visOptions(selectedBy = "enzymeNAME", 
             highlightNearest = TRUE, 
             nodesIdSelection = TRUE) %>%
  visPhysics(stabilization = FALSE)

#detach('package:visNetwork')





library(threejs)
library(htmlwidgets)
library(igraph)

net.js <- g
graph_attr(net.js, "layout") <- NULL 
gjs <- graphjs(net.js, main="Network!", bg="gray10", showLabels=F, stroke=F, 
               curvature=0.1, attraction=0.9, repulsion=0.8, opacity=0.9)
print(gjs)

gjs.an <- graphjs(net.js, bg="gray10", showLabels=T, stroke=F, 
                  
                  layout=list(#layout_randomly(net.js, dim=3)))   # too crowded
                    layout_with_fr(net.js,  dim=3)))  # good
#layout_with_drl(net.js, dim=3)))  # good
#layout_on_sphere(net.js)))    # not good
# vertex.color=list(V(net.js)$color, "gray", "orange", 
#                   V(net.js)$color))
print(gjs.an)



# library(networkD3)
# 
# links.d3 <- data.frame(from=as.numeric(factor(mylinks$from))-1, 
#                        to=as.numeric(factor(mylinks$to))-1 )
# mynodes <- as.data.frame(mynodes)
# colnames(mynodes) <- c("idn","enzymeNAME","enzymeFreq")
# nodes.d3 <- cbind(idn=factor(mynodes$idn, levels=mynodes$idn), mynodes) 
# 
# forceNetwork(Links = links.d3, Nodes = nodes.d3, Source="from", Target="to",
#              NodeID = "idn",linkWidth = 1,  Group = "enzymeNAME",
#              linkColour = "#afafaf", fontSize=12, zoom=T, legend=T,
#              opacity = 0.8, charge=-300, Nodesize=6, <- nodeSize is not node size, it's the number of column for node size 
#              width = 600, height = 400)