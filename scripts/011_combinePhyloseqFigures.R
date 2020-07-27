# Part 1: 
# combine phyloseq ordination and network outputs from 
# CheckM
# dRep
# GTDB

library(ggpubr)

out_dir = "/Users/12705859/Desktop/metapigs_dry/"

top <- ggarrange(cm_ordination_plot,
          dRep_ordination_plot,
          gt_ordination_plot,
          labels=c("A","C","E"),
          common.legend=TRUE,
          ncol=3)

bottom <- ggarrange(cm_network_plot,
          dRep_network_plot,
          gt_network_plot,
          labels=c("B","D","F"),
          common.legend=TRUE,
          ncol=3)

top_bottom <- ggarrange(top,bottom,nrow=2)

pdf(paste0(out_dir,"all_phylo_clustering.pdf"))
top_bottom
dev.off()


