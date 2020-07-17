
# script to merge figures produced for metapigs_dry 

library(ggpubr) 

setwd("/Users/12705859/Desktop/metapigs_dry/figures")

# must have run:
# - checkm_phyloseq.R
# - dRep_phyloseq.R
# - gtdbtk_4_phyloseq.R


######################

# ordination & network

all_ordination <- ggarrange(cm_ordination_plot,
          dRep_ordination_plot,
          gt_ordination_plot,
          nrow=1,
          ncol=3,
          labels=c("A","B","C"),
          common.legend=TRUE)

all_network <- ggarrange(cm_network_plot,
          dRep_network_plot,
          gt_network_plot,
          nrow=1,
          ncol=3,
          labels=c("D","E","F"),
          common.legend=TRUE)

all_ordination_network <- ggarrange(all_ordination,
                                    all_network, 
                                    nrow=2,
                                    common.legend=FALSE)

pdf("all_phylo_clustering.pdf")
all_ordination_network
dev.off()


######################





