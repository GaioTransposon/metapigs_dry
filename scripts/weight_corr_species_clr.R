
sel_frequent <- sel %>% 
  group_by(species,date) %>% 
  dplyr::mutate(num=sum(n)) %>% 
  dplyr::filter(num>20) 


require(plyr)
func <- function(xx)
{
  return(data.frame(COR = cor(xx$weight, xx$all_bins_value, method = "pearson")) %>% 
           dplyr::filter(!between(COR, -0.3, 0.3))) # filter out weak correlations 
}


res_raw <- ddply(sel_frequent, .(date,species), func)

sel_frequent_clr <- sel_frequent %>% 
  group_by(pig,date) %>%
  dplyr::mutate(all_bins_value=as.numeric(clr(all_bins_value))) 
res_clr <- ddply(sel_frequent_clr, .(date,species), func)

NROW(res_raw)
NROW(res_clr)

View(res_raw)
View(res_clr)

NROW(which(res_raw$species %in% res_clr$species == TRUE))

pdf(paste0(out_dir,"gt_corr_weight_species_clr.pdf"))
for (row in 1:nrow(res_clr)) {
  
  speciess <- res_clr[row,2]
  sel0 <- subset(sel, species %in% speciess) %>%
    dplyr::filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|
                    date=="t8"|date=="t10")
  
  print(sel0 %>%
          ggplot(., aes(x=weight,y=log(all_bins_value)))+
          geom_point(size=0.5) +
          facet_grid(~date, scale="free", shrink = TRUE)+
          stat_smooth(method="lm", se=TRUE) +
          stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size=2)+
          ggtitle(speciess))
  
}
dev.off()
library(compositions)
pdf(paste0(out_dir,"gt_corr_weight_species_clr_also_present_in_raw.pdf"))
keep <- res_raw$species
res_clr2 <- dplyr::filter(res_clr, species %in% keep) 

# keep unique species (all time points will be plotted anyway)
#res_clr3 <- res_clr2 %>% dplyr::select(species) %>% distinct()

# or keep only the ones where more than one time point showed a strong correlation
res_clr3 <- res_clr2 %>% group_by(species) %>% filter(n()>1) %>% dplyr::select(species) %>% distinct()

for (row in 1:nrow(res_clr3)) {
  speciess <- res_clr3[row,1]
  sel0 <- subset(sel, species %in% speciess) %>%
    dplyr::filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|
                    date=="t8"|date=="t10")
  
  print(sel0 %>%
          ggplot(., aes(x=weight,y=as.numeric(clr(all_bins_value))))+
          geom_point(size=0.5) +
          facet_grid(~date, scale="free", shrink = TRUE)+
          stat_smooth(method="lm", se=TRUE) +
          stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size=2)+
          ggtitle(speciess))
  
}
dev.off()
