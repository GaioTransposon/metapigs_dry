
NROW(df2)
NROW(taxa_mat)
NROW(sample_df)

head(df2)
head(taxa_mat)
head(sample_df)

df2_taxamat <- inner_join(df2,taxa_mat)
head(df2_taxamat)

all <- inner_join(df2_taxamat,sample_df)


sel <- all
sel <- cSplit(indt = sel, "sample", sep = "_", drop = NA)
sel$date = sel$sample_1
# select those with weight data 
sel %>% 
  dplyr::filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|
                  date=="t8"|date=="t10") %>%
  dplyr::filter(species=="Eubacterium callanderi") %>% 
  ggplot(., aes(x=weight,y=log(all_bins_value)))+
  geom_point(size=0.5) +
  facet_wrap(~date, scale="free", shrink = TRUE)+
  stat_smooth(method="lm", se=TRUE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)




sel <- all
sel <- cSplit(indt = sel, "sample", sep = "_", drop = NA)
sel$date = sel$sample_1
sel$n = 1
# select those with weight data 
sell <- sel %>% 
  dplyr::filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|
                  date=="t8"|date=="t10") %>%
  group_by(species,date) %>% 
  dplyr::mutate(num=sum(n)) %>% 
  dplyr::filter(num>20)

require(plyr)
func <- function(xx)
{
  return(data.frame(COR = cor(xx$weight, xx$all_bins_value, method = "spearman")))
}


res <- ddply(sell, .(date,species), func)

View(ress)

# filter out weak correlations 
ress <- res %>% 
  dplyr::filter(!between(COR, -0.3, 0.3))


# fish out from the ress output 
# selecting species and date 
# and plot 


