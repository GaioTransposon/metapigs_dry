

# understanding clr (centered-log ratio)

library(tidyverse)
library(ggplot2)
library(dplyr)
library(data.table)
library(compositions)


setwd("/Users/12705859/Desktop/metapigs_dry/drafts/")


######################################################################################################
######################################################################################################

# outputs generated: 
# 1.1 raw data, clr, plot - parallel coordinates
# 1.2 raw data, clr, plot - barplot
# 2.1 raw data, mean, clr, plot - parallel coordinates
# 2.2 raw data, mean, clr, plot - barplot
# 3.1 raw data, plot - parallel coordinates
# 3.2 raw data, plot - barplot
# 4.1 raw data, mean, plot - parallel coordinates
# 4.2 raw data, mean, plot - barplot

# the first part will use the data from the library "compositions"
# the second part uses my microbiome (phyla) data df 
# df is generated in checkm_output_analysis.R and looks like this: 

# cohort    pig      bin secondary_cluster date       value      taxa_1           taxa_2             taxa_3
# 1: DScour  14159 bins.100             482_2   t0    297.3006 k__Bacteria             <NA>               <NA>
#   2: DScour  14159 bins.100             482_2   t1    126.2060 k__Bacteria             <NA>               <NA>
#   3: DScour  14159 bins.100             482_2  t10    253.3927 k__Bacteria             <NA>               <NA>
#   4: DScour  14159 bins.100             482_2   t2    351.3063 k__Bacteria             <NA> 


######################################################################################################
######################################################################################################

# 1st part ("compositions" library data)

data(Hydrochem)

Hydrochem2 <- Hydrochem %>%
  dplyr::select(Code,Date,H,Na,K,Mg,Ca,Sr,Ba,NH4,Cl,NO3,PO4,SO4,HCO3,TOC)

Hydrochem2 <- as.data.frame(Hydrochem2)

Hydrochem2 <- reshape2::melt(data = Hydrochem2, id.vars = c("Code","Date"), 
                             measure.vars = c("H","Na","K","Mg","Ca","Sr","Ba","NH4","Cl","NO3","PO4","SO4","HCO3","TOC"))


# 1. raw, clr, plot

Hydrochem3 <- Hydrochem2 %>%
  mutate(value=clr(value))

jpeg("1.1.jpg")
ggplot(Hydrochem3, aes(x=Date, y=value, group=variable, color=variable)) + 
  geom_line() + geom_point(size=0.8)
dev.off()

jpeg("1.2.jpg")
ggplot(Hydrochem3, aes(fill=variable, y=value, x=Date)) + 
  geom_bar(position="fill", stat="identity")
dev.off()

# 2. raw, mean, clr, plot

Hydrochem3 <- Hydrochem2 %>%
  group_by(variable,Date) %>%
  dplyr::summarise(value=mean(value)) 

Hydrochem3 <- as.data.frame(Hydrochem3)

Hydrochem3 <- Hydrochem3 %>%
  mutate(value=clr(value))

Hydrochem3 <- as.data.frame(Hydrochem3)

jpeg("2.1.jpg")
ggplot(Hydrochem3, aes(x=Date, y=value, group=variable, color=variable)) + 
  geom_line() + geom_point(size=0.8)
dev.off()

jpeg("2.2.jpg")
ggplot(Hydrochem3, aes(fill=variable, y=value, x=Date)) + 
  geom_bar(position="fill", stat="identity")
dev.off()

# 3. raw, plot


Hydrochem3 <- Hydrochem2 

Hydrochem3 <- as.data.frame(Hydrochem3)

jpeg("3.1.jpg")
ggplot(Hydrochem3, aes(x=Date, y=value, group=variable, color=variable)) + 
  geom_line() + geom_point(size=0.8)
dev.off()

jpeg("3.2.jpg")
ggplot(Hydrochem3, aes(fill=variable, y=value, x=Date)) + 
  geom_bar(position="fill", stat="identity")
dev.off()

# 4. 

Hydrochem3 <- Hydrochem2 %>%
  group_by(variable,Date) %>%
  dplyr::summarise(value=mean(value)) 

Hydrochem3 <- as.data.frame(Hydrochem3)

Hydrochem3 <- as.data.frame(Hydrochem3)

jpeg("4.1.jpg")
ggplot(Hydrochem3, aes(x=Date, y=value, group=variable, color=variable)) + 
  geom_line() + geom_point(size=0.8)
dev.off()

jpeg("4.2.jpg")
ggplot(Hydrochem3, aes(fill=variable, y=value, x=Date)) + 
  geom_bar(position="fill", stat="identity")
dev.off()



######################################################################################################
######################################################################################################

# 2nd part (my data)

Hydrochem2 <- df %>%
  dplyr::select(date,taxa_2,value)

Hydrochem2 <- na.omit(Hydrochem2)




# 1. raw, clr, plot

Hydrochem3 <- Hydrochem2 %>%
  mutate(value=clr(value))

jpeg("1.1.jpg")
ggplot(Hydrochem3, aes(x=date, y=value, group=taxa_2, color=taxa_2)) + 
  geom_line() + geom_point(size=0.8)
dev.off()

jpeg("1.2.jpg")
ggplot(Hydrochem3, aes(fill=taxa_2, y=value, x=date)) + 
  geom_bar(position="fill", stat="identity")
dev.off()

# 2. raw, mean, clr, plot

Hydrochem3 <- Hydrochem2 %>%
  group_by(taxa_2,date) %>%
  dplyr::summarise(value=mean(value)) 

Hydrochem3 <- as.data.frame(Hydrochem3)

Hydrochem3 <- Hydrochem3 %>%
  mutate(value=clr(value))

Hydrochem3 <- as.data.frame(Hydrochem3)

jpeg("2.1.jpg")
ggplot(Hydrochem3, aes(x=date, y=value, group=taxa_2, color=taxa_2)) + 
  geom_line() + geom_point(size=0.8)
dev.off()

jpeg("2.2.jpg")
ggplot(Hydrochem3, aes(fill=taxa_2, y=value, x=date)) + 
  geom_bar(position="fill", stat="identity")
dev.off()

# 3. raw, plot


Hydrochem3 <- Hydrochem2 

jpeg("3.1.jpg")
ggplot(Hydrochem3, aes(x=date, y=value, group=taxa_2, color=taxa_2)) + 
  geom_line() + geom_point(size=0.8)
dev.off()

jpeg("3.2.jpg")
ggplot(Hydrochem3, aes(fill=taxa_2, y=value, x=date)) + 
  geom_bar(position="fill", stat="identity")
dev.off()

# 4. 

Hydrochem3 <- Hydrochem2 %>%
  group_by(taxa_2,date) %>%
  dplyr::summarise(value=mean(value)) 

Hydrochem3 <- as.data.frame(Hydrochem3)

Hydrochem3 <- as.data.frame(Hydrochem3)

jpeg("4.1.jpg")
ggplot(Hydrochem3, aes(x=date, y=value, group=taxa_2, color=taxa_2)) + 
  geom_line() + geom_point(size=0.8)
dev.off()

jpeg("4.2.jpg")
ggplot(Hydrochem3, aes(fill=taxa_2, y=value, x=date)) + 
  geom_bar(position="fill", stat="identity")
dev.off()






