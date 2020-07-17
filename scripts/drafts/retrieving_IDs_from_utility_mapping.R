library(dplyr)
library(tidyr)
library(data.table)
library(splitstackshape)
library(ggplot2)

map_level4ec_uniref90 <- read_table2("Desktop/map_level4ec_uniref90.txt", 
                                     col_names = FALSE)

map_level4ec_uniref90 <- as.data.frame(map_level4ec_uniref90)
long <- map_level4ec_uniref90 %>%
  pivot_longer(cols=-X1)
long_uniref90IDs <- as.data.frame(unique(long[,3]))
fwrite(x=long_uniref90IDs, file="/Users/12705859/Desktop/uniref90_ec_filtered_IDs.txt", col.names = FALSE)

