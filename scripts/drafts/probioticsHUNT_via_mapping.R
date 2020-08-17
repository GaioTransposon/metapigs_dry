library(readr)
library(dplyr)

new_vs_T12 <- read_delim("Desktop/new_vs_T12.tsv", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)

z <- new_vs_T12 %>% 
  arrange(desc(X3))

View(z)


z <- new_vs_T12 %>%
  group_by(X1) %>%
  dplyr::mutate(diff = X2 - lag(X2, default = X2[1]))

# count the consecutive numbers

z1 <- z %>%
  group_by(X1, grp = with(rle(diff), rep(seq_along(lengths), lengths))) %>%
  mutate(Counter = seq_along(grp)) %>%
  ungroup() %>%
  select(-grp)

View(z1)  
