library(readr)
library(dplyr)
library(tidyr)
library(pheatmap)

X14208_ja31 <- 
  read_csv("~/Downloads/host_finding/14208/mge_contact_maps/14208_ja31.ARG.bin.contig.contacts.csv") %>%
  pivot_longer(
    cols = starts_with("k141"),
    names_to = "contig_name",
    values_to = "Ja31",
    values_drop_na = TRUE
  )

X14208_fe7 <- 
  read_csv("~/Downloads/host_finding/14208/mge_contact_maps/14208_fe7.ARG.bin.contig.contacts.csv") %>%
  pivot_longer(
    cols = starts_with("k141"),
    names_to = "contig_name",
    values_to = "Fe7",
    values_drop_na = TRUE
  )

X14208_fe14 <- 
  read_csv("~/Downloads/host_finding/14208/mge_contact_maps/14208_fe14.ARG.bin.contig.contacts.csv") %>%
  pivot_longer(
    cols = starts_with("k141"),
    names_to = "contig_name",
    values_to = "Fe14",
    values_drop_na = TRUE
  )

X14208_fe21 <- 
  read_csv("~/Downloads/host_finding/14208/mge_contact_maps/14208_fe21.ARG.bin.contig.contacts.csv") %>%
  pivot_longer(
    cols = starts_with("k141"),
    names_to = "contig_name",
    values_to = "Fe21",
    values_drop_na = TRUE
  )

X14208_fe28 <- 
  read_csv("~/Downloads/host_finding/14208/mge_contact_maps/14208_fe28.ARG.bin.contig.contacts.csv") %>%
  pivot_longer(
    cols = starts_with("k141"),
    names_to = "contig_name",
    values_to = "Fe28",
    values_drop_na = TRUE
  )

NROW(X14208_ja31)
NROW(X14208_fe7)


z <- full_join(X14208_ja31,X14208_fe7)
z <- full_join(z, X14208_fe14)
z <- full_join(z, X14208_fe21)
z <- full_join(z, X14208_fe28)
z
test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")

z <- as.data.frame(z)
z$sample <- paste0(z$X1,"_",z$contig_name)
rownames(z) <- z$sample
z$sample <- NULL
z$X1 <- NULL
z$contig_name <- NULL
z[is.na(z)] <- 0

z <- z[c("Ja31", "Fe7", "Fe14","Fe21","Fe28")]

m <- as.matrix(z)

m <- m[rowSums(m==0, na.rm=TRUE)<ncol(m), ]
NROW(m)

heatmap(m)
heatmap(m, scale="row")

