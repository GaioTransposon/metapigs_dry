
# PROCEED to all: 

# for each sample (pig,date), sum up the counts that fall within one species (same species assigned to distinct bins)
df3 <- df1 %>%
  group_by(pig,gOTU,date) %>%
  dplyr::summarize(norm_value = sum(value)) 

# no lib size normalization 

df3 <- as.data.frame(df3)
df3$sample = paste0(df3$date,"_",df3$pig)
head(df3)

# pivot wider
df3 <- df3 %>%
  dplyr::select(sample,gOTU,norm_value) %>%
  pivot_wider(names_from = sample, values_from = norm_value, values_fill = list(norm_value = 0))

feat <- as.data.frame(df3)
which(is.na(feat[,1]))

rownames(feat) <- feat[,1]
feat[,1] <- NULL

head(feat)
dim(feat)

# is the sum of each columns 1? 
colSums(feat)
# yes 

# ready! 

######################################################################

# CREATE METADATA TABLE (like meta.crc.zeller)

theseAREtheSamples <- as.data.frame(colnames(feat))
colnames(theseAREtheSamples) <- "sample"

df1$sample <- paste0(df1$date,"_",df1$pig)

df1 <- df1 %>%
  dplyr::select(sample,cohort,pig,date,breed,birth_day) %>%
  distinct()

# add anther grouping: date+cohort:
df1$group <- paste0(df1$date,"_",df1$cohort)

# add another grouping: date+breed+birth_day:
df1$group2 <- paste0(df1$date,"_",df1$breed,"_",df1$birth_day)

head(theseAREtheSamples)

meta <- left_join(theseAREtheSamples,df1)

rownames(meta) <- meta[,1]
meta[,1] <- NULL

# ready! 

######################################################################
######################################################################

# SIAMCAT starts! 

class(feat)
class(feat.crc.zeller)
class(meta)
class(meta.crc.zeller)

head(feat)
head(feat.crc.zeller)
head(meta)
head(meta.crc.zeller)

colnames(feat)==rownames(meta)





# counts
rnaseq <- read.csv("~/Downloads/CoDaGuide/rnaseq-x.csv", row.names=1)

# metadata
rnaseq.annot <- read.csv("~/Downloads/CoDaGuide/rnaseq-y.csv", row.names=1)
head(rnaseq.annot)
head(meta)



# transpose 
rnaseq <- t(rnaseq)
feat <- t(feat)

colnames(rnaseq)
colnames(feat)


colSums(feat)

# Now we can replace zeros with a small value
#  the ``p-counts'' option has the function return
#  pseudo-counts instead of proportions
rnaseq.no0 <- cmultRepl(rnaseq, output = "p-counts")

feat.no0 <- cmultRepl(feat, output = "p-counts")




set.seed(1)

# propr expects:
#  `counts': the data matrix with rows as samples
#  `metric': the proportionality metric to calculate
#  `ivar': the log-ratio transform reference
library(propr)

pr <- propr(counts = rnaseq.no0,
            metric = "rho",
            ivar = "clr")

pr_mine <- propr(counts = feat.no0,
            metric = "rho",
            ivar = "clr")


# We can select a good cutoff for ‘rho’
# by permuting the FDR at various cutoffs
# Below, we use [0, .05, ..., .95, 1]
# The d3 = TRUE argument makes the plot 3D -- it is optional
getNetwork(pr, 0.9, col1 = up, col2 = down, d3 = TRUE) # <-- Figure 3

sink("1-fdr-propr-no0.txt")
pr@fdr
sink()

# Let’s visualize using a strict cutoff
getNetwork(pr, cutoff = 0.9, col1 = up)
getResults(pr, cutoff = 0.9)


# The propd function tests for events where the proportionality factor 
# (i.e., the magnitude of x/y ) differs between the experimental groups.

# propd expects:
# ‘counts’: the data matrix with rows as samples
# ‘group’: the class labels
set.seed(1)
library(propr)
pd <- propd(counts = rnaseq.no0,
            group = rnaseq.annot$Treatment)
# Calculate an exact p-value
pd <- updateF(pd)
getResults(pd)

pd <- propd(counts = feat,
            group = meta$date)
# Calculate an exact p-value
pd <- updateF(pd)
getResults(pd)

###################################################################
# This section makes the figure that shows the top propd hits
###################################################################

# Get significant propd results that contain
#  Nfkb1 in the numerator or denominator
df <- pd@results
ref.i <- which(colnames(feat) == "Methanobrevibacter_A smithii__758")
index <- (df$Partner == ref.i & df$lrm1 > df$lrm2) |
  (df$Pair == ref.i & df$lrm1 < df$lrm2)
sum(index)
pd.test <- pd
pd.test@results <- pd.test@results[index,]
getResults(pd.test, .8)

# Get the ratios for plotting
df <- getRatios(pd.test, .8, include = NA, or = TRUE, 
                melt = TRUE)

# Flip axes
df$value <- -1 * df$value
pairs <- sapply(as.character(df$variable), strsplit, split = "/")
df$variable <- sapply(pairs, function(i) paste0(i[2], "/", i[1]))

# Set order of axes and colors
df$variable <- factor(df$variable, levels = unique(df$variable))
df$group <- factor(pd.test@group, levels = c("t0", "t2"))


ggplot2::ggplot(df, ggplot2::aes_string(x = "variable", 
                                        y = "value", group = "id", col = "group")) + 
  ggplot2::geom_line() + 
  ggplot2::theme_bw() + ggplot2::scale_colour_brewer(palette = "Set2", 
                                                     name = "Group") + 
  ggplot2::xlab("Feature Pair") + ggplot2::ylab("Log-Ratio Abundance") + 
  ggplot2::ggtitle("Sample-wise Distribution of Log-Ratios Across Pairs") + 
  ggplot2::theme(axis.text.x = element_blank()) + 
  ggplot2::theme(text = ggplot2::element_text(size = 18), 
                 plot.title = ggplot2::element_text(size = 24))




sort(colnames(feat))
