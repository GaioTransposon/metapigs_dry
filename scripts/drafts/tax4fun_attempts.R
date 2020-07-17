

library(Tax4Fun)
library(themetagenomics)
library(Rcpp)

download_ref(
  destination="/Users/12705859/Desktop/metapigs_base/sortmerna/",
  reference = c("silva_ko"),
  overwrite = FALSE,
  verbose = FALSE
)
load("/Users/12705859/Desktop/metapigs_base/sortmerna/t4f_ref_profiles.rds")


tmp <- tempdir()
download_ref(tmp,reference='silva_ko',overwrite=FALSE)
FUNCTIONS <- predict(TOPICS,reference='silva_ko',reference_path=tmp,cn_normalize=TRUE,sample_normalize=FALSE,scalar=25)
})


library(themetagenomics)
GEVERS$OTU[1:5,1:5]
DAVID$TAX

otu_table1 <- otu_table
mode(otu_table1) <- "integer"
head(otu_table1)
NROW(otu_table1)
NROW(taxa_table)

predicted_functions <- t4f(otu_table=otu_table1,rows_are_taxa=TRUE,
                           tax_table=taxa_table,reference_path='/var/folders/zz/zyxvpxvq6csfxvn_n0ddy18c3bgh_3/T//RtmpRs9YLQ',
                           type='uproc',short=TRUE,cn_normalize=TRUE,
                           sample_normalize=FALSE,scalar=NULL,drop=TRUE)
getwd()
library(biomformat)
myBiom <- make_biom(data = otu_table1,observation_metadata = taxa_table,sample_metadata = sample_df)
biom::write_biom(x = myBiom, biom_file = "/Users/12705859/Desktop/metapigs_base/sortmerna/my_biom")


pleeease <- importQIIMEBiomData(inputFiles = "/Users/12705859/Desktop/metapigs_base/sortmerna/my_biom")

head(pleeease)

t4f(otu_table = pleeease,rows_are_taxa=TRUE,reference_path='/var/folders/zz/zyxvpxvq6csfxvn_n0ddy18c3bgh_3/T//RtmpRs9YLQ',
    drop=TRUE)

