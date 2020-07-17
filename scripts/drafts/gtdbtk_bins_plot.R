
############################################################################

setwd("~/Desktop/metapigs_dry/gtdbtk")
basedir = "~/Desktop/metapigs_dry/gtdbtk/"

bac120_ar122_metadata_r89 <- read_delim("~/Desktop/metapigs_dry/bac120_ar122_metadata_r89.tsv", 
                                        "\t", escape_double = FALSE, trim_ws = TRUE)

bac120_ar122_metadata_r89_reduced <- bac120_ar122_metadata_r89 %>%
  select(checkm_marker_lineage,
         gtdb_taxonomy,
         lsu_silva_23s_taxonomy,
         ncbi_organism_name,
         ncbi_species_taxid,
         ncbi_taxid,
         ncbi_taxonomy,
         ncbi_taxonomy_unfiltered,
         ssu_gg_taxonomy,
         ssu_silva_taxonomy)

unique(bac120_ar122_metadata_r89_reduced$gtdb_taxonomy)



taxa_NCBIed <- read_csv(paste0(basedir,"ncbi_retrieved_taxa"),col_names = TRUE)

unique(taxa_NCBIed$taxa_simple) #534
unique(classifiable$name) #570
570-534 #=36

colnames(taxa_NCBIed)[colnames(taxa_NCBIed) == 'taxa_simple'] <- 'name'

taxa_to_search <- unique(classifiable$name)
taxa_to_search <- as.data.frame(taxa_to_search)
colnames(taxa_to_search)[colnames(taxa_to_search) == 'taxa_to_search'] <- 'name'
unique(taxa_to_search$name)

x <- anti_join(taxa_to_search,taxa_NCBIed)
unique(x$name)
# split to get the genus name only 
xx <- cSplit(x, "name"," ")
xxx <- as.data.frame(unique(xx$name_1))
colnames(xxx)[colnames(xxx) == 'unique(xx$name_1)'] <- 'name'


all_taxa_reduced <- cSplit(bac120_ar122_metadata_r89_reduced, "ncbi_organism_name"," ")
all_taxa_reduced <- all_taxa_reduced %>%
  select(ncbi_organism_name_01,ncbi_taxonomy)
all_taxa_reduced <- all_taxa_reduced %>%
  group_by(ncbi_organism_name_01) %>%
  slice(1)    # slice 1 as you don't want multiple hits for the same genus, just 1 per distinct genus ()
z <- merge(xxx,all_taxa_reduced,by.x="name",by.y="ncbi_organism_name_01")
head(z)
colnames(z)

# merge higher taxonomy info (retreived from taxize) to the bins 
most <- inner_join(classifiable,taxa_NCBIed)
NROW(most)
head(most)

# these are the bins for which higher taxonomic levels couldn't be retrieved with taxize
difficult <- anti_join(classifiable,taxa_NCBIed)
NROW(difficult)

# get the genus level only instead (discard species resolution)
difficult <- cSplit(difficult, "name"," ")
difficult <- difficult %>% 
  rename(name = name_1) %>% 
  select(pig,bin,name)

# this time, to infer higher level taxonomy, use the taxonomic info from "bac120_ar122_metadata_r89_reduced" 
rest <- inner_join(difficult,z)
head(rest)
NROW(rest)

# quite good : 
NROW(classifiable)
NROW(most)+NROW(rest)

# now we'll have to merge these, but the column names won't be identical...
# we need some cleaning in orer to make it look like: 
# COLUMNS: pig, bin, name, domain, phylum, class, order, family, genus
head(most)
head(rest)



# cleaning "most"
most$phylum=paste0(most$`3`,"_",most$`4`)
most <- most %>%
  rename(domain = `2`,
         class = `5`,
         order = `6`,
         family = `7`,
         genus = `8`,
         species = `9`) %>%
  select(pig,bin,name,domain,phylum,class,order,family,genus,species)
unique(most$phylum)


# cleaning "rest"
rest <- cSplit(rest, "ncbi_taxonomy",";")
remove_underscores <- function(x) (gsub(".*__","",x))
rest <- rest %>% mutate_at(c("ncbi_taxonomy_1",
                     "ncbi_taxonomy_2",
                     "ncbi_taxonomy_3",
                     "ncbi_taxonomy_4",
                     "ncbi_taxonomy_5",
                     "ncbi_taxonomy_6",
                     "ncbi_taxonomy_7"), remove_underscores)
colnames(rest) <- c("pig","bin","name","domain","phylum","class","order","family","genus","species")


# finally concatenate "most" and "rest"
all_classified <- rbind(most,rest)
NROW(all_classified)

# just making sure the bins unique to each subject are unique (and not duplicated)
NROW(unique(paste0(all_classified$pig,all_classified$bin)))==NROW(all_classified)
# yes, they are. 


domain_phylum <- all_classified %>%
  select(name,domain,phylum) 


domain_phylum$phylum <- gsub("\\<Firmicutes\\>","Terrabacteria group_Firmicutes",domain_phylum$phylum)
domain_phylum$phylum <- gsub("\\<Actinobacteria\\>","Terrabacteria group_Actinobacteria",domain_phylum$phylum)
domain_phylum$phylum <- gsub("\\<Spirochaetes_Spirochaetia\\>","Spirochaetes",domain_phylum$phylum)
domain_phylum$phylum <- gsub("\\<Fusobacteria_Fusobacteriia\\>","Fusobacteria",domain_phylum$phylum)
domain_phylum$phylum <- gsub("\\<Synergistetes_Synergistia\\>","Synergistetes",domain_phylum$phylum)
domain_phylum$phylum <- gsub("\\<Deferribacteres_Deferribacteres\\>","Deferribacteres",domain_phylum$phylum)

domain_phylum$phylum <- gsub("group"," ",domain_phylum$phylum)
domain_phylum$phylum <- gsub("_","\n",domain_phylum$phylum)

unique(domain_phylum$phylum)


# most abundant (>0.45%)
most_ab_counts <- setDT(domain_phylum)[, .(Freq = .N), by = .(phylum)]
# remove lower those that appear less than 3, calculate percentage
most_ab_counts <- most_ab_counts %>%
  filter(!Freq<2) %>%
  filter(!phylum=="NA\nNA") %>%
  mutate(perc=round(Freq/sum(Freq)*100,2)) %>%
  # most abundant 
  filter(perc>1)
most_ab_counts$label <- paste(paste(most_ab_counts$phylum,most_ab_counts$perc,sep = "\n"),"%")

# least abundant (<0.45%)
least_ab_counts <- setDT(domain_phylum)[, .(Freq = .N), by = .(phylum)]
# remove lower those that appear less than 3, calculate percentage
least_ab_counts <- least_ab_counts %>%
  filter(!Freq<2) %>%
  filter(!phylum=="NA\nNA") %>%
  mutate(perc=round(Freq/sum(Freq)*100,2)) %>%
  # most abundant 
  filter(perc<1)
least_ab_counts$label <- paste(paste(least_ab_counts$phylum,least_ab_counts$perc,sep = "\n"),"%")



pdf("treemap_gtdbtk_phyla.pdf")
treemap(most_ab_counts, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Phyla distribution from all MAGs (gtdbtk db) - most abundant", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
treemap(least_ab_counts, #Your data frame object
        index=c("label"),  #A list of your categorical variables
        vSize = "Freq",  #This is your quantitative variable
        type="index", #Type sets the organization and color scheme of your treemap
        title="Phyla distribution from all MAGs (gtdbtk db) - least abundant", #Customize your title
        fontsize.title = 15 #Change the font size of the title
        #fontsize.labels = 8
)
dev.off()

