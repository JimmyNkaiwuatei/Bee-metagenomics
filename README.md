# Bee-metagenomics
read.table("/cloud/project/otu_table.csv" , sep = ",", header = FALSE, row.names = 1)

readRDS("/cloud/project/taxa_image.rds")

View("otu_table.csv")
View("taxa_image.rds")

# Construct a phyloseq object from dada2 outputs

library(phyloseq)
library(ggplot2)

otu <- read.csv("/cloud/project/otu_table.csv")
otu <- t(otu)

tax <- readRDS("/cloud/project/taxa_image.rds")

meta_data <- read.table("/cloud/project/meta_data.csv" , sep = ",", header = TRUE, row.names = 1)
meta <- sample_data(meta_data)


meta <- read.table("/cloud/project/meta_data.csv" , sep = ",", header = TRUE, row.names = 1)

otu <- read.table("/cloud/project/otu_table.csv" , sep = ",", header = TRUE, row.names = 1)

tax <- readRDS("/cloud/project/taxa_image.rds")

otudat = otu_table(as.matrix(otu) ,taxa_are_rows = TRUE)

taxdat = tax_table(as.matrix(tax))


physeq <- phyloseq(transotu, taxdat)
physeq

physeq <-phyloseq(otudat, taxdat)
plot_bar(physeq, fill = NULL, title = TRUE, facet_grid = NULL)
