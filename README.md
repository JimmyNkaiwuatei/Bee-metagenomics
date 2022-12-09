
# bee metagenomic analysis for fungal determination
## 1. step1. dada2 analysis

## loading listed files usinng dada2 package

```{r}
library("dada2")
path <-"/home/jimmy/BeeMetagenomics_dada2_results"
list.files(path)
```
# filter and trim
# Sort ensures forward/reverse reads are in same order

```{r}
fnFs <- sort(list.files(path, pattern="_1.fastq"))
fnFs
fnRs <- sort(list.files(path, pattern="_2.fastq"))
fnRs
```
# Extracting sample names

```{r}
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
sample.names
```

# Specify the full path to the fnFs and fnRs
```{r}
fnFs <- file.path(path, fnFs)
fnFs
fnRs <- file.path(path, fnRs)
fnRs
```

# visualizing the quality profiles of the forward reads
```{r}
plotQualityProfile(fnFs[1:2])
```
 plotQualityProfile(fnFs[1:2])
# visualizing the quality profiles of the reverse reads
```{r}
plotQualityProfile(fnRs[1:2])
```
# Performing filtration and trimming: to truncate the reads and trim the nucleotides where the quality profiles of the reads crashes. 
# Define the filenames for the filtered fastq.gz files:
# Place filtered files in filtered/ subdirectory
```{r}
filt_path <- file.path(path, "filtered") 
filt_path
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtFs
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))
filtRs
```

```{r}
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,180),
maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
compress=TRUE, multithread=TRUE)
out
```
```{r}
errF <- learnErrors(filtFs, multithread=TRUE)
```
```{r}
errR <- learnErrors(filtRs, multithread=TRUE)
```
```{r}
plotErrors(errF, nominalQ=TRUE)
```
```{r}
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```
```{r}
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```


```{r}
dadaFs[[1]]
```
```{r}
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```
```{r}
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```
```{r}
table(nchar(getSequences(seqtab)))
```
```{r}
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
```
```{r}
sum(seqtab.nochim)/sum(seqtab)
```
```{r}
getN <- function(x) sum(getUniques(x))
getN
```
```{r}
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
track
```
```{r}
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
colnames
```
```{r}
rownames(track) <- sample.names
rownames(track)
```
```{r}
head(track)
```
```{r}
taxa <- assignTaxonomy(seqtab.nochim, tryRC=TRUE,  "/home/jimmy/BeeMetagenomics_dada2_results/sh_general_release_10.05.2021/sh_general_release_dynamic_10.05.2021_dev.fasta", multithread=TRUE)
```
# Removing sequence rownames for display only
```{r}
taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)
```

```{r}
save.image(file = "/home/jimmy/result/dada2_assigne_tax_image.RData")
```
```{r}
saveRDS(taxa, "taxa_table.rds")
```
```{r}
saveRDS(seqtab.nochim, "otu_table.rds")
```

```{r}
load("/home/jimmy/result/dada2_assigne_tax_image.RData")
```

##visualization
##The R script 
# metagenomic analysis using dada2

# loading of packages

```{r}
library("phyloseq")
library("vegan")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library("janitor")
library("coin")
library("reshape2")
library("ggnewscale")
library("MicrobiotaProcess")
library("patchwork")
library("dada2")
library("csv")
library("dplyr")
library("ggtree")
library("VennDiagram")
library("UpSetR")
#BiocManager::install("MicrobiotaProcess")
install.packages("MicrobiotaProcess")
library("MicrobiotaProcess")
library("utils")
library("utf8")
library("metagMisc")
```


```{r}
#BiocManager::install("devtools")
library("devtools")
library("tidyverse")
#devtools::install_github(repo= "jbisanz/qiime2R")
library(qiime2R)
asv_table <- read.table("/home/icipe/Downloads/dada2/ASV_table.tsv")
```

# loading the data

```{r}
metadata <- read.table("/home/icipe/Downloads/Geo_metadata.tsv", sep = "\t", row.names = 1, header= TRUE)

ot <- read.table("/home/icipe/Downloads/dada2/ASV_table.tsv", header = TRUE, sep = "\t",  row.names = 1 )
class(ot)
ot <- as.matrix(ot)
class(ot)

tax <- read.table("/home/icipe/Downloads/dada2/ASV_tax_species.tsv", header = TRUE, sep = "\t",  row.names = 1)
#tax <- data.frame(tax)
Tax <- as.matrix(tax)
#tax <- t(tax)

ot <- phyloseq::otu_table(ot, taxa_are_rows = TRUE)
tax <- phyloseq::tax_table(Tax)
metadata <- phyloseq::sample_data(metadata)

```
```{r}
ps <- phyloseq(ot, metadata, tax)
ps
```
```{r}
# we now create a third object called random for merging with the other three object
library("ape")
random_tree <- rtree(ntaxa(physeq), rooted = TRUE, tip.label = taxa_names(physeq))
plot(random_tree)
```

```{r}
#merging the preceeding 3 objects.
ps <- merge_phyloseq(physeq, random_tree)
ps
```

```{r}
# Removing the NAs
tx <- phyloseq_to_df(ps, addtax = T, addtot = F, addmaxrank = F)
tx <-tx[, c(8, 12:150 )]

tx <- tx %>%
  # recode empty strings "" by NAs
 na_if("") %>%
  # remove NAs
  na.omit()
tx
```

```{r}
cumulation <- tx %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation
```
```{r}
#merging the blast taxonomic classification to blast abundance table
merged_data <- tx
write.csv(merged_data, file="/home/icipe/Downloads/dada2/merged_data.csv")
# Group the data
grouped_data <- tx %>% group_by(Genus) %>%summarise_if(is.numeric, sum)
#grouped_data<-grouped_data[-c(1),]
```


```{r}
saveRDS(grouped_data,"/home/icipe/Downloads/dada2/grouped_data.rds")
write.table(tx,"/home/icipe/Downloads/dada2/tx.csv")
```


```{r}
# Grouping data by location
Baringo<-grouped_data[,c(1,2:4)]
Bungoma<-grouped_data[,c(1,5:10)]
Bomet<-grouped_data[,c(1,11:16)]
Busia<-grouped_data[,c(1,17:19)]
Baringo2<-grouped_data[,c(1,20:23)]
Coast<-grouped_data[,c(1,24:26)]
Elgon<-grouped_data[,c(1,27:31)]
Eldama_Ravine<-grouped_data[,c(1,32:34)]
Kakamega<-grouped_data[,c(1,35:46)]
ICIPE<-grouped_data[,c(1,47:54)]
Isiolo<-grouped_data[,c(1,55:62)]
Kasanga_Mwingi<-grouped_data[,c(1,63:67)]
Kasika_Mwingi<-grouped_data[,c(1,68:72)]
Kib_Kakamega<-grouped_data[,c(1,73:74)]
Kisii<-grouped_data[,c(1,75:88)]
Kitui<-grouped_data[,c(1,89:93)]
KNN_Nyamira<-grouped_data[,c(1,94)]
maralal<-grouped_data[,c(1,95:98)]
Maru_Nairobi<-grouped_data[,c(1,99:100)]
Meru<-grouped_data[,c(1,101:106)]
Narok<-grouped_data[,c(1,107:109)]
phi<-grouped_data[,c(1,110)]
saf_Nairobi<-grouped_data[,c(1,111:116)]
Taita<-grouped_data[,c(1,117:123)]
Thog_Nairobi<-grouped_data[,c(1,124:126)]
Uasin_Gishu<-grouped_data[,c(1,127:132)]
Ukasic_Mwingi<-grouped_data[,c(1,133:135)]
Vin_Kakamega<-grouped_data[,c(1,136:139)]
Voic_Taita <- grouped_data[1,140]
```

```{r}
all<-grouped_data[,c(1,2:140)]
all
```

```{r}
Baringo_total <-Baringo%>% adorn_totals(c("col"))
Baringo_total <- mutate(Baringo_total, Baringo=rowSums(Baringo_total[5])/4)
Baringo_total <- Baringo_total[,c(1,6)]

Bungoma_total <-Bungoma%>% adorn_totals(c("col"))
Bungoma_total <- mutate(Bungoma_total, Bungoma=rowSums(Bungoma_total[8])/7)
Bungoma_total <- Bungoma_total[,c(1,9)]

Bomet_total <-Bomet%>% adorn_totals(c("col"))
Bomet_total <- mutate(Bomet_total, Bomet=rowSums(Bomet_total[8])/7)
Bomet_total <- Bomet_total[,c(1,9)]

Busia_total<-Busia%>% adorn_totals(c("col"))
Busia_total <- mutate(Busia_total, Busia=rowSums(Busia_total[5])/4)
Busia_total <- Busia_total[,c(1,6)]

Baringo2_total<-Baringo2%>% adorn_totals(c("col"))
Baringo2_total <- mutate(Baringo2_total, Baringo2=rowSums(Baringo2_total[6])/5)
Baringo2_total <- Baringo2_total[,c(1,7)]

Coast_total<-Coast%>% adorn_totals(c("col"))
Coast_total <- mutate(Coast_total, Coast=rowSums(Coast_total[5])/4)
Coast_total <- Coast_total[,c(1,6)]

Elgon_total<-Elgon%>% adorn_totals(c("col"))
Elgon_total <- mutate(Elgon_total, Elgon=rowSums(Elgon_total[7])/6)
Elgon_total <- Elgon_total[,c(1,8)]

Eldama_Ravine_total<-Eldama_Ravine%>% adorn_totals(c("col"))
Eldama_Ravine_total <- mutate(Eldama_Ravine_total, Eldama_Ravine=rowSums(Eldama_Ravine_total[5])/4)
Eldama_Ravine_total <- Eldama_Ravine_total[,c(1,6)]

Kakamega_total<-Kakamega%>% adorn_totals(c("col"))
Kakamega_total <- mutate(Kakamega_total, Kakamega=rowSums(Kakamega_total[14])/13)
Kakamega_total <- Kakamega_total[,c(1,15)]

ICIPE_total <-ICIPE%>% adorn_totals(c("col"))
ICIPE_total <- mutate(ICIPE_total, ICIPE=rowSums(ICIPE_total[10])/9)
ICIPE_total <- ICIPE_total[,c(1,11)]

Isiolo_total<-Isiolo%>% adorn_totals(c("col"))
Isiolo_total <- mutate(Isiolo_total, Isiolo=rowSums(Isiolo_total[10])/9)
Isiolo_total <- Isiolo_total[,c(1,11)]

Kasanga_Mwingi_total <- Kasanga_Mwingi%>%adorn_totals(c("col"))
Kasanga_Mwingi_total<- mutate(Kasanga_Mwingi_total, Kasanga_Mwingi=rowSums(Kasanga_Mwingi_total[7])/6)
Kasanga_Mwingi_total <- Kasanga_Mwingi_total[,c(1,8)]

Kasika_Mwingi_total<-Kasika_Mwingi%>% adorn_totals(c("col"))
Kasika_Mwingi_total <- mutate(Kasika_Mwingi_total, Kasika_Mwingi=rowSums(Kasika_Mwingi_total[7])/6)
Kasika_Mwingi_total <- Kasika_Mwingi_total[,c(1,8)]

Kib_Kakamega_total<-Kib_Kakamega%>% adorn_totals(c("col"))
Kib_Kakamega_total <- mutate(Kib_Kakamega_total, Kib_Kakamega=rowSums(Kib_Kakamega_total[4])/3)
Kib_Kakamega_total <- Kib_Kakamega_total[,c(1,5)]

Kisii_total<-Kisii%>% adorn_totals(c("col"))
Kisii_total <- mutate(Kisii_total, Kisii=rowSums(Kisii_total[16])/15)
Kisii_total <- Kisii_total[,c(1,17)]

Kitui_total<-Kitui%>% adorn_totals(c("col"))
Kitui_total <- mutate(Kitui_total, Kitui=rowSums(Kitui_total[7])/6)
Kitui_total <- Kitui_total[,c(1,8)]

KNN_Nyamira_total<-KNN_Nyamira%>% adorn_totals(c("col"))
KNN_Nyamira_total <- mutate(KNN_Nyamira_total, KNN_Nyamira=rowSums(KNN_Nyamira_total[3])/2)
KNN_Nyamira_total <- KNN_Nyamira_total[,c(1,4)]

maralal_total<-maralal%>% adorn_totals(c("col"))
maralal_total <- mutate(maralal_total, maralal=rowSums(maralal_total[6])/5)
maralal_total <- maralal_total[,c(1,7)]

Maru_Nairobi_total<-Maru_Nairobi%>% adorn_totals(c("col"))
Maru_Nairobi_total <- mutate(Maru_Nairobi_total, Maru_Nairobi=rowSums(Maru_Nairobi_total[4])/3)
Maru_Nairobi_total <- Maru_Nairobi_total[,c(1,5)]

Meru_total<-Meru%>% adorn_totals(c("col"))
Meru_total <- mutate(Meru_total, Meru=rowSums(Meru_total[8])/7)
Meru_total <- Meru_total[,c(1,9)]

Narok_total<-Narok%>% adorn_totals(c("col"))
Narok_total <- mutate(Narok_total, Narok=rowSums(Narok_total[5])/4)
Narok_total <- Narok_total[,c(1,6)]

phi_total<-phi%>% adorn_totals(c("col"))
phi_total <- mutate(phi_total, phi=rowSums(phi_total[3])/2)
phi_total <- phi_total[,c(1,4)]

saf_Nairobi_total<-saf_Nairobi%>% adorn_totals(c("col"))
saf_Nairobi_total <- mutate(saf_Nairobi_total, saf_Nairobi=rowSums(saf_Nairobi_total[8])/7)
saf_Nairobi_total <- saf_Nairobi_total[,c(1,9)]

Taita_total<-Taita%>% adorn_totals(c("col"))
Taita_total <- mutate(Taita_total, Taita=rowSums(Taita_total[9])/8)
Taita_total <- Taita_total[,c(1,10)]

Thog_Nairobi_total<-Thog_Nairobi%>% adorn_totals(c("col"))
Thog_Nairobi_total <- mutate(Thog_Nairobi_total, Thog_Nairobi=rowSums(Thog_Nairobi_total[5])/4)
Thog_Nairobi_total <- Thog_Nairobi_total[,c(1,6)]

Uasin_Gishu_total<-Uasin_Gishu%>% adorn_totals(c("col"))
Uasin_Gishu_total <- mutate(Uasin_Gishu_total, Uasin_Gishu=rowSums(Uasin_Gishu_total[8])/7)
Uasin_Gishu_total <- Uasin_Gishu_total[,c(1,9)]

Ukasic_Mwingi_total<-Ukasic_Mwingi%>% adorn_totals(c("col"))
Ukasic_Mwingi_total <- mutate(Ukasic_Mwingi_total, Ukasic_Mwingi=rowSums(Ukasic_Mwingi_total[5])/4)
Ukasic_Mwingi_total <- Ukasic_Mwingi_total[,c(1,6)]

Vin_Kakamega_total<-Vin_Kakamega%>% adorn_totals(c("col"))
Vin_Kakamega_total <- mutate(Vin_Kakamega_total, Vin_Kakamega=rowSums(Vin_Kakamega_total[6])/5)
Vin_Kakamega_total <- Vin_Kakamega_total[,c(1,7)]

Voic_Taita_total <-Voic_Taita%>% adorn_totals(c("col"))
Voic_Taita_total <- mutate(Voic_Taita_total, Voic_Taita=rowSums(Voic_Taita_total[3])/2)
Voic_Taita_total <- Voic_Taita_total[,c(1,4)]
```

```{r}
all_total <- all %>% adorn_totals(c("col"))
all_total <- mutate(all_total, all=rowSums(all_total[141])/140)
all_total <- all_total[,c(1,142)]
```


```{r}
# Merging
merged <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE),
list(Baringo_total, Bungoma_total, Bomet_total, Busia_total, Baringo2_total, Coast_total, Elgon_total, Eldama_Ravine_total,Kakamega_total, ICIPE_total,Isiolo_total,Kasanga_Mwingi_total,Kasika_Mwingi_total,Kib_Kakamega_total,Kisii_total,Kitui_total,KNN_Nyamira,maralal_total,Maru_Nairobi_total,Meru_total,Narok_total,phi_total,saf_Nairobi_total,Taita_total,Thog_Nairobi_total, Uasin_Gishu_total,Ukasic_Mwingi_total,Vin_Kakamega_total,all_total))
```

```{r}
names(merged)<-c('Genus','Baringo','Bungoma','Bomet','Busia','Baringo2','Coast','Elgon','Eldama_Ravine','Kakamega','ICIPE','Isiolo','Kasanga_Mwingi','Kasika_Mwingi','Kib_Kakamega','Kisii','Kitui','KNN_Nyamira','maralal','Maru_Nairobi','Meru','Narok','phi','saf_Nairobi','Taita','Thog_Nairobi','Uasin_Gishu','Ukasic_Mwingi','Vin_Kakamega','All')
```

```{r}
#calculating the total abundance per genus and ordering from the most abundant to the lowest
cumulation <- merged %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100
```

```{r}
head(cumulation, n=20)
head(cumulation$Genus, n=51)
```

```{r}
genus_Rep <- c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus")
```


```{r}
group <- merged
```


```{r}
Baringo<-group[,c(1:2)]
Bungoma<- group[,c(1:3)] 
Bomet<- group[,c(1:4)]
Busia<- group[,c(1:5)]
Baringo2<- group[,c(1:6)]
Coast<- group[,c(1:7)]
Elgon<- group[,c(1:8)]
Eldama_Ravine<- group[,c(1:9)]
Kakamega<- group[,c(1:10)]
ICIPE<- group[,c(1:11)]
Isiolo<- group[,c(1:12)]
Kasanga_Mwingi<- group[,c(1:13)]
Kasika_Mwingi<- group[,c(1:14)]
Kib_Kakamega<- group[,c(1:15)]
Kisii<- group[,c(1:16)]
Kitui<- group[,c(1:17)]
KNN_Nyamira<- group[,c(1:18)]
maralal<- group[,c(1:19)]
Maru_Nairobi<- group[,c(1:20)]
Meru<- group[,c(1:21)]
Narok<- group[,c(1:22)]
phi<- group[,c(1:23)]
saf_Nairobi<- group[,c(1:24)]
Taita<- group[,c(1:25)]
Thog_Nairobi<- group[,c(1:26)]
Uasin_Gishu<- group[,c(1:27)]
Ukasic_Mwingi<- group[,c(1:28)]
Vin_Kakamega<- group[,c(1:29)]
All<- group[,c(1:30)]
```


```{r}
# Viewing Sample Diversity
#install.packages("janitor")
library(janitor)
```


```{r}
#converting the abudances into percentage
bar_all <- adorn_percentages(All, denominator = "col", na.rm = T)
bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all<-bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all
write.csv(dist_all, "/home/icipe/Downloads/dada2/Abundance.csv")
```
```{r}
#gathering the data
bar_all <- bar_all %>%
  gather(value = "abundance", key = "Site", -Genus)
bar_all <- as.data.frame(gsub("\\(", " (", as.matrix(bar_all)))
```

```{r}
# coerce the dataframe columns into respective data type
bar_all$Genus <- as.factor(bar_all$Genus)
bar_all$Site <- as.character(bar_all$Site)
bar_all$abundance <- as.numeric(bar_all$abundance)
```


```{r}
#ordering the data for plotting
bar_all$Genus <- reorder(bar_all$Genus, bar_all$abundance)
bar_all$Genus <- factor(bar_all$Genus, levels=rev(levels(bar_all$Genus)))
bar_all$Genus <- factor(bar_all$Genus, 
                        levels=c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus"))
```

```{r}
# Defining the color pallete
myPalette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#99000D", "#E6AB02", "#A6761D", "#666666","#FDCDAC", "#1F78B4", "#B2DF8A", "#33A02C", "#CBD5E8", "#E31A1C", "#FDBF6F", "#FF7F00","#4A1486","#C0C0C0","#B3E2CD","#FFFF33", "#5172b2","#F4CAE4", "#E6F5C9", "#FCFBFD","#139BF1","#09FF00","#065535", "#1D91C0", "#C0FFEE","#B35806","#0C2C84","#D0ED0E","#092617","#8A17FD","#4D5D53","#E48400","#6082B6","#316689","#EEA47FFF","#738678","#645452","#EEA47FFF", "#00539CFF", "#FC766AFF", "#42EADDFF", "#00A4CCFF", "#69B3BB", "#B589D6","#D1DFB7","#97BC62FF","#D198C5FF","#000000", "#CBCE91FF", "#616247FF", "#D64161FF","#435E55FF", "#DD4132FF","#CE4A7EFF", "#BD7F37FF","#FFA351FF","#185E57")
```


```{r}
# Definig the names in italics
guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 10, 
face = "italic", colour = "Black", angle = 0)))
```

# Plotting barplot

```{r}
#library(Cairo)
library(forcats)

#plotting the barplot 
p_all <- ggplot(bar_all,aes(x = fct_inorder(Site), y = abundance), labs(fill= Genus), group=row.names(bar_all))+ xlab("Site")+ ylab("abundance") + geom_col(aes(fill = Genus),position = position_stack(reverse = FALSE))+
  theme(axis.text.x = element_text(angle = 72,size = 15, hjust = 1, face = "italic", family = "Arial"))+
   scale_fill_manual(values = myPalette)+
  #guides(fill = guide_legend(reverse = FALSE))+
  guide_italics+
  theme(legend.text = element_text(size = 8, colour = "black", face = "italic", family = "Arial"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 8, family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 8, family = "Arial"))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.15, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.position = "right", legend.justification = "top", legend.direction = "vertical", legend.text = element_text(size = 10))+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 8, family = "Arial"))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 8, family = "Arial"))
```


```{r}
#install.packages('extrafont')
library(extrafont)
#font_import
#BiocManager::install("ggpubr")
```


```{r}
#change the misspelled names
ps@sam_data[ps@sam_data == 'liotrigona'] <- 'Liotrigona'
```

```{r}
#check if it worked
ps@sam_data@.Data[[2]]
```

```{r}
# check the number of taxa
ntaxa(ps)
```

```{r}
# Per county
grouped_data


Baringo<-grouped_data[,c(1,2:4,20:23)]
Bungoma<-grouped_data[,c(1,5:10)]
Bomet<-grouped_data[,c(1,11:16)]
Busia<-grouped_data[,c(1,17:19)]
Coast<-grouped_data[,c(1,24:26)]
Elgon<-grouped_data[,c(1,27:31)]
Eldama_Ravine<-grouped_data[,c(1,32:34)]
Kakamega<-grouped_data[,c(1,35:46,73:74,136:139)]
Nairobi<-grouped_data[,c(1,47:54,99:100,110:116,124:126)]
Isiolo<-grouped_data[,c(1,55:62)]
Mwingi<-grouped_data[,c(1,63:67,68:72,133:135)]
Kisii<-grouped_data[,c(1,75:88)]
Kitui<-grouped_data[,c(1,89:93)]
Nyamira<-grouped_data[,c(1,94)]
Maralal<-grouped_data[,c(1,95:98)]
Meru<-grouped_data[,c(1,101:106)]
Narok<-grouped_data[,c(1,107:109)]
Taita<-grouped_data[,c(1,117:123,140)]
Uasin_Gishu<-grouped_data[,c(1,127:132)]
Voic_Taita<-grouped_data[,c(1,140)]
all<-grouped_data[,c(1,2:140)]

Baringo_total <-Baringo%>% adorn_totals(c("col"))
Baringo_total <- mutate(Baringo_total, Baringo=rowSums(Baringo_total[9])/8)
Baringo_total <- Baringo_total[,c(1,10)]

Bungoma_total <-Bungoma%>% adorn_totals(c("col"))
Bungoma_total <- mutate(Bungoma_total, Bungoma=rowSums(Bungoma_total[8])/7)
Bungoma_total <- Bungoma_total[,c(1,9)]

Bomet_total <-Bomet%>% adorn_totals(c("col"))
Bomet_total <- mutate(Bomet_total, Bomet=rowSums(Bomet_total[8])/7)
Bomet_total <- Bomet_total[,c(1,9)]

Busia_total<-Busia%>% adorn_totals(c("col"))
Busia_total <- mutate(Busia_total, Busia=rowSums(Busia_total[5])/4)
Busia_total <- Busia_total[,c(1,6)]


Coast_total<-Coast%>% adorn_totals(c("col"))
Coast_total <- mutate(Coast_total, Coast=rowSums(Coast_total[5])/4)
Coast_total <- Coast_total[,c(1,6)]

Elgon_total<-Elgon%>% adorn_totals(c("col"))
Elgon_total <- mutate(Elgon_total, Elgon=rowSums(Elgon_total[7])/6)
Elgon_total <- Elgon_total[,c(1,8)]

Eldama_Ravine_total<-Eldama_Ravine%>% adorn_totals(c("col"))
Eldama_Ravine_total <- mutate(Eldama_Ravine_total, Eldama_Ravine=rowSums(Eldama_Ravine_total[5])/4)
Eldama_Ravine_total <- Eldama_Ravine_total[,c(1,6)]

Kakamega_total<-Kakamega%>% adorn_totals(c("col"))
Kakamega_total <- mutate(Kakamega_total, Kakamega=rowSums(Kakamega_total[19])/20)
Kakamega_total <- Kakamega_total[,c(1,21)]

Nairobi_total <-Nairobi%>% adorn_totals(c("col"))
Nairobi_total <- mutate(Nairobi_total, Nairobi=rowSums(Nairobi_total[22])/21)
Nairobi_total <- Nairobi_total[,c(1,23)]

Isiolo_total<-Isiolo%>% adorn_totals(c("col"))
Isiolo_total <- mutate(Isiolo_total, Isiolo=rowSums(Isiolo_total[10])/9)
Isiolo_total <- Isiolo_total[,c(1,11)]

Mwingi_total <- Mwingi%>%adorn_totals(c("col"))
Mwingi_total<- mutate(Mwingi_total, Mwingi=rowSums(Mwingi_total[15])/14)
Mwingi_total <- Mwingi_total[,c(1,16)]


Kisii_total<-Kisii%>% adorn_totals(c("col"))
Kisii_total <- mutate(Kisii_total, Kisii=rowSums(Kisii_total[16])/15)
Kisii_total <- Kisii_total[,c(1,17)]

Kitui_total<-Kitui%>% adorn_totals(c("col"))
Kitui_total <- mutate(Kitui_total, Kitui=rowSums(Kitui_total[7])/6)
Kitui_total <- Kitui_total[,c(1,8)]

Nyamira_total<-Nyamira%>% adorn_totals(c("col"))
Nyamira_total <- mutate(Nyamira_total, Nyamira=rowSums(Nyamira_total[3])/2)
Nyamira_total <- Nyamira_total[,c(1,4)]

Maralal_total<-Maralal%>% adorn_totals(c("col"))
Maralal_total <- mutate(Maralal_total, Maralal=rowSums(Maralal_total[6])/5)
Maralal_total <- Maralal_total[,c(1,7)]

Meru_total<-Meru%>% adorn_totals(c("col"))
Meru_total <- mutate(Meru_total, Meru=rowSums(Meru_total[8])/7)
Meru_total <- Meru_total[,c(1,9)]

Narok_total<-Narok%>% adorn_totals(c("col"))
Narok_total <- mutate(Narok_total, Narok=rowSums(Narok_total[5])/4)
Narok_total <- Narok_total[,c(1,6)]

Taita_total<-Taita%>% adorn_totals(c("col"))
Taita_total <- mutate(Taita_total, Taita=rowSums(Taita_total[10])/9)
Taita_total <- Taita_total[,c(1,11)]

Uasin_Gishu_total<-Uasin_Gishu%>% adorn_totals(c("col"))
Uasin_Gishu_total <- mutate(Uasin_Gishu_total, Uasin_Gishu=rowSums(Uasin_Gishu_total[8])/7)
Uasin_Gishu_total <- Uasin_Gishu_total[,c(1,9)]

all_total <- all %>% adorn_totals(c("col"))
all_total <- mutate(all_total, all=rowSums(all_total[141])/140)
all_total <- all_total[,c(1,142)]

merged1 <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE),
                 list(Baringo_total, Bungoma_total, Bomet_total, Busia_total,Coast_total, Elgon_total, Eldama_Ravine_total,Kakamega_total, Nairobi_total,Isiolo_total,Mwingi_total,Kisii_total,Kitui_total,Nyamira_total,Maralal_total,Meru_total,Narok_total,Taita_total,Uasin_Gishu_total,all_total))

cumulation <- merged1 %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100

genus_Rep <-c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus")

group <- aggregate(merged1[-1], list(Genus = replace(merged1$Genus,!(merged1$Genus %in% genus_Rep), "Others")), sum)

Baringo<-group[,c(1:2)]
Bungoma<- group[,c(1:3)] 
Bomet<- group[,c(1:4)]
Busia<- group[,c(1:5)]
Coast<- group[,c(1:6)]
Elgon<- group[,c(1:7)]
Eldama_Ravine<- group[,c(1:8)]
Kakamega<- group[,c(1:9)]
Nairobi<- group[,c(1:10)]
Isiolo<- group[,c(1:11)]
Mwingi<- group[,c(1:12)]
Kisii<- group[,c(1:13)]
Kitui<- group[,c(1:14)]
Nyamira<- group[,c(1:15)]
Maralal<- group[,c(1:16)]
Meru<- group[,c(1:17)]
Narok<- group[,c(1:18)]
Taita<- group[,c(1:19)]
Uasin_Gishu<- group[,c(1:20)]
All<- group[,c(1:20)]
# metagenomic analysis using dada2

# loading of packages

```{r}
library("phyloseq")
library("vegan")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library("janitor")
library("coin")
library("reshape2")
library("ggnewscale")
library("MicrobiotaProcess")
library("patchwork")
library("dada2")
library("csv")
library("dplyr")
library("ggtree")
library("VennDiagram")
library("UpSetR")
#BiocManager::install("MicrobiotaProcess")
install.packages("MicrobiotaProcess")
library("MicrobiotaProcess")
library("utils")
library("utf8")
library("metagMisc")
```


```{r}
#BiocManager::install("devtools")
library("devtools")
library("tidyverse")
#devtools::install_github(repo= "jbisanz/qiime2R")
library(qiime2R)
asv_table <- read.table("/home/icipe/Downloads/dada2/ASV_table.tsv")
```

# loading the data

```{r}
metadata <- read.table("/home/icipe/Downloads/Geo_metadata.tsv", sep = "\t", row.names = 1, header= TRUE)

ot <- read.table("/home/icipe/Downloads/dada2/ASV_table.tsv", header = TRUE, sep = "\t",  row.names = 1 )
class(ot)
ot <- as.matrix(ot)
class(ot)

tax <- read.table("/home/icipe/Downloads/dada2/ASV_tax_species.tsv", header = TRUE, sep = "\t",  row.names = 1)
#tax <- data.frame(tax)
Tax <- as.matrix(tax)
#tax <- t(tax)

ot <- phyloseq::otu_table(ot, taxa_are_rows = TRUE)
tax <- phyloseq::tax_table(Tax)
metadata <- phyloseq::sample_data(metadata)

```
```{r}
ps <- phyloseq(ot, metadata, tax)
ps
```
```{r}
# we now create a third object called random for merging with the other three object
library("ape")
random_tree <- rtree(ntaxa(physeq), rooted = TRUE, tip.label = taxa_names(physeq))
plot(random_tree)
```

```{r}
#merging the preceeding 3 objects.
ps <- merge_phyloseq(physeq, random_tree)
ps
```

```{r}
# Removing the NAs
tx <- phyloseq_to_df(ps, addtax = T, addtot = F, addmaxrank = F)
tx <-tx[, c(8, 12:150 )]

tx <- tx %>%
  # recode empty strings "" by NAs
 na_if("") %>%
  # remove NAs
  na.omit()
tx
```

```{r}
cumulation <- tx %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation
```
```{r}
#merging the blast taxonomic classification to blast abundance table
merged_data <- tx
write.csv(merged_data, file="/home/icipe/Downloads/dada2/merged_data.csv")
# Group the data
grouped_data <- tx %>% group_by(Genus) %>%summarise_if(is.numeric, sum)
#grouped_data<-grouped_data[-c(1),]
```


```{r}
saveRDS(grouped_data,"/home/icipe/Downloads/dada2/grouped_data.rds")
write.table(tx,"/home/icipe/Downloads/dada2/tx.csv")
```


```{r}
# Grouping data by location
Baringo<-grouped_data[,c(1,2:4)]
Bungoma<-grouped_data[,c(1,5:10)]
Bomet<-grouped_data[,c(1,11:16)]
Busia<-grouped_data[,c(1,17:19)]
Baringo2<-grouped_data[,c(1,20:23)]
Coast<-grouped_data[,c(1,24:26)]
Elgon<-grouped_data[,c(1,27:31)]
Eldama_Ravine<-grouped_data[,c(1,32:34)]
Kakamega<-grouped_data[,c(1,35:46)]
ICIPE<-grouped_data[,c(1,47:54)]
Isiolo<-grouped_data[,c(1,55:62)]
Kasanga_Mwingi<-grouped_data[,c(1,63:67)]
Kasika_Mwingi<-grouped_data[,c(1,68:72)]
Kib_Kakamega<-grouped_data[,c(1,73:74)]
Kisii<-grouped_data[,c(1,75:88)]
Kitui<-grouped_data[,c(1,89:93)]
KNN_Nyamira<-grouped_data[,c(1,94)]
maralal<-grouped_data[,c(1,95:98)]
Maru_Nairobi<-grouped_data[,c(1,99:100)]
Meru<-grouped_data[,c(1,101:106)]
Narok<-grouped_data[,c(1,107:109)]
phi<-grouped_data[,c(1,110)]
saf_Nairobi<-grouped_data[,c(1,111:116)]
Taita<-grouped_data[,c(1,117:123)]
Thog_Nairobi<-grouped_data[,c(1,124:126)]
Uasin_Gishu<-grouped_data[,c(1,127:132)]
Ukasic_Mwingi<-grouped_data[,c(1,133:135)]
Vin_Kakamega<-grouped_data[,c(1,136:139)]
Voic_Taita <- grouped_data[1,140]
```

```{r}
all<-grouped_data[,c(1,2:140)]
all
```

```{r}
Baringo_total <-Baringo%>% adorn_totals(c("col"))
Baringo_total <- mutate(Baringo_total, Baringo=rowSums(Baringo_total[5])/4)
Baringo_total <- Baringo_total[,c(1,6)]

Bungoma_total <-Bungoma%>% adorn_totals(c("col"))
Bungoma_total <- mutate(Bungoma_total, Bungoma=rowSums(Bungoma_total[8])/7)
Bungoma_total <- Bungoma_total[,c(1,9)]

Bomet_total <-Bomet%>% adorn_totals(c("col"))
Bomet_total <- mutate(Bomet_total, Bomet=rowSums(Bomet_total[8])/7)
Bomet_total <- Bomet_total[,c(1,9)]

Busia_total<-Busia%>% adorn_totals(c("col"))
Busia_total <- mutate(Busia_total, Busia=rowSums(Busia_total[5])/4)
Busia_total <- Busia_total[,c(1,6)]

Baringo2_total<-Baringo2%>% adorn_totals(c("col"))
Baringo2_total <- mutate(Baringo2_total, Baringo2=rowSums(Baringo2_total[6])/5)
Baringo2_total <- Baringo2_total[,c(1,7)]

Coast_total<-Coast%>% adorn_totals(c("col"))
Coast_total <- mutate(Coast_total, Coast=rowSums(Coast_total[5])/4)
Coast_total <- Coast_total[,c(1,6)]

Elgon_total<-Elgon%>% adorn_totals(c("col"))
Elgon_total <- mutate(Elgon_total, Elgon=rowSums(Elgon_total[7])/6)
Elgon_total <- Elgon_total[,c(1,8)]

Eldama_Ravine_total<-Eldama_Ravine%>% adorn_totals(c("col"))
Eldama_Ravine_total <- mutate(Eldama_Ravine_total, Eldama_Ravine=rowSums(Eldama_Ravine_total[5])/4)
Eldama_Ravine_total <- Eldama_Ravine_total[,c(1,6)]

Kakamega_total<-Kakamega%>% adorn_totals(c("col"))
Kakamega_total <- mutate(Kakamega_total, Kakamega=rowSums(Kakamega_total[14])/13)
Kakamega_total <- Kakamega_total[,c(1,15)]

ICIPE_total <-ICIPE%>% adorn_totals(c("col"))
ICIPE_total <- mutate(ICIPE_total, ICIPE=rowSums(ICIPE_total[10])/9)
ICIPE_total <- ICIPE_total[,c(1,11)]

Isiolo_total<-Isiolo%>% adorn_totals(c("col"))
Isiolo_total <- mutate(Isiolo_total, Isiolo=rowSums(Isiolo_total[10])/9)
Isiolo_total <- Isiolo_total[,c(1,11)]

Kasanga_Mwingi_total <- Kasanga_Mwingi%>%adorn_totals(c("col"))
Kasanga_Mwingi_total<- mutate(Kasanga_Mwingi_total, Kasanga_Mwingi=rowSums(Kasanga_Mwingi_total[7])/6)
Kasanga_Mwingi_total <- Kasanga_Mwingi_total[,c(1,8)]

Kasika_Mwingi_total<-Kasika_Mwingi%>% adorn_totals(c("col"))
Kasika_Mwingi_total <- mutate(Kasika_Mwingi_total, Kasika_Mwingi=rowSums(Kasika_Mwingi_total[7])/6)
Kasika_Mwingi_total <- Kasika_Mwingi_total[,c(1,8)]

Kib_Kakamega_total<-Kib_Kakamega%>% adorn_totals(c("col"))
Kib_Kakamega_total <- mutate(Kib_Kakamega_total, Kib_Kakamega=rowSums(Kib_Kakamega_total[4])/3)
Kib_Kakamega_total <- Kib_Kakamega_total[,c(1,5)]

Kisii_total<-Kisii%>% adorn_totals(c("col"))
Kisii_total <- mutate(Kisii_total, Kisii=rowSums(Kisii_total[16])/15)
Kisii_total <- Kisii_total[,c(1,17)]

Kitui_total<-Kitui%>% adorn_totals(c("col"))
Kitui_total <- mutate(Kitui_total, Kitui=rowSums(Kitui_total[7])/6)
Kitui_total <- Kitui_total[,c(1,8)]

KNN_Nyamira_total<-KNN_Nyamira%>% adorn_totals(c("col"))
KNN_Nyamira_total <- mutate(KNN_Nyamira_total, KNN_Nyamira=rowSums(KNN_Nyamira_total[3])/2)
KNN_Nyamira_total <- KNN_Nyamira_total[,c(1,4)]

maralal_total<-maralal%>% adorn_totals(c("col"))
maralal_total <- mutate(maralal_total, maralal=rowSums(maralal_total[6])/5)
maralal_total <- maralal_total[,c(1,7)]

Maru_Nairobi_total<-Maru_Nairobi%>% adorn_totals(c("col"))
Maru_Nairobi_total <- mutate(Maru_Nairobi_total, Maru_Nairobi=rowSums(Maru_Nairobi_total[4])/3)
Maru_Nairobi_total <- Maru_Nairobi_total[,c(1,5)]

Meru_total<-Meru%>% adorn_totals(c("col"))
Meru_total <- mutate(Meru_total, Meru=rowSums(Meru_total[8])/7)
Meru_total <- Meru_total[,c(1,9)]

Narok_total<-Narok%>% adorn_totals(c("col"))
Narok_total <- mutate(Narok_total, Narok=rowSums(Narok_total[5])/4)
Narok_total <- Narok_total[,c(1,6)]

phi_total<-phi%>% adorn_totals(c("col"))
phi_total <- mutate(phi_total, phi=rowSums(phi_total[3])/2)
phi_total <- phi_total[,c(1,4)]

saf_Nairobi_total<-saf_Nairobi%>% adorn_totals(c("col"))
saf_Nairobi_total <- mutate(saf_Nairobi_total, saf_Nairobi=rowSums(saf_Nairobi_total[8])/7)
saf_Nairobi_total <- saf_Nairobi_total[,c(1,9)]

Taita_total<-Taita%>% adorn_totals(c("col"))
Taita_total <- mutate(Taita_total, Taita=rowSums(Taita_total[9])/8)
Taita_total <- Taita_total[,c(1,10)]

Thog_Nairobi_total<-Thog_Nairobi%>% adorn_totals(c("col"))
Thog_Nairobi_total <- mutate(Thog_Nairobi_total, Thog_Nairobi=rowSums(Thog_Nairobi_total[5])/4)
Thog_Nairobi_total <- Thog_Nairobi_total[,c(1,6)]

Uasin_Gishu_total<-Uasin_Gishu%>% adorn_totals(c("col"))
Uasin_Gishu_total <- mutate(Uasin_Gishu_total, Uasin_Gishu=rowSums(Uasin_Gishu_total[8])/7)
Uasin_Gishu_total <- Uasin_Gishu_total[,c(1,9)]

Ukasic_Mwingi_total<-Ukasic_Mwingi%>% adorn_totals(c("col"))
Ukasic_Mwingi_total <- mutate(Ukasic_Mwingi_total, Ukasic_Mwingi=rowSums(Ukasic_Mwingi_total[5])/4)
Ukasic_Mwingi_total <- Ukasic_Mwingi_total[,c(1,6)]

Vin_Kakamega_total<-Vin_Kakamega%>% adorn_totals(c("col"))
Vin_Kakamega_total <- mutate(Vin_Kakamega_total, Vin_Kakamega=rowSums(Vin_Kakamega_total[6])/5)
Vin_Kakamega_total <- Vin_Kakamega_total[,c(1,7)]

Voic_Taita_total <-Voic_Taita%>% adorn_totals(c("col"))
Voic_Taita_total <- mutate(Voic_Taita_total, Voic_Taita=rowSums(Voic_Taita_total[3])/2)
Voic_Taita_total <- Voic_Taita_total[,c(1,4)]
```

```{r}
all_total <- all %>% adorn_totals(c("col"))
all_total <- mutate(all_total, all=rowSums(all_total[141])/140)
all_total <- all_total[,c(1,142)]
```


```{r}
# Merging
merged <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE),
list(Baringo_total, Bungoma_total, Bomet_total, Busia_total, Baringo2_total, Coast_total, Elgon_total, Eldama_Ravine_total,Kakamega_total, ICIPE_total,Isiolo_total,Kasanga_Mwingi_total,Kasika_Mwingi_total,Kib_Kakamega_total,Kisii_total,Kitui_total,KNN_Nyamira,maralal_total,Maru_Nairobi_total,Meru_total,Narok_total,phi_total,saf_Nairobi_total,Taita_total,Thog_Nairobi_total, Uasin_Gishu_total,Ukasic_Mwingi_total,Vin_Kakamega_total,all_total))
```

```{r}
names(merged)<-c('Genus','Baringo','Bungoma','Bomet','Busia','Baringo2','Coast','Elgon','Eldama_Ravine','Kakamega','ICIPE','Isiolo','Kasanga_Mwingi','Kasika_Mwingi','Kib_Kakamega','Kisii','Kitui','KNN_Nyamira','maralal','Maru_Nairobi','Meru','Narok','phi','saf_Nairobi','Taita','Thog_Nairobi','Uasin_Gishu','Ukasic_Mwingi','Vin_Kakamega','All')
```

```{r}
#calculating the total abundance per genus and ordering from the most abundant to the lowest
cumulation <- merged %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100
```

```{r}
head(cumulation, n=20)
head(cumulation$Genus, n=51)
```

```{r}
genus_Rep <- c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus")
```


```{r}
group <- merged
```


```{r}
Baringo<-group[,c(1:2)]
Bungoma<- group[,c(1:3)] 
Bomet<- group[,c(1:4)]
Busia<- group[,c(1:5)]
Baringo2<- group[,c(1:6)]
Coast<- group[,c(1:7)]
Elgon<- group[,c(1:8)]
Eldama_Ravine<- group[,c(1:9)]
Kakamega<- group[,c(1:10)]
ICIPE<- group[,c(1:11)]
Isiolo<- group[,c(1:12)]
Kasanga_Mwingi<- group[,c(1:13)]
Kasika_Mwingi<- group[,c(1:14)]
Kib_Kakamega<- group[,c(1:15)]
Kisii<- group[,c(1:16)]
Kitui<- group[,c(1:17)]
KNN_Nyamira<- group[,c(1:18)]
maralal<- group[,c(1:19)]
Maru_Nairobi<- group[,c(1:20)]
Meru<- group[,c(1:21)]
Narok<- group[,c(1:22)]
phi<- group[,c(1:23)]
saf_Nairobi<- group[,c(1:24)]
Taita<- group[,c(1:25)]
Thog_Nairobi<- group[,c(1:26)]
Uasin_Gishu<- group[,c(1:27)]
Ukasic_Mwingi<- group[,c(1:28)]
Vin_Kakamega<- group[,c(1:29)]
All<- group[,c(1:30)]
```


```{r}
# Viewing Sample Diversity
#install.packages("janitor")
library(janitor)
```


```{r}
#converting the abudances into percentage
bar_all <- adorn_percentages(All, denominator = "col", na.rm = T)
bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all<-bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all
write.csv(dist_all, "/home/icipe/Downloads/dada2/Abundance.csv")
```
```{r}
#gathering the data
bar_all <- bar_all %>%
  gather(value = "abundance", key = "Site", -Genus)
bar_all <- as.data.frame(gsub("\\(", " (", as.matrix(bar_all)))
```

```{r}
# coerce the dataframe columns into respective data type
bar_all$Genus <- as.factor(bar_all$Genus)
bar_all$Site <- as.character(bar_all$Site)
bar_all$abundance <- as.numeric(bar_all$abundance)
```


```{r}
#ordering the data for plotting
bar_all$Genus <- reorder(bar_all$Genus, bar_all$abundance)
bar_all$Genus <- factor(bar_all$Genus, levels=rev(levels(bar_all$Genus)))
bar_all$Genus <- factor(bar_all$Genus, 
                        levels=c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus"))
```

```{r}
# Defining the color pallete
myPalette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#99000D", "#E6AB02", "#A6761D", "#666666","#FDCDAC", "#1F78B4", "#B2DF8A", "#33A02C", "#CBD5E8", "#E31A1C", "#FDBF6F", "#FF7F00","#4A1486","#C0C0C0","#B3E2CD","#FFFF33", "#5172b2","#F4CAE4", "#E6F5C9", "#FCFBFD","#139BF1","#09FF00","#065535", "#1D91C0", "#C0FFEE","#B35806","#0C2C84","#D0ED0E","#092617","#8A17FD","#4D5D53","#E48400","#6082B6","#316689","#EEA47FFF","#738678","#645452","#EEA47FFF", "#00539CFF", "#FC766AFF", "#42EADDFF", "#00A4CCFF", "#69B3BB", "#B589D6","#D1DFB7","#97BC62FF","#D198C5FF","#000000", "#CBCE91FF", "#616247FF", "#D64161FF","#435E55FF", "#DD4132FF","#CE4A7EFF", "#BD7F37FF","#FFA351FF","#185E57")
```


```{r}
# Definig the names in italics
guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 10, 
face = "italic", colour = "Black", angle = 0)))
```

# Plotting barplot

```{r}
#library(Cairo)
library(forcats)

#plotting the barplot 
p_all <- ggplot(bar_all,aes(x = fct_inorder(Site), y = abundance), labs(fill= Genus), group=row.names(bar_all))+ xlab("Site")+ ylab("abundance") + geom_col(aes(fill = Genus),position = position_stack(reverse = FALSE))+
  theme(axis.text.x = element_text(angle = 72,size = 15, hjust = 1, face = "italic", family = "Arial"))+
   scale_fill_manual(values = myPalette)+
  #guides(fill = guide_legend(reverse = FALSE))+
  guide_italics+
  theme(legend.text = element_text(size = 8, colour = "black", face = "italic", family = "Arial"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 8, family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 8, family = "Arial"))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.15, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.position = "right", legend.justification = "top", legend.direction = "vertical", legend.text = element_text(size = 10))+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 8, family = "Arial"))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 8, family = "Arial"))
```


```{r}
#install.packages('extrafont')
library(extrafont)
#font_import
#BiocManager::install("ggpubr")
```


```{r}
#change the misspelled names
ps@sam_data[ps@sam_data == 'liotrigona'] <- 'Liotrigona'
```

```{r}
#check if it worked
ps@sam_data@.Data[[2]]
```

```{r}
# check the number of taxa
ntaxa(ps)
```

```{r}
# Per county
grouped_data


Baringo<-grouped_data[,c(1,2:4,20:23)]
Bungoma<-grouped_data[,c(1,5:10)]
Bomet<-grouped_data[,c(1,11:16)]
Busia<-grouped_data[,c(1,17:19)]
Coast<-grouped_data[,c(1,24:26)]
Elgon<-grouped_data[,c(1,27:31)]
Eldama_Ravine<-grouped_data[,c(1,32:34)]
Kakamega<-grouped_data[,c(1,35:46,73:74,136:139)]
Nairobi<-grouped_data[,c(1,47:54,99:100,110:116,124:126)]
Isiolo<-grouped_data[,c(1,55:62)]
Mwingi<-grouped_data[,c(1,63:67,68:72,133:135)]
Kisii<-grouped_data[,c(1,75:88)]
Kitui<-grouped_data[,c(1,89:93)]
Nyamira<-grouped_data[,c(1,94)]
Maralal<-grouped_data[,c(1,95:98)]
Meru<-grouped_data[,c(1,101:106)]
Narok<-grouped_data[,c(1,107:109)]
Taita<-grouped_data[,c(1,117:123,140)]
Uasin_Gishu<-grouped_data[,c(1,127:132)]
Voic_Taita<-grouped_data[,c(1,140)]
all<-grouped_data[,c(1,2:140)]

Baringo_total <-Baringo%>% adorn_totals(c("col"))
Baringo_total <- mutate(Baringo_total, Baringo=rowSums(Baringo_total[9])/8)
Baringo_total <- Baringo_total[,c(1,10)]

Bungoma_total <-Bungoma%>% adorn_totals(c("col"))
Bungoma_total <- mutate(Bungoma_total, Bungoma=rowSums(Bungoma_total[8])/7)
Bungoma_total <- Bungoma_total[,c(1,9)]

Bomet_total <-Bomet%>% adorn_totals(c("col"))
Bomet_total <- mutate(Bomet_total, Bomet=rowSums(Bomet_total[8])/7)
Bomet_total <- Bomet_total[,c(1,9)]

Busia_total<-Busia%>% adorn_totals(c("col"))
Busia_total <- mutate(Busia_total, Busia=rowSums(Busia_total[5])/4)
Busia_total <- Busia_total[,c(1,6)]


Coast_total<-Coast%>% adorn_totals(c("col"))
Coast_total <- mutate(Coast_total, Coast=rowSums(Coast_total[5])/4)
Coast_total <- Coast_total[,c(1,6)]

Elgon_total<-Elgon%>% adorn_totals(c("col"))
Elgon_total <- mutate(Elgon_total, Elgon=rowSums(Elgon_total[7])/6)
Elgon_total <- Elgon_total[,c(1,8)]

Eldama_Ravine_total<-Eldama_Ravine%>% adorn_totals(c("col"))
Eldama_Ravine_total <- mutate(Eldama_Ravine_total, Eldama_Ravine=rowSums(Eldama_Ravine_total[5])/4)
Eldama_Ravine_total <- Eldama_Ravine_total[,c(1,6)]

Kakamega_total<-Kakamega%>% adorn_totals(c("col"))
Kakamega_total <- mutate(Kakamega_total, Kakamega=rowSums(Kakamega_total[19])/20)
Kakamega_total <- Kakamega_total[,c(1,21)]

Nairobi_total <-Nairobi%>% adorn_totals(c("col"))
Nairobi_total <- mutate(Nairobi_total, Nairobi=rowSums(Nairobi_total[22])/21)
Nairobi_total <- Nairobi_total[,c(1,23)]

Isiolo_total<-Isiolo%>% adorn_totals(c("col"))
Isiolo_total <- mutate(Isiolo_total, Isiolo=rowSums(Isiolo_total[10])/9)
Isiolo_total <- Isiolo_total[,c(1,11)]

Mwingi_total <- Mwingi%>%adorn_totals(c("col"))
Mwingi_total<- mutate(Mwingi_total, Mwingi=rowSums(Mwingi_total[15])/14)
Mwingi_total <- Mwingi_total[,c(1,16)]


Kisii_total<-Kisii%>% adorn_totals(c("col"))
Kisii_total <- mutate(Kisii_total, Kisii=rowSums(Kisii_total[16])/15)
Kisii_total <- Kisii_total[,c(1,17)]

Kitui_total<-Kitui%>% adorn_totals(c("col"))
Kitui_total <- mutate(Kitui_total, Kitui=rowSums(Kitui_total[7])/6)
Kitui_total <- Kitui_total[,c(1,8)]

Nyamira_total<-Nyamira%>% adorn_totals(c("col"))
Nyamira_total <- mutate(Nyamira_total, Nyamira=rowSums(Nyamira_total[3])/2)
Nyamira_total <- Nyamira_total[,c(1,4)]

Maralal_total<-Maralal%>% adorn_totals(c("col"))
Maralal_total <- mutate(Maralal_total, Maralal=rowSums(Maralal_total[6])/5)
Maralal_total <- Maralal_total[,c(1,7)]

Meru_total<-Meru%>% adorn_totals(c("col"))
Meru_total <- mutate(Meru_total, Meru=rowSums(Meru_total[8])/7)
Meru_total <- Meru_total[,c(1,9)]

Narok_total<-Narok%>% adorn_totals(c("col"))
Narok_total <- mutate(Narok_total, Narok=rowSums(Narok_total[5])/4)
Narok_total <- Narok_total[,c(1,6)]

Taita_total<-Taita%>% adorn_totals(c("col"))
Taita_total <- mutate(Taita_total, Taita=rowSums(Taita_total[10])/9)
Taita_total <- Taita_total[,c(1,11)]

Uasin_Gishu_total<-Uasin_Gishu%>% adorn_totals(c("col"))
Uasin_Gishu_total <- mutate(Uasin_Gishu_total, Uasin_Gishu=rowSums(Uasin_Gishu_total[8])/7)
Uasin_Gishu_total <- Uasin_Gishu_total[,c(1,9)]

all_total <- all %>% adorn_totals(c("col"))
all_total <- mutate(all_total, all=rowSums(all_total[141])/140)
all_total <- all_total[,c(1,142)]

merged1 <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE),
                 list(Baringo_total, Bungoma_total, Bomet_total, Busia_total,Coast_total, Elgon_total, Eldama_Ravine_total,Kakamega_total, Nairobi_total,Isiolo_total,Mwingi_total,Kisii_total,Kitui_total,Nyamira_total,Maralal_total,Meru_total,Narok_total,Taita_total,Uasin_Gishu_total,all_total))

cumulation <- merged1 %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100

genus_Rep <-c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus")

group <- aggregate(merged1[-1], list(Genus = replace(merged1$Genus,!(merged1$Genus %in% genus_Rep), "Others")), sum)

Baringo<-group[,c(1:2)]
Bungoma<- group[,c(1:3)] 
Bomet<- group[,c(1:4)]
Busia<- group[,c(1:5)]
Coast<- group[,c(1:6)]
Elgon<- group[,c(1:7)]
Eldama_Ravine<- group[,c(1:8)]
Kakamega<- group[,c(1:9)]
Nairobi<- group[,c(1:10)]
Isiolo<- group[,c(1:11)]
Mwingi<- group[,c(1:12)]
Kisii<- group[,c(1:13)]
Kitui<- group[,c(1:14)]
Nyamira<- group[,c(1:15)]
Maralal<- group[,c(1:16)]
Meru<- group[,c(1:17)]
Narok<- group[,c(1:18)]
Taita<- group[,c(1:19)]
Uasin_Gishu<- group[,c(1:20)]
All<- group[,c(1:20)]

bar_all <- adorn_percentages(All, denominator = "col", na.rm = T)
bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all<-bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all

bar_all <- bar_all %>%
  gather(value = "abundance", key = "Site", -Genus)
bar_all <- as.data.frame(gsub("\\(", " (", as.matrix(bar_all)))

bar_all$Genus <- as.factor(bar_all$Genus)
bar_all$Site <- as.character(bar_all$Site)
bar_all$abundance <- as.numeric(bar_all$abundance)

bar_all$Genus <- reorder(bar_all$Genus, bar_all$abundance)
bar_all$Genus <- factor(bar_all$Genus, levels=rev(levels(bar_all$Genus)))
bar_all$Genus <- factor(bar_all$Genus, 
                        levels=genus_Rep)

myPalette <- c("#D95F02", "#7570B3", "#E7298A", "#99000D", "#E6AB02", "#A6761D", "#666666","#FDCDAC", "#1F78B4", "#B2DF8A", "#33A02C", "#CBD5E8", "#E31A1C", "#FDBF6F", "#FF7F00","#4A1486","#C0C0C0","#B3E2CD","#FFFF33", "#5172b2","#F4CAE4", "#E6F5C9", "#FCFBFD","#139BF1","#09FF00","#065535", "#1D91C0", "#C0FFEE","#B35806","#0C2C84","#D0ED0E","#092617","#499976","#4D5D53","#E48400","#6082B6","#316689","#CEFB02","#738678","#645452","#EEA47FFF", "#00539CFF", "#FC766AFF", "#42EADDFF", "#00A4CCFF", "#69B3BB", "#B589D6","#D1DFB7","#97BC62FF","#D198C5FF","#000000")
length(myPalette)

guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 15, 
                                                                       face = "italic", colour = "Black", angle = 0)))

p_all <- ggplot(bar_all,aes(x = fct_inorder(Site), y = abundance), labs(fill= Genus), group=row.names(bar_all))+ xlab("County")+ ylab("abundance") + geom_col(aes(fill = Genus),position = position_stack(reverse = FALSE))+
  theme(axis.text.x = element_text(angle = 72,size = 15, hjust = 1, face = "italic", family = "Arial"))+
  scale_fill_manual(values = myPalette)+
  #guides(fill = guide_legend(reverse = FALSE))+
  guide_italics+
  theme(legend.text = element_text(size = 8, colour = "black", face = "italic", family = "Arial"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 8, family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 25, family = "Arial"))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.15, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.position = "right", legend.justification = "top", legend.direction = "vertical", legend.text = element_text(size = 10))+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 15, family = "Arial"))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 15, family = "Arial"))

plot4_60<-p_all + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + xlab("County")
plot4_60

ggsave("/home/icipe/Downloads/dada2/plot4_60.jpeg", width=12, height=12, dpi=600)
```

```{r}
merged <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE),
list(Baringo_total, Bungoma_total, Bomet_total, Busia_total,Coast_total, Elgon_total, Eldama_Ravine_total,Kakamega_total, Nairobi_total,Isiolo_total,Mwingi_total,Kisii_total,Kitui_total,Nyamira_total,Maralal_total,Meru_total,Narok_total,Taita_total,Uasin_Gishu_total,all_total))
```

```{r}
#calculating the total abundance per genus and ordering from the most abundant to the lowest
cumulation <- merged %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100

head(cumulation$Genus, n=51)
```

```{r}
genus_Rep <- c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus")
```
# metagenomic analysis using dada2

# loading of packages

```{r}
library("phyloseq")
library("vegan")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library("janitor")
library("coin")
library("reshape2")
library("ggnewscale")
library("MicrobiotaProcess")
library("patchwork")
library("dada2")
library("csv")
library("dplyr")
library("ggtree")
library("VennDiagram")
library("UpSetR")
#BiocManager::install("MicrobiotaProcess")
install.packages("MicrobiotaProcess")
library("MicrobiotaProcess")
library("utils")
library("utf8")
library("metagMisc")
```


```{r}
#BiocManager::install("devtools")
library("devtools")
library("tidyverse")
#devtools::install_github(repo= "jbisanz/qiime2R")
library(qiime2R)
asv_table <- read.table("/home/icipe/Downloads/dada2/ASV_table.tsv")
```

# loading the data

```{r}
metadata <- read.table("/home/icipe/Downloads/Geo_metadata.tsv", sep = "\t", row.names = 1, header= TRUE)

ot <- read.table("/home/icipe/Downloads/dada2/ASV_table.tsv", header = TRUE, sep = "\t",  row.names = 1 )
class(ot)
ot <- as.matrix(ot)
class(ot)

tax <- read.table("/home/icipe/Downloads/dada2/ASV_tax_species.tsv", header = TRUE, sep = "\t",  row.names = 1)
#tax <- data.frame(tax)
Tax <- as.matrix(tax)
#tax <- t(tax)

ot <- phyloseq::otu_table(ot, taxa_are_rows = TRUE)
tax <- phyloseq::tax_table(Tax)
metadata <- phyloseq::sample_data(metadata)

```
```{r}
ps <- phyloseq(ot, metadata, tax)
ps
```
```{r}
# we now create a third object called random for merging with the other three object
library("ape")
random_tree <- rtree(ntaxa(physeq), rooted = TRUE, tip.label = taxa_names(physeq))
plot(random_tree)
```

```{r}
#merging the preceeding 3 objects.
ps <- merge_phyloseq(physeq, random_tree)
ps
```

```{r}
# Removing the NAs
tx <- phyloseq_to_df(ps, addtax = T, addtot = F, addmaxrank = F)
tx <-tx[, c(8, 12:150 )]

tx <- tx %>%
  # recode empty strings "" by NAs
 na_if("") %>%
  # remove NAs
  na.omit()
tx
```

```{r}
cumulation <- tx %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation
```
```{r}
#merging the blast taxonomic classification to blast abundance table
merged_data <- tx
write.csv(merged_data, file="/home/icipe/Downloads/dada2/merged_data.csv")
# Group the data
grouped_data <- tx %>% group_by(Genus) %>%summarise_if(is.numeric, sum)
#grouped_data<-grouped_data[-c(1),]
```


```{r}
saveRDS(grouped_data,"/home/icipe/Downloads/dada2/grouped_data.rds")
write.table(tx,"/home/icipe/Downloads/dada2/tx.csv")
```


```{r}
# Grouping data by location
Baringo<-grouped_data[,c(1,2:4)]
Bungoma<-grouped_data[,c(1,5:10)]
Bomet<-grouped_data[,c(1,11:16)]
Busia<-grouped_data[,c(1,17:19)]
Baringo2<-grouped_data[,c(1,20:23)]
Coast<-grouped_data[,c(1,24:26)]
Elgon<-grouped_data[,c(1,27:31)]
Eldama_Ravine<-grouped_data[,c(1,32:34)]
Kakamega<-grouped_data[,c(1,35:46)]
ICIPE<-grouped_data[,c(1,47:54)]
Isiolo<-grouped_data[,c(1,55:62)]
Kasanga_Mwingi<-grouped_data[,c(1,63:67)]
Kasika_Mwingi<-grouped_data[,c(1,68:72)]
Kib_Kakamega<-grouped_data[,c(1,73:74)]
Kisii<-grouped_data[,c(1,75:88)]
Kitui<-grouped_data[,c(1,89:93)]
KNN_Nyamira<-grouped_data[,c(1,94)]
maralal<-grouped_data[,c(1,95:98)]
Maru_Nairobi<-grouped_data[,c(1,99:100)]
Meru<-grouped_data[,c(1,101:106)]
Narok<-grouped_data[,c(1,107:109)]
phi<-grouped_data[,c(1,110)]
saf_Nairobi<-grouped_data[,c(1,111:116)]
Taita<-grouped_data[,c(1,117:123)]
Thog_Nairobi<-grouped_data[,c(1,124:126)]
Uasin_Gishu<-grouped_data[,c(1,127:132)]
Ukasic_Mwingi<-grouped_data[,c(1,133:135)]
Vin_Kakamega<-grouped_data[,c(1,136:139)]
Voic_Taita <- grouped_data[1,140]
```

```{r}
all<-grouped_data[,c(1,2:140)]
all
```

```{r}
Baringo_total <-Baringo%>% adorn_totals(c("col"))
Baringo_total <- mutate(Baringo_total, Baringo=rowSums(Baringo_total[5])/4)
Baringo_total <- Baringo_total[,c(1,6)]

Bungoma_total <-Bungoma%>% adorn_totals(c("col"))
Bungoma_total <- mutate(Bungoma_total, Bungoma=rowSums(Bungoma_total[8])/7)
Bungoma_total <- Bungoma_total[,c(1,9)]

Bomet_total <-Bomet%>% adorn_totals(c("col"))
Bomet_total <- mutate(Bomet_total, Bomet=rowSums(Bomet_total[8])/7)
Bomet_total <- Bomet_total[,c(1,9)]

Busia_total<-Busia%>% adorn_totals(c("col"))
Busia_total <- mutate(Busia_total, Busia=rowSums(Busia_total[5])/4)
Busia_total <- Busia_total[,c(1,6)]

Baringo2_total<-Baringo2%>% adorn_totals(c("col"))
Baringo2_total <- mutate(Baringo2_total, Baringo2=rowSums(Baringo2_total[6])/5)
Baringo2_total <- Baringo2_total[,c(1,7)]

Coast_total<-Coast%>% adorn_totals(c("col"))
Coast_total <- mutate(Coast_total, Coast=rowSums(Coast_total[5])/4)
Coast_total <- Coast_total[,c(1,6)]

Elgon_total<-Elgon%>% adorn_totals(c("col"))
Elgon_total <- mutate(Elgon_total, Elgon=rowSums(Elgon_total[7])/6)
Elgon_total <- Elgon_total[,c(1,8)]

Eldama_Ravine_total<-Eldama_Ravine%>% adorn_totals(c("col"))
Eldama_Ravine_total <- mutate(Eldama_Ravine_total, Eldama_Ravine=rowSums(Eldama_Ravine_total[5])/4)
Eldama_Ravine_total <- Eldama_Ravine_total[,c(1,6)]

Kakamega_total<-Kakamega%>% adorn_totals(c("col"))
Kakamega_total <- mutate(Kakamega_total, Kakamega=rowSums(Kakamega_total[14])/13)
Kakamega_total <- Kakamega_total[,c(1,15)]

ICIPE_total <-ICIPE%>% adorn_totals(c("col"))
ICIPE_total <- mutate(ICIPE_total, ICIPE=rowSums(ICIPE_total[10])/9)
ICIPE_total <- ICIPE_total[,c(1,11)]

Isiolo_total<-Isiolo%>% adorn_totals(c("col"))
Isiolo_total <- mutate(Isiolo_total, Isiolo=rowSums(Isiolo_total[10])/9)
Isiolo_total <- Isiolo_total[,c(1,11)]

Kasanga_Mwingi_total <- Kasanga_Mwingi%>%adorn_totals(c("col"))
Kasanga_Mwingi_total<- mutate(Kasanga_Mwingi_total, Kasanga_Mwingi=rowSums(Kasanga_Mwingi_total[7])/6)
Kasanga_Mwingi_total <- Kasanga_Mwingi_total[,c(1,8)]

Kasika_Mwingi_total<-Kasika_Mwingi%>% adorn_totals(c("col"))
Kasika_Mwingi_total <- mutate(Kasika_Mwingi_total, Kasika_Mwingi=rowSums(Kasika_Mwingi_total[7])/6)
Kasika_Mwingi_total <- Kasika_Mwingi_total[,c(1,8)]

Kib_Kakamega_total<-Kib_Kakamega%>% adorn_totals(c("col"))
Kib_Kakamega_total <- mutate(Kib_Kakamega_total, Kib_Kakamega=rowSums(Kib_Kakamega_total[4])/3)
Kib_Kakamega_total <- Kib_Kakamega_total[,c(1,5)]

Kisii_total<-Kisii%>% adorn_totals(c("col"))
Kisii_total <- mutate(Kisii_total, Kisii=rowSums(Kisii_total[16])/15)
Kisii_total <- Kisii_total[,c(1,17)]

Kitui_total<-Kitui%>% adorn_totals(c("col"))
Kitui_total <- mutate(Kitui_total, Kitui=rowSums(Kitui_total[7])/6)
Kitui_total <- Kitui_total[,c(1,8)]

KNN_Nyamira_total<-KNN_Nyamira%>% adorn_totals(c("col"))
KNN_Nyamira_total <- mutate(KNN_Nyamira_total, KNN_Nyamira=rowSums(KNN_Nyamira_total[3])/2)
KNN_Nyamira_total <- KNN_Nyamira_total[,c(1,4)]

maralal_total<-maralal%>% adorn_totals(c("col"))
maralal_total <- mutate(maralal_total, maralal=rowSums(maralal_total[6])/5)
maralal_total <- maralal_total[,c(1,7)]

Maru_Nairobi_total<-Maru_Nairobi%>% adorn_totals(c("col"))
Maru_Nairobi_total <- mutate(Maru_Nairobi_total, Maru_Nairobi=rowSums(Maru_Nairobi_total[4])/3)
Maru_Nairobi_total <- Maru_Nairobi_total[,c(1,5)]

Meru_total<-Meru%>% adorn_totals(c("col"))
Meru_total <- mutate(Meru_total, Meru=rowSums(Meru_total[8])/7)
Meru_total <- Meru_total[,c(1,9)]

Narok_total<-Narok%>% adorn_totals(c("col"))
Narok_total <- mutate(Narok_total, Narok=rowSums(Narok_total[5])/4)
Narok_total <- Narok_total[,c(1,6)]

phi_total<-phi%>% adorn_totals(c("col"))
phi_total <- mutate(phi_total, phi=rowSums(phi_total[3])/2)
phi_total <- phi_total[,c(1,4)]

saf_Nairobi_total<-saf_Nairobi%>% adorn_totals(c("col"))
saf_Nairobi_total <- mutate(saf_Nairobi_total, saf_Nairobi=rowSums(saf_Nairobi_total[8])/7)
saf_Nairobi_total <- saf_Nairobi_total[,c(1,9)]

Taita_total<-Taita%>% adorn_totals(c("col"))
Taita_total <- mutate(Taita_total, Taita=rowSums(Taita_total[9])/8)
Taita_total <- Taita_total[,c(1,10)]

Thog_Nairobi_total<-Thog_Nairobi%>% adorn_totals(c("col"))
Thog_Nairobi_total <- mutate(Thog_Nairobi_total, Thog_Nairobi=rowSums(Thog_Nairobi_total[5])/4)
Thog_Nairobi_total <- Thog_Nairobi_total[,c(1,6)]

Uasin_Gishu_total<-Uasin_Gishu%>% adorn_totals(c("col"))
Uasin_Gishu_total <- mutate(Uasin_Gishu_total, Uasin_Gishu=rowSums(Uasin_Gishu_total[8])/7)
Uasin_Gishu_total <- Uasin_Gishu_total[,c(1,9)]

Ukasic_Mwingi_total<-Ukasic_Mwingi%>% adorn_totals(c("col"))
Ukasic_Mwingi_total <- mutate(Ukasic_Mwingi_total, Ukasic_Mwingi=rowSums(Ukasic_Mwingi_total[5])/4)
Ukasic_Mwingi_total <- Ukasic_Mwingi_total[,c(1,6)]

Vin_Kakamega_total<-Vin_Kakamega%>% adorn_totals(c("col"))
Vin_Kakamega_total <- mutate(Vin_Kakamega_total, Vin_Kakamega=rowSums(Vin_Kakamega_total[6])/5)
Vin_Kakamega_total <- Vin_Kakamega_total[,c(1,7)]

Voic_Taita_total <-Voic_Taita%>% adorn_totals(c("col"))
Voic_Taita_total <- mutate(Voic_Taita_total, Voic_Taita=rowSums(Voic_Taita_total[3])/2)
Voic_Taita_total <- Voic_Taita_total[,c(1,4)]
```

```{r}
all_total <- all %>% adorn_totals(c("col"))
all_total <- mutate(all_total, all=rowSums(all_total[141])/140)
all_total <- all_total[,c(1,142)]
```


```{r}
# Merging
merged <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE),
list(Baringo_total, Bungoma_total, Bomet_total, Busia_total, Baringo2_total, Coast_total, Elgon_total, Eldama_Ravine_total,Kakamega_total, ICIPE_total,Isiolo_total,Kasanga_Mwingi_total,Kasika_Mwingi_total,Kib_Kakamega_total,Kisii_total,Kitui_total,KNN_Nyamira,maralal_total,Maru_Nairobi_total,Meru_total,Narok_total,phi_total,saf_Nairobi_total,Taita_total,Thog_Nairobi_total, Uasin_Gishu_total,Ukasic_Mwingi_total,Vin_Kakamega_total,all_total))
```

```{r}
names(merged)<-c('Genus','Baringo','Bungoma','Bomet','Busia','Baringo2','Coast','Elgon','Eldama_Ravine','Kakamega','ICIPE','Isiolo','Kasanga_Mwingi','Kasika_Mwingi','Kib_Kakamega','Kisii','Kitui','KNN_Nyamira','maralal','Maru_Nairobi','Meru','Narok','phi','saf_Nairobi','Taita','Thog_Nairobi','Uasin_Gishu','Ukasic_Mwingi','Vin_Kakamega','All')
```

```{r}
#calculating the total abundance per genus and ordering from the most abundant to the lowest
cumulation <- merged %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100
```

```{r}
head(cumulation, n=20)
head(cumulation$Genus, n=51)
```

```{r}
genus_Rep <- c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus")
```


```{r}
group <- merged
```


```{r}
Baringo<-group[,c(1:2)]
Bungoma<- group[,c(1:3)] 
Bomet<- group[,c(1:4)]
Busia<- group[,c(1:5)]
Baringo2<- group[,c(1:6)]
Coast<- group[,c(1:7)]
Elgon<- group[,c(1:8)]
Eldama_Ravine<- group[,c(1:9)]
Kakamega<- group[,c(1:10)]
ICIPE<- group[,c(1:11)]
Isiolo<- group[,c(1:12)]
Kasanga_Mwingi<- group[,c(1:13)]
Kasika_Mwingi<- group[,c(1:14)]
Kib_Kakamega<- group[,c(1:15)]
Kisii<- group[,c(1:16)]
Kitui<- group[,c(1:17)]
KNN_Nyamira<- group[,c(1:18)]
maralal<- group[,c(1:19)]
Maru_Nairobi<- group[,c(1:20)]
Meru<- group[,c(1:21)]
Narok<- group[,c(1:22)]
phi<- group[,c(1:23)]
saf_Nairobi<- group[,c(1:24)]
Taita<- group[,c(1:25)]
Thog_Nairobi<- group[,c(1:26)]
Uasin_Gishu<- group[,c(1:27)]
Ukasic_Mwingi<- group[,c(1:28)]
Vin_Kakamega<- group[,c(1:29)]
All<- group[,c(1:30)]
```


```{r}
# Viewing Sample Diversity
#install.packages("janitor")
library(janitor)
```


```{r}
#converting the abudances into percentage
bar_all <- adorn_percentages(All, denominator = "col", na.rm = T)
bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all<-bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all
write.csv(dist_all, "/home/icipe/Downloads/dada2/Abundance.csv")
```
```{r}
#gathering the data
bar_all <- bar_all %>%
  gather(value = "abundance", key = "Site", -Genus)
bar_all <- as.data.frame(gsub("\\(", " (", as.matrix(bar_all)))
```

```{r}
# coerce the dataframe columns into respective data type
bar_all$Genus <- as.factor(bar_all$Genus)
bar_all$Site <- as.character(bar_all$Site)
bar_all$abundance <- as.numeric(bar_all$abundance)
```


```{r}
#ordering the data for plotting
bar_all$Genus <- reorder(bar_all$Genus, bar_all$abundance)
bar_all$Genus <- factor(bar_all$Genus, levels=rev(levels(bar_all$Genus)))
bar_all$Genus <- factor(bar_all$Genus, 
                        levels=c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus"))
```

```{r}
# Defining the color pallete
myPalette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#99000D", "#E6AB02", "#A6761D", "#666666","#FDCDAC", "#1F78B4", "#B2DF8A", "#33A02C", "#CBD5E8", "#E31A1C", "#FDBF6F", "#FF7F00","#4A1486","#C0C0C0","#B3E2CD","#FFFF33", "#5172b2","#F4CAE4", "#E6F5C9", "#FCFBFD","#139BF1","#09FF00","#065535", "#1D91C0", "#C0FFEE","#B35806","#0C2C84","#D0ED0E","#092617","#8A17FD","#4D5D53","#E48400","#6082B6","#316689","#EEA47FFF","#738678","#645452","#EEA47FFF", "#00539CFF", "#FC766AFF", "#42EADDFF", "#00A4CCFF", "#69B3BB", "#B589D6","#D1DFB7","#97BC62FF","#D198C5FF","#000000", "#CBCE91FF", "#616247FF", "#D64161FF","#435E55FF", "#DD4132FF","#CE4A7EFF", "#BD7F37FF","#FFA351FF","#185E57")
```


```{r}
# Definig the names in italics
guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 10, 
face = "italic", colour = "Black", angle = 0)))
```

# Plotting barplot

```{r}
#library(Cairo)
library(forcats)

#plotting the barplot 
p_all <- ggplot(bar_all,aes(x = fct_inorder(Site), y = abundance), labs(fill= Genus), group=row.names(bar_all))+ xlab("Site")+ ylab("abundance") + geom_col(aes(fill = Genus),position = position_stack(reverse = FALSE))+
  theme(axis.text.x = element_text(angle = 72,size = 15, hjust = 1, face = "italic", family = "Arial"))+
   scale_fill_manual(values = myPalette)+
  #guides(fill = guide_legend(reverse = FALSE))+
  guide_italics+
  theme(legend.text = element_text(size = 8, colour = "black", face = "italic", family = "Arial"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 8, family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 8, family = "Arial"))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.15, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.position = "right", legend.justification = "top", legend.direction = "vertical", legend.text = element_text(size = 10))+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 8, family = "Arial"))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 8, family = "Arial"))
```


```{r}
#install.packages('extrafont')
library(extrafont)
#font_import
#BiocManager::install("ggpubr")
```


```{r}
#change the misspelled names
ps@sam_data[ps@sam_data == 'liotrigona'] <- 'Liotrigona'
```

```{r}
#check if it worked
ps@sam_data@.Data[[2]]
```

```{r}
# check the number of taxa
ntaxa(ps)
```

```{r}
# Per county
grouped_data


Baringo<-grouped_data[,c(1,2:4,20:23)]
Bungoma<-grouped_data[,c(1,5:10)]
Bomet<-grouped_data[,c(1,11:16)]
Busia<-grouped_data[,c(1,17:19)]
Coast<-grouped_data[,c(1,24:26)]
Elgon<-grouped_data[,c(1,27:31)]
Eldama_Ravine<-grouped_data[,c(1,32:34)]
Kakamega<-grouped_data[,c(1,35:46,73:74,136:139)]
Nairobi<-grouped_data[,c(1,47:54,99:100,110:116,124:126)]
Isiolo<-grouped_data[,c(1,55:62)]
Mwingi<-grouped_data[,c(1,63:67,68:72,133:135)]
Kisii<-grouped_data[,c(1,75:88)]
Kitui<-grouped_data[,c(1,89:93)]
Nyamira<-grouped_data[,c(1,94)]
Maralal<-grouped_data[,c(1,95:98)]
Meru<-grouped_data[,c(1,101:106)]
Narok<-grouped_data[,c(1,107:109)]
Taita<-grouped_data[,c(1,117:123,140)]
Uasin_Gishu<-grouped_data[,c(1,127:132)]
Voic_Taita<-grouped_data[,c(1,140)]
all<-grouped_data[,c(1,2:140)]

Baringo_total <-Baringo%>% adorn_totals(c("col"))
Baringo_total <- mutate(Baringo_total, Baringo=rowSums(Baringo_total[9])/8)
Baringo_total <- Baringo_total[,c(1,10)]

Bungoma_total <-Bungoma%>% adorn_totals(c("col"))
Bungoma_total <- mutate(Bungoma_total, Bungoma=rowSums(Bungoma_total[8])/7)
Bungoma_total <- Bungoma_total[,c(1,9)]

Bomet_total <-Bomet%>% adorn_totals(c("col"))
Bomet_total <- mutate(Bomet_total, Bomet=rowSums(Bomet_total[8])/7)
Bomet_total <- Bomet_total[,c(1,9)]

Busia_total<-Busia%>% adorn_totals(c("col"))
Busia_total <- mutate(Busia_total, Busia=rowSums(Busia_total[5])/4)
Busia_total <- Busia_total[,c(1,6)]


Coast_total<-Coast%>% adorn_totals(c("col"))
Coast_total <- mutate(Coast_total, Coast=rowSums(Coast_total[5])/4)
Coast_total <- Coast_total[,c(1,6)]

Elgon_total<-Elgon%>% adorn_totals(c("col"))
Elgon_total <- mutate(Elgon_total, Elgon=rowSums(Elgon_total[7])/6)
Elgon_total <- Elgon_total[,c(1,8)]

Eldama_Ravine_total<-Eldama_Ravine%>% adorn_totals(c("col"))
Eldama_Ravine_total <- mutate(Eldama_Ravine_total, Eldama_Ravine=rowSums(Eldama_Ravine_total[5])/4)
Eldama_Ravine_total <- Eldama_Ravine_total[,c(1,6)]

Kakamega_total<-Kakamega%>% adorn_totals(c("col"))
Kakamega_total <- mutate(Kakamega_total, Kakamega=rowSums(Kakamega_total[19])/20)
Kakamega_total <- Kakamega_total[,c(1,21)]

Nairobi_total <-Nairobi%>% adorn_totals(c("col"))
Nairobi_total <- mutate(Nairobi_total, Nairobi=rowSums(Nairobi_total[22])/21)
Nairobi_total <- Nairobi_total[,c(1,23)]

Isiolo_total<-Isiolo%>% adorn_totals(c("col"))
Isiolo_total <- mutate(Isiolo_total, Isiolo=rowSums(Isiolo_total[10])/9)
Isiolo_total <- Isiolo_total[,c(1,11)]

Mwingi_total <- Mwingi%>%adorn_totals(c("col"))
Mwingi_total<- mutate(Mwingi_total, Mwingi=rowSums(Mwingi_total[15])/14)
Mwingi_total <- Mwingi_total[,c(1,16)]


Kisii_total<-Kisii%>% adorn_totals(c("col"))
Kisii_total <- mutate(Kisii_total, Kisii=rowSums(Kisii_total[16])/15)
Kisii_total <- Kisii_total[,c(1,17)]

Kitui_total<-Kitui%>% adorn_totals(c("col"))
Kitui_total <- mutate(Kitui_total, Kitui=rowSums(Kitui_total[7])/6)
Kitui_total <- Kitui_total[,c(1,8)]

Nyamira_total<-Nyamira%>% adorn_totals(c("col"))
Nyamira_total <- mutate(Nyamira_total, Nyamira=rowSums(Nyamira_total[3])/2)
Nyamira_total <- Nyamira_total[,c(1,4)]

Maralal_total<-Maralal%>% adorn_totals(c("col"))
Maralal_total <- mutate(Maralal_total, Maralal=rowSums(Maralal_total[6])/5)
Maralal_total <- Maralal_total[,c(1,7)]

Meru_total<-Meru%>% adorn_totals(c("col"))
Meru_total <- mutate(Meru_total, Meru=rowSums(Meru_total[8])/7)
Meru_total <- Meru_total[,c(1,9)]

Narok_total<-Narok%>% adorn_totals(c("col"))
Narok_total <- mutate(Narok_total, Narok=rowSums(Narok_total[5])/4)
Narok_total <- Narok_total[,c(1,6)]

Taita_total<-Taita%>% adorn_totals(c("col"))
Taita_total <- mutate(Taita_total, Taita=rowSums(Taita_total[10])/9)
Taita_total <- Taita_total[,c(1,11)]

Uasin_Gishu_total<-Uasin_Gishu%>% adorn_totals(c("col"))
Uasin_Gishu_total <- mutate(Uasin_Gishu_total, Uasin_Gishu=rowSums(Uasin_Gishu_total[8])/7)
Uasin_Gishu_total <- Uasin_Gishu_total[,c(1,9)]

all_total <- all %>% adorn_totals(c("col"))
all_total <- mutate(all_total, all=rowSums(all_total[141])/140)
all_total <- all_total[,c(1,142)]

merged1 <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE),
                 list(Baringo_total, Bungoma_total, Bomet_total, Busia_total,Coast_total, Elgon_total, Eldama_Ravine_total,Kakamega_total, Nairobi_total,Isiolo_total,Mwingi_total,Kisii_total,Kitui_total,Nyamira_total,Maralal_total,Meru_total,Narok_total,Taita_total,Uasin_Gishu_total,all_total))

cumulation <- merged1 %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100

genus_Rep <-c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus")

group <- aggregate(merged1[-1], list(Genus = replace(merged1$Genus,!(merged1$Genus %in% genus_Rep), "Others")), sum)

Baringo<-group[,c(1:2)]
Bungoma<- group[,c(1:3)] 
Bomet<- group[,c(1:4)]
Busia<- group[,c(1:5)]
Coast<- group[,c(1:6)]
Elgon<- group[,c(1:7)]
Eldama_Ravine<- group[,c(1:8)]
Kakamega<- group[,c(1:9)]
Nairobi<- group[,c(1:10)]
Isiolo<- group[,c(1:11)]
Mwingi<- group[,c(1:12)]
Kisii<- group[,c(1:13)]
Kitui<- group[,c(1:14)]
Nyamira<- group[,c(1:15)]
Maralal<- group[,c(1:16)]
Meru<- group[,c(1:17)]
Narok<- group[,c(1:18)]
Taita<- group[,c(1:19)]
Uasin_Gishu<- group[,c(1:20)]
All<- group[,c(1:20)]

bar_all <- adorn_percentages(All, denominator = "col", na.rm = T)
bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all<-bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all

bar_all <- bar_all %>%
  gather(value = "abundance", key = "Site", -Genus)
bar_all <- as.data.frame(gsub("\\(", " (", as.matrix(bar_all)))

bar_all$Genus <- as.factor(bar_all$Genus)
bar_all$Site <- as.character(bar_all$Site)
bar_all$abundance <- as.numeric(bar_all$abundance)

bar_all$Genus <- reorder(bar_all$Genus, bar_all$abundance)
bar_all$Genus <- factor(bar_all$Genus, levels=rev(levels(bar_all$Genus)))
bar_all$Genus <- factor(bar_all$Genus, 
                        levels=genus_Rep)

myPalette <- c("#D95F02", "#7570B3", "#E7298A", "#99000D", "#E6AB02", "#A6761D", "#666666","#FDCDAC", "#1F78B4", "#B2DF8A", "#33A02C", "#CBD5E8", "#E31A1C", "#FDBF6F", "#FF7F00","#4A1486","#C0C0C0","#B3E2CD","#FFFF33", "#5172b2","#F4CAE4", "#E6F5C9", "#FCFBFD","#139BF1","#09FF00","#065535", "#1D91C0", "#C0FFEE","#B35806","#0C2C84","#D0ED0E","#092617","#499976","#4D5D53","#E48400","#6082B6","#316689","#CEFB02","#738678","#645452","#EEA47FFF", "#00539CFF", "#FC766AFF", "#42EADDFF", "#00A4CCFF", "#69B3BB", "#B589D6","#D1DFB7","#97BC62FF","#D198C5FF","#000000")
length(myPalette)

guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 15, 
                                                                       face = "italic", colour = "Black", angle = 0)))

p_all <- ggplot(bar_all,aes(x = fct_inorder(Site), y = abundance), labs(fill= Genus), group=row.names(bar_all))+ xlab("County")+ ylab("abundance") + geom_col(aes(fill = Genus),position = position_stack(reverse = FALSE))+
  theme(axis.text.x = element_text(angle = 72,size = 15, hjust = 1, face = "italic", family = "Arial"))+
  scale_fill_manual(values = myPalette)+
  #guides(fill = guide_legend(reverse = FALSE))+
  guide_italics+
  theme(legend.text = element_text(size = 8, colour = "black", face = "italic", family = "Arial"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 8, family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 25, family = "Arial"))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.15, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.position = "right", legend.justification = "top", legend.direction = "vertical", legend.text = element_text(size = 10))+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 15, family = "Arial"))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 15, family = "Arial"))

plot4_60<-p_all + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + xlab("County")
plot4_60

ggsave("/home/icipe/Downloads/dada2/plot4_60.jpeg", width=12, height=12, dpi=600)
```

```{r}
merged <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE),
list(Baringo_total, Bungoma_total, Bomet_total, Busia_total,Coast_total, Elgon_total, Eldama_Ravine_total,Kakamega_total, Nairobi_total,Isiolo_total,Mwingi_total,Kisii_total,Kitui_total,Nyamira_total,Maralal_total,Meru_total,Narok_total,Taita_total,Uasin_Gishu_total,all_total))
```

```{r}
#calculating the total abundance per genus and ordering from the most abundant to the lowest
cumulation <- merged %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100

head(cumulation$Genus, n=51)
```

```{r}
genus_Rep <- c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus")
```

```{r}
group <- aggregate(merged[-1], list(Genus = replace(merged$Genus,!(merged$Genus %in% genus_Rep), "Others")), sum)

Baringo<-group[,c(1:2)]
Bungoma<- group[,c(1:3)] 
Bomet<- group[,c(1:4)]
Busia<- group[,c(1:5)]
Coast<- group[,c(1:6)]
Elgon<- group[,c(1:7)]
Eldama_Ravine<- group[,c(1:8)]
Kakamega<- group[,c(1:9)]
Nairobi<- group[,c(1:10)]
Isiolo<- group[,c(1:11)]
Mwingi<- group[,c(1:12)]
Kisii<- group[,c(1:13)]
Kitui<- group[,c(1:14)]
Nyamira<- group[,c(1:15)]
Maralal<- group[,c(1:16)]
Meru<- group[,c(1:17)]
Narok<- group[,c(1:18)]
Taita<- group[,c(1:19)]
Uasin_Gishu<- group[,c(1:20)]
All<- group[,c(1:20)]

```

```{r}
#converting the abudances into percentage
bar_all <- adorn_percentages(All, denominator = "col", na.rm = T)
bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all<-bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all
write.csv(dist_all, "/home/icipe/Downloads/dada2/Abundance.csv")

#gathering the data
bar_all <- bar_all %>%
  gather(value = "abundance", key = "Site", -Genus)
bar_all <- as.data.frame(gsub("\\(", " (", as.matrix(bar_all)))

# coerce the dataframe columns into respective data type
bar_all$Genus <- as.factor(bar_all$Genus)
bar_all$Site <- as.character(bar_all$Site)
bar_all$abundance <- as.numeric(bar_all$abundance)

#ordering the data for plotting
bar_all$Genus <- reorder(bar_all$Genus, bar_all$abundance)
bar_all$Genus <- factor(bar_all$Genus, levels=rev(levels(bar_all$Genus)))
bar_all$Genus <- factor(bar_all$Genus, 
                        levels=c("Emilia","Erechtites","Raphanus","Bidens","Ageratum","Dovyalis","Macadamia","Brassica","Baccharoides","Teclea","Vigna","Syzygium","Senecio","Artemisia","Salvadora","Flueggea","Amaranthus","Mikania","Parthenium","Blepharispermum","Xanthium","Clitoria","Cleome","Euryops","Sambucus","Cyathula","Ligustrum","Sersalisia","Poupartia","Achyranthes","Calopyxis","Scalesia","Phytolacca","Antiaris","Tropaeolum","Lablab","Zornia","Oldenlandia","Conyza","Combretum","Crassocephalum","Bougainvillea","Croton","Digera","Artocarpus", "Ormocarpum","Tridax","Fraxinus","Acilepis","Senna","Erucastrum","Erlangea","Allium","Helianthus","Biancaea","Ajuga","Hirpicium","Micrococca","Indigofera","Melaleuca","Others"))

# Defining the color pallete
myPalette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#99000D", "#E6AB02", "#A6761D", "#666666","#FDCDAC", "#1F78B4", "#B2DF8A", "#33A02C", "#CBD5E8", "#E31A1C", "#FDBF6F", "#FF7F00","#4A1486","#C0C0C0","#B3E2CD","#FFFF33", "#5172b2","#F4CAE4", "#E6F5C9", "#FCFBFD","#139BF1","#09FF00","#065535", "#1D91C0", "#C0FFEE","#B35806","#0C2C84","#D0ED0E","#092617","#499976","#4D5D53","#E48400","#6082B6","#316689","#CEFB02","#738678","#645452","#EEA47FFF", "#00539CFF", "#FC766AFF", "#42EADDFF", "#00A4CCFF", "#69B3BB", "#B589D6","#D1DFB7","#97BC62FF","#D198C5FF","#000000", "#CBCE91FF", "#616247FF", "#D64161FF","#435E55FF", "#DD4132FF","#CE4A7EFF", "#BD7F37FF","#FFA351FF","#185E57")
length(myPalette)

# Definig the names in italics
guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 15, 
face = "italic", colour = "Black", angle = 0)))
```

```{r}
group <- aggregate(merged[-1], list(Genus = replace(merged$Genus,!(merged$Genus %in% genus_Rep), "Others")), sum)

Baringo<-group[,c(1:2)]
Bungoma<- group[,c(1:3)] 
Bomet<- group[,c(1:4)]
Busia<- group[,c(1:5)]
Coast<- group[,c(1:6)]
Elgon<- group[,c(1:7)]
Eldama_Ravine<- group[,c(1:8)]
Kakamega<- group[,c(1:9)]
Nairobi<- group[,c(1:10)]
Isiolo<- group[,c(1:11)]
Mwingi<- group[,c(1:12)]
Kisii<- group[,c(1:13)]
Kitui<- group[,c(1:14)]
Nyamira<- group[,c(1:15)]
Maralal<- group[,c(1:16)]
Meru<- group[,c(1:17)]
Narok<- group[,c(1:18)]
Taita<- group[,c(1:19)]
Uasin_Gishu<- group[,c(1:20)]
All<- group[,c(1:20)]
```

```{r}
merged <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE),
list(Baringo_total, Bungoma_total, Bomet_total, Busia_total,Coast_total, Elgon_total, Eldama_Ravine_total,Kakamega_total, Nairobi_total,Isiolo_total,Mwingi_total,Kisii_total,Kitui_total,Nyamira_total,Maralal_total,Meru_total,Narok_total,Taita_total,Uasin_Gishu_total,all_total))
```

```{r}
#converting the abudances into percentage
bar_all <- adorn_percentages(All, denominator = "col", na.rm = T)
bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all<-bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all
write.csv(dist_all, "/home/icipe/Downloads/dada2/Abundance.csv")
```

```{r}
#gathering the data
bar_all <- bar_all %>%
  gather(value = "abundance", key = "Site", -Genus)
bar_all <- as.data.frame(gsub("\\(", " (", as.matrix(bar_all)))

# coerce the dataframe columns into respective data type
bar_all$Genus <- as.factor(bar_all$Genus)
bar_all$Site <- as.character(bar_all$Site)
bar_all$abundance <- as.numeric(bar_all$abundance)
```

```{r}
#ordering the data for plotting
bar_all$Genus <- reorder(bar_all$Genus, bar_all$abundance)
bar_all$Genus <- factor(bar_all$Genus, levels=rev(levels(bar_all$Genus)))
bar_all$Genus <- factor(bar_all$Genus, 
                        levels=c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus"))
```


```{r}
# Defining the color pallete
myPalette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#99000D", "#E6AB02", "#A6761D", "#666666","#FDCDAC", "#1F78B4", "#B2DF8A", "#33A02C", "#CBD5E8", "#E31A1C", "#FDBF6F", "#FF7F00","#4A1486","#C0C0C0","#B3E2CD","#FFFF33", "#5172b2","#F4CAE4", "#E6F5C9", "#FCFBFD","#139BF1","#09FF00","#065535", "#1D91C0", "#C0FFEE","#B35806","#0C2C84","#D0ED0E","#092617","#8A17FD","#4D5D53","#E48400","#6082B6","#316689","#EEA47FFF","#738678","#645452","#EEA47FFF", "#00539CFF", "#FC766AFF", "#42EADDFF", "#00A4CCFF", "#69B3BB", "#B589D6","#D1DFB7","#97BC62FF","#D198C5FF","#000000", "#CBCE91FF", "#616247FF", "#D64161FF","#435E55FF", "#DD4132FF","#CE4A7EFF", "#BD7F37FF","#FFA351FF","#185E57")
```


```{r}
#plotting the barplot 
p_all <- ggplot(bar_all,aes(x = fct_inorder(Site), y = abundance), labs(fill= Genus), group=row.names(bar_all))+ xlab("Site")+ ylab("abundance") + geom_col(aes(fill = Genus),position = position_stack(reverse = FALSE))+
  theme(axis.text.x = element_text(angle = 72,size = 15, hjust = 1, face = "italic", family = "Arial"))+
   scale_fill_manual(values = myPalette)+
  #guides(fill = guide_legend(reverse = FALSE))+
  guide_italics+
  theme(legend.text = element_text(size = 8, colour = "black", face = "italic", family = "Arial"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 8, family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 25, family = "Arial"))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.15, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.position = "right", legend.justification = "top", legend.direction = "vertical", legend.text = element_text(size = 10))+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 15, family = "Arial"))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 15, family = "Arial"))
```


```{r}
# Definig the names in italics
guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 10, 
face = "italic", colour = "Black", angle = 0)))
```


```{r}
#install.packages('extrafont')
library(extrafont)
#font_import()

plot4_60<-p_all + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + xlab("Site")
plot4_60
```

```{r}
save.image("./Abundance.RData")

#alpha diversity

#Shannon

plotA = plot_richness(ps, x="site",color="genus", measures=("Shannon") ) +
  geom_boxplot(color = "blue")+
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -72))
plotA
```
```{r}
ggsave("/home/icipe/Downloads/dada2/Shanon_graph.jpeg", width=12, height = 12, dpi = 600 )
```


```{r}
#Simpson
plotB <-plot_richness(ps, x="site",color="genus", measures=c("Simpson") ) +
  geom_boxplot(color = "red")+
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))
plotB
```


```{r}
#Chao1
plotC <-plot_richness(ps, x="site",color="genus", measures=c("Chao1") ) +
geom_boxplot(color = "red")+
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))
plotC
```


```{r}
#ACE
plotD = plot_richness(ps, x="site",color="genus", measures=c("ACE") ) +
  geom_boxplot(color = "red") +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))
plotD
```

```{r}
#arrange the plots
mega_plot = ggarrange(plotA,plotB,plotC,plotD, ncol = 2, nrow = 2)
mega_plot
```


```{r}
#save the plots into different formats(for genus)
ggsave("alphaplot.jpeg", plot = mega_plot, width = 12, height = 12, dpi = 600)
ggsave("alphaplot.png", plot = mega_plot, width = 12 , height = 12, dpi = 600)
ggsave("alphaplot.svg", plot = mega_plot, width = 12, height = 12 , dpi = 600)
ggsave("alphaplot.tiff", plot = mega_plot, width = 12, height = 12, dpi = 600)
```

#Beta diversity
#PCA
```{r}
Site <- c('Baringo','Bungoma','Bomet','Busia','Baringo2','Coast','Elgon','Eldama_Ravine','Kakamega','ICIPE','Isiolo','Kasanga_Mwingi','Kasika_Mwingi','Kib_Kakamega','Kisii','Kitui','KNN_Nyamira','maralal','Maru_Nairobi','Meru','Narok','phi','saf_Nairobi','Taita','Thog_Nairobi','Uasin_Gishu','Ukasic_Mwingi','Vin_Kakamega')
```

```{r}
p1 <- ggplot(ps, aes(x=Site, y=genus_Rep)) +
  geom_boxplot(fill=c("blue","red")) +
  labs(title= 'Richness', x= ' ', y= '', tag = "A") +
  geom_point()
p1
```


```{r}
#Plot the PCA
pcaplot1 = ggordpoint(obj = pca, biplot = T, speciesannot = T,
                      factorNames = ("genus"), ellipse = T) +
  scale_color_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#FFA351FF")) +
  scale_fill_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA","#FFA351FF"))
pcaplot1
```

```{r}
#PCOA
pcoa = get_pcoa(obj = ps, distmethod = "euclidean", method = "hellinger")
pcoa
```


```{r}
#Plot the PCOA
pcoaplot1 = ggordpoint(obj = pcoa, biplot = T, speciesannot = F,showsample = T,
                        factorNames=c("genus"), ellipse = T) +
  scale_color_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#FFA351FF")) +
  scale_fill_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#FFA351FF"))
pcoaplot1
```


```{r}
#save the plots## Alpha Diversity

ps3
set.seed(1024)
rareres <- get_rarecurve(obj=ps3, chunks=400)

p_rare <- ggrarecurve(obj=rareres,
                      indexNames=c("Observe","Chao1","ACE","Shannon"),
                      ) +
          theme(legend.spacing.y=unit(0.01,"cm"),
                legend.text=element_text(size=4))

prare1 <- ggrarecurve(obj=rareres, factorNames="Animal_Type",
                # metagenomic analysis using dada2

# loading of packages

```{r}
library("phyloseq")
library("vegan")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library("janitor")
library("coin")
library("reshape2")
library("ggnewscale")
library("MicrobiotaProcess")
library("patchwork")
library("dada2")
library("csv")
library("dplyr")
library("ggtree")
library("VennDiagram")
library("UpSetR")
#BiocManager::install("MicrobiotaProcess")
install.packages("MicrobiotaProcess")
library("MicrobiotaProcess")
library("utils")
library("utf8")
library("metagMisc")
```


```{r}
#BiocManager::install("devtools")
library("devtools")
library("tidyverse")
#devtools::install_github(repo= "jbisanz/qiime2R")
library(qiime2R)
asv_table <- read.table("/home/icipe/Downloads/dada2/ASV_table.tsv")
```

# loading the data

```{r}
metadata <- read.table("/home/icipe/Downloads/Geo_metadata.tsv", sep = "\t", row.names = 1, header= TRUE)

ot <- read.table("/home/icipe/Downloads/dada2/ASV_table.tsv", header = TRUE, sep = "\t",  row.names = 1 )
class(ot)
ot <- as.matrix(ot)
class(ot)

tax <- read.table("/home/icipe/Downloads/dada2/ASV_tax_species.tsv", header = TRUE, sep = "\t",  row.names = 1)
#tax <- data.frame(tax)
Tax <- as.matrix(tax)
#tax <- t(tax)

ot <- phyloseq::otu_table(ot, taxa_are_rows = TRUE)
tax <- phyloseq::tax_table(Tax)
metadata <- phyloseq::sample_data(metadata)

```
```{r}
ps <- phyloseq(ot, metadata, tax)
ps
```
```{r}
# we now create a third object called random for merging with the other three object
library("ape")
random_tree <- rtree(ntaxa(physeq), rooted = TRUE, tip.label = taxa_names(physeq))
plot(random_tree)
```

```{r}
#merging the preceeding 3 objects.
ps <- merge_phyloseq(physeq, random_tree)
ps
```

```{r}
# Removing the NAs
tx <- phyloseq_to_df(ps, addtax = T, addtot = F, addmaxrank = F)
tx <-tx[, c(8, 12:150 )]

tx <- tx %>%
  # recode empty strings "" by NAs
 na_if("") %>%
  # remove NAs
  na.omit()
tx
```

```{r}
cumulation <- tx %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation
```
```{r}
#merging the blast taxonomic classification to blast abundance table
merged_data <- tx
write.csv(merged_data, file="/home/icipe/Downloads/dada2/merged_data.csv")
# Group the data
grouped_data <- tx %>% group_by(Genus) %>%summarise_if(is.numeric, sum)
#grouped_data<-grouped_data[-c(1),]
```


```{r}
saveRDS(grouped_data,"/home/icipe/Downloads/dada2/grouped_data.rds")
write.table(tx,"/home/icipe/Downloads/dada2/tx.csv")
```


```{r}
# Grouping data by location
Baringo<-grouped_data[,c(1,2:4)]
Bungoma<-grouped_data[,c(1,5:10)]
Bomet<-grouped_data[,c(1,11:16)]
Busia<-grouped_data[,c(1,17:19)]
Baringo2<-grouped_data[,c(1,20:23)]
Coast<-grouped_data[,c(1,24:26)]
Elgon<-grouped_data[,c(1,27:31)]
Eldama_Ravine<-grouped_data[,c(1,32:34)]
Kakamega<-grouped_data[,c(1,35:46)]
ICIPE<-grouped_data[,c(1,47:54)]
Isiolo<-grouped_data[,c(1,55:62)]
Kasanga_Mwingi<-grouped_data[,c(1,63:67)]
Kasika_Mwingi<-grouped_data[,c(1,68:72)]
Kib_Kakamega<-grouped_data[,c(1,73:74)]
Kisii<-grouped_data[,c(1,75:88)]
Kitui<-grouped_data[,c(1,89:93)]
KNN_Nyamira<-grouped_data[,c(1,94)]
maralal<-grouped_data[,c(1,95:98)]
Maru_Nairobi<-grouped_data[,c(1,99:100)]
Meru<-grouped_data[,c(1,101:106)]
Narok<-grouped_data[,c(1,107:109)]
phi<-grouped_data[,c(1,110)]
saf_Nairobi<-grouped_data[,c(1,111:116)]
Taita<-grouped_data[,c(1,117:123)]
Thog_Nairobi<-grouped_data[,c(1,124:126)]
Uasin_Gishu<-grouped_data[,c(1,127:132)]
Ukasic_Mwingi<-grouped_data[,c(1,133:135)]
Vin_Kakamega<-grouped_data[,c(1,136:139)]
Voic_Taita <- grouped_data[1,140]
```

```{r}
all<-grouped_data[,c(1,2:140)]
all
```

```{r}
Baringo_total <-Baringo%>% adorn_totals(c("col"))
Baringo_total <- mutate(Baringo_total, Baringo=rowSums(Baringo_total[5])/4)
Baringo_total <- Baringo_total[,c(1,6)]

Bungoma_total <-Bungoma%>% adorn_totals(c("col"))
Bungoma_total <- mutate(Bungoma_total, Bungoma=rowSums(Bungoma_total[8])/7)
Bungoma_total <- Bungoma_total[,c(1,9)]

Bomet_total <-Bomet%>% adorn_totals(c("col"))
Bomet_total <- mutate(Bomet_total, Bomet=rowSums(Bomet_total[8])/7)
Bomet_total <- Bomet_total[,c(1,9)]

Busia_total<-Busia%>% adorn_totals(c("col"))
Busia_total <- mutate(Busia_total, Busia=rowSums(Busia_total[5])/4)
Busia_total <- Busia_total[,c(1,6)]

Baringo2_total<-Baringo2%>% adorn_totals(c("col"))
Baringo2_total <- mutate(Baringo2_total, Baringo2=rowSums(Baringo2_total[6])/5)
Baringo2_total <- Baringo2_total[,c(1,7)]

Coast_total<-Coast%>% adorn_totals(c("col"))
Coast_total <- mutate(Coast_total, Coast=rowSums(Coast_total[5])/4)
Coast_total <- Coast_total[,c(1,6)]

Elgon_total<-Elgon%>% adorn_totals(c("col"))
Elgon_total <- mutate(Elgon_total, Elgon=rowSums(Elgon_total[7])/6)
Elgon_total <- Elgon_total[,c(1,8)]

Eldama_Ravine_total<-Eldama_Ravine%>% adorn_totals(c("col"))
Eldama_Ravine_total <- mutate(Eldama_Ravine_total, Eldama_Ravine=rowSums(Eldama_Ravine_total[5])/4)
Eldama_Ravine_total <- Eldama_Ravine_total[,c(1,6)]

Kakamega_total<-Kakamega%>% adorn_totals(c("col"))
Kakamega_total <- mutate(Kakamega_total, Kakamega=rowSums(Kakamega_total[14])/13)
Kakamega_total <- Kakamega_total[,c(1,15)]

ICIPE_total <-ICIPE%>% adorn_totals(c("col"))
ICIPE_total <- mutate(ICIPE_total, ICIPE=rowSums(ICIPE_total[10])/9)
ICIPE_total <- ICIPE_total[,c(1,11)]

Isiolo_total<-Isiolo%>% adorn_totals(c("col"))
Isiolo_total <- mutate(Isiolo_total, Isiolo=rowSums(Isiolo_total[10])/9)
Isiolo_total <- Isiolo_total[,c(1,11)]

Kasanga_Mwingi_total <- Kasanga_Mwingi%>%adorn_totals(c("col"))
Kasanga_Mwingi_total<- mutate(Kasanga_Mwingi_total, Kasanga_Mwingi=rowSums(Kasanga_Mwingi_total[7])/6)
Kasanga_Mwingi_total <- Kasanga_Mwingi_total[,c(1,8)]

Kasika_Mwingi_total<-Kasika_Mwingi%>% adorn_totals(c("col"))
Kasika_Mwingi_total <- mutate(Kasika_Mwingi_total, Kasika_Mwingi=rowSums(Kasika_Mwingi_total[7])/6)
Kasika_Mwingi_total <- Kasika_Mwingi_total[,c(1,8)]

Kib_Kakamega_total<-Kib_Kakamega%>% adorn_totals(c("col"))
Kib_Kakamega_total <- mutate(Kib_Kakamega_total, Kib_Kakamega=rowSums(Kib_Kakamega_total[4])/3)
Kib_Kakamega_total <- Kib_Kakamega_total[,c(1,5)]

Kisii_total<-Kisii%>% adorn_totals(c("col"))
Kisii_total <- mutate(Kisii_total, Kisii=rowSums(Kisii_total[16])/15)
Kisii_total <- Kisii_total[,c(1,17)]

Kitui_total<-Kitui%>% adorn_totals(c("col"))
Kitui_total <- mutate(Kitui_total, Kitui=rowSums(Kitui_total[7])/6)
Kitui_total <- Kitui_total[,c(1,8)]

KNN_Nyamira_total<-KNN_Nyamira%>% adorn_totals(c("col"))
KNN_Nyamira_total <- mutate(KNN_Nyamira_total, KNN_Nyamira=rowSums(KNN_Nyamira_total[3])/2)
KNN_Nyamira_total <- KNN_Nyamira_total[,c(1,4)]

maralal_total<-maralal%>% adorn_totals(c("col"))
maralal_total <- mutate(maralal_total, maralal=rowSums(maralal_total[6])/5)
maralal_total <- maralal_total[,c(1,7)]

Maru_Nairobi_total<-Maru_Nairobi%>% adorn_totals(c("col"))
Maru_Nairobi_total <- mutate(Maru_Nairobi_total, Maru_Nairobi=rowSums(Maru_Nairobi_total[4])/3)
Maru_Nairobi_total <- Maru_Nairobi_total[,c(1,5)]

Meru_total<-Meru%>% adorn_totals(c("col"))
Meru_total <- mutate(Meru_total, Meru=rowSums(Meru_total[8])/7)
Meru_total <- Meru_total[,c(1,9)]

Narok_total<-Narok%>% adorn_totals(c("col"))
Narok_total <- mutate(Narok_total, Narok=rowSums(Narok_total[5])/4)
Narok_total <- Narok_total[,c(1,6)]

phi_total<-phi%>% adorn_totals(c("col"))
phi_total <- mutate(phi_total, phi=rowSums(phi_total[3])/2)
phi_total <- phi_total[,c(1,4)]

saf_Nairobi_total<-saf_Nairobi%>% adorn_totals(c("col"))
saf_Nairobi_total <- mutate(saf_Nairobi_total, saf_Nairobi=rowSums(saf_Nairobi_total[8])/7)
saf_Nairobi_total <- saf_Nairobi_total[,c(1,9)]

Taita_total<-Taita%>% adorn_totals(c("col"))
Taita_total <- mutate(Taita_total, Taita=rowSums(Taita_total[9])/8)
Taita_total <- Taita_total[,c(1,10)]

Thog_Nairobi_total<-Thog_Nairobi%>% adorn_totals(c("col"))
Thog_Nairobi_total <- mutate(Thog_Nairobi_total, Thog_Nairobi=rowSums(Thog_Nairobi_total[5])/4)
Thog_Nairobi_total <- Thog_Nairobi_total[,c(1,6)]

Uasin_Gishu_total<-Uasin_Gishu%>% adorn_totals(c("col"))
Uasin_Gishu_total <- mutate(Uasin_Gishu_total, Uasin_Gishu=rowSums(Uasin_Gishu_total[8])/7)
Uasin_Gishu_total <- Uasin_Gishu_total[,c(1,9)]

Ukasic_Mwingi_total<-Ukasic_Mwingi%>% adorn_totals(c("col"))
Ukasic_Mwingi_total <- mutate(Ukasic_Mwingi_total, Ukasic_Mwingi=rowSums(Ukasic_Mwingi_total[5])/4)
Ukasic_Mwingi_total <- Ukasic_Mwingi_total[,c(1,6)]

Vin_Kakamega_total<-Vin_Kakamega%>% adorn_totals(c("col"))
Vin_Kakamega_total <- mutate(Vin_Kakamega_total, Vin_Kakamega=rowSums(Vin_Kakamega_total[6])/5)
Vin_Kakamega_total <- Vin_Kakamega_total[,c(1,7)]

Voic_Taita_total <-Voic_Taita%>% adorn_totals(c("col"))
Voic_Taita_total <- mutate(Voic_Taita_total, Voic_Taita=rowSums(Voic_Taita_total[3])/2)
Voic_Taita_total <- Voic_Taita_total[,c(1,4)]
```

```{r}
all_total <- all %>% adorn_totals(c("col"))
all_total <- mutate(all_total, all=rowSums(all_total[141])/140)
all_total <- all_total[,c(1,142)]
```


```{r}
# Merging
merged <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE),
list(Baringo_total, Bungoma_total, Bomet_total, Busia_total, Baringo2_total, Coast_total, Elgon_total, Eldama_Ravine_total,Kakamega_total, ICIPE_total,Isiolo_total,Kasanga_Mwingi_total,Kasika_Mwingi_total,Kib_Kakamega_total,Kisii_total,Kitui_total,KNN_Nyamira,maralal_total,Maru_Nairobi_total,Meru_total,Narok_total,phi_total,saf_Nairobi_total,Taita_total,Thog_Nairobi_total, Uasin_Gishu_total,Ukasic_Mwingi_total,Vin_Kakamega_total,all_total))
```

```{r}
names(merged)<-c('Genus','Baringo','Bungoma','Bomet','Busia','Baringo2','Coast','Elgon','Eldama_Ravine','Kakamega','ICIPE','Isiolo','Kasanga_Mwingi','Kasika_Mwingi','Kib_Kakamega','Kisii','Kitui','KNN_Nyamira','maralal','Maru_Nairobi','Meru','Narok','phi','saf_Nairobi','Taita','Thog_Nairobi','Uasin_Gishu','Ukasic_Mwingi','Vin_Kakamega','All')
```

```{r}
#calculating the total abundance per genus and ordering from the most abundant to the lowest
cumulation <- merged %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100
```

```{r}
head(cumulation, n=20)
head(cumulation$Genus, n=51)
```

```{r}
genus_Rep <- c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus")
```


```{r}
group <- merged
```


```{r}
Baringo<-group[,c(1:2)]
Bungoma<- group[,c(1:3)] 
Bomet<- group[,c(1:4)]
Busia<- group[,c(1:5)]
Baringo2<- group[,c(1:6)]
Coast<- group[,c(1:7)]
Elgon<- group[,c(1:8)]
Eldama_Ravine<- group[,c(1:9)]
Kakamega<- group[,c(1:10)]
ICIPE<- group[,c(1:11)]
Isiolo<- group[,c(1:12)]
Kasanga_Mwingi<- group[,c(1:13)]
Kasika_Mwingi<- group[,c(1:14)]
Kib_Kakamega<- group[,c(1:15)]
Kisii<- group[,c(1:16)]
Kitui<- group[,c(1:17)]
KNN_Nyamira<- group[,c(1:18)]
maralal<- group[,c(1:19)]
Maru_Nairobi<- group[,c(1:20)]
Meru<- group[,c(1:21)]
Narok<- group[,c(1:22)]
phi<- group[,c(1:23)]
saf_Nairobi<- group[,c(1:24)]
Taita<- group[,c(1:25)]
Thog_Nairobi<- group[,c(1:26)]
Uasin_Gishu<- group[,c(1:27)]
Ukasic_Mwingi<- group[,c(1:28)]
Vin_Kakamega<- group[,c(1:29)]
All<- group[,c(1:30)]
```


```{r}
# Viewing Sample Diversity
#install.packages("janitor")
library(janitor)
```


```{r}
#converting the abudances into percentage
bar_all <- adorn_percentages(All, denominator = "col", na.rm = T)
bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all<-bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all
write.csv(dist_all, "/home/icipe/Downloads/dada2/Abundance.csv")
```
```{r}
#gathering the data
bar_all <- bar_all %>%
  gather(value = "abundance", key = "Site", -Genus)
bar_all <- as.data.frame(gsub("\\(", " (", as.matrix(bar_all)))
```

```{r}
# coerce the dataframe columns into respective data type
bar_all$Genus <- as.factor(bar_all$Genus)
bar_all$Site <- as.character(bar_all$Site)
bar_all$abundance <- as.numeric(bar_all$abundance)
```


```{r}
#ordering the data for plotting
bar_all$Genus <- reorder(bar_all$Genus, bar_all$abundance)
bar_all$Genus <- factor(bar_all$Genus, levels=rev(levels(bar_all$Genus)))
bar_all$Genus <- factor(bar_all$Genus, 
                        levels=c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus"))
```

```{r}
# Defining the color pallete
myPalette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#99000D", "#E6AB02", "#A6761D", "#666666","#FDCDAC", "#1F78B4", "#B2DF8A", "#33A02C", "#CBD5E8", "#E31A1C", "#FDBF6F", "#FF7F00","#4A1486","#C0C0C0","#B3E2CD","#FFFF33", "#5172b2","#F4CAE4", "#E6F5C9", "#FCFBFD","#139BF1","#09FF00","#065535", "#1D91C0", "#C0FFEE","#B35806","#0C2C84","#D0ED0E","#092617","#8A17FD","#4D5D53","#E48400","#6082B6","#316689","#EEA47FFF","#738678","#645452","#EEA47FFF", "#00539CFF", "#FC766AFF", "#42EADDFF", "#00A4CCFF", "#69B3BB", "#B589D6","#D1DFB7","#97BC62FF","#D198C5FF","#000000", "#CBCE91FF", "#616247FF", "#D64161FF","#435E55FF", "#DD4132FF","#CE4A7EFF", "#BD7F37FF","#FFA351FF","#185E57")
```


```{r}
# Definig the names in italics
guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 10, 
face = "italic", colour = "Black", angle = 0)))
```

# Plotting barplot

```{r}
#library(Cairo)
library(forcats)

#plotting the barplot 
p_all <- ggplot(bar_all,aes(x = fct_inorder(Site), y = abundance), labs(fill= Genus), group=row.names(bar_all))+ xlab("Site")+ ylab("abundance") + geom_col(aes(fill = Genus),position = position_stack(reverse = FALSE))+
  theme(axis.text.x = element_text(angle = 72,size = 15, hjust = 1, face = "italic", family = "Arial"))+
   scale_fill_manual(values = myPalette)+
  #guides(fill = guide_legend(reverse = FALSE))+
  guide_italics+
  theme(legend.text = element_text(size = 8, colour = "black", face = "italic", family = "Arial"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 8, family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 8, family = "Arial"))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.15, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.position = "right", legend.justification = "top", legend.direction = "vertical", legend.text = element_text(size = 10))+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 8, family = "Arial"))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 8, family = "Arial"))
```


```{r}
#install.packages('extrafont')
library(extrafont)
#font_import
#BiocManager::install("ggpubr")
```


```{r}
#change the misspelled names
ps@sam_data[ps@sam_data == 'liotrigona'] <- 'Liotrigona'
```

```{r}
#check if it worked
ps@sam_data@.Data[[2]]
```

```{r}
# check the number of taxa
ntaxa(ps)
```

```{r}
# Per county
grouped_data


Baringo<-grouped_data[,c(1,2:4,20:23)]
Bungoma<-grouped_data[,c(1,5:10)]
Bomet<-grouped_data[,c(1,11:16)]
Busia<-grouped_data[,c(1,17:19)]
Coast<-grouped_data[,c(1,24:26)]
Elgon<-grouped_data[,c(1,27:31)]
Eldama_Ravine<-grouped_data[,c(1,32:34)]
Kakamega<-grouped_data[,c(1,35:46,73:74,136:139)]
Nairobi<-grouped_data[,c(1,47:54,99:100,110:116,124:126)]
Isiolo<-grouped_data[,c(1,55:62)]
Mwingi<-grouped_data[,c(1,63:67,68:72,133:135)]
Kisii<-grouped_data[,c(1,75:88)]
Kitui<-grouped_data[,c(1,89:93)]
Nyamira<-grouped_data[,c(1,94)]
Maralal<-grouped_data[,c(1,95:98)]
Meru<-grouped_data[,c(1,101:106)]
Narok<-grouped_data[,c(1,107:109)]
Taita<-grouped_data[,c(1,117:123,140)]
Uasin_Gishu<-grouped_data[,c(1,127:132)]
Voic_Taita<-grouped_data[,c(1,140)]
all<-grouped_data[,c(1,2:140)]

Baringo_total <-Baringo%>% adorn_totals(c("col"))
Baringo_total <- mutate(Baringo_total, Baringo=rowSums(Baringo_total[9])/8)
Baringo_total <- Baringo_total[,c(1,10)]

Bungoma_total <-Bungoma%>% adorn_totals(c("col"))
Bungoma_total <- mutate(Bungoma_total, Bungoma=rowSums(Bungoma_total[8])/7)
Bungoma_total <- Bungoma_total[,c(1,9)]

Bomet_total <-Bomet%>% adorn_totals(c("col"))
Bomet_total <- mutate(Bomet_total, Bomet=rowSums(Bomet_total[8])/7)
Bomet_total <- Bomet_total[,c(1,9)]

Busia_total<-Busia%>% adorn_totals(c("col"))
Busia_total <- mutate(Busia_total, Busia=rowSums(Busia_total[5])/4)
Busia_total <- Busia_total[,c(1,6)]


Coast_total<-Coast%>% adorn_totals(c("col"))
Coast_total <- mutate(Coast_total, Coast=rowSums(Coast_total[5])/4)
Coast_total <- Coast_total[,c(1,6)]

Elgon_total<-Elgon%>% adorn_totals(c("col"))
Elgon_total <- mutate(Elgon_total, Elgon=rowSums(Elgon_total[7])/6)
Elgon_total <- Elgon_total[,c(1,8)]

Eldama_Ravine_total<-Eldama_Ravine%>% adorn_totals(c("col"))
Eldama_Ravine_total <- mutate(Eldama_Ravine_total, Eldama_Ravine=rowSums(Eldama_Ravine_total[5])/4)
Eldama_Ravine_total <- Eldama_Ravine_total[,c(1,6)]

Kakamega_total<-Kakamega%>% adorn_totals(c("col"))
Kakamega_total <- mutate(Kakamega_total, Kakamega=rowSums(Kakamega_total[19])/20)
Kakamega_total <- Kakamega_total[,c(1,21)]

Nairobi_total <-Nairobi%>% adorn_totals(c("col"))
Nairobi_total <- mutate(Nairobi_total, Nairobi=rowSums(Nairobi_total[22])/21)
Nairobi_total <- Nairobi_total[,c(1,23)]

Isiolo_total<-Isiolo%>% adorn_totals(c("col"))
Isiolo_total <- mutate(Isiolo_total, Isiolo=rowSums(Isiolo_total[10])/9)
Isiolo_total <- Isiolo_total[,c(1,11)]

Mwingi_total <- Mwingi%>%adorn_totals(c("col"))
Mwingi_total<- mutate(Mwingi_total, Mwingi=rowSums(Mwingi_total[15])/14)
Mwingi_total <- Mwingi_total[,c(1,16)]


Kisii_total<-Kisii%>% adorn_totals(c("col"))
Kisii_total <- mutate(Kisii_total, Kisii=rowSums(Kisii_total[16])/15)
Kisii_total <- Kisii_total[,c(1,17)]

Kitui_total<-Kitui%>% adorn_totals(c("col"))
Kitui_total <- mutate(Kitui_total, Kitui=rowSums(Kitui_total[7])/6)
Kitui_total <- Kitui_total[,c(1,8)]

Nyamira_total<-Nyamira%>% adorn_totals(c("col"))
Nyamira_total <- mutate(Nyamira_total, Nyamira=rowSums(Nyamira_total[3])/2)
Nyamira_total <- Nyamira_total[,c(1,4)]

Maralal_total<-Maralal%>% adorn_totals(c("col"))
Maralal_total <- mutate(Maralal_total, Maralal=rowSums(Maralal_total[6])/5)
Maralal_total <- Maralal_total[,c(1,7)]

Meru_total<-Meru%>% adorn_totals(c("col"))
Meru_total <- mutate(Meru_total, Meru=rowSums(Meru_total[8])/7)
Meru_total <- Meru_total[,c(1,9)]

Narok_total<-Narok%>% adorn_totals(c("col"))
Narok_total <- mutate(Narok_total, Narok=rowSums(Narok_total[5])/4)
Narok_total <- Narok_total[,c(1,6)]

Taita_total<-Taita%>% adorn_totals(c("col"))
Taita_total <- mutate(Taita_total, Taita=rowSums(Taita_total[10])/9)
Taita_total <- Taita_total[,c(1,11)]

Uasin_Gishu_total<-Uasin_Gishu%>% adorn_totals(c("col"))
Uasin_Gishu_total <- mutate(Uasin_Gishu_total, Uasin_Gishu=rowSums(Uasin_Gishu_total[8])/7)
Uasin_Gishu_total <- Uasin_Gishu_total[,c(1,9)]

all_total <- all %>% adorn_totals(c("col"))
all_total <- mutate(all_total, all=rowSums(all_total[141])/140)
all_total <- all_total[,c(1,142)]

merged1 <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE),
                 list(Baringo_total, Bungoma_total, Bomet_total, Busia_total,Coast_total, Elgon_total, Eldama_Ravine_total,Kakamega_total, Nairobi_total,Isiolo_total,Mwingi_total,Kisii_total,Kitui_total,Nyamira_total,Maralal_total,Meru_total,Narok_total,Taita_total,Uasin_Gishu_total,all_total))

cumulation <- merged1 %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100

genus_Rep <-c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus")

group <- aggregate(merged1[-1], list(Genus = replace(merged1$Genus,!(merged1$Genus %in% genus_Rep), "Others")), sum)

Baringo<-group[,c(1:2)]
Bungoma<- group[,c(1:3)] 
Bomet<- group[,c(1:4)]
Busia<- group[,c(1:5)]
Coast<- group[,c(1:6)]
Elgon<- group[,c(1:7)]
Eldama_Ravine<- group[,c(1:8)]
Kakamega<- group[,c(1:9)]
Nairobi<- group[,c(1:10)]
Isiolo<- group[,c(1:11)]
Mwingi<- group[,c(1:12)]
Kisii<- group[,c(1:13)]
Kitui<- group[,c(1:14)]
Nyamira<- group[,c(1:15)]
Maralal<- group[,c(1:16)]
Meru<- group[,c(1:17)]
Narok<- group[,c(1:18)]
Taita<- group[,c(1:19)]
Uasin_Gishu<- group[,c(1:20)]
All<- group[,c(1:20)]

bar_all <- adorn_percentages(All, denominator = "col", na.rm = T)
bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all<-bar_al# metagenomic analysis using dada2


#the R script
#loading of packages

```{r}
library("phyloseq")
library("vegan")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library("janitor")
library("coin")
library("reshape2")
library("ggnewscale")
library("MicrobiotaProcess")
library("patchwork")
library("dada2")
library("csv")
library("dplyr")
library("ggtree")
library("VennDiagram")
library("UpSetR")
#BiocManager::install("MicrobiotaProcess")
install.packages("MicrobiotaProcess")
library("MicrobiotaProcess")
library("utils")
library("utf8")
library("metagMisc")
```

```{r}
#BiocManager::install("devtools")
library("devtools")
library("tidyverse")
#devtools::install_github(repo= "jbisanz/qiime2R")
library(qiime2R)
asv_table <- read.table("/home/icipe/Downloads/dada2/ASV_table.tsv")
```

# loading the data

```{r}
metadata <- read.table("/home/icipe/Downloads/Geo_metadata.tsv", sep = "\t", row.names = 1, header= TRUE)

ot <- read.table("/home/icipe/Downloads/dada2/ASV_table.tsv", header = TRUE, sep = "\t",  row.names = 1 )
class(ot)
ot <- as.matrix(ot)
class(ot)

tax <- read.table("/home/icipe/Downloads/dada2/ASV_tax_species.tsv", header = TRUE, sep = "\t",  row.names = 1)
#tax <- data.frame(tax)
Tax <- as.matrix(tax)
#tax <- t(tax)

ot <- phyloseq::otu_table(ot, taxa_are_rows = TRUE)
tax <- phyloseq::tax_table(Tax)
metadata <- phyloseq::sample_data(metadata)
```

```{r}
ps <- phyloseq(ot, metadata, tax)
ps
```

```{r}
# we now create a third object called random for merging with the other three object
library("ape")
random_tree <- rtree(ntaxa(physeq), rooted = TRUE, tip.label = taxa_names(physeq))
plot(random_tree)
```

```{r}
#merging the preceeding 3 objects.
ps <- merge_phyloseq(physeq, random_tree)
ps
```

```{r}
# Removing the NAs
tx <- phyloseq_to_df(ps, addtax = T, addtot = F, addmaxrank = F)
tx <-tx[, c(8, 12:150 )]

tx <- tx %>%
  # recode empty strings "" by NAs
 na_if("") %>%
  # remove NAs
  na.omit()
tx
```

```{r}
cumulation <- tx %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation
```

```{r}
#merging the blast taxonomic classification to blast abundance table
merged_data <- tx
write.csv(merged_data, file="/home/icipe/Downloads/dada2/merged_data.csv")
# Group the data
grouped_data <- tx %>% group_by(Genus) %>%summarise_if(is.numeric, sum)
#grouped_data<-grouped_data[-c(1),]
```

```{r}
saveRDS(grouped_data,"/home/icipe/Downloads/dada2/grouped_data.rds")
write.table(tx,"/home/icipe/Downloads/dada2/tx.csv")
```


```{r}
# Grouping data by location

Baringo<-grouped_data[,c(1,2:4)]
Bungoma<-grouped_data[,c(1,5:10)]
Bomet<-grouped_data[,c(1,11:16)]
Busia<-grouped_data[,c(1,17:19)]
Baringo2<-grouped_data[,c(1,20:23)]
Coast<-grouped_data[,c(1,24:26)]
Elgon<-grouped_data[,c(1,27:31)]
Eldama_Ravine<-grouped_data[,c(1,32:34)]
Kakamega<-grouped_data[,c(1,35:46)]
ICIPE<-grouped_data[,c(1,47:54)]
Isiolo<-grouped_data[,c(1,55:62)]
Kasanga_Mwingi<-grouped_data[,c(1,63:67)]
Kasika_Mwingi<-grouped_data[,c(1,68:72)]
Kib_Kakamega<-grouped_data[,c(1,73:74)]
Kisii<-grouped_data[,c(1,75:88)]
Kitui<-grouped_data[,c(1,89:93)]
KNN_Nyamira<-grouped_data[,c(1,94)]
maralal<-grouped_data[,c(1,95:98)]
Maru_Nairobi<-grouped_data[,c(1,99:100)]
Meru<-grouped_data[,c(1,101:106)]
Narok<-grouped_data[,c(1,107:109)]
phi<-grouped_data[,c(1,110)]
saf_Nairobi<-grouped_data[,c(1,111:116)]
Taita<-grouped_data[,c(1,117:123)]
Thog_Nairobi<-grouped_data[,c(1,124:126)]
Uasin_Gishu<-grouped_data[,c(1,127:132)]
Ukasic_Mwingi<-grouped_data[,c(1,133:135)]
Vin_Kakamega<-grouped_data[,c(1,136:139)]
Voic_Taita <- grouped_data[1,140]
```

```{r}
all<-grouped_data[,c(1,2:140)]
all
```

```{r}
Baringo_total <-Baringo%>% adorn_totals(c("col"))
Baringo_total <- mutate(Baringo_total, Baringo=rowSums(Baringo_total[5])/4)
Baringo_total <- Baringo_total[,c(1,6)]

Bungoma_total <-Bungoma%>% adorn_totals(c("col"))
Bungoma_total <- mutate(Bungoma_total, Bungoma=rowSums(Bungoma_total[8])/7)
Bungoma_total <- Bungoma_total[,c(1,9)]

Bomet_total <-Bomet%>% adorn_totals(c("col"))
Bomet_total <- mutate(Bomet_total, Bomet=rowSums(Bomet_total[8])/7)
Bomet_total <- Bomet_total[,c(1,9)]

Busia_total<-Busia%>% adorn_totals(c("col"))
Busia_total <- mutate(Busia_total, Busia=rowSums(Busia_total[5])/4)
Busia_total <- Busia_total[,c(1,6)]

Baringo2_total<-Baringo2%>% adorn_totals(c("col"))
Baringo2_total <- mutate(Baringo2_total, Baringo2=rowSums(Baringo2_total[6])/5)
Baringo2_total <- Baringo2_total[,c(1,7)]

Coast_total<-Coast%>% adorn_totals(c("col"))
Coast_total <- mutate(Coast_total, Coast=rowSums(Coast_total[5])/4)
Coast_total <- Coast_total[,c(1,6)]

Elgon_total<-Elgon%>% adorn_totals(c("col"))
Elgon_total <- mutate(Elgon_total, Elgon=rowSums(Elgon_total[7])/6)
Elgon_total <- Elgon_total[,c(1,8)]

Eldama_Ravine_total<-Eldama_Ravine%>% adorn_totals(c("col"))
Eldama_Ravine_total <- mutate(Eldama_Ravine_total, Eldama_Ravine=rowSums(Eldama_Ravine_total[5])/4)
Eldama_Ravine_total <- Eldama_Ravine_total[,c(1,6)]

Kakamega_total<-Kakamega%>% adorn_totals(c("col"))
Kakamega_total <- mutate(Kakamega_total, Kakamega=rowSums(Kakamega_total[14])/13)
Kakamega_total <- Kakamega_total[,c(1,15)]

ICIPE_total <-ICIPE%>% adorn_totals(c("col"))
ICIPE_total <- mutate(ICIPE_total, ICIPE=rowSums(ICIPE_total[10])/9)
ICIPE_total <- ICIPE_total[,c(1,11)]

Isiolo_total<-Isiolo%>% adorn_totals(c("col"))
Isiolo_total <- mutate(Isiolo_total, Isiolo=rowSums(Isiolo_total[10])/9)
Isiolo_total <- Isiolo_total[,c(1,11)]

Kasanga_Mwingi_total <- Kasanga_Mwingi%>%adorn_totals(c("col"))
Kasanga_Mwingi_total<- mutate(Kasanga_Mwingi_total, Kasanga_Mwingi=rowSums(Kasanga_Mwingi_total[7])/6)
Kasanga_Mwingi_total <- Kasanga_Mwingi_total[,c(1,8)]

Kasika_Mwingi_total<-Kasika_Mwingi%>% adorn_totals(c("col"))
Kasika_Mwingi_total <- mutate(Kasika_Mwingi_total, Kasika_Mwingi=rowSums(Kasika_Mwingi_total[7])/6)
Kasika_Mwingi_total <- Kasika_Mwingi_total[,c(1,8)]

Kib_Kakamega_total<-Kib_Kakamega%>% adorn_totals(c("col"))
Kib_Kakamega_total <- mutate(Kib_Kakamega_total, Kib_Kakamega=rowSums(Kib_Kakamega_total[4])/3)
Kib_Kakamega_total <- Kib_Kakamega_total[,c(1,5)]

Kisii_total<-Kisii%>% adorn_totals(c("col"))
Kisii_total <- mutate(Kisii_total, Kisii=rowSums(Kisii_total[16])/15)
Kisii_total <- Kisii_total[,c(1,17)]

Kitui_total<-Kitui%>% adorn_totals(c("col"))
Kitui_total <- mutate(Kitui_total, Kitui=rowSums(Kitui_total[7])/6)
Kitui_total <- Kitui_total[,c(1,8)]

KNN_Nyamira_total<-KNN_Nyamira%>% adorn_totals(c("col"))
KNN_Nyamira_total <- mutate(KNN_Nyamira_total, KNN_Nyamira=rowSums(KNN_Nyamira_total[3])/2)
KNN_Nyamira_total <- KNN_Nyamira_total[,c(1,4)]

maralal_total <- maralal%>% adorn_totals(c("col"))
maralal_total <- mutate(maralal_total, maralal=rowSums(maralal_total[6])/5)
maralal_total <- maralal_total[,c(1,7)]

Maru_Nairobi_total <- Maru_Nairobi%>% adorn_totals(c("col"))
Maru_Nairobi_total <- mutate(Maru_Nairobi_total, Maru_Nairobi=rowSums(Maru_Nairobi_total[4])/3)
Maru_Nairobi_total <- Maru_Nairobi_total[,c(1,5)]

Meru_total<-Meru%>% adorn_totals(c("col"))
Meru_total <- mutate(Meru_total, Meru=rowSums(Meru_total[8])/7)
Meru_total <- Meru_total[,c(1,9)]

Narok_total<-Narok%>% adorn_totals(c("col"))
Narok_total <- mutate(Narok_total, Narok=rowSums(Narok_total[5])/4)
Narok_total <- Narok_total[,c(1,6)]

phi_total<-phi%>% adorn_totals(c("col"))
phi_total <- mutate(phi_total, phi=rowSums(phi_total[3])/2)
phi_total <- phi_total[,c(1,4)]

saf_Nairobi_total<-saf_Nairobi%>% adorn_totals(c("col"))
saf_Nairobi_total <- mutate(saf_Nairobi_total, saf_Nairobi=rowSums(saf_Nairobi_total[8])/7)
saf_Nairobi_total <- saf_Nairobi_total[,c(1,9)]

Taita_total<-Taita%>% adorn_totals(c("col"))
Taita_total <- mutate(Taita_total, Taita=rowSums(Taita_total[9])/8)
Taita_total <- Taita_total[,c(1,10)]

Thog_Nairobi_total<-Thog_Nairobi%>% adorn_totals(c("col"))
Thog_Nairobi_total <- mutate(Thog_Nairobi_total, Thog_Nairobi=rowSums(Thog_Nairobi_total[5])/4)
Thog_Nairobi_total <- Thog_Nairobi_total[,c(1,6)]

Uasin_Gishu_total<-Uasin_Gishu%>% adorn_totals(c("col"))
Uasin_Gishu_total <- mutate(Uasin_Gishu_total, Uasin_Gishu=rowSums(Uasin_Gishu_total[8])/7)
Uasin_Gishu_total <- Uasin_Gishu_total[,c(1,9)]

Ukasic_Mwingi_total<-Ukasic_Mwingi%>% adorn_totals(c("col"))
Ukasic_Mwingi_total <- mutate(Ukasic_Mwingi_total, Ukasic_Mwingi=rowSums(Ukasic_Mwingi_total[5])/4)
Ukasic_Mwingi_total <- Ukasic_Mwingi_total[,c(1,6)]

Vin_Kakamega_total<-Vin_Kakamega%>% adorn_totals(c("col"))
Vin_Kakamega_total <- mutate(Vin_Kakamega_total, Vin_Kakamega=rowSums(Vin_Kakamega_total[6])/5)
Vin_Kakamega_total <- Vin_Kakamega_total[,c(1,7)]

Voic_Taita_total <-Voic_Taita%>% adorn_totals(c("col"))
Voic_Taita_total <- mutate(Voic_Taita_total, Voic_Taita=rowSums(Voic_Taita_total[3])/2)
Voic_Taita_total <- Voic_Taita_total[,c(1,4)]
```

```{r}
all_total <- all %>% adorn_totals(c("col"))
all_total <- mutate(all_total, all=rowSums(all_total[141])/140)
all_total <- all_total[,c(1,142)]
```

```{r}
# Merging
merged <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE),
list(Baringo_total, Bungoma_total, Bomet_total, Busia_total, Baringo2_total, Coast_total, Elgon_total, Eldama_Ravine_total,Kakamega_total, ICIPE_total,Isiolo_total,Kasanga_Mwingi_total,Kasika_Mwingi_total,Kib_Kakamega_total,Kisii_total,Kitui_total,KNN_Nyamira,maralal_total,Maru_Nairobi_total,Meru_total,Narok_total,phi_total,saf_Nairobi_total,Taita_total,Thog_Nairobi_total, Uasin_Gishu_total,Ukasic_Mwingi_total,Vin_Kakamega_total,all_total))
```

```{r}
names(merged)<-c('Genus','Baringo','Bungoma','Bomet','Busia','Baringo2','Coast','Elgon','Eldama_Ravine','Kakamega','ICIPE','Isiolo','Kasanga_Mwingi','Kasika_Mwingi','Kib_Kakamega','Kisii','Kitui','KNN_Nyamira','maralal','Maru_Nairobi','Meru','Narok','phi','saf_Nairobi','Taita','Thog_Nairobi','Uasin_Gishu','Ukasic_Mwingi','Vin_Kakamega','All')
```

```{r}
#calculating the total abundance per genus and ordering from the most abundant to the lowest
cumulation <- merged %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100
```

```{r}
head(cumulation, n=20)
head(cumulation$Genus, n=51)
```

```{r}
genus_Rep <- c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus")
```

```{r}
group <- merged
```

```{r}
Baringo<-group[,c(1:2)]
Bungoma<- group[,c(1:3)] 
Bomet<- group[,c(1:4)]
Busia<- group[,c(1:5)]
Baringo2<- group[,c(1:6)]
Coast<- group[,c(1:7)]
Elgon<- group[,c(1:8)]
Eldama_Ravine<- group[,c(1:9)]
Kakamega<- group[,c(1:10)]
ICIPE<- group[,c(1:11)]
Isiolo<- group[,c(1:12)]
Kasanga_Mwingi<- group[,c(1:13)]
Kasika_Mwingi<- group[,c(1:14)]
Kib_Kakamega<- group[,c(1:15)]
Kisii<- group[,c(1:16)]
Kitui<- group[,c(1:17)]
KNN_Nyamira<- group[,c(1:18)]
maralal<- group[,c(1:19)]
Maru_Nairobi<- group[,c(1:20)]
Meru<- group[,c(1:21)]
Narok<- group[,c(1:22)]
phi<- group[,c(1:23)]
saf_Nairobi<- group[,c(1:24)]
Taita<- group[,c(1:25)]
Thog_Nairobi<- group[,c(1:26)]
Uasin_Gishu<- group[,c(1:27)]
Ukasic_Mwingi<- group[,c(1:28)]
Vin_Kakamega<- group[,c(1:29)]
All<- group[,c(1:30)]
```

```{r}
# Viewing Sample Diversity
#install.packages("janitor")
library(janitor)
```

```{r}
#converting the abudances into percentage
bar_all <- adorn_percentages(All, denominator = "col", na.rm = T)
bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all<-bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all
write.csv(dist_all, "/home/icipe/Downloads/dada2/Abundance.csv")
```
```{r}
#gathering the data
bar_all <- bar_all %>%
  gather(value = "abundance", key = "Site", -Genus)
bar_all <- as.data.frame(gsub("\\(", " (", as.matrix(bar_all)))
```

```{r}
# coerce the dataframe columns into respective data type
bar_all$Genus <- as.factor(bar_all$Genus)
bar_all$Site <- as.character(bar_all$Site)
bar_all$abundance <- as.numeric(bar_all$abundance)
```

```{r}
#ordering the data for plotting
bar_all$Genus <- reorder(bar_all$Genus, bar_all$abundance)
bar_all$Genus <- factor(bar_all$Genus, levels=rev(levels(bar_all$Genus)))
bar_all$Genus <- factor(bar_all$Genus, 
                        levels=c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus"))
```

```{r}
# Defining the color pallete
myPalette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#99000D", "#E6AB02", "#A6761D", "#666666","#FDCDAC", "#1F78B4", "#B2DF8A", "#33A02C", "#CBD5E8", "#E31A1C", "#FDBF6F", "#FF7F00","#4A1486","#C0C0C0","#B3E2CD","#FFFF33", "#5172b2","#F4CAE4", "#E6F5C9", "#FCFBFD","#139BF1","#09FF00","#065535", "#1D91C0", "#C0FFEE","#B35806","#0C2C84","#D0ED0E","#092617","#8A17FD","#4D5D53","#E48400","#6082B6","#316689","#EEA47FFF","#738678","#645452","#EEA47FFF", "#00539CFF", "#FC766AFF", "#42EADDFF", "#00A4CCFF", "#69B3BB", "#B589D6","#D1DFB7","#97BC62FF","#D198C5FF","#000000", "#CBCE91FF", "#616247FF", "#D64161FF","#435E55FF", "#DD4132FF","#CE4A7EFF", "#BD7F37FF","#FFA351FF","#185E57")
```

```{r}
# Definig the names in italics
guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 10, 
face = "italic", colour = "Black", angle = 0)))
```

# Plotting barplot

```{r}
#library(Cairo)
library(forcats)
 
p_all <- ggplot(bar_all,aes(x = fct_inorder(Site), y = abundance), labs(fill= Genus), group=row.names(bar_all))+ xlab("Site")+ ylab("abundance") + geom_col(aes(fill = Genus),position = position_stack(reverse = FALSE))+
  theme(axis.text.x = element_text(angle = 72,size = 15, hjust = 1, face = "italic", family = "Arial"))+
   scale_fill_manual(values = myPalette)+
  #guides(fill = guide_legend(reverse = FALSE))+
  guide_italics+
  theme(legend.text = element_text(size = 8, colour = "black", face = "italic", family = "Arial"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 8, family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 8, family = "Arial"))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.15, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.position = "right", legend.justification = "top", legend.direction = "vertical", legend.text = element_text(size = 10))+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 8, family = "Arial"))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 8, family = "Arial"))
```

```{r}
#install.packages('extrafont')
library(extrafont)
#font_import
#BiocManager::install("ggpubr")
```


```{r}
#change the misspelled names
ps@sam_data[ps@sam_data == 'liotrigona'] <- 'Liotrigona'
```

```{r}
#check if it worked
ps@sam_data@.Data[[2]]
```

```{r}
# check the number of taxa
ntaxa(ps)
```

```{r}
# group data Per county
grouped_data


Baringo<-grouped_data[,c(1,2:4,20:23)]
Bungoma<-grouped_data[,c(1,5:10)]
Bomet<-grouped_data[,c(1,11:16)]
Busia<-grouped_data[,c(1,17:19)]
Coast<-grouped_data[,c(1,24:26)]
Elgon<-grouped_data[,c(1,27:31)]
Eldama_Ravine<-grouped_data[,c(1,32:34)]
Kakamega<-grouped_data[,c(1,35:46,73:74,136:139)]
Nairobi<-grouped_data[,c(1,47:54,99:100,110:116,124:126)]
Isiolo<-grouped_data[,c(1,55:62)]
Mwingi<-grouped_data[,c(1,63:67,68:72,133:135)]
Kisii<-grouped_data[,c(1,75:88)]
Kitui<-grouped_data[,c(1,89:93)]
Nyamira<-grouped_data[,c(1,94)]
Maralal<-grouped_data[,c(1,95:98)]
Meru<-grouped_data[,c(1,101:106)]
Narok<-grouped_data[,c(1,107:109)]
Taita<-grouped_data[,c(1,117:123,140)]
Uasin_Gishu<-grouped_data[,c(1,127:132)]
Voic_Taita<-grouped_data[,c(1,140)]
all<-grouped_data[,c(1,2:140)]

Baringo_total <-Baringo%>% adorn_totals(c("col"))
Baringo_total <- mutate(Baringo_total, Baringo=rowSums(Baringo_total[9])/8)
Baringo_total <- Baringo_total[,c(1,10)]

Bungoma_total <-Bungoma%>% adorn_totals(c("col"))
Bungoma_total <- mutate(Bungoma_total, Bungoma=rowSums(Bungoma_total[8])/7)
Bungoma_total <- Bungoma_total[,c(1,9)]

Bomet_total <-Bomet%>% adorn_totals(c("col"))
Bomet_total <- mutate(Bomet_total, Bomet=rowSums(Bomet_total[8])/7)
Bomet_total <- Bomet_total[,c(1,9)]

Busia_total<-Busia%>% adorn_totals(c("col"))
Busia_total <- mutate(Busia_total, Busia=rowSums(Busia_total[5])/4)
Busia_total <- Busia_total[,c(1,6)]


Coast_total<-Coast%>% adorn_totals(c("col"))
Coast_total <- mutate(Coast_total, Coast=rowSums(Coast_total[5])/4)
Coast_total <- Coast_total[,c(1,6)]

Elgon_total<-Elgon%>% adorn_totals(c("col"))
Elgon_total <- mutate(Elgon_total, Elgon=rowSums(Elgon_total[7])/6)
Elgon_total <- Elgon_total[,c(1,8)]

Eldama_Ravine_total<-Eldama_Ravine%>% adorn_totals(c("col"))
Eldama_Ravine_total <- mutate(Eldama_Ravine_total, Eldama_Ravine=rowSums(Eldama_Ravine_total[5])/4)
Eldama_Ravine_total <- Eldama_Ravine_total[,c(1,6)]

Kakamega_total<-Kakamega%>% adorn_totals(c("col"))
Kakamega_total <- mutate(Kakamega_total, Kakamega=rowSums(Kakamega_total[19])/20)
Kakamega_total <- Kakamega_total[,c(1,21)]

Nairobi_total <-Nairobi%>% adorn_totals(c("col"))
Nairobi_total <- mutate(Nairobi_total, Nairobi=rowSums(Nairobi_total[22])/21)
Nairobi_total <- Nairobi_total[,c(1,23)]

Isiolo_total<-Isiolo%>% adorn_totals(c("col"))
Isiolo_total <- mutate(Isiolo_total, Isiolo=rowSums(Isiolo_total[10])/9)
Isiolo_total <- Isiolo_total[,c(1,11)]

Mwingi_total <- Mwingi%>%adorn_totals(c("col"))
Mwingi_total<- mutate(Mwingi_total, Mwingi=rowSums(Mwingi_total[15])/14)
Mwingi_total <- Mwingi_total[,c(1,16)]

Kisii_total<-Kisii%>% adorn_totals(c("col"))
Kisii_total <- mutate(Kisii_total, Kisii=rowSums(Kisii_total[16])/15)
Kisii_total <- Kisii_total[,c(1,17)]

Kitui_total<-Kitui%>% adorn_totals(c("col"))
Kitui_total <- mutate(Kitui_total, Kitui=rowSums(Kitui_total[7])/6)
Kitui_total <- Kitui_total[,c(1,8)]

Nyamira_total<-Nyamira%>% adorn_totals(c("col"))
Nyamira_total <- mutate(Nyamira_total, Nyamira=rowSums(Nyamira_total[3])/2)
Nyamira_total <- Nyamira_total[,c(1,4)]

Maralal_total<-Maralal%>% adorn_totals(c("col"))
Maralal_total <- mutate(Maralal_total, Maralal=rowSums(Maralal_total[6])/5)
Maralal_total <- Maralal_total[,c(1,7)]

Meru_total<-Meru%>% adorn_totals(c("col"))
Meru_total <- mutate(Meru_total, Meru=rowSums(Meru_total[8])/7)
Meru_total <- Meru_total[,c(1,9)]

Narok_total<-Narok%>% adorn_totals(c("col"))
Narok_total <- mutate(Narok_total, Narok=rowSums(Narok_total[5])/4)
Narok_total <- Narok_total[,c(1,6)]

Taita_total<-Taita%>% adorn_totals(c("col"))
Taita_total <- mutate(Taita_total, Taita=rowSums(Taita_total[10])/9)
Taita_total <- Taita_total[,c(1,11)]

Uasin_Gishu_total<-Uasin_Gishu%>% adorn_totals(c("col"))
Uasin_Gishu_total <- mutate(Uasin_Gishu_total, Uasin_Gishu=rowSums(Uasin_Gishu_total[8])/7)
Uasin_Gishu_total <- Uasin_Gishu_total[,c(1,9)]

all_total <- all %>% adorn_totals(c("col"))
all_total <- mutate(all_total, all=rowSums(all_total[141])/140)
all_total <- all_total[,c(1,142)]

merged1 <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE),
                 list(Baringo_total, Bungoma_total, Bomet_total, Busia_total,Coast_total, Elgon_total, Eldama_Ravine_total,Kakamega_total, Nairobi_total,Isiolo_total,Mwingi_total,Kisii_total,Kitui_total,Nyamira_total,Maralal_total,Meru_total,Narok_total,Taita_total,Uasin_Gishu_total,all_total))

cumulation <- merged1 %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100

genus_Rep <-c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus")

group <- aggregate(merged1[-1], list(Genus = replace(merged1$Genus,!(merged1$Genus %in% genus_Rep), "Others")), sum)

Baringo<-group[,c(1:2)]
Bungoma<- group[,c(1:3)] 
Bomet<- group[,c(1:4)]
Busia<- group[,c(1:5)]
Coast<- group[,c(1:6)]
Elgon<- group[,c(1:7)]
Eldama_Ravine<- group[,c(1:8)]
Kakamega<- group[,c(1:9)]
Nairobi<- group[,c(1:10)]
Isiolo<- group[,c(1:11)]
Mwingi<- group[,c(1:12)]
Kisii<- group[,c(1:13)]
Kitui<- group[,c(1:14)]
Nyamira<- group[,c(1:15)]
Maralal<- group[,c(1:16)]
Meru<- group[,c(1:17)]
Narok<- group[,c(1:18)]
Taita<- group[,c(1:19)]
Uasin_Gishu<- group[,c(1:20)]
All<- group[,c(1:20)]

bar_all <- adorn_percentages(All, denominator = "col", na.rm = T)
bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all<-bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all

bar_all <- bar_all %>%
  gather(value = "abundance", key = "Site", -Genus)
bar_all <- as.data.frame(gsub("\\(", " (", as.matrix(bar_all)))

bar_all$Genus <- as.factor(bar_all$Genus)
bar_all$Site <- as.character(bar_all$Site)
bar_all$abundance <- as.numeric(bar_all$abundance)

bar_all$Genus <- reorder(bar_all$Genus, bar_all$abundance)
bar_all$Genus <- factor(bar_all$Genus, levels=rev(levels(bar_all$Genus)))
bar_all$Genus <- factor(bar_all$Genus, 
                        levels=genus_Rep)

myPalette <- c("#D95F02", "#7570B3", "#E7298A", "#99000D", "#E6AB02", "#A6761D", "#666666","#FDCDAC", "#1F78B4", "#B2DF8A", "#33A02C", "#CBD5E8", "#E31A1C", "#FDBF6F", "#FF7F00","#4A1486","#C0C0C0","#B3E2CD","#FFFF33", "#5172b2","#F4CAE4", "#E6F5C9", "#FCFBFD","#139BF1","#09FF00","#065535", "#1D91C0", "#C0FFEE","#B35806","#0C2C84","#D0ED0E","#092617","#499976","#4D5D53","#E48400","#6082B6","#316689","#CEFB02","#738678","#645452","#EEA47FFF", "#00539CFF", "#FC766AFF", "#42EADDFF", "#00A4CCFF", "#69B3BB", "#B589D6","#D1DFB7","#97BC62FF","#D198C5FF","#000000")
length(myPalette)

guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 15, 
                                                                       face = "italic", colour = "Black", angle = 0)))

p_all <- ggplot(bar_all,aes(x = fct_inorder(Site), y = abundance), labs(fill= Genus), group=row.names(bar_all))+ xlab("County")+ ylab("abundance") + geom_col(aes(fill = Genus),position = position_stack(reverse = FALSE))+
  theme(axis.text.x = element_text(angle = 72,size = 15, hjust = 1, face = "italic", family = "Arial"))+
  scale_fill_manual(values = myPalette)+
  #guides(fill = guide_legend(reverse = FALSE))+
  guide_italics+
  theme(legend.text = element_text(size = 8, colour = "black", face = "italic", family = "Arial"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 8, family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 25, family = "Arial"))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.15, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.position = "right", legend.justification = "top", legend.direction = "vertical", legend.text = element_text(size = 10))+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 15, family = "Arial"))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 15, family = "Arial"))

plot4_60<-p_all + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + xlab("County")
plot4_60

ggsave("/home/icipe/Downloads/dada2/plot4_60.jpeg", width=12, height=12, dpi=600)
```

```{r}
merged <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE),
list(Baringo_total, Bungoma_total, Bomet_total, Busia_total,Coast_total, Elgon_total, Eldama_Ravine_total,Kakamega_total, Nairobi_total,Isiolo_total,Mwingi_total,Kisii_total,Kitui_total,Nyamira_total,Maralal_total,Meru_total,Narok_total,Taita_total,Uasin_Gishu_total,all_total))
```

```{r}
#calculating the total abundance per genus and ordering from the most abundant to the lowest
cumulation <- merged %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100

head(cumulation$Genus, n=51)
```

```{r}
genus_Rep <- c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus")
```

```{r}
group <- aggregate(merged[-1], list(Genus = replace(merged$Genus,!(merged$Genus %in% genus_Rep), "Others")), sum)

Baringo<-group[,c(1:2)]
Bungoma<- group[,c(1:3)] 
Bomet<- group[,c(1:4)]
Busia<- group[,c(1:5)]
Coast<- group[,c(1:6)]
Elgon<- group[,c(1:7)]
Eldama_Ravine<- group[,c(1:8)]
Kakamega<- group[,c(1:9)]
Nairobi<- group[,c(1:10)]
Isiolo<- group[,c(1:11)]
Mwingi<- group[,c(1:12)]
Kisii<- group[,c(1:13)]
Kitui<- group[,c(1:14)]
Nyamira<- group[,c(1:15)]
Maralal<- group[,c(1:16)]
Meru<- group[,c(1:17)]
Narok<- group[,c(1:18)]
Taita<- group[,c(1:19)]
Uasin_Gishu<- group[,c(1:20)]
All<- group[,c(1:20)]

```

```{r}
#converting the abudances into percentage
bar_all <- adorn_percentages(All, denominator = "col", na.rm = T)
bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all<-bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all
write.csv(dist_all, "/home/icipe/Downloads/dada2/Abundance.csv")

#gathering the data
bar_all <- bar_all %>%
  gather(value = "abundance", key = "Site", -Genus)
bar_all <- as.data.frame(gsub("\\(", " (", as.matrix(bar_all)))

# coerce the dataframe columns into respective data type
bar_all$Genus <- as.factor(bar_all$Genus)
bar_all$Site <- as.character(bar_all$Site)
bar_all$abundance <- as.numeric(bar_all$abundance)

#ordering the data for plotting
bar_all$Genus <- reorder(bar_all$Genus, bar_all$abundance)
bar_all$Genus <- factor(bar_all$Genus, levels=rev(levels(bar_all$Genus)))
bar_all$Genus <- factor(bar_all$Genus, 
                        

# Defining the color pallete
myPalette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#99000D", "#E6AB02", "#A6761D", "#666666","#FDCDAC", "#1F78B4", "#B2DF8A", "#33A02C", "#CBD5E8", "#E31A1C", "#FDBF6F", "#FF7F00","#4A1486","#C0C0C0","#B3E2CD","#FFFF33", "#5172b2","#F4CAE4", "#E6F5C9", "#FCFBFD","#139BF1","#09FF00","#065535", "#1D91C0", "#C0FFEE","#B35806","#0C2C84","#D0ED0E","#092617","#499976","#4D5D53","#E48400","#6082B6","#316689","#CEFB02","#738678","#645452","#EEA47FFF", "#00539CFF", "#FC766AFF", "#42EADDFF", "#00A4CCFF", "#69B3BB", "#B589D6","#D1DFB7","#97BC62FF","#D198C5FF","#000000", "#CBCE91FF", "#616247FF", "#D64161FF","#435E55FF", "#DD4132FF","#CE4A7EFF", "#BD7F37FF","#FFA351FF","#185E57")
length(myPalette)

# Definig the names in italics
guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 15, 
face = "italic", colour = "Black", angle = 0)))
```

```{r}
group <- aggregate(merged[-1], list(Genus = replace(merged$Genus,!(merged$Genus %in% genus_Rep), "Others")), sum)

Baringo<-group[,c(1:2)]
Bungoma<- group[,c(1:3)] 
Bomet<- group[,c(1:4)]
Busia<- group[,c(1:5)]
Coast<- group[,c(1:6)]
Elgon<- group[,c(1:7)]
Eldama_Ravine<- group[,c(1:8)]
Kakamega<- group[,c(1:9)]
Nairobi<- group[,c(1:10)]
Isiolo<- group[,c(1:11)]
Mwingi<- group[,c(1:12)]
Kisii<- group[,c(1:13)]
Kitui<- group[,c(1:14)]
Nyamira<- group[,c(1:15)]
Maralal<- group[,c(1:16)]
Meru<- group[,c(1:17)]
Narok<- group[,c(1:18)]
Taita<- group[,c(1:19)]
Uasin_Gishu<- group[,c(1:20)]
All<- group[,c(1:20)]
```

```{r}
merged <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE),
list(Baringo_total, Bungoma_total, Bomet_total, Busia_total,Coast_total, Elgon_total, Eldama_Ravine_total,Kakamega_total, Nairobi_total,Isiolo_total,Mwingi_total,Kisii_total,Kitui_total,Nyamira_total,Maralal_total,Meru_total,Narok_total,Taita_total,Uasin_Gishu_total,all_total))
```

```{r}
#converting the abudances into percentage
bar_all <- adorn_percentages(All, denominator = "col", na.rm = T)
bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all<-bar_all %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all
write.csv(dist_all, "/home/icipe/Downloads/dada2/Abundance.csv")
```

```{r}
#gathering the data
bar_all <- bar_all %>%
  gather(value = "abundance", key = "Site", -Genus)
bar_all <- as.data.frame(gsub("\\(", " (", as.matrix(bar_all)))

# coerce the dataframe columns into respective data type
bar_all$Genus <- as.factor(bar_all$Genus)
bar_all$Site <- as.character(bar_all$Site)
bar_all$abundance <- as.numeric(bar_all$abundance)
```

```{r}
#ordering the data for plotting
bar_all$Genus <- reorder(bar_all$Genus, bar_all$abundance)
bar_all$Genus <- factor(bar_all$Genus, levels=rev(levels(bar_all$Genus)))
bar_all$Genus <- factor(bar_all$Genus, 
                        levels=c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus"))
```

```{r}
# Defining the color pallete
myPalette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#99000D", "#E6AB02", "#A6761D", "#666666","#FDCDAC", "#1F78B4", "#B2DF8A", "#33A02C", "#CBD5E8", "#E31A1C", "#FDBF6F", "#FF7F00","#4A1486","#C0C0C0","#B3E2CD","#FFFF33", "#5172b2","#F4CAE4", "#E6F5C9", "#FCFBFD","#139BF1","#09FF00","#065535", "#1D91C0", "#C0FFEE","#B35806","#0C2C84","#D0ED0E","#092617","#8A17FD","#4D5D53","#E48400","#6082B6","#316689","#EEA47FFF","#738678","#645452","#EEA47FFF", "#00539CFF", "#FC766AFF", "#42EADDFF", "#00A4CCFF", "#69B3BB", "#B589D6","#D1DFB7","#97BC62FF","#D198C5FF","#000000", "#CBCE91FF", "#616247FF", "#D64161FF","#435E55FF", "#DD4132FF","#CE4A7EFF", "#BD7F37FF","#FFA351FF","#185E57")
```


```{r}
#plotting the barplot 
p_all <- ggplot(bar_all,aes(x = fct_inorder(Site), y = abundance), labs(fill= Genus), group=row.names(bar_all))+ xlab("Site")+ ylab("abundance") + geom_col(aes(fill = Genus),position = position_stack(reverse = FALSE))+
  theme(axis.text.x = element_text(angle = 72,size = 15, hjust = 1, face = "italic", family = "Arial"))+
   scale_fill_manual(values = myPalette)+
  #guides(fill = guide_legend(reverse = FALSE))+
  guide_italics+
  theme(legend.text = element_text(size = 8, colour = "black", face = "italic", family = "Arial"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 8, family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 25, family = "Arial"))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.15, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.position = "right", legend.justification = "top", legend.direction = "vertical", legend.text = element_text(size = 10))+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 15, family = "Arial"))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 15, family = "Arial"))
```


```{r}
# Definig the names in italics
guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 10, 
face = "italic", colour = "Black", angle = 0)))
```


```{r}
#install.packages('extrafont')
library(extrafont)
#font_import()

plot4_60<-p_all + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + xlab("Site")
plot4_60
```

```{r}
save.image("./Abundance.RData")

#alpha diversity

#Shannon

plotA = plot_richness(ps, x="site",color="genus", measures=("Shannon") ) +
  geom_boxplot(color = "blue")+
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -72))
plotA
```

```{r}
ggsave("/home/icipe/Downloads/dada2/Shanon_graph.jpeg", width=12, height = 12, dpi = 600 )
```

```{r}
#Simpson
plotB <-plot_richness(ps, x="site",color="genus", measures=c("Simpson") ) +
  geom_boxplot(color = "red")+
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))
plotB
```


```{r}
#Chao1
plotC <-plot_richness(ps, x="site",color="genus", measures=c("Chao1") ) +
geom_boxplot(color = "red")+
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))
plotC
```


```{r}
#ACE
plotD = plot_richness(ps, x="site",color="genus", measures=c("ACE") ) +
  geom_boxplot(color = "red") +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))
plotD
```

```{r}
#arrange the plots
mega_plot = ggarrange(plotA,plotB,plotC,plotD, ncol = 2, nrow = 2)
mega_plot
```


```{r}
#save the plots into different formats(for genus)
ggsave("alphaplot.jpeg", plot = mega_plot, width = 12, height = 12, dpi = 600)
ggsave("alphaplot.png", plot = mega_plot, width = 12 , height = 12, dpi = 600)
ggsave("alphaplot.svg", plot = mega_plot, width = 12, height = 12 , dpi = 600)
ggsave("alphaplot.tiff", plot = mega_plot, width = 12, height = 12, dpi = 600)
```

#Beta diversity
#PCA
```{r}
Site <- c('Baringo','Bungoma','Bomet','Busia','Baringo2','Coast','Elgon','Eldama_Ravine','Kakamega','ICIPE','Isiolo','Kasanga_Mwingi','Kasika_Mwingi','Kib_Kakamega','Kisii','Kitui','KNN_Nyamira','maralal','Maru_Nairobi','Meru','Narok','phi','saf_Nairobi','Taita','Thog_Nairobi','Uasin_Gishu','Ukasic_Mwingi','Vin_Kakamega')
```

```{r}
p1 <- ggplot(ps, aes(x=Site, y=genus_Rep)) +
  geom_boxplot(fill=c("blue","red")) +
  labs(title= 'Richness', x= ' ', y= '', tag = "A") +
  geom_point()
p1
```

```{r}
#Plot the PCA
pcaplot1 = ggordpoint(obj = pca, biplot = T, speciesannot = T,
                      factorNames = ("genus"), ellipse = T) +
  scale_color_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#FFA351FF")) +
  scale_fill_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA","#FFA351FF"))
pcaplot1
```

```{r}
#PCOA
pcoa = get_pcoa(obj = ps, distmethod = "euclidean", method = "hellinger")
pcoa
```


```{r}
#Plot the PCOA
pcoaplot1 = ggordpoint(obj = pcoa, biplot = T, speciesannot = F,showsample = T,
                        factorNames=c("genus"), ellipse = T) +
  scale_color_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#FFA351FF")) +
  scale_fill_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#FFA351FF"))
pcoaplot1
```


```{r}
#save the plots## Alpha Diversity

ps3
set.seed(1024)
rareres <- get_rarecurve(obj=ps3, chunks=400)

p_rare <- ggrarecurve(obj=rareres,
                      indexNames=c("Observe","Chao1","ACE","Shannon"),
                      ) +
          theme(legend.spacing.y=unit(0.01,"cm"),
                legend.text=element_text(size=4))

prare1 <- ggrarecurve(obj=rareres, factorNames="Animal_Type",
                      indexNames=c("Observe","Chao1","ACE","Shannon")
                      ) +
          scale_fill_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA"))+
          scale_color_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA"))+
          theme_bw()+
          theme(axis.text=element_text(size=8), panel.grid=element_blank(),
                strip.background = element_rect(colour=NA,fill="grey"),
                strip.text.x = element_text(face="bold"))          

prare2 <- ggrarecurve(obj=rareres,
                      factorNames="Animal_Type",
                      shadow=FALSE,
                      indexNames=c("Observe", "Chao1", "ACE","Shannon")
                      ) +
          scale_color_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA"))+
          theme_bw()+
          theme(axis.text=element_text(size=8), panel.grid=element_blank(),
                strip.background = element_rect(colour=NA,fill="grey"),
                strip.text.x = element_text(face="bold"))
p_rare / prare1 / prare2
```

```{r}
ggsave("betaplot_PCA.jpeg", plot = pcaplot1, width = 12, height = 12, dpi = 600)
ggsave("betaplot_PCA.png", plot = pcaplot1, width = 12, height = 12, dpi = 600)
ggsave("betaplot_pcoa.png", plot = pcoaplot1, width = 12, height = 12, dpi = 600)
ggsave("betaplot_pcoa.jpeg", plot = pcoaplot1, width = 12, height = 12, dpi = 600)
```


l %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()
dist_all

bar_all <- bar_all %>%
  gather(value = "abundance", key = "Site", -Genus)
bar_all <- as.data.frame(gsub("\\(", " (", as.matrix(bar_all)))

bar_all$Genus <- as.factor(bar_all$Genus)
bar_all$Site <- as.character(bar_all$Site)
bar_all$abundance <- as.numeric(bar_all$abundance)

bar_all$Genus <- reorder(bar_all$Genus, bar_all$abundance)
bar_all$Genus <- factor(bar_all$Genus, levels=rev(levels(bar_all$Genus)))
bar_all$Genus <- factor(bar_all$Genus, 
                        levels=genus_Rep)

myPalette <- c("#D95F02", "#7570B3", "#E7298A", "#99000D", "#E6AB02", "#A6761D", "#666666","#FDCDAC", "#1F78B4", "#B2DF8A", "#33A02C", "#CBD5E8", "#E31A1C", "#FDBF6F", "#FF7F00","#4A1486","#C0C0C0","#B3E2CD","#FFFF33", "#5172b2","#F4CAE4", "#E6F5C9", "#FCFBFD","#139BF1","#09FF00","#065535", "#1D91C0", "#C0FFEE","#B35806","#0C2C84","#D0ED0E","#092617","#499976","#4D5D53","#E48400","#6082B6","#316689","#CEFB02","#738678","#645452","#EEA47FFF", "#00539CFF", "#FC766AFF", "#42EADDFF", "#00A4CCFF", "#69B3BB", "#B589D6","#D1DFB7","#97BC62FF","#D198C5FF","#000000")
length(myPalette)

guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 15, 
                                                                       face = "italic", colour = "Black", angle = 0)))

p_all <- ggplot(bar_all,aes(x = fct_inorder(Site), y = abundance), labs(fill= Genus), group=row.names(bar_all))+ xlab("County")+ ylab("abundance") + geom_col(aes(fill = Genus),position = position_stack(reverse = FALSE))+
  theme(axis.text.x = element_text(angle = 72,size = 15, hjust = 1, face = "italic", family = "Arial"))+
  scale_fill_manual(values = myPalette)+
  #guides(fill = guide_legend(reverse = FALSE))+
  guide_italics+
  theme(legend.text = element_text(size = 8, colour = "black", face = "italic", family = "Arial"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 8, family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 25, family = "Arial"))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.15, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.position = "right", legend.justification = "top", legend.direction = "vertical", legend.text = element_text(size = 10))+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 15, family = "Arial"))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 15, family = "Arial"))

plot4_60<-p_all + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + xlab("County")
plot4_60

ggsave("/home/icipe/Downloads/dada2/plot4_60.jpeg", width=12, height=12, dpi=600)
```

```{r}
merged <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE),
list(Baringo_total, Bungoma_total, Bomet_total, Busia_total,Coast_total, Elgon_total, Eldama_Ravine_total,Kakamega_total, Nairobi_total,Isiolo_total,Mwingi_total,Kisii_total,Kitui_total,Nyamira_total,Maralal_total,Meru_total,Narok_total,Taita_total,Uasin_Gishu_total,all_total))
```

```{r}
#calculating the total abundance per genus and ordering from the most abundant to the lowest
cumulation <- merged %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100

head(cumulation$Genus, n=51)
```

```{r}
genus_Rep <- c("Zygosaccharomyces","Starmerella","Cercospora","Cladosporium","Unguiculariopsis","Epicoccum","Alternaria","Metschnikowia","Trichosporonoides","Curvularia","Acidea","Wallemia","Aureobasidium","Hannaella","Papiliotrema","Coniothyrium","Hormonema","Melampsora","Fusarium","Wickerhamiella","Wickerhamomyces","Periconia","Hanseniaspora","Erysiphe","Aspergillus","Antealophiotrema","Choanephora","Exserohilum","Phaeosphaeria","Quambalaria","Tilletia","Yamadazyma","Pycnopeziza","Nomuraea","Bipolaris","Neosetophoma","Paraphaeosphaeria","Podosphaera","Penicillium","Pyrenochaeta","Cryptococcus","Phlebiopsis","Amphobotrys","Nigrospora","Parapyrenochaeta","Leptospora","Setophoma","Symmetrospora","Saccharomyces","Neopestalotiopsis","Gymnopilus")
```

```{r}
group <- aggregate(merged[-1], list(Genus = replace(merged$Genus,!(merged$Genus %in% genus_Rep), "Others")), sum)

Baringo<-group[,c(1:2)]
Bungoma<- group[,c(1:3)] 
Bomet<- group[,c(1:4)]
Busia<- group[,c(1:5)]
Coast<- group[,c(1:6)]
Elgon<- group[,c(1:7)]
Eldama_Ravine<- group[,c(1:8)]
Kakamega<- group[,c(1:9)]
Nairobi<- group[,c(1:10)]
Isiolo<- group[,c(1:11)]
Mwingi<- group[,c(1:12)]
Kisii<- group[,c(1:13)]
Kitui<- group[,c(1:14)]
Nyamira<- group[,c(1:15)]
Maralal<- group[,c(1:16)]
Meru<- group[,c(1:17)]
Narok<- group[,c(1:18)]
Taita<- group[,c(1:19)]
Uasin_Gishu<- group[,c(1:20)]
All<- group[,c(1:20)]
```

```{r}
#save the plots into different formats(for genus)
ggsave("alphaplot.jpeg", plot = mega_plot, width = 12, height = 12, dpi = 600)
ggsave("alphaplot.png", plot = mega_plot, width = 12 , height = 12, dpi = 600)
ggsave("alphaplot.svg", plot = mega_plot, width = 12, height = 12 , dpi = 600)
ggsave("alphaplot.tiff", plot = mega_plot, width = 12, height = 12, dpi = 600)
```
