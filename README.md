
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

 
