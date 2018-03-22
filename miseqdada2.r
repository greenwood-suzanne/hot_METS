---
title: "Untitled"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r}
path <- "~/miseq/MiSeq_SOP"
list.files(path)
#sets the path to the folder containing the miseq_sop files
```


```{r}
# formats the forward reads to SAMPLENAME_R1_001.fastq and the reverse reads to SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# gets the file names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```


```{r}
#plots the first two reads forward reads to see if they need trimming
plotQualityProfile(fnFs[1:2])
```
```{r}
#plots the first two reads reverse reads to see if they need trimming
plotQualityProfile(fnRs[1:2])
```

```{r}
#creates a new subdirectory called filteres
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
```

```{r}
#filter and trims the seqeunces, the forward are cut after 245 bps and the reverse are cut after 200 bps
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(245,200),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) 
head(out)
```

```{r}
#estimates the error rates of the forward reads
errF <- learnErrors(filtFs, multithread=TRUE)
```

```{r}
#estimates the error rates of the reverse reads
errR <- learnErrors(filtRs, multithread=TRUE)
```

```{r}
#plots the error rates
plotErrors(errF, nominalQ=TRUE)
```

```{r}
#combines identical sequence reads into unique sequences
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

```{r}
#infers the sequence variants in each forward sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
```

```{r}
#infers the sequence variants in each reverse sample
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

```{r}
#inspects the class object returned by dada
dadaFs[[1]]
```

```{r}
#merges the forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
```

```{r}
head(mergers[[1]])
```
```{r}
#constructs a sequenec table
seqtab <- makeSequenceTable(mergers)
```

```{r}
dim(seqtab)
```

```{r}
#view the distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```

```{r}
#removes chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
```

```{r}
#calculates the fraction of chimeric sequnces
sum(seqtab.nochim)/sum(seqtab)
```

```{r}
#view the number of reads that were maintained throughout the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
```

```{r}
#assigns the taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/homes/cprintzis/miseq/MiSeq_SOP/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
```

```{r}
#asigns the species level taxonomy
taxa <- addSpecies(taxa, "/homes/cprintzis/miseq/MiSeq_SOP/silva_species_assignment_v132.fa.gz")
#doesnt work since just concatenate leaves N's in the sequence
```

```{r}
#view the taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```
```{r}
#evaluate dada2's accuracy
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
```

```{r}
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```

```{r}
#a data.frame is created to read the data
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
```

```{r}
#a phylo seq object is created
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps
```

```{r}
#visualizes the alpha-diversity
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When") + theme_bw()
```

```{r}
#ordination picks out distinctions between early and late samples
ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
```
```{r}
plot_ordination(ps, ord.nmds.bray, color="When", title="Bray NMDS")
```

```{r}
#plots the taxonomic distribution
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")
```

