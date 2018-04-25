library(dada2); packageVersion("dada2")

now <- Sys.time()
now

args = commandArgs(trailingOnly = TRUE)
path <- args[1]
setwd(path)

#path <- "~/raw" # CHANGE ME to the directory containing the fastq files after unzipping.
#path <- "~/crosshack/crossteam" 
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#png(filename = "fqual.png")
#plotQualityProfile(fnFs[1:2])
#dev.off()

#png(filename = "rqual.png")
#plotQualityProfile(fnRs[1:2])
#dev.off()

filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(320,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, trimLeft=c(15,15), matchIDs=TRUE) 
head(out)

now <- Sys.time()
now

#png(filename = "filtFqual.png")
#plotQualityProfile(filtFs[1:2])
#dev.off()

#png(filename = "filtRqual.png")
#plotQualityProfile(filtRs[1:2])
#dev.off()

errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

now <- Sys.time()
now

#plotErrors(errF, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=FALSE)

dadaRs <- dada(derepRs, err=errR, multithread=FALSE)

dadaFs[[1]]
now <- Sys.time()
now 

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

now <- Sys.time()
now

seqtab <- makeSequenceTable(mergers)

dim(seqtab)

table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, "~/raw/silva_nr_v132_train_set.fa.gz", multithread=FALSE, tryRC=TRUE)

taxa <- addSpecies(taxa, "~/raw/silva_species_assignment_v132.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

now <- Sys.time()
now

#phyloseq
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

now <- Sys.time()
now

#gets the location from the sample string
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "_S"), `[`, 1)
location <- substr(subject,21,26)
location <- gsub('-', '',location)
location <- gsub('O', 'H20',location)
location <- gsub('1USA', 'USA', location)


#gets the sample number from the sample string
sample_num <-substr(subject, 14, 16)
sample_num <- gsub('-', '',sample_num)
samdf <- data.frame(Location = location, SampleName = sample_num)
rownames(samdf) <- samples.out

#a phylo seq object is created
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps

png(filename = "alpha_diversity.png", width = 2000, height = 1200)
plot_richness(ps, x= "SampleName", measures=c("Shannon", "Simpson"), color="Location") + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))
dev.off()

#ords.nmds.bray <- ordinate(ps, method = "NMDS", distance = "bray")

#png(filename = "ordination.png")
#plot_ordination(ps, ord.nmds.bray, color="Location", title="Bray NMDS")
#dev.off()

now <- Sys.time()
now

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

#creates a bar plot of the family
png(filename = "familybarplot.png", width = 2000, height = 1200)
plot_bar(ps.top20, x = "SampleName", fill="Family") + facet_wrap(~Location, scales="free_x") + theme(axis.text.x=element_text(angle=90,hjust=1))
dev.off()

#creates a bar plot of the species
#png(filename = "speciesbarplot.png",width = 2000, height = 1200)
#plot_bar(ps.top20, x= "SampleName", fill="Species") + facet_wrap(~Location, scales="free_x") + theme(axis.text.x=element_text(angle=90,hjust=1))
#dev.off()

now <- Sys.time()
now
