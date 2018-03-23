library(dada2)
packageVersion("dada2") 

  ####Here we grab our files and set up lists of forward and reverse reads####
path = "/homes/sgreenwood1/raw"
list.files(path)
fnFs = sort(list.files(path, pattern= "R1_001.fastq", full.names = TRUE))
fnRs = sort(list.files(path, pattern= "R2_001.fastq", full.names = TRUE))
sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1)

  ###Next we take a look at quality###
plotQualityProfile(fnFs[1:2]) #some trash reads but most are alright
#this is looking only at two files
#we need to do some quality trimming here but it isnt 
#done in the tutorial so we'll have to figure out how

plotQualityProfile(fnRs[1:2])
#the reverse read quality isnt great. more than half the reads have phred scores below 25/20
#again, only looking at first two

  ####Assign the filenames for filtered fastq files####
filt_path <- file.path(path, "filtered") 
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
#filtered forwards
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))
#filtered reverses

  ####Now we filter####
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncQ= 2, truncLen = c(240, 160), maxN=0, maxEE=c(2,2), rm.phix=TRUE,compress=TRUE, multithread=TRUE, trimLeft=c(8,9), matchIDs = TRUE) 
#out is a table displaying the file name, number of reads to start with
#and number of reads after truncating at the given positions (210 on forward and 120 on reverse to keep scores above 20)
head(out)
#prints the first few rows showing the original read counts next to the filtered counts


  ####Examine errors####
errF <- learnErrors(filtFs, multithread=TRUE)
#tells us how amny unique reads we have in each sample
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
#oints are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence.
#The red line shows the error rates expected under the nominal definition of the Q-value. 
plotErrors(errR, nominalQ=TRUE)

  ####Dereplication Process####
derepFs <- derepFastq(filtFs, verbose=TRUE)
# Name the derep-class objects by the sample names
#Dereplication combines all identical sequencing reads into into
#“unique sequences” with a corresponding “abundance”: the number of 
#reads with that unique sequence
names(derepFs) <- sample.names
#gives us a table of all the unique sequences from all the forwards and how many times they occur
#when we go back, there will be reverse steps too
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepRs) <- sample.names

####denoising####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
#get sequence variants from each sample
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#do the same for the reverse set
dadaFs[[1]]


####Merge Paired Reads####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
#had to add justConcatenate to get this working
head(mergers[[1]])


####Construct a Sequence Table####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#matrix rows are the samples and columns are the variants so 
#were only getting sample #170 wiht 2518 variants

####Remove Chimeras####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) #slightly high
sum(seqtab.nochim)/sum(seqtab) #slightly low



####Track Reads Through Pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
#show us how many reads we have for each file as we filter down

####Assign Taxonomy####
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=TRUE)
#classify sequence variants taxonomically
taxa <- addSpecies(taxa, "silva_species_assignment_v132.fa.gz")
#add species level assignments by exact matching

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
#display the taxon info


####Evaluate Accuracy####
unqs.mock <- seqtab.nochim["GSF856-Dugas-101-QH2O",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

#lets investigate the mock community to set up something to compare to
#4 sample sequences found in mock community

mock.ref <- getSequences(file.path(path, "raw/GSF856-Dugas-101-QH2O_S159_L001_R1_001.fastq"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
#0 matches 

####PHYLOSEQ PORTION####
library(phyloseq) 
packageVersion("phyloseq")
library(ggplot2) 
packageVersion("ggplot2")

#import tble into phyloseq
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "_L001"), `[`, 1)
origin <- substr(subject, 18, length(subject))
origin <- substr(origin, 4,10)
print(origin)
samdf <- data.frame(Subject=subject, Where=origin)
samdf$Where <- origin
#samdf$Where[samdf$Where!="USA"] <- "Ghana"
rownames(samdf) <- samples.out

#phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samdf),tax_table(taxa))
ps

#taxonomy distribution plot:
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="origin", fill="Family") + facet_wrap(~Where, scales="free_x")

