library(dada2)
packageVersion("dada2") 

args = commandArgs(trailingOnly=TRUE)
path<-args[1]
#the user should put their working directory in the command to run the script
setwd(path)

####Here we grab our files and set up lists of forward and reverse reads####
#path = "/homes/sgreenwood1/crossteam"
list.files(path)
#hardcoded path will be removed but for simplicity's sake we have it here for this activity
fnFs = sort(list.files(path, pattern= "R1_001.fastq", full.names = TRUE))
fnRs = sort(list.files(path, pattern= "R2_001.fastq", full.names = TRUE))
sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1)

###Next we take a look at quality###
png(filename = "fqual.png")
plotQualityProfile(fnFs[1:2]) #some trash reads but most are alright
#this is looking only at two files
#we need to do some quality trimming here but it isnt 
#done in the tutorial so we'll have to figure out how
dev.off()
png(filename = "rqual.png")
plotQualityProfile(fnRs[1:2])
#the reverse read quality isnt great. more than half the reads have phred scores below 25/20
#again, only looking at first two
dev.off()

####Assign the filenames for filtered fastq files####
filt_path <- file.path(path, "filtered") 
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
#filtered forwards
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))
#filtered reverses

####Now we filter####
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(320,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, trimLeft=c(15,15), matchIDs=TRUE) 
#out is a table displaying the file name, number of reads to start with
#and number of reads after truncating at the given positions
head(out)
                     
#prints the first few rows showing the original read counts next to the filtered counts
png(filename = "filtFqual.png")
plotQualityProfile(filtFs[1:2]) 
dev.off()

png(filename = "filtRqual.png")
plotQualityProfile(filtRs[1:2]) 
dev.off()

####Examine errors####
errF <- learnErrors(filtFs, multithread=TRUE)
#tells us how amny unique reads we have in each sample
errR <- learnErrors(filtRs, multithread=TRUE)
#this is needed for the dada function

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


####Remove Chimeras####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) 
sum(seqtab.nochim)/sum(seqtab)



####Track Reads Through Pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nochim")
rownames(track) <- sample.names
head(track)
#show us how many reads we have for each file as we filter down

####Assign Taxonomy####
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=FALSE)
#classify sequence variants taxonomically
taxa <- addSpecies(taxa, "silva_species_assignment_v132.fa.gz")
#add species level assignments by exact matching

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
#display the taxon info





