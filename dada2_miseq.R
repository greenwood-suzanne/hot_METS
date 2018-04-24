library(dada2)
library(optparse)


option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Path to working directory folder", metavar="character"),
  make_option(c("-F", "--For"), type="character", default=320, 
              help="Forward read length", metavar="character"),
  make_option(c("-R", "--Rev"), type="character", default=220, 
              help="Reverse read length", metavar="character"),
  make_option(c("-T", "--TrimL"), type="character", default=15, 
              help="Primer nucleotides trimmed", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

now <- Sys.time()
now

####Process the rest of the arguments####
print("Parameter list:")
sprintf("Forward read length: %s", opt$For)
sprintf("Reverse read length: %s", opt$Rev)
sprintf("Primer nucleotides trimmed: %s", opt$trimL)

#the arguments are:
  #1: working directory path; WHOLE PATH
  #2: desired forward read length
  #3: desired reverse read length
  #4: how many nucls to trim from the L. ATM that is the same nmber on fwd and rev

####Here we grab our files and set up lists of forward and reverse reads####
#path = "/homes/sgreenwood1/MiSeq_SOP"
path <- opt$file
setwd(path)
list.files(path)
fnFs = sort(list.files(path, pattern= "R1_001.fastq", full.names = TRUE))
fnRs = sort(list.files(path, pattern= "R2_001.fastq", full.names = TRUE))
sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1)

###Next we take a look at raw quality###
png(filename = "Output/Fqual.png")
plotQualityProfile(fnFs)
dev.off()

png(filename = "Output/Rqual.png")
plotQualityProfile(fnRs)
dev.off()
 
####Assign the filenames for filtered fastq.gz files####
filt_path <- file.path(path, "filtered") 
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
#filtered forwards
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
#filtered reverses

print("Filtered files will appear in 'Filtered' subdirectory.")

print("Now filtering reads...")

####Now we filter####
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(as.integer(opt$For),as.integer(opt$Rev)),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, trimLeft=c(as.integer(opt$TrimL),as.integer(opt$TrimL)), matchIDs=TRUE)
#out is a table displaying the file name, number of reads to start with
#and number of reads after truncating at the given positions (240 on forward and 160 on reverse to keep scores above 20)
print("Filtered. Here's the first few files' reads: ")
head(out)
#prints the first few rows showing the original read counts next to the filtered counts

#now we look at filtered quality
print("Now generating filtered read quality plots...")
png(filename = "Output/filtFqual.png")
plotQualityProfile(filtFs) 
dev.off()
print("Filtered forward read quality profile now available in 'Output' subdirectory.")

png(filename = "Output/filtRqual.png")
plotQualityProfile(filtRs) 
dev.off()
print("Filtered reverse read quality profile now available in 'Output' subdirectory.")

####Examine errors####
errF <- learnErrors(filtFs, multithread=TRUE)
#tells us how amny unique reads we have in each sample
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
#oints are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence.
#The red line shows the error rates expected under the nominal definition of the Q-value. 
#there are reverse steps too
plotErrors(errF, nominalQ=TRUE)

####Dereplication Process####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
#Dereplication combines all identical sequencing reads into into
#“unique sequences” with a corresponding “abundance”: the number of 
#reads with that unique sequence

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
#gives us a table of all the unique sequences from all the forwards and how many times they occur
#when we go back, there will be reverse steps too

####denoising####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
#get sequence variants from each sample
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#do the same for the reverse set
dadaFs[[1]]


####Merge Paired Reads####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate= TRUE, verbose=TRUE)
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
print("This is the number of reads after each step  so far:")
head(track)
#show us how many reads we have for each file as we filter down

####Assign Taxonomy####
setwd(path)
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=TRUE)
#classify sequence variants taxonomically
taxa <- addSpecies(taxa, "silva_species_assignment_v132.fa.gz")
#add species level assignments by exact matching
#error encountered


taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.table("Output/taxonomy.txt", sep="\t") #save the taxonomy table to output folder
#display the taxon info

####Evaluate Accuracy####
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

#lets investigate the mock community 

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")


####PHYLOSEQ PORTION####
library(phyloseq) 
packageVersion("phyloseq")
library(ggplot2) 
packageVersion("ggplot2")

#import tble into phyloseq
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

#phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps

#alpha diversity
print("Plotting alpha diversity...")
png(filename= "Output/alphadiversity.png")                        
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When") + theme_bw()
dev.off()
                        
ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
png(filename= "Output/ordination.png")                        
plot_ordination(ps, ord.nmds.bray, color="When", title="Bray NMDS")
dev.off()                        

#taxonomy distribution plot:
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
print("Plotting family level taxonomy: ")     
png(filename= "Output/family.png")                                    
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")
dev.off()                                    
                                    
print("Plotting species level taxonomy: ") 
png(filename= "Output/species.png")                                    
plot_bar(ps.top20, x="Day", fill="Species") + facet_wrap(~When, scales="free_x")   
dev.off()
                                    
print("All done!")                                  
now <- Sys.time()
now
