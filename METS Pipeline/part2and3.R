#Pipeline Part II: This part begins with the filtering and trimming of raw data and runs through plotting.
library(optparse)
library(dada2)

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

####Here we have to again define fwd and reverse reads####
path <- opt$file
#the user should put their working directory in the command to run the script
setwd(path)

####Process the rest of the arguments####
#print("Parameter list:")
#sprintf("Forward read length: %s", opt$For)
#sprintf("Reverse read length: %s", opt$Rev)
#sprintf("Primer nucleotides trimmed: %s", opt$trimL)

#the arguments are:
#1: working directory path; WHOLE PATH
#2: desired forward read length
#3: desired reverse read length
#4: how many nucls to trim from the L. ATM that is the same nmber on fwd and rev

#raw files
fnFs = sort(list.files(path, pattern= "R1_001.fastq", full.names = TRUE))
fnRs = sort(list.files(path, pattern= "R2_001.fastq", full.names = TRUE))
sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1)

####Assign the filenames for filtered fastq files####
filt_path <- file.path(path, "filtered") 
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
#filtered forwards
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))
#filtered reverses

print("Filtered files will appear in 'Filtered' subdirectory.")

print("Now filtering reads...")


####Now we filter####
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(as.integer(opt$For),as.integer(opt$Rev)),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, trimLeft=c(as.integer(opt$TrimL),as.integer(opt$TrimL)), matchIDs=TRUE)  
#out is a table displaying the file name, number of reads to start with
#and number of reads after truncating at the given positions

print("Finished filtering. Here's what the filtered files look like:  ")

head(out)
#prints the first few rows showing the original read counts next to the filtered counts
#write.table(out, "Output/filteredandtrimmed.txt", sep="\t")
#out table saved to txt file. little ugly but does the trick

####Output filtered read quality####
print("Now generating filtered read quality plots...")
png(filename = "Output/filtFqual.png")
plotQualityProfile(filtFs) 
dev.off()
print("Filtered forward read quality profile now available in 'Output' subdirectory.")

png(filename = "Output/filtRqual.png")
plotQualityProfile(filtRs) 
dev.off()
print("Filtered reverse read quality profile now available in 'Output' subdirectory.")

print("If you are happy with the filtering so far, move on to Part III. Part III will work with the unique sequences to produce taxonomy table and comparison plots.")


#Pipeline Part III: This is the part of the pipeline that takes filtered reads, finds the unique
#sequences, merges paired end reads, assigns taxonomy, creates the taxonomy table and phyloseq
#plots. 

now <- Sys.time()
now


####Examine errors####
print("Examining errors...")
errF <- learnErrors(filtFs, multithread=FALSE)
#tells us how amny unique reads we have in each sample
errR <- learnErrors(filtRs, multithread=FALSE)
#this is needed for the dada function

setwd(path)
png(filename= "Output/errF.png")
plotErrors(errF, nominalQ=TRUE)
#oints are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence.
#The red line shows the error rates expected under the nominal definition of the Q-value. 
dev.off()

png(filename= "Output/errR.png")
plotErrors(errR, nominalQ=TRUE)
dev.off()
print("Error information is now in Output.")

setwd(filt_path)
####Dereplication Process####
print("Now dereplicating...")
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
print("Now denoising...")
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
#get sequence variants from each sample
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#do the same for the reverse set
dadaFs[[1]]


####Merge Paired Reads####
print("Now merging paired reads...")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
#had to add justConcatenate to get this working
head(mergers[[1]])


####Construct a Sequence Table####
print("Creating sequence table...")
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#matrix rows are the samples and columns are the variants so 


####Remove Chimeras####
print("Removing chimeras...")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim) 
sum(seqtab.nochim)/sum(seqtab)


####Track Reads Through Pipeline
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nochim")
rownames(track) <- sample.names
print( "Let's track the reads left after each step:" )
head(track)

####Assign Taxonomy####
print("Assigning taxonomy...this will take a while.")
setwd(path)
getwd()
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=FALSE)
#classify sequence variants taxonomically
taxa <- addSpecies(taxa, "silva_species_assignment_v132.fa.gz")
#add species level assignments by exact matching

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
print("Here is the taxonomy table:")
taxa.print
#display the taxon info
write.table(taxa.print, "Output/taxonomy.txt", sep="\t")

####Phyloseq####
library(phyloseq)
library(ggplot2)

now <- Sys.time()
now
print("Right version")
print("Beginning comparison plotting...")

#the variables for the sample data and the METS data are different. use the right chunk.

####gets the location from the filenames
#samples.out gets all of the filenames
#parse the name of the sample with strsplit and only look at the left side
#location creates a substring of the indices 21 to 26 (in this case USA and GHANA)
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "_S"), `[`, 1)
location <- substr(subject,21,26)
#remove the "-", extra 0, and the 1USA from location
location <- gsub('-', '',location)
location <- gsub('O', 'H20',location)
location <- gsub('1USA', 'USA', location)

#gets the sample number from the sample string
#sample_num creates a substring of the numbers from the filenames
sample_num <-substr(subject, 14, 16)
#removes the "-" from sample_num
sample_num <- gsub('-', '',sample_num)
#creates the data frame with a column called location and a column called SampleName
samdf <- data.frame(Location = location, SampleName = sample_num)
rownames(samdf) <- samples.out

#a phylo seq object is created
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps

print("Calculating alpha diversity...")
png(filename = "Output/alpha_diversity.png", width = 2000, height = 1200)
plot_richness(ps, x= "SampleName", measures=c("Shannon", "Simpson"), color="Location") + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))
dev.off()
print("alpha diversity plot saved to Output.")

#ords.nmds.bray <- ordinate(ps, method = "NMDS", distance = "bray")
#print("Creating ordination plot...")
#png(filename = "Output/ordination.png")
#plot_ordination(ps, ord.nmds.bray, color="Sample", title="Bray NMDS")
#dev.off()
#print("Ordination plot saved to Output.")

now <- Sys.time()
now

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

print("Creating family-level bar plot...")
png(filename = "Output/familybarplot.png", width = 2000, height = 1200)
plot_bar(ps.top20, x = "SampleName", fill="Family") + facet_wrap(~Location, scales="free_x") + theme(axis.text.x=element_text(angle=90,hjust=1))
dev.off()
print("Family-level bar plot saved to Output.")

print("Creating species-level bar plot...")
png(filename = "Output/speciesbarplot.png")
plot_bar(ps.top20, x= "SampleName", fill="Species") + facet_wrap(~Location, scales="free_x") + theme(axis.text.x=element_text(angle=90,hjust=1))
dev.off()
print("Species-level bar plot saved to Output.")

now <- Sys.time()
now

print("This is the end of the pipeline. All tables are tab-delimited and all plots are saved as .png files in the Output subdirectory within the working directory provided.")


now <- Sys.time()
now  