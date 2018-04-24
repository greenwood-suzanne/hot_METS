#Pipeline Part III: This is the part of the pipeline that takes filtered reads, finds the unique
#sequences, merges paired end reads, assigns taxonomy, creates the taxonomy table and phyloseq
#plots. 

library(optparse)
library(dada2)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Path to working directory folder", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

now <- Sys.time()
now

#run this with the same directory for an arg
path<-opt$file
#the user should put their working directory in the command to run the script
setwd(path)
getwd()
filtpath = paste(path, "filtered", sep = "/", collapse = NULL)
print(filtpath)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs = sort(list.files(path, pattern= "R1_001.fastq", full.names = TRUE))
fnRs = sort(list.files(path, pattern= "R2_001.fastq", full.names = TRUE))
sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#samples.names is needed to name the dereplicated ones later. that's why we pass the main directory as the wd to start

# Forward and reverse filtered files
setwd(filtpath) #now we move to filtered to actually derep and denoise
getwd()
filtFs = sort(list.files(filtpath, pattern= "F_filt.fastq", full.names = TRUE))
filtRs = sort(list.files(filtpath, pattern= "R_filt.fastq", full.names = TRUE))
print("We start with: ")
list.files(filtpath)
#assign vars for the filtered fwds and revs

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

setwd(filtpath)
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
#getN <- function(x) sum(getUniques(x))

# track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
# colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nochim")
# rownames(track) <- sample.names
# print( "Let's track the reads left after each step:" )
# head(track)
# #show us how many reads we have for each file as we filter down
# write.table(track, "Output/tracker.txt", sep="\t")
# print("Tracker table is now in Output.")

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

#the variables for the sample data and the METS data are different. use the right file.

####gets the location from the sample string FOR MISEQ SAMPLE DATA####
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out


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
####Only applies to MiSeq samples; no mock group in Mets data####
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
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
png(filename = "Output/familybarplot.png", width = 2000, height = 1200))
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x") + theme(axis.text.x=element_text(angle=90,hjust=1))
dev.off()
print("Family-level bar plot saved to Output.")

print("Creating species-level bar plot...")
png(filename = "Output/speciesbarplot.png")
plot_bar(ps.top20, x="Day", fill="Species") + facet_wrap(~When, scales="free_x") + theme(axis.text.x=element_text(angle=90,hjust=1))
dev.off()
print("Species-level bar plot saved to Output.")

now <- Sys.time()
now

print("This is the end of the pipeline. All tables are tab-delimited and all plots are saved as .png files in the Output subdirectory within the working directory provided.")

                                    
now <- Sys.time()
now  
