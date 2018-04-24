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

####Here we have to again define fwd and reverse reads####
path <- opt$file
#the user should put their working directory in the command to run the script
setwd(path)

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
write.table(out, "Output/filteredandtrimmed.txt", sep="\t")
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

