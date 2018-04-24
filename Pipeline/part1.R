#Pipeline Part I: This part takes in the raw data and outputs the quality profile. The user should
#use this information to decide on filtering parameters. Filtering is the first step of Part II.

library(optparse)
library(dada2)

####This is the flag function####
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
#used to determine run time^

####Getting Set Up####
args = commandArgs(trailingOnly = TRUE)
path <- opt$file
#the user should put their working directory in the command to run the script
#path = "/homes/sgreenwood1/crossteam"


print("Here is what we're working with: ")
list.files(path)
#this lets us see the contents of the working directory to make sure we're looking at the right
#files and they're all reading in

dir.create("Output")

print("An Output folder has been created in the working directory.")

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs = sort(list.files(path, pattern= "R1_001.fastq", full.names = TRUE))
fnRs = sort(list.files(path, pattern= "R2_001.fastq", full.names = TRUE))
sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#now we have all our forward and reverse files grouped

###Next we take a look at quality###
print("Assessing read quality.")

png(filename = "Output/fqual.png")
plotQualityProfile(fnFs) #some trash reads but most are alright
#plots the quality profiles for all forward reads and outputs to Output
dev.off()

print("Forward read quality assessed.")

png(filename = "Output/rqual.png")
plotQualityProfile(fnRs)
#plots the quality profiles for all reverse reads. Outputs to Output
dev.off()
print("Reverse read quality assessed.")

print("The quality profiles for the forward and reverse reads have been saved in Output. Use this quality information to choose parameters for filtering and trimming.")

now = Sys.time()
now
#lets us know the run time for part I
