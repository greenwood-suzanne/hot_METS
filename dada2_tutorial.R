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
#we need to do some quality trimming here but it isnt 
#done in the tutorial so we'll have to figure out how

plotQualityProfile(fnRs[1:2])
#the reverse read quality isnt great. more than half the reads have phred scores below 25/20

  ####Assign the filenames for filtered fastq.gz files####
filt_path <- file.path(path, "filtered") 
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
#filtered forwards
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
#filtered reverses

  ####Now we filter####
out = filterAndTrim(fnFs, filtFs, truncQ= 20, truncLen = 210, maxN=0, maxEE=c(4,4), rm.phix=TRUE,compress=TRUE, multithread=TRUE) 
#out is a table displaying the file name, number of reads to start with
#and number of reads after truncating at the given positions (210 on forward and 120 on reverse to keep scores above 20)
head(out)
#prints the first few rows showing the original read counts next to the filtered counts


  ####Examine errors####
errF <- learnErrors(filtFs, multithread=TRUE)
#tells us how amny unique reads we have in each sample
#skipped the errR step since we arent looking at our reverse reads at the moment

plotErrors(errF, nominalQ=TRUE)
#oints are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence.
#The red line shows the error rates expected under the nominal definition of the Q-value. 
#there are reverse steps too

  ####Dereplication Process####
derepFs <- derepFastq(filtFs, verbose=TRUE)
# Name the derep-class objects by the sample names
#Dereplication combines all identical sequencing reads into into
#“unique sequences” with a corresponding “abundance”: the number of 
#reads with that unique sequence
names(derepFs) <- sample.names
#gives us a table of all the unique sequences from all the forwards and how many times they occur
#when we go back, there will be reverse steps too

  ####denoising####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
#infer sequence variants in each sample
#there would also be reverse steps 

dadaFs[[1]]
#261 real sequence variants from 4221 unique sequences in the first sample

  ####Merge Paired Reads####


