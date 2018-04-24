# hot_METS
Group Repo for CompBio METS project

# Introduction:
  This is the METS project, a pipeline for microbiome analysis that can be run in the command line. It comes in 3 parts, all of which are contained in the Pipeline folder of this repository. You will need sample data to ensure that you have correctly installed all dependencies as well as the Silva v132 database files. These are also available in this repository. 
  
          
# Software Requirements:
This pipeline can be run as a series of Rscripts through command line. The required packages are:

 dada2: https://benjjneb.github.io/dada2/dada-installation.html
	
phyloseq: https://joey711.github.io/phyloseq/install.html 

ggplot2: http://ggplot2.org 
	

# To Import from GitHub:
In addition to the crossteam directory (see Introduction) you will need the crossteam.R file to run the pipeline. 
Import this file from GitHub using:

	git clone https://github.com/greenwood-suzanne/hot_METS.git

	git fork https://github.com/greenwood-suzanne/hot_METS.git

  
# Implementation:
   The pipeline will take in this directory for input and runs through the following major steps:
   
   -Filter and trim the raw data: Remove PCR/Sequencing primers, exclude all N bases, exclude all reads with low PHRED scores.
      
   -Evaluate errors: this is a step neessary for getting taxonomy information later on. Outputs plots of error x quality. 
      (May take a few minutes)
      
   -Pull out just the unique sequences: Remove duplicates from our filtered data set
      
   -Merge the paired reads: Here we combine the forward and reverse reads so that we can ascertain taxonomy information 
      in subsequent steps
      
   -Remove chimeras: We want to have only properly matched forward and reverse reads 
      
   -Plotting with phyloseq:  The output bar graphs are created for the taxonomic groups present in our
      two sample groups: Ghana and US for comparison .

In order to run the pipeline through the command line your unix commands should follow this format:
 	
	Rscript part1.R -f <your_directory_here>

 	Rscript part2.R -f <your_directory> -F <forward_filter_length> -R <reverse_filter_length> -T a<#nucls_to_trim_from_left>

	 Rscript part3.R -f <your_directory_here>

Part 2 must come before part 3, but part 1 is optional if the user already knows their desired filtering parameters. The flags (-char) are required. For help, try:
	
	Rscript part#.R -h



# Results:																						
The first few lines of the tables that are produced throughout the running of the pipeline can be viewed in the nohup.out file that
will appear after running is complete.

	Tables:
		-out: out is a table that summarizes the number of reads in each file before and after filtering by PHRED score, 
		removing Ns and removing primers.
		-seqtab: probably not useful to look at for this activity; lists occurances of reads in each file
		-seqtab.nochim: same as seqtab but with chimeras removed. 
		-track: track is a table that shows you how many reads are in each file after each step of the pipeline. 
		Recall that each of the steps that make up the columns of the table are removing "bad" reads whether they be
		full of N bases, duplicates, or contain chimeras.
		-taxa: preliminary taxonomy table. Will likely take several minutes to run.
		-taxa.print: taxonomy table adjusted for display. This contains all the taxonomy information we are able to obtain
		from checking against the silva version 132 database.
			
	PNG Files:
		*note: If there is trouble viewing the image files, RStudio is a nice option. Simply navigate to your home directory 
		under the files tab on the right hand side of the screen. There you should see your copy of the crossteam directory and
		all the .png files among your documents. There is also the option of SCP-ing the images to your own computer for viewing.
		-filtFqual.png: qualilty plot of forward reads after filtering
		-filtRqual.png: qualilty plot of reverse reads after filtering
		-fqual.png: qualilty plot of forward reads before filtering
		-Rqual.png: qualilty plot of reverse reads before filtering
			
You will notice a new subdirectory "filtered" within your crossteam directory. This contains the filtered files. 
All of your raw data files are unchanged; all filtering/denoising/dereplicating, etc. was done within the filtered subdirectory.
