# hot_METS
Group Repo for CompBio METS project

# Introduction:
  The following file provides information regarding implementation and requirements for the METS project, a pipeline for terminal-based microbiome analysis. 
  
  
# Requirements:

 * Sample data 
 	* located in the MiSeq_SOP folder on this repository
 * Silva v132 database files. 
  	* DADA2-formatted Silva training files can be found at: https://zenodo.org/record/1172783#.WtvfsWbMyu5
 * Scripts 'part1.R' and 'part2_miseq.R'
  	* located in 'Sample_Pipeline' folder
	
 *The METS Pipeline folder contains the same files as Sample_Pipeline, however input parameters have been specifically ammended for the METS data (not published in repository due to size constraints).
  
          
# Software Requirements:
The pipeline can be run as a series of R scripts through the command line. The required packages include:

dada2: https://benjjneb.github.io/dada2/dada-installation.html
	
phyloseq: https://joey711.github.io/phyloseq/install.html 

ggplot2: http://ggplot2.org 
	

# To Import from GitHub:
In addition to the crossteam directory (see Introduction) you will need the crossteam.R file to run the pipeline. 
Import this file from GitHub using:

	git clone https://github.com/greenwood-suzanne/hot_METS.git

	git fork https://github.com/greenwood-suzanne/hot_METS.git

  
# Implementation:
   The pipeline will input in the specified folder directory containing your FASTQ reads and will perform the following major steps:
   
   Part I:
   -Evaluate quality of raw data
   
   Part II:
   -Filter and trim the raw data: Removes PCR/Sequencing primers, excluding all N bases and all reads with PHRED scores below specified threshold.
   
   -Evaluate errors: neessary for the acquisition of taxonomy information in latter steps. Outputs plots of error x quality. 
      (May take a few minutes) 
      
   -Sort out only unique sequences: removes duplicates from filtered data set
   
   -Merge paired reads: combines the forward and reverse reads to prepare for the ascertainment of taxonomy information 
      in subsequent steps
      
   -Remove chimeras: isolates properly-paired forward and reverse reads 
   
   Part III:
   -Plot with Phyloseq: creates bar-graph output for the taxonomic groups present in the
      two-sample groups (Ghana and US for comparison)

In order to run the pipeline through the command line, the Unix commands should follow this format:
 	
	Rscript part1.R 
	-f <working_directory_folder>

 	Rscript part2and3.R 
	-f <working_directory_folder> 
	-F <forward_filter_length> 
	-R <reverse_filter_length> 
	-T <#nucls_to_trim_from_left>
	
*For MiSeq_SOP sample data, we recommend using -F 240 -R 160 -T 0

For assistance, please consult the 'help' prompt:
	
	Rscript part#.R -h


# Results:																						
Only the first few lines of the tables will be displayed in the terminal. All plots are saved as .png files in a folder titled 'Output' within the specified working directory. The graphical output is as follows:

	Tables:
		-out: table summarizing the number of reads in each file before and after filtering by PHRED score and after 		the removal of Ns and primers
		-seqtab: list containing reads occurances in each file
		-seqtab.nochim: seqtab without chimeras 
		-track: table displaying frequency of reads in each file after each operation of the pipeline 
		(Each step should display a decrease in reads, as low quality data is truncated from the original dataset.)
		-taxa: preliminary taxonomy table - (will likely require several minutes to run)
		-taxa.print: taxonomy table from alignment with Silva v. 132 database - (adjusted for display)
					
	PNG Files:
		-filtFqual.png: qualilty plot of forward reads after filtering
		-filtRqual.png: qualilty plot of reverse reads after filtering
		-fqual.png: qualilty plot of forward reads before filtering
		-Rqual.png: qualilty plot of reverse reads before filtering
		-familybarplot.png: taxonomy bar plot on family-level
		
		*note: If issues arise while attempting to view image files, RStudio provides a user-friendly graphical interface. Simply navigate to the home directory under the files tab on the right-hand side of the screen. There one can view a copy of the provided directory and all .png images associated with the input files. Images can also be copied to one's personal terminal using the Unix command 'scp'.
			
A new subdirectory "filtered" will be created within the specified working directory, containing the filtered reads. 
All of original raw data files remain unaltered; all filtering/denoising/dereplicating, etc. was performed within the '~/filtered' subdirectory.
