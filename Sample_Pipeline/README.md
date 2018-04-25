# METS Pipeline

Suzanne Greenwood, Christina Printzis, Nathaniel Polley

This pipeline is designed for microbiome comparison between two groups using paired-end 16S rRNA sequencing reads in FASTQ format.
Aftering filtering, trimming and further cleaning the raw data, taxonomy information is determined and outputted in a taxonomy table and comparison plots. Also included are alpha diversity plots and ordination plots.
In order to use this pipeline for FASTQ files not containing 16S rRNA sequencing reads, the database used to determine taxonomy must be changed. The default is Silva v132. 
Both the training set file and the database file must be in the working directory for taxonomic assignment.
The pipeline requires that input FASTQ files be named in the following format:

			Forward Reads: <sample information>_R1_001.fastq
				ex: GSF856-Dugas-99-1297-USA_S157_L001_R1_001.fastq
			Reverse Reads: <sample information>_R2_001.fastq
				ex: GSF856-Dugas-99-1297-USA_S157_L001_R2_001.fastq
				
The pipeline will not recognize which are forward and which are reverse reads unless the files bear the correct suffixes. 
The group name should be part of the <sample information> portion of the file name. (i.e. USA, Ghana, etc.)

# Before running this pipeline:
You must install the R packages dada2, phyloseq and ggplot2.
The installation instructions can be found at:

dada2: https://benjjneb.github.io/dada2/dada-installation.html
	 
phyloseq: https://joey711.github.io/phyloseq/install.html 
	
ggplot2: http://ggplot2.org 

	
The pipeline is broken into 2 R scripts that can be run in command line.

# Part I:
Analyzes the quality of the raw sequencing reads and outputs quality profiles for each appropriately-named file in the working directory.
	output files are .png format, and appear in a newly created folder within the provided working directory, Output.
	Output file names: fqual.png -- forward read quality profile
			   rqual.png -- reverse read quality profile

To run Part I: 
	
		Rscript part1.R -f <path_to_working_directory>
	
The -f is required for the program to recognize your path.
for help, type: 
	
		Rscript part1.R -h

If raw sequencing read quality is already known, and desired filtering parameters are known,  Part I can be skipped.


# Part II: 
Part II filters and trims the reads in all files, outputs quality profiles for the filtered data, and creates a new folder within the working directory called Filtered which contains the filtered reads.
	Output files are in .png or .txt formats and will again appear in the Output folder.
	Output file names: filtFqual.png -- quality profile of filtered forward files
			   filtRqual.png -- quality profile of filtered reverse files
			   filteredandtrimmed.txt -- tab-delimited table displaying file name	#reads in input	#reads in output

 
To run Parts II and III: 
	
	Rscript part2_miseq.R -f arg1 -F arg2 -R arg3  -T arg4>

where arg#:

 #1: working directory path; WHOLE PATH
			  
 #2: desired forward read length (recommended 240 for MiSeq_SOP data)
			  
 #3: desired reverse read length (recommended 160 for MiSeq_SOP data)
			  
 #4: how many nucleotides to trim from the left end. (recommended 0 for MiSeq_SOP data)
			  
			  
The flags can be also typed out as: --file, --For, --Rev, and --TrimL, respectively. The -h/--help flag is also available on parts 2 and 3.
 
 
Part III is the part of the pipeline that takes filtered reads, finds the unique sequences, merges paired end reads, assigns taxonomy, creates the taxonomy table and phyloseq plots.
	This is the most time-consuming part of the pipeline. 
	This is where the Silva v132 database is used.
	Output files: .png and .txt formats and will be located in the Output sub-directory
	Output file names: errF.png -- error plots for forward reads
			   errR.png -- error plots for reverse reads
			   taxonomy.txt -- tab-delimited taxonomy table
			   alpha_diversity.png -- alpha diversity plots using Shannon and Simpson indices.
			   familybarplot.png -- taxonomy comparison plot of the two groups at the family level
			   speciesbarplot.png -- taxonomy comparison plot of the two groups at the species level 


The command looks very much like that for Part I. The various commands within Part III will change the working directory for the appropriate sub-directory as needed (i.e. filtered, Output)
