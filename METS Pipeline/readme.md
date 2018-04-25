# METS Pipeline

This pipeline is broken into 2 R scripts for microbiome analysis. The version of the pipeline in this folder is designed specifically for the METS dataset.
To adapt this pipeline to other data sets, the phyloseq portion would need to be edited with the appropriate variables for plotting.

In order to run the pipeline through the command line your unix commands should follow this format:

`Rscript part1.R -f <your_directory_here>`

`Rscript part2and3.R -f <your_directory> -F <forward_filter_length> -R <reverse_filter_length> -T a<#nucls_to_trim_from_left>`

*For METS data, we recommend using -F 320 -R 220 -T 15

For help, try:

`Rscript part#.R -h`
