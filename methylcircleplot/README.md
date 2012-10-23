methylcircleplot.R
======================

This guide will go through the process to make methylation figures (and also NOME-seq figures) similar to the one below:
http://codingenes.files.wordpress.com/2012/08/fig0.png

Prerequisites
-------------
Software:
[R](http://www.r-project.org/) must be installed
  If you are new to R, I would suggest also installing [R-studio](http://www.rstudio.com/)
[Bioconductor](http://bioconductor.org/install/) must be installed
[Biostrings](http://bioconductor.org/packages/release/bioc/html/Biostrings.html) package in bioconductor
[Script file](https://raw.github.com/ying-w/bioinformatics-figures/methylcircleplot/master/methylcircleplot.R) Needs to be saved in the current working directory

Sequence:
Reference to be used to identify potential CpG and GpCs sites
Bisulfate converted sequence (files must end in .txt or .fasta)
1. In a multi-fasta file  

	>clone1
	atggatgttttaggttttttagaggatggttgagtgggttgtaaggatag
	gtcgagagggtgtagtgttaataggttttgtggtgcgatggggtattcga
	>clone2
	atggatgttttaggttttttagaggatggttgagtgggttgtaaggatag
	gtcgagagggtgtagtgttaataggttttgtggtgcgatggggtattcga
	>clone3
	atggatgttttaggttttttagaggatggttgagtgggttgtaaggatag
	gtcgagagggtgtagtgttaataggttttgtggtgcgatggggtattcga
	
2. In a folder with each file being a single-fasta file
3. In a file with each line being a different set

	atggatgttttaggttttttagaggatggttgagtgggttgtaaggataggtcgagagggtgtagtgttaataggttttgtggtgcgatggggtattcga
	atggatgttttaggttttttagaggatggttgagtgggttgtaaggataggtcgagagggtgtagtgttaataggttttgtggtgcgatggggtattcga
	atggatgttttaggttttttagaggatggttgagtgggttgtaaggataggtcgagagggtgtagtgttaataggttttgtggtgcgatggggtattcga

Primer sequence (optional) in 5'->3' orientation

Introduction
------------
Talk about how the program works here

Directions
----------
start up an R (or Rstudio) session and load the script using 

	source("methylation_figure.R") 

If the command does not work, make sure the R session is in the same folder as the location of the file you downloaded (use `getwd()` to check and `setwd()` to change)

Next, specify the reference sequence and file/folder with bisulfite sequence

	ref.seq = "???"
	bis.seq = "clone.fasta"

Lastly, specify primers or any other options you want

	fwd.primer = "???"
	rev.primer = "???"
	
To generate the figure run the following command

	methylcircleplot(ref.seq, bis.seq, fwd.primer, rev.primer)
	
Optional parameters
-------------------
rev.comp -- default = FALSE, reverse complements bisulfate sequence
size -- default = 2, size of circles
scaling -- default = 1, enter number between 0 and 1 with 1 being max distance
reference -- default = FALSE, show reference sequence at top?
NOME -- default = 0, GpC detection: 0 disabled (default) | 1 plot together | 2 plot apart
noaxis -- default = FALSE, Turn off x/y-axis labels
col.um -- default = "white", color of unmethylated CpG
col.me -- default = "black", color of methylated CpG
col.gme -- default = "lightgreen", color of methylated GpC
col.gum -- default = "lightblue", color of unmethylated GpC
verbose -- default = TRUE, Display diagnostic messages
cloneName -- default = NULL, Specify sample names (Y-axis labels)
getAln -- default = FALSE, return alignment between reference and samples

