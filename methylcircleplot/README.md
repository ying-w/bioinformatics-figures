methylcircleplot.R
======================

This guide will go through the process to make methylation figures (and also NOME-seq figures) similar to the one below:
http://codingenes.files.wordpress.com/2012/08/fig0.png

Prerequisites
-------------
Software:
* [R](http://www.r-project.org/) must be installed
  If you are new to R, I would suggest also installing [R-studio](http://www.rstudio.com/)
* [Bioconductor](http://bioconductor.org/install/) must be installed
* [Biostrings](http://bioconductor.org/packages/release/bioc/html/Biostrings.html) package in bioconductor
* [Script file](https://raw.github.com/ying-w/bioinformatics-figures/master/methylcircleplot/methylcircleplot.R) Needs to be saved in the current working directory

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

Overview
--------
Talk about how the program works here

Example
-------
start up an R (or Rstudio) session and load the script using 

	source("methylcircleplot.R") 

If the command does not work, make sure the R session is in the same folder as the location of the file you downloaded (use `getwd()` to check and `setwd()` to change)

Next, specify the reference sequence and file/folder with bisulfite sequence

~~~~
ref.seq = "GGAGTGAAGGCGGGACTTGTGCGGTTACCAGCGGAAATGCCTCGGGGTCAGAAGTCGCAGGAGAG
	ATAGACAGCTGCTGAACCAATGGGACCAGCGGATGGGGCGGATGTTATCTACCATTGGTGAACGTTAGAAAC
	GAATAGCAGCCAATGAATCAGCTGGGGGGGGCGGAGCAGTGACGTTTATTGCGGAGGGGGCCGCTTCGAATC
	GGCGGCGGCCAGCTTGGTGGCCTGGGCCAATGAACGGCCTCCAACGAGCAGGGCCTTCACCAATCGGCGGCC
	TCCACGACGGGGCTGGGGGAGGGTATAT"
bis.seq = "clone.fasta"
~~~~
Lastly, specify primers or any other options you want

	fwd.primer = "GAGAAGAAAAAGTTTAGATTTTATAG"
	rev.primer = "AAACACCCCAATAAATCAATC"

To generate the figure run the following command

	methylcircleplot(ref.seq, bis.seq, fwd.primer, rev.primer, reference=TRUE, NOME=2)
	
After you have created a suitable figure you can save it in PNG format by using the following command:
	
	png(filename="clone.png", width = 480, height = 480, units = "px")
	methylcircleplot(ref.seq, bis.seq, fwd.primer, rev.primer, reference=reference, NOME=NOME)
	dev.off()
	
It is possible to adjust the width and height to be larger numbers.
For more options on output format (such as tiff):
http://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/png.html
you can also export to pdf this way

Whole script below:

~~~~ R
#load script
source("methylcircleplot.R") 
#specify reference
ref.seq = "GGAGTGAAGGCGGGACTTGTGCGGTTACCAGCGGAAATGCCTCGGGGTCAGAAGTCGCAGGAGAG
	ATAGACAGCTGCTGAACCAATGGGACCAGCGGATGGGGCGGATGTTATCTACCATTGGTGAACGTTAGAAAC
	GAATAGCAGCCAATGAATCAGCTGGGGGGGGCGGAGCAGTGACGTTTATTGCGGAGGGGGCCGCTTCGAATC
	GGCGGCGGCCAGCTTGGTGGCCTGGGCCAATGAACGGCCTCCAACGAGCAGGGCCTTCACCAATCGGCGGCC
	TCCACGACGGGGCTGGGGGAGGGTATAT"
#specify bisulfite converted sequence file
bis.seq = "clone.fasta"
#specify primers
fwd.primer = "GAGAAGAAAAAGTTTAGATTTTATAG"
rev.primer = "AAACACCCCAATAAATCAATC"
#specify image name + size
png(filename="myclones.png", width=600, height=600)
#generate image
methylcircleplot(ref.seq, bis.seq, fwd.primer, rev.primer, reference=TRUE)
#save image
dev.off()
~~~~


Optional parameters
-------------------
Example of difference can be found [here](http://codingenes.wordpress.com/2012/08/23/script-methylation-figure-generation/#more-57)
rev.comp -- default = FALSE, reverse complements bisulfate sequence

size -- default = 2, size of circles

scaling -- default = 1, enter number between 0 and 1 with anything less than 1 moving the bubbles closer together

reference -- default = FALSE, show reference sequence at top?

NOME -- default = 0, GpC detection: 0 disabled (default) | 1 plot together | 2 plot apart

noaxis -- default = FALSE, Turn off x/y-axis labels

col.um -- default = "white", color of unmethylated CpG

col.me -- default = "black", color of methylated CpG

col.gme -- default = "lightgreen", color of methylated GpC

col.gum -- default = "aliceblue", color of unmethylated GpC

verbose -- default = TRUE, Display diagnostic messages

sampleName -- default = NULL, Specify sample names (Y-axis labels)

getAln -- default = FALSE, Return alignment between reference and samples of class PairwiseAlignedFixedSubject 

Additional Notes
----------------
Talk about issues with resizing and scaling

List of colors that R allows can be found [here](http://research.stowers-institute.org/efg/R/Color/Chart/)