methylcircleplot.R
==================

This guide will go through the process to make methylation figures (and also NOME-seq figures) similar to the one below:
![image](http://codingenes.files.wordpress.com/2012/10/fig-title.png)

Prerequisites
-------------
Software:
* [R](http://www.r-project.org/) must be installed. 
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

###Details
start up an R (or Rstudio) session and load the script using 

	source("methylcircleplot.R") 

If the command does not work, make sure the R session is in the same folder as the location of the file you downloaded (use `getwd()` to check and `setwd()` to change)

Next, specify the reference sequence and file/folder with bisulfite sequence. Remember to surround them with quotes.

~~~~
ref.seq = "GGAGTGAAGGCGGGACTTGTGCGGTTACCAGCGGAAATGCCTCGGGGTCAGAAGTCGCAGGAGAG
	ATAGACAGCTGCTGAACCAATGGGACCAGCGGATGGGGCGGATGTTATCTACCATTGGTGAACGTTAGAAAC
	GAATAGCAGCCAATGAATCAGCTGGGGGGGGCGGAGCAGTGACGTTTATTGCGGAGGGGGCCGCTTCGAATC
	GGCGGCGGCCAGCTTGGTGGCCTGGGCCAATGAACGGCCTCCAACGAGCAGGGCCTTCACCAATCGGCGGCC
	TCCACGACGGGGCTGGGGGAGGGTATAT"
bis.seq = "clone.fasta"
~~~~
Lastly, specify primers. (Optional step)

	fwd.primer = "GAGAAGAAAAAGTTTAGATTTTATAG"
	rev.primer = "AAACACCCCAATAAATCAATC"

To generate the figure run the following command (notice some options specified at the end)

	methylcircleplot(ref.seq, bis.seq, fwd.primer, rev.primer, reference=TRUE, NOME=2)
	
After you have created a suitable figure you can save it in PNG format by using the following command:
	
	png(filename="clone.png", width = 480, height = 480, units = "px")
	methylcircleplot(ref.seq, bis.seq, fwd.primer, rev.primer, reference=reference, NOME=NOME)
	dev.off()
	
It is possible to adjust the width and height to be larger numbers.
For more options on output format (such as tiff):

http://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/png.html

It is also possible to export to pdf this way

Optional parameters
-------------------
Example of what the parameters look like can be found [here](http://codingenes.wordpress.com/2012/08/23/script-methylation-figure-generation/#more-57)

rev.comp -- default = FALSE, reverse complements bisulfate sequence

size -- default = 2, size of circles

scaling -- default = 1, distance between circles, accepts value between 0 and 1 with < 1 moving the bubbles closer together

reference -- default = FALSE, show reference sequence at top

NOME -- default = 0, GpC detection: 0 - disabled (default) | 1 - plot together | 2 - plot GpC below CpG

noaxis -- default = FALSE, Turn off x/y-axis and labels

col.um -- default = "white", color of unmethylated CpG

col.me -- default = "black", color of methylated CpG

col.gme -- default = "lightgreen", color of methylated GpC

col.gum -- default = "aliceblue", color of unmethylated GpC

verbose -- default = TRUE, Display diagnostic messages

sampleName -- default = NULL, Specify sample names (Y-axis labels)

getAln -- default = FALSE, Return alignment between reference and samples (PairwiseAlignedFixedSubject object, see Biostrings for more details)

Additional Notes
----------------
By default, R will set the size of the figure and space out the rows as much as possible. 
It is possible to specify how much distance should be between the rows but this only works for a certain figure size.
If figure dimensions is increased (the png() command), you will need to adjust size and scaling values to keep the look the same
I have yet to find a way to 'autoscale' size and spacing to account for changes in figure size

List of colors in R can be found [here](http://research.stowers-institute.org/efg/R/Color/Chart/)
