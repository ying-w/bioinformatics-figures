methylcircleplot.R
==================

This guide will go through the process to make methylation figures (and also NOME-seq figures) similar to the one below:
![image](http://codingenes.files.wordpress.com/2012/10/fig-title.png)

Prerequisites
-------------
**Software:**
* [R](http://www.r-project.org/) must be installed. 
  If you are new to R, I would suggest also installing [R-studio](http://www.rstudio.com/)
* [Bioconductor](http://bioconductor.org/install/) must be installed
* [Biostrings](http://bioconductor.org/packages/release/bioc/html/Biostrings.html) package in bioconductor
* [Script file](https://raw.github.com/ying-w/bioinformatics-figures/master/methylcircleplot/methylcircleplot.R) Needs to be saved in the current working directory

**Sequence:**

Reference to be used to identify potential CpG and GpCs sites

Bisulfate converted sequence (files must end in .txt or .fasta)

1. In a multi-fasta file  
	
        >methyl
        tgggctgaaatactgggttcacccatatttaggattttaggcgggtgggcaagtaagaattga
        ggagtggttttagaaataattggtatacgaatatttaatggatgttttaggttttttagagga
        tggttgagtgggttgtaaggataggtcgagaatggctggacacctggcttcag
        >unmethyl
        tgggctgaaatactgggttcacccatatttaggattttaggcgggtgggtaagcaagaattga
        ggagtggctttagaaataattggcatatgaatatttaatggatgttttaggctttttagagga
        tggctgagtgggctgtaaggataggctgagaatggctggacacctggcttcag
        >mix-error
        tgggctgaaatactgggttcacccatatttaggattttaggcgggtgggtaagtaagaattga
        ggagtggttttagaaataattggcatatgaatatttaatggatgttttaggctttttagagga
        tggctgagtggggtgtaaggataggtcgagaatggctggacacctggcttcag

2. In a folder with each file being a single-fasta file
3. In a file with each line being a different set

        tgggctgaaatactgggttcacccatatttaggattttaggcgggtgggcaagtaagaattgaggagtggttttagaaataattggtatacgaatatttaatggatgttttaggttttttagaggatggttgagtgggttgtaaggataggtcgagaatggctggacacctggcttcag
        tgggctgaaatactgggttcacccatatttaggattttaggcgggtgggtaagcaagaattgaggagtggctttagaaataattggcatatgaatatttaatggatgttttaggctttttagaggatggctgagtgggctgtaaggataggctgagaatggctggacacctggcttcag
        tgggctgaaatactgggttcacccatatttaggattttaggcgggtgggtaagtaagaattgaggagtggttttagaaataattggcatatgaatatttaatggatgttttaggctttttagaggatggctgagtggggtgtaaggataggtcgagaatggctggacacctggcttcag

Primer sequence (optional) in 5'->3' orientation

Overview
--------
This program uses reference sequence to find CpG and GpC sites. Primers are used to isolate the sequence of interest. 
Sequence of interest is then aligned to reference sequence to find corresponding CpG and GpC sites. 
These sites are checked for methylation and result is displayed with circles.

Example
-------
![image](http://codingenes.files.wordpress.com/2012/11/fig-example.png)
~~~~ R
#load script
source("methylcircleplot.R") 
#specify reference
ref.seq = "atatctaggactctaggcgggtgggtaa
	gcaagaactgaggagtggccccagaaataattggcacacga
	acattcaatggatgttttaggctctccagaggatggctgag
	tgggctgtaaggacaggccgaga"
#specify bisulfite converted sequence file
bis.seq = "clone.fasta" #multi-fasta example file above
#specify primers
fwd.primer = "tgggctgaaatactgggttcaccc"
rev.primer = "atggctggacacctggcttcag"
#specify image name + size
png("myclones.png", width=550, height=400)
#generate image
methylcircleplot(ref.seq, bis.seq, fwd.primer, 
  rev.primer, reference=TRUE, NOME=2, col.gme = "lightgrey", col.gum = "white")
#save image
dev.off()
~~~~

###Details
start up an R (or Rstudio) session and load the script using 

	source("methylcircleplot.R") 

If the command does not work, make sure the R session is in the same folder as the location of the file you downloaded (use `getwd()` to check and `setwd()` to change)

Next, specify the reference sequence and file/folder with bisulfite sequence. Remember to surround them with quotes.

~~~~
ref.seq = "atatctaggactctaggcgggtgggtaa
	gcaagaactgaggagtggccccagaaataattggcacacga
	acattcaatggatgttttaggctctccagaggatggctgag
	tgggctgtaaggacaggccgaga"
bis.seq = "clone.fasta"
~~~~
Lastly, specify primers. (Optional step)

	fwd.primer = "tgggctgaaatactgggttcaccc"
	rev.primer = "atggctggacacctggcttcag"

To generate the figure run the following command (notice some options specified at the end)

	methylcircleplot(ref.seq, bis.seq, fwd.primer, rev.primer, reference=TRUE, NOME=2)
	
After you have created a suitable figure you can save it in PNG format by using the following command:
	
	png(filename="clone.png", width = 600, height = 600, units = "px")
	methylcircleplot(ref.seq, bis.seq, fwd.primer, rev.primer, reference=TRUE, NOME=2)
	dev.off()
	
It is possible to adjust the width and height to be larger numbers.
For more options on output format (such as tiff): see http://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/png.html

It is also possible to export to pdf this way (see additional notes below)

Optional parameters
-------------------
While reference (`ref.seq`) and bisulfite sequence (`bis.seq`) are required, there are other parameters that can be set.

Some examples of what the parameters look like can be found [here](http://codingenes.wordpress.com/2012/08/23/script-methylation-figure-generation/#more-57)

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

showNumUnconverted -- default = FALSE, Show number of unconverted Cs

cloneName -- default = NULL, Specify clone names (Y-axis labels)

cloneOrder -- default = NULL, Specify the order in which to rearrange the rows (clones)

getAln -- default = FALSE, Return alignment between reference and samples (PairwiseAlignedFixedSubject object, see Biostrings for more details)

Additional Notes
----------------
By default, R will set the size (resolution) of the figure and space out the rows as much as possible. 
It is possible to specify how much distance should be between the rows but this only works for a certain figure size/resolution.
If figure dimensions is increased (see png() command), you will need to adjust size and scaling values to keep the look the same
I have yet to find a way to 'autoscale' size and spacing to account for changes in figure size/resolution.

By having a `png()` command before and a `dev.off()` command after, R will save the last plot generated before `dev.off()` in the .png file specified.
It is also possible to export `pdf()` and `tiff()` by replacing `png()`. `pdf()` and `svg()` require width and height to be specified in inches such as `svg("clone.svg", width=5.5, height=4)` (default is 7x7)
Occasionally you may encounter some odd banding patterns while using `png()` to export. 
To work around this, I typically export with `pdf()` and use the linux `convert` utility to change the file from .pdf to .png. 
Another workaround would be to export using `svg()` and then convert to png (using Inkscape). 
I have found that svg files produced by R sometimes have font incompatibility issues with Adobe Illustrator, to 
work around this and keep everything vectorized I first open the .svg file using Inkscape, then export a .emf file for Illustrator to open.

List of colors in R can be found [here](http://research.stowers-institute.org/efg/R/Color/Chart/)

