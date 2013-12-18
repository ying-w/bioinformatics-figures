bioinformatics-figures
======================

Collection of scripts and examples that will help visualize sequencing data. 
Most of the scripts are for visualizing sequencing overlaps from ChIP-seq / DNase-seq / FAIRE-seq experiments using venn diagrams

* [Venn Diagram examples](Venn Diagram examples/) -- Making venn diagrams in R using limma / Vennerable / VennDiagram libraries (input format conversion)
* [makeVenn](makeVenn/makeVenn.md)/makeVenn.R -- will create 2-5 way venn diagrams from GRanges objects using Vennerable library
* [methylcircleplot](methylcircleplot)/methylcircleplot.R -- will create circles that either filled or not filled depending on methylation status
* [igv.R](igv.R) -- R function to take genomic coordinates and get (pre-exisiting) [igv](http://www.broadinstitute.org/igv/) session to output snapshots at thoes locations 
* [imageplot](myImagePlot.R) -- Is a simple example of making a plot in R with colors (corresponding to expression level) [example](http://www.phaget4.org/R/image_matrix.html)

A tip for others using R or knitR on github: to highlight code on in markdown documents use 

    ```S

instead of

    ```R
    
