# Example of how to use makeVenn.R {WIP}

## Intro

I created these function because I got frustrated with how slow the ChIPseqAnno package was and 
the preceieved lack of utilities for making good venn diagrams with GRanges. 

This document is still a work in progress is currently only a brief outline of functions and an example.

## Functions

*peak2GRanges*(bedfile, type="macs", skip=0)
INCOMPLETE. The goal of this function is to convert peak caller output to GRanges
- bedfile is the name file that is peak caller output (typically a bed file)

*createResultMatrix*(typ, fo)
Generate a 'result matrix' (required for every venn diagram)
- results matrix has n columns where n is the number of sets being compared and nrow(fo) rows
- typ is a vector differentiating between the different GRanges
- fo is result of `findOverlaps()`

*extractOverlap*(..., res, typ)
Given res and typ, it will return a true/false vector indexed by whatever column name(s) it is passed in '...'

*printOverlap*(..., res, typ)
Print out the overlaps in a compressed way, need to explain this better when I remember what I did.

*readinGRanges*(...)
Put a bunch of GRanges in, returns result and type matrix `c(res, typ)`

*createOverlapMatrix*(res, typ) 
This function will create an overlap matrix. An overlap matrix is a human readable matrix that enumerates 
all possible overlaps. See source for more details and examples.

*createVenn*(res, typ, overlap = NA, doWeights = FALSE, ...)
Create venn diagram using Vennerable library (can draw up to 5-way venn diagrams, if there are libraries that draw better
I would be happy to stick it in).

*makeVennRunall*(...)
Put in GRanges, get a venn diagram (type and result matrix are returned invisibily)

*makeVennExample*()
INCOMPLETE This example only works if you have named `high1_peaks.bed` `high2_peaks.bed` and `high_peaks.bed` 
in the current directory with score in the 5th column.

## Example

Make a venn diagram from 3 GRanges objects named small/medium/large (in the future `peak2GRanges()` will take over this step)


```r
source("makeVenn.R")

tmpbed = read.table("small.narrowPeak")
    tmpgrg = GRanges(seqnames = Rle(tmpbed[,1]), ranges = IRanges(start=as.numeric(tmpbed[,2]),
end=as.numeric(tmpbed[,3]), names=tmpbed[,4]), score=tmpbed[,9])
small = tmpgrg
tmpbed = read.table("medium.narrowPeak")
    tmpgrg = GRanges(seqnames = Rle(tmpbed[,1]), ranges = IRanges(start=as.numeric(tmpbed[,2]),
end=as.numeric(tmpbed[,3]), names=tmpbed[,4]), score=tmpbed[,9])
medium = tmpgrg
tmpbed = read.table("large.narrowPeak")
tmpgrg = GRanges(seqnames = Rle(tmpbed[,1]), ranges = IRanges(start=as.numeric(tmpbed[,2]),
    end=as.numeric(tmpbed[,3]), names=tmpbed[,4]), score=tmpbed[,9])
large = tmpgrg

head(small)
```

```
## GRanges with 6 ranges and 1 metadata column:
##     seqnames                 ranges strand |            score
##        <Rle>              <IRanges>  <Rle> |        <numeric>
##   .    chr17 [  8057661,   8057930]      * | 2.29885307640994
##   .     chr5 [172296142, 172296382]      * | 2.29885307640994
##   .     chr5 [133802146, 133802339]      * | 2.29885307640994
##   .     chr6 [ 26285669,  26285953]      * | 2.29885307640994
##   .     chr3 [188665820, 188666023]      * | 2.29885307640994
##   .     chr6 [ 26189252,  26189536]      * | 2.29885307640994
##   ---
##   seqlengths:
##     chr1 chr10 chr11 chr12 chr13 chr14 ...  chr6  chr7  chr8  chr9  chrX
##       NA    NA    NA    NA    NA    NA ...    NA    NA    NA    NA    NA
```

```r
head(medium)
```

```
## GRanges with 6 ranges and 1 metadata column:
##     seqnames                 ranges strand |            score
##        <Rle>              <IRanges>  <Rle> |        <numeric>
##   .    chr17 [ 57902900,  57903239]      * | 3.53617953213724
##   .     chr5 [159247148, 159247438]      * | 3.53617953213724
##   .    chr10 [ 74008539,  74008891]      * | 3.53617953213724
##   .     chr5 [172296109, 172296456]      * | 3.53617953213724
##   .    chr20 [ 45962885,  45963193]      * | 3.53617953213724
##   .     chr8 [129169584, 129169867]      * | 3.53617953213724
##   ---
##   seqlengths:
##     chr1 chr10 chr11 chr12 chr13 chr14 ...  chr7  chr8  chr9  chrX  chrY
##       NA    NA    NA    NA    NA    NA ...    NA    NA    NA    NA    NA
```

```r
head(large)
```

```
## GRanges with 6 ranges and 1 metadata column:
##     seqnames                 ranges strand |            score
##        <Rle>              <IRanges>  <Rle> |        <numeric>
##   .     chr5 [159247124, 159247479]      * | 4.12205196263325
##   .    chr10 [ 74008476,  74008936]      * | 4.12205196263325
##   .    chr17 [ 57902860,  57903277]      * | 4.12205196263325
##   .     chr1 [  8250012,   8250391]      * | 4.12205196263325
##   .     chr2 [101736958, 101737460]      * | 4.12205196263325
##   .    chr12 [ 52625529,  52625906]      * | 4.12205196263325
##   ---
##   seqlengths:
##     chr1 chr10 chr11 chr12 chr13 chr14 ...  chr7  chr8  chr9  chrX  chrY
##       NA    NA    NA    NA    NA    NA ...    NA    NA    NA    NA    NA
```


You can then put it into `makeVennRunall(small, medium, large)` or go through each of the steps individually.
Skip using `readinGRanges()` since this function is not complete and might change later.


```r
glg = GRangesList(small, medium, large)  #GRanges list
typ = rep(as.character(substitute(list(small, medium, large)))[-1L], as.numeric(lapply(glg, 
    length)))
fo = findOverlaps(unlist(glg), ignoreSelf = T)  #find overlaps
res = createResultMatrix(typ, fo)  #results matrix
overlap = createOverlapMatrix(res, typ)
createVenn(res, typ, overlap)  #will return overlap matrix
```

![plot of chunk unnamed-chunk-2](https://raw.github.com/ying-w/bioinformatics-figures/master/makeVenn/figure/unnamed-chunk-2.png) 

```
##        large medium small
## all      890    887   865
## large     81   8969    98
## medium  9021      4     6
## small    103      7     0
## unique 13310    416    40
```


Make weighted venn diagram. When the differences in intersections are so great, the venn digram does not look good weighted.

```r
createVenn(res, typ, overlap, weighted = TRUE)  #will print out overlap matrix
```

![plot of chunk unnamed-chunk-3](https://raw.github.com/ying-w/bioinformatics-figures/master/makeVenn/figure/unnamed-chunk-3.png) 

```
##        large medium small
## all      890    887   865
## large     81   8969    98
## medium  9021      4     6
## small    103      7     0
## unique 13310    416    40
```


Isolate regions in `large` that overlap all regions

```r
# large must be specified first in extractOverlap
large_intersect_all = large[extractOverlap("large", "medium", "small", res = res, 
    typ = typ)]
length(large_intersect_all)
```

```
## [1] 890
```

```r
head(large_intersect_all)
```

```
## GRanges with 6 ranges and 1 metadata column:
##     seqnames                 ranges strand |            score
##        <Rle>              <IRanges>  <Rle> |        <numeric>
##   .     chr5 [159247124, 159247479]      * | 4.12205196263325
##   .    chr10 [ 74008476,  74008936]      * | 4.12205196263325
##   .    chr17 [ 57902860,  57903277]      * | 4.12205196263325
##   .     chr2 [101736958, 101737460]      * | 4.12205196263325
##   .    chr12 [ 52625529,  52625906]      * | 4.12205196263325
##   .     chr5 [172296068, 172296532]      * | 4.12205196263325
##   ---
##   seqlengths:
##     chr1 chr10 chr11 chr12 chr13 chr14 ...  chr7  chr8  chr9  chrX  chrY
##       NA    NA    NA    NA    NA    NA ...    NA    NA    NA    NA    NA
```


Isolate regions that are unique to `medium`

```r
medium_unique = medium[extractOverlap("medium", res = res, typ = typ)]
length(medium_unique)
```

```
## [1] 416
```

```r
head(medium_unique)
```

```
## GRanges with 6 ranges and 1 metadata column:
##     seqnames               ranges strand |            score
##        <Rle>            <IRanges>  <Rle> |        <numeric>
##   .    chr19 [39174799, 39174872]      * | 3.53617953213724
##   .    chr15 [59492386, 59492538]      * | 3.53617953213724
##   .    chr14 [51471910, 51472010]      * | 2.93459943821543
##   .    chr20 [48922435, 48922719]      * |  2.6998377258693
##   .    chr11 [29141291, 29141575]      * |  2.6268875962802
##   .    chr11 [66654038, 66654322]      * | 2.54203878013516
##   ---
##   seqlengths:
##     chr1 chr10 chr11 chr12 chr13 chr14 ...  chr7  chr8  chr9  chrX  chrY
##       NA    NA    NA    NA    NA    NA ...    NA    NA    NA    NA    NA
```


Isolate regions that are self overlaps in `medium`

```r
medium_self_overlap = medium[extractOverlap("medium", "medium", res = res, typ = typ)]
length(medium_self_overlap)
```

```
## [1] 4
```

```r
head(medium_self_overlap)
```

```
## GRanges with 4 ranges and 1 metadata column:
##     seqnames                 ranges strand |             score
##        <Rle>              <IRanges>  <Rle> |         <numeric>
##   .     chr8 [  8521987,   8522271]      * |  1.80817124057177
##   .     chr7 [130538426, 130538710]      * |  1.03774313022021
##   .    chr12 [ 92996740,  92997024]      * | 0.981889396176367
##   .    chr18 [  3711022,   3711306]      * | 0.963924634554715
##   ---
##   seqlengths:
##     chr1 chr10 chr11 chr12 chr13 chr14 ...  chr7  chr8  chr9  chrX  chrY
##       NA    NA    NA    NA    NA    NA ...    NA    NA    NA    NA    NA
```


Note to self: I followed the example shown here: https://github.com/yihui/knitr-examples
