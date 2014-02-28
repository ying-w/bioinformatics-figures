# GREAT to Entrez ID
[Great](http://great.stanford.edu/) is a program that can associate ChIP-seq signal regions with genes. Using their helpful web interface will output a list of genes that are associated with input regions. However, when using human (hg19) the list of associated genes are listed by gene symbol. The gene symbols used are quite outdated and gene symbols are not good identifiers. Thus, I have attempted to map all all genes used by GREAT's interface for hg19 to Entrez gene IDs.

## Approach
Genes for GREAT's hg19 are [based off](http://bejerano.stanford.edu/help/display/GREAT/Genes) an old knownGene table. Some of these gene symbols have changed. Using UCSC table browser, I ouput chr/start/stop/sym/locusID(entrez gene)/refseq/ensembl from knownCanonical and knownGenes table. I will use chr+TSS and gene symbols to map the gene symbols used by great to a entrez gene ID.

## Files
```
$ head Hg19.great2.0.genes.txt 
5	chr1	69090	+	OR4F5
11	chr1	367658	+	OR4F16
17	chr1	622034	-	OR4F16
30	chr1	861120	+	SAMD11
31	chr1	894679	-	NOC2L
32	chr1	895966	+	KLHL17
35	chr1	935552	-	HES4
36	chr1	948846	+	ISG15
37	chr1	955502	+	AGRN
39	chr1	1009687	-	RNF223
$ wc -l Hg19.great2.0.genes.txt #17744 enteries but only 17524 unique symbols
17744 Hg19.great2.0.genes.txt
$ zcat canonical.gz | head
#hg19.knownCanonical.chrom	hg19.knownCanonical.chromStart	hg19.knownCanonical.chromEnd	hg19.knownCanonical.clusterId	hg19.kgXref.geneSymbol	hg19.knownGene.strand	hg19.knownToEnsembl.value	hg19.knownToLocusLink.value	hg19.knownToRefSeq.value
chr1	11873	14409	1	DDX11L1	+	ENST00000518655	100287102	NR_046018
chr1	14361	19759	2	WASH7P	-	ENST00000438504	653635	NR_024540
chr1	14406	29370	3	WASH7P	-	ENST00000438504	653635	NR_024540
chr1	34610	36081	4	FAM138F	-	ENST00000417324	641702	NR_026820
chr1	69090	70008	5	OR4F5	+	ENST00000335137	79501	NM_001005484
chr1	134772	140566	6	LOC729737	-	ENST00000493797	729737	NR_039983
chr1	321083	321115	7	DQ597235	+	n/a	n/a	n/a
chr1	321145	321207	8	DQ599768	+	n/a	n/a	n/a
chr1	322036	326938	9	LOC100133331	+	ENST00000440038	100133331	NR_028327
$ zcat knowngeneall.gz | head
#hg19.knownGene.chrom	hg19.knownGene.strand	hg19.knownGene.txStart	hg19.knownGene.txEnd	hg19.kgXref.geneSymbol	hg19.knownToEnsembl.value	hg19.knownToLocusLink.value	hg19.knownToRefSeq.value
chr1	+	11873	14409	DDX11L1	ENST00000456328	100287102	NR_046018
chr1	+	11873	14409	DDX11L1	ENST00000456328	100287102	NR_046018
chr1	+	11873	14409	DDX11L1	ENST00000518655	100287102	NR_046018
chr1	-	14361	16765	WASH7P	ENST00000423562	653635	NR_024540
chr1	-	16857	17751	WASH7P	ENST00000541675	653635	NR_024540
chr1	-	15795	18061	WASH7P	ENST00000488147	653635	NR_024540
chr1	-	14361	19759	WASH7P	ENST00000438504	653635	NR_024540
chr1	-	14361	19759	WASH7P	ENST00000438504	653635	NR_024540
chr1	-	14361	19759	WASH7P	ENST00000438504	653635	NR_024540
```

Columns for GREAT are:

	<ucscClusterId>	<tssChrom>	<tssCoord>	<tssStrand>	<geneSymbol>

*caution* First column for GREAT is <ucscClusterId> which is ID for transcript in knownCanonical, however this ID changes 

*caution* Even though GREAT uses knownCanonical, the gene symbols are not unique

Columns for canonical are:

    <chr>   <start> <stop>  <clusterID> <symbol>    <strand>    <ensembl>   <entrez>    <refseq>

Columns for knowngeneall are:

    <chr>   <strand>    <start> <stop>  <symbol>    <ensembl>   <entrez>    <refseq>

All columns are tab seperated

## Matching

```r
# read in all GREAT genes
greatall = read.table("Hg19.great2.0.genes.txt", stringsAsFactors = F)
greatout = cbind(greatall, entrez = NA, ensembl = NA, refseq = NA)
greatall = cbind(greatall, chrtss = paste(greatall[, 2], greatall[, 3], sep = ":"), 
    stringsAsFactors = F)
head(greatall)
```

```
##   V1   V2     V3 V4     V5      chrtss
## 1  5 chr1  69090  +  OR4F5  chr1:69090
## 2 11 chr1 367658  + OR4F16 chr1:367658
## 3 17 chr1 622034  - OR4F16 chr1:622034
## 4 30 chr1 861120  + SAMD11 chr1:861120
## 5 31 chr1 894679  -  NOC2L chr1:894679
## 6 32 chr1 895966  + KLHL17 chr1:895966
```

```r

# read in current knownGenes from UCSC (retrieved 2/2014)
ucscall = read.table("knowngeneall.gz", sep = "\t", stringsAsFactors = F)
ucscall = cbind(ucscall, chrtss = paste(ucscall[, 1], ucscall[, 3], sep = ":"), 
    stringsAsFactors = F)
rev = ucscall[, 2] == "-"
ucscall[rev, 9] = paste(ucscall[rev, 1], ucscall[rev, 4], sep = ":")
head(ucscall)
```

```
##     V1 V2    V3    V4      V5              V6        V7        V8
## 1 chr1  + 11873 14409 DDX11L1 ENST00000456328 100287102 NR_046018
## 2 chr1  + 11873 14409 DDX11L1 ENST00000456328 100287102 NR_046018
## 3 chr1  + 11873 14409 DDX11L1 ENST00000518655 100287102 NR_046018
## 4 chr1  - 14361 16765  WASH7P ENST00000423562    653635 NR_024540
## 5 chr1  - 16857 17751  WASH7P ENST00000541675    653635 NR_024540
## 6 chr1  - 15795 18061  WASH7P ENST00000488147    653635 NR_024540
##       chrtss
## 1 chr1:11873
## 2 chr1:11873
## 3 chr1:11873
## 4 chr1:16765
## 5 chr1:17751
## 6 chr1:18061
```

```r

# match using TSS
fixme = is.na(match(greatall[, 6], ucscall[, 9]))
table(fixme)  #992
```

```
## fixme
## FALSE  TRUE 
## 16752   992
```

```r
# greatall[head(which(fixme)),]
greatout[!fixme, 6] = ucscall[match(greatall[, 6], ucscall[, 9])[!fixme], 7]  # entrez
greatout[!fixme, 7] = ucscall[match(greatall[, 6], ucscall[, 9])[!fixme], 6]  # ensembl
greatout[!fixme, 8] = ucscall[match(greatall[, 6], ucscall[, 9])[!fixme], 8]  # refseq
head(greatout)
```

```
##   V1   V2     V3 V4     V5 entrez         ensembl       refseq
## 1  5 chr1  69090  +  OR4F5  79501 ENST00000335137 NM_001005484
## 2 11 chr1 367658  + OR4F16 729759 ENST00000426406 NM_001005277
## 3 17 chr1 622034  - OR4F16 729759 ENST00000332831 NM_001005277
## 4 30 chr1 861120  + SAMD11 148398 ENST00000342066    NM_152486
## 5 31 chr1 894679  -  NOC2L  26155 ENST00000327044    NM_015658
## 6 32 chr1 895966  + KLHL17 339451 ENST00000338591    NM_198317
```

```r

# match using gene symbol from ucsc
canonical = read.table("canonical.gz", sep = "\t", stringsAsFactors = F)
head(canonical)
```

```
##     V1     V2     V3 V4        V5 V6              V7        V8
## 1 chr1  11873  14409  1   DDX11L1  + ENST00000518655 100287102
## 2 chr1  14361  19759  2    WASH7P  - ENST00000438504    653635
## 3 chr1  14406  29370  3    WASH7P  - ENST00000438504    653635
## 4 chr1  34610  36081  4   FAM138F  - ENST00000417324    641702
## 5 chr1  69090  70008  5     OR4F5  + ENST00000335137     79501
## 6 chr1 134772 140566  6 LOC729737  - ENST00000493797    729737
##             V9
## 1    NR_046018
## 2    NR_024540
## 3    NR_024540
## 4    NR_026820
## 5 NM_001005484
## 6    NR_039983
```

```r
fixme2 = fixme & is.na(match(greatall[, 5], canonical[, 5]))
table(fixme2)  #51
```

```
## fixme2
## FALSE  TRUE 
## 17693    51
```

```r
# greatall[head(which(fixme2)),]
greatout[!fixme2, 6] = ucscall[match(greatall[, 6], ucscall[, 9])[!fixme2], 
    7]  # entrez
greatout[!fixme2, 7] = ucscall[match(greatall[, 6], ucscall[, 9])[!fixme2], 
    6]  # ensembl
greatout[!fixme2, 8] = ucscall[match(greatall[, 6], ucscall[, 9])[!fixme2], 
    8]  # refseq

# match using gene symbol against alias/depreciated symbols
library(illuminaHumanv4.db)
```

```
## Loading required package: AnnotationDbi
## Loading required package: BiocGenerics
## 
## Attaching package: 'BiocGenerics'
## 
## The following object(s) are masked from 'package:stats':
## 
##     xtabs
## 
## The following object(s) are masked from 'package:base':
## 
##     anyDuplicated, cbind, colnames, duplicated, eval, Filter,
##     Find, get, intersect, lapply, Map, mapply, mget, order, paste,
##     pmax, pmax.int, pmin, pmin.int, Position, rbind, Reduce,
##     rep.int, rownames, sapply, setdiff, table, tapply, union,
##     unique
## 
## Loading required package: Biobase
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Loading required package: org.Hs.eg.db
## Loading required package: DBI
```

```r
alias2ez = mget(as.character(greatall[fixme2, 5]), org.Hs.egALIAS2EG, ifnotfound = NA)
alias2ez[lapply(alias2ez, length) > 1] = NA  #nonunique 
fixme3 = is.na(alias2ez)
table(fixme3)  #6
```

```
## fixme3
## FALSE  TRUE 
##    45     6
```

```r

greatall[fixme2, ][fixme3, ]
```

```
##          V1    V2       V3 V4           V5         chrtss
## 3130   5017 chr11 64863682  +      C11orf2 chr11:64863682
## 3157   5055 chr11 65547822  - DKFZp761E198 chr11:65547822
## 3366   5357 chr11 89586009  +    LOC399940 chr11:89586009
## 3371   5362 chr11 89724838  -    LOC399940 chr11:89724838
## 9120  15399 chr19 53784810  +    LOC646508 chr19:53784810
## 14459 24268  chr6 66011353  +     AK302514  chr6:66011353
```

```r
# C11orf2 -> VPS51 DKFZp761E198 -> AP5B1 LOC399940 -> TRIM53AP (closest
# match) LOC399940 -> TRIM53AP (closest match) LOC646508 -> FAM90A27P
# (closest match) AK302514 -> LOC441155
alias2ez[[which(fixme3)[1]]] = unique(ucscall[ucscall[, 5] == "VPS51", 7])
alias2ez[[which(fixme3)[2]]] = unique(ucscall[ucscall[, 5] == "AP5B1", 7])
alias2ez[[which(fixme3)[3]]] = unique(ucscall[ucscall[, 5] == "TRIM53AP", 7])
alias2ez[[which(fixme3)[4]]] = unique(ucscall[ucscall[, 5] == "TRIM53AP", 7])
alias2ez[[which(fixme3)[5]]] = unique(ucscall[ucscall[, 5] == "FAM90A27P", 7])
alias2ez[[which(fixme3)[6]]] = unique(ucscall[ucscall[, 5] == "LOC441155", 7])

greatout[fixme2, 6] = unlist(alias2ez)  # entrez

# replace NA with 'n/a'
greatout[is.na(greatout[, 6]), 6] = "n/a"
greatout[is.na(greatout[, 7]), 7] = "n/a"
greatout[is.na(greatout[, 8]), 8] = "n/a"

head(greatout)
```

```
##   V1   V2     V3 V4     V5 entrez         ensembl       refseq
## 1  5 chr1  69090  +  OR4F5  79501 ENST00000335137 NM_001005484
## 2 11 chr1 367658  + OR4F16 729759 ENST00000426406 NM_001005277
## 3 17 chr1 622034  - OR4F16 729759 ENST00000332831 NM_001005277
## 4 30 chr1 861120  + SAMD11 148398 ENST00000342066    NM_152486
## 5 31 chr1 894679  -  NOC2L  26155 ENST00000327044    NM_015658
## 6 32 chr1 895966  + KLHL17 339451 ENST00000338591    NM_198317
```

```r

# write.table(greatout,
# file=paste('Hg19.great2.0.genes.entrez.',format(Sys.time(),
# '%Y-%m-%d'),'.txt',sep=''),sep='\t',row.names = FALSE, quote=F)

sessionInfo()
```

```
## R version 2.15.3 (2013-03-01)
## Platform: x86_64-pc-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_US.utf8       LC_NUMERIC=C             
##  [3] LC_TIME=en_US.utf8        LC_COLLATE=en_US.utf8    
##  [5] LC_MONETARY=en_US.utf8    LC_MESSAGES=en_US.utf8   
##  [7] LC_PAPER=C                LC_NAME=C                
##  [9] LC_ADDRESS=C              LC_TELEPHONE=C           
## [11] LC_MEASUREMENT=en_US.utf8 LC_IDENTIFICATION=C      
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] illuminaHumanv4.db_1.16.0 org.Hs.eg.db_2.8.0       
## [3] RSQLite_0.11.4            DBI_0.2-7                
## [5] AnnotationDbi_1.20.7      Biobase_2.18.0           
## [7] BiocGenerics_0.4.0        knitr_1.5                
## 
## loaded via a namespace (and not attached):
## [1] AnnotationForge_1.0.3 evaluate_0.5.1        formatR_0.10         
## [4] IRanges_1.16.6        parallel_2.15.3       stats4_2.15.3        
## [7] stringr_0.6.2         tools_2.15.3
```


# Running great with custom set of TSS
[Download](http://bejerano.stanford.edu/help/display/GREAT/Download) and follow the instructions in the README

Also follow install directions for [kent tools](http://genome-source.cse.ucsc.edu/gitweb/?p=kent.git;a=blob;f=src/product/README.building.source) up to 3a

[intersectBed](http://bedtools.readthedocs.org/en/latest/content/tools/intersect.html?highlight=intersectbed) can be used to find ChIP-seq regions overlapping gene regions.

The last annotation release from NCBI using hg19 is v105 and can be found [here](ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/)
