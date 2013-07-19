#Make venn diagram from GenomicRanges
#Author: Ying Wu daiyingw@gmail.com
#Current Update: March 2013
#License: ___ (still deciding)
#Version: TESTING

require(GenomicRanges)

peak2GRanges = function(bedfile, type="macs", skip=0) {

	#The goal of this function is to convert peak caller output to GRanges
	#bedfile is the name file that is peak caller output (typically a bed file)
	
	#type == "macs xls" #start pos +1, score -> -log?(fdr)
	#type == "macs" #see below, BED6?
	#type == "zinba" #read w/header=T, must convert score to >1, has self overlaps 
	#type == "homer" #
	# cname = c("PeakID", "chr", "start", "end", "strand", "Normalized.Tag.Count", "focus.ratio",
	# "findPeaks.Score", "Total.Tags", "Control.Tags", "Fold.Change.vs.Control", "p.val.vs.ctrl",
	# "Fold.Change.vs.Local", "p.val.vs.Local", "Clonal.Fold.Change")
	# g1.r = GRanges(seqnames = g1$chr, ranges = IRanges(start=as.numeric(g1$start),
		# end=as.numeric(g1$end), names=g1$PeakID),
		# strand=g1$strand, score=g1$Normalized.Tag.Count)
	#careful of keeping elementmetadata the same
	
	#TODO: implement other types, keep in mind some have header
	#	macs 2 and macs 1.4 might have different output
	grgb = read.table(bedfile, sep="\t", skip=skip)
	#type==macs
	grg = GRanges(seqnames = Rle(grgb[,1]), ranges = IRanges(start=as.numeric(grgb[,2]),
		end=as.numeric(grgb[,3]), names=grgb[,4]),
		strand=Rle("+"), score=grgb[,5])
	grg
	#self overlaps can be removed using reduce() but it would also get rid of metadata
}
	
createResultMatrix = function(typ, fo) {
 
	#generate result matrix (required for every venn diagram)
	#results matrix has n columns where n is the number of sets being compared
	#results matrix has nrow(fo) rows
	
	typlvl = levels(factor(typ)) #type (which GRanges) -> factor -> level
	res = matrix(0, nrow=length(typ), ncol=length(typlvl)) #result matrix
	colnames(res) = typlvl
	for(i in 1:ncol(res))
	{ #since num of sets being compared should be relatively low I did not apply this
		tmp = table(queryHits(fo)[typ[subjectHits(fo)] == typlvl[i]])
		res[as.numeric(names(tmp)),i] = tmp
	}
	res
}

#self overlaps will be double counted
#self overlaps can be odd numbers:
#<- 1 ->
#<--- 1 --->                   #self overlap
#       <--- 1 --->            #2 overlap
#              <--- 2 --->     #1 overlap

#TODO: examples and what happens when there is only one parameter
extractOverlap = function(..., res, typ) {  
    #http://stackoverflow.com/questions/3057341/how-to-use-rs-ellipsis-feature-when-writing-your-own-function
    
	#This is to read in case where 1st argument is a list (probably a better way to do this)
	if(length(list(...)) == 1) { argv = as.list(...) }
	else { argv = list(...) }
	#cat(paste("DEBUG: ",paste(argv, collapse=" "),"\n"))
	
	#TODO: input checking
	if(length(argv) == 0) { stop("Must specify at least one set\n\tUsage: extractOverlap(1,2,res=res,typ=typ)") }
	#if(length(argv) > ncol(res)) { stop("number of overlaps greater than sets in venn diagram") }
	argv = as.character(argv)
	if(!all(argv %in% colnames(res))) { 
		stop("Invalid set: ", paste(argv[!argv %in% colnames(res)], collapase=""), 
		"\n\tpossible sets are: ", paste(colnames(res), collapase="")) 
	}
	
	basecol = colnames(res) == argv[1] 
	curtyp = typ == argv[1]
	ret = NULL
	
	#must wrap w/as.matrix to account for case of 1 row
	if(length(argv) == 1) { ret = apply(as.matrix(res[curtyp, !basecol] == 0), 1, all) } #unique to base
	if(length(argv) == 2 && argv[1] == argv[2]) { #self overlap
		ret = res[curtyp, basecol] > 0 & apply(as.matrix(res[curtyp, !basecol] == 0), 1, all) 
	}
	#if(length(unique(argv)) != length(argv)) { warning("Duplicate set removed") }
	argv = unique(argv) #get rid of duplicates (should not exist anyways)
	
    #TODO: does not catch case where 3 elements all the same length(argv) == 3 && length(unique(argv)) == 1
    # Use == 0 instead of >= 1 because harder to deal with multiple overlaps with the latter
	if(length(argv) >= 2 && argv[1] != argv[2]) { #everything else
		othercol = as.matrix(res[curtyp, !(colnames(res) %in% argv)])
		if(ncol(othercol) == 0) { 
			othercol = rep(TRUE, nrow(othercol))
		} else {
			othercol = apply(othercol == 0, 1, all) #all others == 0
		}
		ret = apply(as.matrix(res[curtyp, colnames(res) %in% argv[-1]]) > 0, 1, all) & #all specified sets are >0
			 othercol #all others == 0
	}
	ret
}

printOverlap = function(..., res, typ) {
	#called by createOverlapMatrix()
	#cat(paste("DEBUG: Processing",paste(as.character(...), collapse=" "),"\n"))
	n = ncol(res)
	tf = factor(typ)
	arow = matrix(0,1,n)
	if(n == length(...)) {
		for(i in 1:n) {
			arow[,i] = sum(extractOverlap(as.list(levels(tf)[c(i:n,1:i)[1:n]]),res=res,typ=typ))
		}
	} else {
		arow = sapply(lapply(strsplit(paste(levels(tf)[c(1:n)], paste(as.character(...),collapse=" ")),split=" "), 
            extractOverlap, res=res, typ=typ),sum)
	}
	arow
}

readinGRanges = function(...) {
	#This is kind of ugly but neccessary to drop elementMetadata mismatches
	tmp = list(...)
	elmMd = lapply(tmp , function(x) names(elementMetadata(x)))
	if(length(unique(lapply(elmMd, length))) == 1) {
		elmMdKeep = lapply(apply(as.data.frame(elmMd),1, unique),length) == 1
		warning("Metadata column(s) ", paste(which(!elmMdKeep), collapse=","), " mismatch and will be dropped")
	} else { 
		elmMdKeep = FALSE
		warning("All metadata column(s) will be dropped")
	}
	for(i in 1:length(tmp)) {
		elementMetadata(tmp[[i]]) = elementMetadata(tmp[[i]])[elmMdKeep]
	}
	glg = GRangesList(tmp)
	#glg = GRangesList(...) #gives error on elementMetadata colname mismatch
	
	if(length(glg) < 2) { stop("Need more ranges to compare") }
	if(length(glg) > 4) { warning("only tested up to 5 ranges") } #this somewhere else
	#typ = rep(2^(0:(length(glg)-1)),as.numeric(lapply(glg, length))) 
	#fo = findOverlaps(c(g1.r, g2.r, g12.r, .ignoreElementMetadata=TRUE), ignoreSelf=T)
	typ = rep(as.character(substitute(list(...)))[-1L], as.numeric(lapply(glg, length))) #since lapplay returns list
	fo = findOverlaps(unlist(glg), ignoreSelf=T)
	
	res = createResultMatrix(typ, fo)
	#cat(paste(c(paste(colnames(res), as.character(substitute(list(...)))[-1L],sep=" = "),""),collapse="\n"))	
	cbind(typ,res)
}

createOverlapMatrix = function(res, typ) {
	#This function will create an overlap matrix
	#An overlap matrix is a human readable matrix that enumerates all possible overlaps
	#An example of this for a 3-way venn diagram is as follows:
	#       A   B	C
	#all	3	3	3
	# A		2	7	11
	# B		6	8	12
	# C		11	12	13
	#unique	5	0	14
	#
	#The 'all' row is overlap of ABC where the number shown is the number of
	#  elements that contribute to the overlap from each of the 3 GRanges
	#  ex. all-A means that A contributes 3 ranges to overlap set
	#The A/B/C rows shows the overlap of A/B/C with each of the columns with
	#  self overlaps included (A-A, B-B, C-C), typically these self overlaps are 0
	#Lastly are the GRanges unique to each
	#
	#In the case of a 4-way venn diagram, an overlaps matrix would look like:
	#		A	B	C	D
	#all
	# AB
	# AC
	# AD
	# BC
	# BD
	# CD
	#self
	#unique
	#
	#This pattern can be generalized so that for an n-way comparison,
	#  the following rows will be needed in the overlaps matrix
	#  (since column is always the same as n)
	#  all (n), n-2 (n-2), n-4 (n-4), ... unique(n-k) where k is the first value > n
	#  all the ones in between (n-1, n-3) will be filled implicitly
	#  ex. in n=3 you will have: All (3), A/B/C (1), unique(-1)
	#  ex. in n=4: All(4), AB/AC/AD/BC/BD/CD (2), self (0), unique(-1)
	#  self is a special case when counter is 0
	#  you can see that in the case of n=3, self is included in (1) as A-A
    #
    #This function makes heavy use of printOverlap which is described above
	
	n = ncol(res)
	tf = factor(typ)
	last_printed = n
	overlap = matrix(0, sum(choose(n,seq(n-2,0,by=-2)))+2,n)
	rownames(overlap) = 1:nrow(overlap) #needs initalization before assignment
	colnames(overlap) = levels(tf)
	current_row = 1
	
	for(i in n:0) {
		#cat(paste("DEBUG: i=", i, "current_row =", current_row, "\n"))
		if(i == n-1) { #special case: all overlap
			overlap[current_row,] = printOverlap(levels(tf),res=res, typ=typ)
			rownames(overlap)[current_row] = "all"
			current_row = current_row+1 
		} else if(i == last_printed-2) { 
			if(i == 0) { #special case: self overlap (only explicitly shown in some)
				overlap[current_row,] = sapply(apply(rbind(levels(tf), levels(tf)), 2, extractOverlap, res=res, typ=typ),sum)
				rownames(overlap)[current_row] = "self"
				current_row = current_row+1 
			} else { # base case, combn in apply is the key part
				overlap[current_row:(current_row+choose(n,i)-1),] = t(apply(combn(levels(tf), i), 2, printOverlap, res=res, typ=typ)) 
				last_printed = i
				rownames(overlap)[current_row:(current_row+choose(n,i)-1)] = apply(combn(levels(tf), i),2,paste, collapse=" ")
				current_row = current_row + choose(n,i)
			}
		}
	}
	overlap[nrow(overlap),] = t(apply(combn(levels(tf), i), 2, printOverlap, res=res, typ=typ))
	rownames(overlap)[nrow(overlap)] = "unique"
	overlap
}

createVenn = function(res, typ, overlap = NULL, name = NULL, weighted = FALSE, main=NULL) {
    # This function is a bit complicated since I've tried to generalize the # of columns that res can have. 
    # I might switch over to using VennDiagram library instead of Vennerable since the latter
    # is not well documented and maintainer does not really respond to inquries
    # VennDiagram can be found here: http://cran.r-project.org/web/packages/VennDiagram/ http://pubmed.gov/21269502
    # both make use of the grid library to draw, see end of function for more details
    # more options here: http://www.biostars.org/p/7713/
    
    if(!require("Vennerable", quietly = TRUE)) { stop("Missing Vennerable library. Please install via:
	install.packages(\"Vennerable\", repos=\"http://R-Forge.R-project.org\", dependencies=TRUE) ")}

	if(any(is.null(overlap))) { overlap = createOverlapMatrix(res, typ) }
	if(length(name) != ncol(res)) { name = colnames(res) } 
    
	#see createOverlapMatrix() for details
	n = ncol(res)
	tf = factor(typ)
	last_printed = n
	counter = matrix(0, sum(choose(n,seq(n-2,0,by=-2)))+2,n)
	rownames(counter) = 1:nrow(counter) #needs initalization before assignment
	colnames(counter) = levels(tf)
	current_row = 1
	weight = 2^(0:(n-1))
	
	for(i in n:0) {
		if(i == n-1) { 
			counter[current_row,] = rep(2^n-1, n)
			rownames(counter)[current_row] = "all"
			current_row = current_row+1 
		} else if(i == last_printed-2) { 
			if(i == 0) {
				counter[current_row,] = rep(0,n)
				rownames(counter)[current_row] = "self"
				current_row = current_row+1 
			} else {
				#good luck -
				#combn() will generate the rows, weight will hold the column 
				#sum(unique()) will give you value for that cell in matrix
				#combn() should be added to every weight giving you each row
				#there might be a more readable way of doing this but should still need 2x apply
				counter[current_row:(current_row+choose(n,i)-1),] = t(apply(combn(weight, i), 2, function(z) { 
					sapply(weight, function(x,y) { sum(unique(c(x,y))) }, z) }))

				last_printed = i
				rownames(counter)[current_row:(current_row+choose(n,i)-1)] = apply(combn(levels(tf), i),2,paste, collapse=" ")
				current_row = current_row + choose(n,i)
			}
		}
	}
	counter[counter %in% weight] = 0 #remove self overlaps
	counter[nrow(counter),] = weight
	rownames(counter)[nrow(counter)] = "unique"

	#create venn diagram using Vennerable
    plot.new() #this is needed for text
	vc = sapply(1:(2^(n)-1), function(x) { round(median(overlap[counter %in% x])) })
	plotVenn(Venn(SetNames=name, Weight=c(0,vc)), doWeights=weighted)
    # cannot pass in ... to above function
    #text(-0.1, 0.01, main) #the coordinates change
    # mtext(main, side=1, at=0) #this is an ugly hack since cannot pass to plotVenn()
    # above mtext command causes "invalid graphics state" errors
    
    # VennDiagram library takes list of overlaps

    # ballpark method
    # convert results matrix (res) into this list
    # ll = as.list(as.data.frame(res))
    # ll2 = lapply(ll, identical, 0)
    # grid.draw(venn.diagram(x=ll2, filename=NULL, scaled=TRUE))
    # this method doesnt work very well since scaling wont work and numbers are off
    
    # convert from Venn() input
	# counter = matrix("", sum(vc)+2,n) #different from above
	# colnames(counter) = levels(tf)
    ## fill up matrix with characters that are unique, VennDiagram uses intercept()
    # vennlist = as.list(as.data.frame(counter))
    ##names(vennlist) = colnames(res)
    # grid.draw(venn.diagram(x=vennlist, filename=NULL, scaled=TRUE))
    ## scaling does not work for 3-way (email on 5/1/2013 says that they will look into it)
    # draw.triple.venn(100, 240, 85, 21, 46, 17, 6, scaled=TRUE, euler.d=TRUE) #example code
    ## once you have 2 non-zero intersections, the scaling code will not work
}

makeVennRunall = function(...) {
	#suppressMessages(gc(verbose=FALSE)) #it likes to talkings
	tmp = readinGRanges(...) #split tmp into res/typ
	typ = tmp[,1]
	res = apply(tmp[,2:ncol(tmp)], 2, as.numeric)
	overlap = createOverlapMatrix(res,typ)
	createVenn(res, typ, overlap)
    invisible(tmp)
}

makeVennExample = function() {
	g1 = read.table("./high1_peaks.bed", sep="\t")
	g2 = read.table("./high2_peaks.bed", sep="\t")
	g12 = read.table("./high_peaks.bed", sep="\t")
	#g1.r = BED2RangedData(g1)
	#g1.r = as(g1.r, "GRanges")

	#note score change
	g1.r = GRanges(seqnames = Rle(g1[,1]), ranges = IRanges(start=as.numeric(g1[,2]),
		end=as.numeric(g1[,3]), names=g1[,4]),
		strand=Rle("+"), score=g1[,5])
	g2.r = GRanges(seqnames = Rle(g2[,1]), ranges = IRanges(start=as.numeric(g2[,2]),
		end=as.numeric(g2[,3]), names=g2[,4]),
		strand=Rle("+"), score=g2[,5])
	g12.r = GRanges(seqnames = Rle(g12[,1]), ranges = IRanges(start=as.numeric(g12[,2]),
		end=as.numeric(g12[,3]), names=g12[,4]),
		strand=Rle("+"), score=g12[,5])

	runall(g1.r, g2.r, g12.r)
}