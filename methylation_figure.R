#DNA methylation figure
#Author: Ying Wu daiyingw@gmail.com
#Date: Aug 2012
#License: ___ for private use only
#
#Step by step examples and directions can be found here:
# http://codingenes.wordpress.com/2012/08/23/script-methylation-figure-generation/
#
#Directions for use:
#Change ref.seq/fwd.primer/rev.primer/bis.seq 
#Run the methylcircleplot function

examplerun = function()
{
#reference sequence (non-bisulfite converted)
ref.seq = "tgggttgaaatattgggtttatttatatttaggattttagacgggtgggt
aagtaagaattgaggagtggttttagaaataattggtatacgaatattta
atggatgttttaggttttttagaggatggttgagtgggttgtaaggatag
gtcgagagggtgtagtgttaataggttttgtggtgcgatggggtattcga
gtaattggtttgtgaggtgttcggtgatttaaggtaggggtgagaggatt
ttgaaggttgaaaatgaaggttttttggggtttcgttttaagggttgttt
tgtttagacgtttttaattttcgtttggaagatataggtagatagcgttc
gttttagttttttttatttttatagttttgtttttttatttatttagggg
gcggggttagaggttaaggttagagggtgggattggggagggagaggtga
AATCGTTTTTAGGTGAGTCGTTTTTTTATTAGGTTTTCGGTTCGGGGTGT
TTATTTTTTTTATGGTTGGATATTTGGTTTTAG"
rev.comp = TRUE #reverse complement bisulfite sequence?
size = 2 #size of methylation bubbles/
scaling = 1 #proportion of figure to use for plotting (useful for adjusting closeness)
fwd.primer = "tgggttgaaatattgggtttattt"
rev.primer = "TATGGTTGGATATTTGGTTTTAG"
#if no primers, leave it as blank, ie ""
#primer should be before bisulfite treatment

#bis.seq can take either txt/fasta file, character vector, NULL (default)
#	if NULL then begin interactive prompt for sequence
#vector:
bis.seq = c(
"NNNNNNNNNNNNGAGCTCGGATCCACTAGTAACGGCCGCCAGTGTGCTGGAATTCGCCCTTCTAAAACCAAATATCCAACCATAAAAAAAATAAACACCCCAAACCAAAAACCTAATAAAAAAACAACTCACCTAAAAACAATTTCACCTCTCCCTCCCCAATCCCACCCTCTAACCTTAACCTCTAACCCCACCCCCTAAATAAATAAAAAAACAAAACTATAAAAATAAAAAAAACTAAAACAAACACTATCTACCTATATCTTCCAAACAAAAATTAAAAACATCTAAACAAAACAACCCTTAAAACAAAACCCCAAAAAACCTTCATTTTCAACCTTCAAAATCCTCTCACCCCTACCTTAAATCACCAAACACCTCACAAACCAATTACTCAAATACCCCATCACACCACAAAACCTATTAACACTACACCCTCTCAACCTATCCTTACAACCCACTCAACCATCCTCTAAAAAACCTAAAACATCCATTAAATATTCATATACCAATTATTTCTAAAACCACTCCTCAATTCTTACTTACCCACCCATCTAAAATCCTAAATATAAATAAACCCAATATTTCAACCCAAAGGGCGAATTCTGCAGATATCCATCACACTGGCGGCCGCTCGAGCATGCATCTAGAGGGCCCAATTCGCCCTATAGTGAGTCGTATTACAATTCACTGGCCGTCGTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGNAATAGCGAAGAGGCNGCACCGATCGCCTTCCCAACAGTTGCGCAGCCTGATGGCGATGGACGCGCCCTGNAGCGCGCATAAGCGCGCGGGTGTGTGTTACGCGCAGCGTGACCGCTACACTGCAGCNCCTAGCGCCGCTCTTTCGCTTCTTCCTTCTTCTCGCACGTCGCGGNTTCCCNTCAGCTCTAATCGGGGCTCCTTNNGGTCCNATTAGNNNTTANNNACTCGACCCAAAACTGATAGGTGANGTCNCGTANTNGGCATCGCCTGANAGACGTTTCGCCNTGACGNGNGTCNNNTNNNNNANANNGANNCTNGTCAANNGGACANNNCACCNANNNNNNTCTTNNTNNANGATTNCGATNGCNNTGTAANNGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
,"NNNNNNNNNNNNGAGCTCGGATCCACTAGTAACGGCCGCCAGTGTGCTGGAATTCGCCCTTCTAAAACCAAATATCCAACCATAAAAAAAATAAACACCCCAAACCAAAAACCTAATAAAAAAACAACTCACCTAAAAACAATTTCACCTCTCCCTCCCCAATCCCACCCTCTAACCTTAACCTCTAACCCCACCCCCTAAATAAATAAAAAAACAAAACTATAAAAATAAAAAAAACTAAAACAAACACTATCTACCTATATCTTCCAAACAAAAATTAAAAACATCTAAACAAAACAACCCTTAAAACAAAACCCCAAAAAACCTTCATTTTCAACCTTCAAAATCCTCTCACCCCTACCTTAAATCACCAAACACCTCACAAACCAATTACTCAAATACCCCATCACACCACAAAACCTATTAACACTACACCCTCTCAACCTATCCTTACAACCCACTCAACCATCCTCTAAAAAACCTAAAACATCCATTAAATATTCATATACCAATTATTTCTAAAACCACTCCTCAATTCTTACTTACCCACCCATCTAAAATCCTAAATATAAATAAACCCAATATTTCAACCCAAAGGGCGAATTCTGCAGATATCCATCACACTGGCGGCCGCTCGAGCATGCATCTAGAGGGCCCAATTCGCCCTATAGTGAGTCGTATTACAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGACGCNNNNNNGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTACACTTGCAGCGCCNANNNNNGCNCTTCGTTCTCTCTNNNCNNNNCNNNNNNNNNANNNNNNNNACANTANNNNANNNNNNNANAANANNNNNCNANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGNNNNNNNNNGNTNNNNNNNNCTCNNNNNNNNNCNNNNNNNNNNTNNNNNNNNNNNNNNTNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
)
#file: (must end in .txt or .fasta)
#bis.seq = "./Desktop/fullset.txt" 
#tip: to change directory go to: Tools | Change Working Dir of RStudio

#RUN THIS FUNCTION
methylcircleplot(ref.seq, bis.seq, fwd.primer, rev.primer, rev.comp, size, scaling)
dev.copy(png, "mefigure.png")
dev.off()
#dev.copy2pdf will keep resized figure dimensions

#There are 3 options for NOME
# 0 disabled (default)
# 1 plot together
# 2 plot apart

#colors can be set, below are defaults
#col.um = "white" (CG unmethylated)
#col.me = "black" (CG methylated)
#col.gme = "darkgrey" (GC methylated)
#col.gum = "lightgrey" (GC unmethylated)
}

################################################################################
##  BEGIN CODE
################################################################################
rc = function(x) { as.character(reverseComplement(as(x, "DNAString"))) }
drawmeCircles = function(x, y, cex=2, colormatrix="white")
{   #only takes first value of y and first length(x) values of colormatrix
	#this is internal function so not much input sanitation done
	#x should be vector of integers
	#y should be a single integer
	abline(h=y[1], lwd=2)
	xyc = cbind(x, y[1], colormatrix[1:length(x)])
	xyc = xyc[order(x),]
	points(as.numeric(xyc[,1]),as.numeric(xyc[,2]), cex = cex, pch=21, bg=as.character(xyc[,3]))
}
methylcircleplot = function(ref.seq, bis.seq = NULL, fwd.primer, rev.primer, 
	rev.comp = FALSE, size = 2, scaling = 1, reference=FALSE, NOME=0, noaxis=FALSE,
	col.um = "white", col.me = "black", col.gme = "darkgrey", col.gum = "lightgrey")
{
	#originally wanted to add reference sequence in monospace characters and have a closeness variable for plotting
	#text(0,1,"AACCCTTTTGGGGG", family="mono", adj=c(0,0))
	#	could not get font= or vfont= to work
	#	adj statement makes it left justified
	#	doesnt seem to be a way to 'stretch' font out so will not scale with plot dim changes
	#closeness scaling w/cex, 0 = max distance apart, 1 = closest together (scale with cex)
	#	under default resolution, cex=2 the bubbles are about 0.5 apart
	#	however, closeness is related to plot dim, for workaround use scaling variable
	#bis.seq can take either txt/fasta file, character vector, NULL (default)
	#	if NULL then begin interactive prompt for sequence
	
	if(!require("Biostrings", quietly = TRUE)) { stop("Missing Biostrings library. Please install via:
	source(\"http://bioconductor.org/biocLite.R\")
	biocLite(\"Biostrings\")\n\n") }
	if(!(NOME == 0 || NOME == 1 || NOME == 2)) { 
		stop("Invalid NOME flag, options are: 0 (disabled), 1 (plot together), 2 (plot seperate) ") 
	} #TRUE is 1
	
	#first try and open bis.seq as file
	#TODO: catch fail open
	if(length(bis.seq) == 1 && (grepl("\\.txt$", bis.seq) || grepl("\\.fasta$", bis.seq))) { 
		fasta = scan(bis.seq, what="character")
		if(is.na(table(grepl("^>", fasta))[2])) { bis.seq = fasta } #no '>' in file (assume all lines are diff clones)
		else{ #quick and dirty parser, can be optimized 
			bi = 0
			tmp=NULL;
			for(i in 1:length(fasta))
			{
				if(grepl("^>", fasta)[i]) { bis.seq[bi]=tmp; bi=bi+1; tmp=NULL; next;} #0 is thrown away
				tmp = paste(tmp, fasta[i], sep="")
			}
			bis.seq[bi] = tmp
		}
	}
	#interactive
	if(!length(bis.seq))
	{
		print("Paste in sequence, seperate different sequences with newlines, two newlines to end") #, press ctrl+z (on windows) or ctrl+d (on unix)
		bis.seq = scan(what="character")
		if(length(bis.seq) < 1) { stop("No bisulfite converted sequence entered") }
	}
	#cleanup: remove all \n, make capital letters (for pattern matching)
	bis.seq = toupper(gsub("[\n ]","",bis.seq)) #best way may be to remove \s+ w/PERL=TRUE
	ref.seq = toupper(gsub("[\n ]","",ref.seq))
	fwd.primer = toupper(fwd.primer) #should not have spaces or newlines
	rev.primer = toupper(rev.primer)
	if(rev.comp) { bis.seq = apply(as.matrix(bis.seq),1,rc) }
	
	#TODO check for [ACTGN] in bis.seq and ref.seq
	
	#strip primers at beginning and end of reference (if any)
	#supress warning and max() to catch case of blank primers
	fr = suppressWarnings(pairwiseAlignment(substr(ref.seq,1,max(1,nchar(fwd.primer))), fwd.primer, type="overlap"))
	rr = suppressWarnings(pairwiseAlignment(substr(ref.seq,nchar(ref.seq)-max(1,nchar(rev.primer)),nchar(ref.seq)), rev.primer, type="overlap"))
	if(nmatch(fr) > nchar(fwd.primer)*0.9 & nmatch(rr) > nchar(rev.primer)*0.9) 
	{
		fr = pairwiseAlignment(ref.seq, fwd.primer, type="overlap")
		rr = pairwiseAlignment(ref.seq, rev.primer, type="overlap")
		ref.seq = substr(ref.seq, end(pattern(fr))+1, start(pattern(rr))-1)
		#for testing
		#ref.seq.interest = substr(ref.seq, end(pattern(fr))+1, start(pattern(rr))-1)
		#fr2 = pairwiseAlignment(substr(ref.seq.interest,1,nchar(fwd.primer)), fwd.primer, type="overlap")
		#rr2 = pairwiseAlignment(substr(ref.seq.interest,nchar(ref.seq)-nchar(rev.primer),nchar(ref.seq)), rev.primer, type="overlap")
	} #could use score but score depends on sequence length
	
	#pairwiseAlignment only returns one match
	fpsa = suppressWarnings(pairwiseAlignment(bis.seq, fwd.primer, type="local-global")) #patternOverlap works better w/Ns
	rpsa = suppressWarnings(pairwiseAlignment(bis.seq, rev.primer, type="local-global")) #patternOverlap works better w/Ns
	bis.seq.interest = substr(bis.seq, end(pattern(fpsa))+1, start(pattern(rpsa))-1)
	if(fwd.primer=="" && rev.primer=="") { bis.seq.interest = bis.seq } #primers already removed
	
	#check lengths
	difflength = nchar(bis.seq.interest) - nchar(ref.seq)
	if(length(which(difflength != 0))) { 
		warning(paste("Sequence length difference detected in the following samples(diff): \n", 
		paste(paste(which(difflength !=0), difflength[which(difflength !=0)], sep="("), collapse=") "), ")", sep=""))
	}
	#difflength used later for plotting (not anymore!)
	
	pwa = pairwiseAlignment(bis.seq.interest, ref.seq, type="overlap") #maybe should be global here
	
	#identify CpGs in reference
	gcpos = "" #length("") is 1, length(NULL) is 0. Useful since plot(x,y) requires length(x) == length(y)
	cgsite = start(matchPattern("CG", ref.seq))
	if(NOME) { 
		gcsite = end(matchPattern("GC", ref.seq))
		ambiguousite = intersect(cgsite, gcsite)
		cgsite = cgsite[!cgsite %in% ambiguousite]
		gcsite = gcsite[!gcsite %in% ambiguousite]
		gcpos = gcsite/nchar(ref.seq)
	}
	mepos = cgsite/nchar(ref.seq) #Methylation positions
	
	#plotting
	plot.new()
	spacer = 0 #number of extra lines, used for spacing reference
	ylength = length(pwa)
	if(NOME == 2) { ylength = ylength*2 + 1 }
	if(reference) 
	{ 
		spacer = 2
		drawmeCircles(c(mepos,gcpos),1, cex = size, colormatrix=c(rep(col.um, length(mepos)),rep(col.gum,length(gcpos))))
	} 
	ypos = 1-seq(0,1, length.out=(ylength+spacer))*scaling
	if(!is.na(ypos[2]) && ypos[1]-ypos[2] < 0.005) { 
		ypos = 1-seq(0,by=0.005, length.out=(ylength+spacer))
		warning("Vertical distance between samples too low")
	}
	
	mesum = rep(0,length(pwa))
	nosum = rep(0,length(pwa))
	for(i in 1:length(pwa))
	{ #plot from top to bottom
		#curheight = ypos[i+spacer]
		abline(h=ypos[i+spacer], lty=3, col="grey10") #only see this line if sample failed
		
		#check if current sample failed 
		#if(i %in% difflength) { next; } #length must match (very conservative)
		if(nchar(bis.seq.interest[i]) == 0) { next; } #primer sequence not found

		mms = mismatchSummary(pwa[i])$subject #this is what holds the interesting data
		mms = mms[mms$Subject == "C" & mms$Pattern == "T",] #should also check for C-N, currently assume everything else methylated
		mestat = cgsite %in% mms$SubjectPosition #methylation status of CpG sites
		mesum[i] = sum(mestat)
		
		if(NOME) {
			nostat = gcsite %in% mms$SubjectPosition #methylation status of GpC sites
			nosum[i] = sum(nostat)
		} 
		if(NOME == 1) {
			drawmeCircles(c(mepos[mestat],mepos[!mestat],gcpos[nostat],gcpos[!nostat]),ypos[i+spacer], cex = size, 
				colormatrix=c(rep(col.um,length(mepos[mestat])), rep(col.me,length(mepos[!mestat])),
							  rep(col.gum,length(gcpos[nostat])), rep(col.gme,length(gcpos[!nostat]))))
		} else {
		drawmeCircles(c(mepos[mestat],mepos[!mestat]),ypos[i+spacer], cex = size, 
			colormatrix=c(rep(col.um, length(mepos[mestat])),rep(col.me,length(mepos[!mestat]))))
		}
	}
	
	if(NOME == 2) { 
		for(i in 1:length(pwa))
		{ #plot from top to bottom
			curheight = ypos[length(pwa)+i+spacer+1] #THIS LINE IS DIFFERENT
			abline(h=curheight, lty=3, col="grey10") #only see this line if sample failed
			
			#check if current experiment failed 
			if(nchar(bis.seq.interest[i]) == 0) { next; } #primer sequence not found
			
			mms = mismatchSummary(pwa[i])$subject #this is what holds the interesting data
			mms = mms[mms$Subject == "C" & mms$Pattern == "T",] #should also check for C-N
			nostat = gcsite %in% mms$SubjectPosition #methylation status of CpG sites

			drawmeCircles(c(gcpos[nostat],gcpos[!nostat]),curheight, cex = size, 
				colormatrix=c(rep(col.gum, length(gcpos[nostat])),rep(col.gme,length(gcpos[!nostat]))))			
		}
	}
	
	if(!noaxis) { 		#http://www.statmethods.net/advgraphs/axes.html
		xaxis.loc = seq(0,nchar(ref.seq),by=100)/nchar(ref.seq)
		xaxis.val = seq(0,nchar(ref.seq),by=100)[xaxis.loc < 0.97] #remove too close
		xaxis.loc = xaxis.loc[xaxis.loc < 0.97] #remove too close
		#can also use mtext()
		axis(1, at = c(xaxis.loc, 1), labels = c(xaxis.val, nchar(ref.seq))) 

		yaxis.val = 1:length(pwa)
		#yaxis.loc = NULL
		if(reference) {
			yaxis.val = c("ref", "", yaxis.val)
			
		} 
		if(NOME == 2) { yaxis.val = c(yaxis.val, "", 1:length(pwa)) }
		#yaxis.loc = seq(1,0,length.out=length(yaxis.val))*scaling
		axis(2, at=ypos, labels=yaxis.val, cex.axis=1)
	}
	if(NOME) { 
		title(main=paste("Bisulfite sequencing", 
		round((1-sum(mesum)/(length(cgsite)*length(pwa)))*100,2), "% methylated",
		"| NOME", round((1-sum(nosum)/(length(gcsite)*length(pwa)))*100,2), "% methylated"))
	}
	else { title(main=paste("Bisulfite sequencing", round((1-sum(mesum)/(length(cgsite)*length(pwa)))*100,2), "% methylated")) }
}

#plotting old way
#points(mepos[mestat],rep(curheight,length(mepos[mestat])), 
#	cex = size, pch=21, col="black", bg=col.um) #unmethylated
#points(mepos[!mestat],rep(curheight,length(mepos[!mestat])),
#	cex = size, pch=21, col="black", bg=col.me) #methylated, could replace w/pch=19
