#DNA methylation figure creator
#Author: Ying Wu daiyingw@gmail.com
#Current Update: October 2012
#License: GPLv3
#
#Step by step examples and directions can be found here:
# https://github.com/ying-w/bioinformatics-figures/methylcricleplot
# note to self, you cannot source https, find another way to source() from github
#
#Quickstart directions:
#Change ref.seq/fwd.primer/rev.primer/bis.seq in examplerun()
#Run the methylcircleplot function by typing examplerun()

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
size = 2 #size of methylation bubbles
scaling = 1 #proportion of figure to use for plotting (useful for adjusting closeness)
fwd.primer = "tgggttgaaatattgggtttattt"
rev.primer = "TATGGTTGGATATTTGGTTTTAG"
#if no primers, leave it as blank, ie ""

#bis.seq can take either txt/fasta file, folder, character vector, NULL (default)
#	if NULL then begin interactive prompt for sequence
#vector:
bis.seq = c(
"NNNNNNNNNNNNGAGCTCGGATCCACTAGTAACGGCCGCCAGTGTGCTGGAATTCGCCCTTCTAAAACCAAATATCCAACCATAAAAAAAATAAACACCCCAAACCAAAAACCTAATAAAAAAACAACTCACCTAAAAACAATTTCACCTCTCCCTCCCCAATCCCACCCTCTAACCTTAACCTCTAACCCCACCCCCTAAATAAATAAAAAAACAAAACTATAAAAATAAAAAAAACTAAAACAAACACTATCTACCTATATCTTCCAAACAAAAATTAAAAACATCTAAACAAAACAACCCTTAAAACAAAACCCCAAAAAACCTTCATTTTCAACCTTCAAAATCCTCTCACCCCTACCTTAAATCACCAAACACCTCACAAACCAATTACTCAAATACCCCATCACACCACAAAACCTATTAACACTACACCCTCTCAACCTATCCTTACAACCCACTCAACCATCCTCTAAAAAACCTAAAACATCCATTAAATATTCATATACCAATTATTTCTAAAACCACTCCTCAATTCTTACTTACCCACCCATCTAAAATCCTAAATATAAATAAACCCAATATTTCAACCCAAAGGGCGAATTCTGCAGATATCCATCACACTGGCGGCCGCTCGAGCATGCATCTAGAGGGCCCAATTCGCCCTATAGTGAGTCGTATTACAATTCACTGGCCGTCGTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGNAATAGCGAAGAGGCNGCACCGATCGCCTTCCCAACAGTTGCGCAGCCTGATGGCGATGGACGCGCCCTGNAGCGCGCATAAGCGCGCGGGTGTGTGTTACGCGCAGCGTGACCGCTACACTGCAGCNCCTAGCGCCGCTCTTTCGCTTCTTCCTTCTTCTCGCACGTCGCGGNTTCCCNTCAGCTCTAATCGGGGCTCCTTNNGGTCCNATTAGNNNTTANNNACTCGACCCAAAACTGATAGGTGANGTCNCGTANTNGGCATCGCCTGANAGACGTTTCGCCNTGACGNGNGTCNNNTNNNNNANANNGANNCTNGTCAANNGGACANNNCACCNANNNNNNTCTTNNTNNANGATTNCGATNGCNNTGTAANNGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
,"NNNNNNNNNNNNGAGCTCGGATCCACTAGTAACGGCCGCCAGTGTGCTGGAATTCGCCCTTCTAAAACCAAATATCCAACCATAAAAAAAATAAACACCCCAAACCAAAAACCTAATAAAAAAACAACTCACCTAAAAACAATTTCACCTCTCCCTCCCCAATCCCACCCTCTAACCTTAACCTCTAACCCCACCCCCTAAATAAATAAAAAAACAAAACTATAAAAATAAAAAAAACTAAAACAAACACTATCTACCTATATCTTCCAAACAAAAATTAAAAACATCTAAACAAAACAACCCTTAAAACAAAACCCCAAAAAACCTTCATTTTCAACCTTCAAAATCCTCTCACCCCTACCTTAAATCACCAAACACCTCACAAACCAATTACTCAAATACCCCATCACACCACAAAACCTATTAACACTACACCCTCTCAACCTATCCTTACAACCCACTCAACCATCCTCTAAAAAACCTAAAACATCCATTAAATATTCATATACCAATTATTTCTAAAACCACTCCTCAATTCTTACTTACCCACCCATCTAAAATCCTAAATATAAATAAACCCAATATTTCAACCCAAAGGGCGAATTCTGCAGATATCCATCACACTGGCGGCCGCTCGAGCATGCATCTAGAGGGCCCAATTCGCCCTATAGTGAGTCGTATTACAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGACGCNNNNNNGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTACACTTGCAGCGCCNANNNNNGCNCTTCGTTCTCTCTNNNCNNNNCNNNNNNNNNANNNNNNNNACANTANNNNANNNNNNNANAANANNNNNCNANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGNNNNNNNNNGNTNNNNNNNNCTCNNNNNNNNNCNNNNNNNNNNTNNNNNNNNNNNNNNTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
)
#file: (must end in .txt or .fasta)
#tip: to change directory go to: Tools | Change Working Dir of RStudio

#RUN THIS FUNCTION
methylcircleplot(ref.seq, bis.seq, fwd.primer, rev.primer, rev.comp, size, scaling)
dev.copy(png, "example-methylation-figure.png")
dev.off()
#dev.copy2pdf will keep resized figure dimensions

#There are 3 options for NOME
# 0 disabled (default)
# 1 plot together
# 2 plot apart

#colors can be set, below are defaults
#col.um = "white" (CG unmethylated)
#col.me = "black" (CG methylated)
#col.gum = "lightgrey" (GC unmethylated)
#col.gme = "darkgrey" (GC methylated)
}

################################################################################
##  BEGIN CODE
################################################################################
rc = function(x) { as.character(reverseComplement(as(x, "DNAString"))) }
drawmeCircles = function(x, y, cex=2, colormatrix="white")
{   #only takes first value of y and first length(x) values of colormatrix
	#expected to only be called by methylcircleplot() so no input sanitation done
	#x should be vector of integers
	#y should be a single integer
	abline(h=y[1], lwd=2)
	xyc = cbind(x, y[1], colormatrix[1:length(x)]) #catch when colormatrix length wrong
	xyc = xyc[order(x),]
	points(as.numeric(xyc[,1]),as.numeric(xyc[,2]), cex = cex, pch=21, bg=as.character(xyc[,3]))
}
methylcircleplot = function(ref.seq, bis.seq = NULL, fwd.primer = "", rev.primer = "", 
	rev.comp = FALSE, size = 2, scaling = 1, reference=FALSE, NOME=0, noaxis=FALSE,
	col.um = "white", col.me = "black", col.gme = "lightgreen", col.gum = "aliceblue",
	verbose = TRUE, sampleName=NULL, getAln = FALSE)
{
	#originally wanted to add reference sequence in monospace characters and have a closeness variable for plotting
	#text(0,1,"AACCCTTTTGGGGG", family="mono", adj=c(0,0))
	#	could not get font= or vfont= to work
	#	adj statement makes it left justified
	#	doesnt seem to be a way to 'stretch' font out so will not scale with plot dim changes
	#closeness scaling w/cex, 0 = max distance apart, 1 = closest together (scale with cex)
	#	under default resolution, cex=2 the bubbles are about 0.5 apart
	#	however, closeness is related to plot dim, for workaround use scaling variable
	#bis.seq can take either txt/fasta file, folder, character vector, NULL (default)
	#	if NULL then begin interactive prompt for sequence
	
	if(!require("Biostrings", quietly = TRUE)) { stop("Missing Biostrings library. Please install by entering:
	source(\"http://bioconductor.org/biocLite.R\")
	biocLite(\"Biostrings\")\n\n") }
	if(!is.numeric(NOME) || !(NOME == 0 || NOME == 1 || NOME == 2)) { 
		stop("Invalid NOME flag, options are: 0 (disabled), 1 (plot together), 2 (plot seperate) ") 
	} #TRUE is 1
	
	#first try and open bis.seq as file
	if(length(bis.seq) == 1) { #possibilities: sequence/file/folder
		if(grepl("\\.txt$", bis.seq) || grepl("\\.fasta$", bis.seq)) { #file
			if(verbose) { message("[info] Reading from file: ", bis.seq) }
			fasta = scan(bis.seq, what="character") #will throw error if file missing
			fa_header = grepl("^>", fasta)
			if(sum(fa_header) == 0) { bis.seq = fasta } #no '>' in file (assume all lines are diff clones)
			else{ #quick and dirty parser, can be optimized 
				bindex = 0
				tmp=NULL;
				for(i in 1:length(fasta)) {
					if(fa_header[i]) { bis.seq[bindex]=tmp; bindex=bindex+1; tmp=NULL; next;} #initial 0 is thrown away
					tmp = paste(tmp, fasta[i], sep="")
				}
				bis.seq[bindex] = tmp
				if(length(sampleName) < 1) { 
					sampleName = fasta[fa_header]
				}
			}
		} else if(length(list.files(bis.seq)) > 0) { #folder
			readlist = list.files(bis.seq)
			readlist = readlist[grepl("\\.txt$", readlist) || grepl("\\.fasta$", readlist)]
			if(length(sampleName) < 1) {
				sampleName = unlist(lapply(strsplit(readlist, "\\."), "[", 1)) #split by '.' take 1st element
			}
			readlist = paste(bis.seq, readlist, sep="/")
			bindex = 0
			
			for(f in readlist) {
				tmp = NULL
				bindex = bindex + 1
				if(verbose) { message("[info] Reading from file: ", f) }
				fasta = scan(f, what="character") #will throw error if file missing
				fa_header = grepl("^>", fasta)
				if(sum(fa_header) > 1) { message("[ATTN] More than one set of sequences found, only using first set") }
				if(sum(fa_header) == 0) { 
					for(i in 1:length(fasta)) { tmp = paste(tmp, fasta[i], sep="") }
				} else { #takes into account case where there is 1+ set of sequence
					for(i in 2:(c(which(fa_header),length(fasta)+1)[2]-1)) { tmp = paste(tmp, fasta[i], sep="") }
				}
				bis.seq[bindex] = tmp
			}
		}
	}
	#interactive
	if(!length(bis.seq)) {
		cat("Paste in sequence, seperate different sequences with newlines, two newlines to end\n") 
		#press ctrl+z (on windows) or ctrl+d (on unix)
		bis.seq = scan(what="character")
		if(length(bis.seq) < 1) { stop("No bisulfite converted sequence entered") }
	}
	
	#cleanup: remove all \n, make capital letters (for pattern matching)
	bis.seq = toupper(gsub("[\n ]","",bis.seq)) #best way may be to remove \s+ w/PERL=TRUE
	ref.seq = toupper(gsub("[\n ]","",ref.seq))
	fwd.primer = toupper(gsub("[\n ]","",fwd.primer))
	rev.primer = toupper(gsub("[\n ]","",rev.primer)) #remember, 5'->3' orientation
	if(rev.comp) { bis.seq = apply(as.matrix(bis.seq),1,rc)
	if(verbose) { message("Sequence has been reverse complemented\n") } }
	
	#check: [ACTGN] in bis.seq and ref.seq
	if(grepl("[^ACTGN]", ref.seq)) { warning("Non-nucleotide character found in reference") }
	if(all(grepl("[^ACTGN]", bis.seq))) { warning("Non-nucleotide character found in bisulfite sequence ", 
		paste(which(grepl("[^ACTGN]", bis.seq)), collapse=" ")) }
	if(grepl("[^ACTGN]", fwd.primer)) { warning("Non-nucleotide character found in foward primer") }
	if(grepl("[^ACTGN]", rev.primer)) { warning("Non-nucleotide character found in reverse primer") }
	
	#cleanup: strip primers at beginning and end of reference (if any)
	#supress warning and max() to catch case of blank primers
	fr = suppressWarnings(pairwiseAlignment(substr(ref.seq,1,max(1,nchar(fwd.primer))), fwd.primer, type="overlap"))
	rr = suppressWarnings(pairwiseAlignment(substr(ref.seq,nchar(ref.seq)-max(1,nchar(rev.primer)),nchar(ref.seq)), rev.primer, type="overlap"))
	
	#reverse primer match not found, try again
	if(nmatch(fr) > nchar(fwd.primer)*0.9 & !(nmatch(rr) > nchar(rev.primer)*0.9)) { 
		if(verbose) { message("[info] Forward primer found in reference but reverse primer not found\n",
		"-Trying reverse complement of reverse primer") }
		rr2 = suppressWarnings(pairwiseAlignment(substr(ref.seq,
			nchar(ref.seq)-max(1,nchar(rev.primer)),nchar(ref.seq)),
			rc(rev.primer), type="overlap"))
		if(nmatch(rr2) > nmatch(rr)) { #reverse complement match better
			if(verbose) { message("[info] Match found using reverse complement of rev.primer") }
			rev.primer = rc(rev.primer)
			if(verbose) { message("[ATTN] Reverse complementing reverse primer, please adjust rev.primer next time") }
		}
		if(nmatch(rr2) > nchar(rev.primer)*0.9) { rr = rr2
		} else { warning("Forward primer found in reference but reverse primer not found") }
	}
	if(!(nmatch(fr) > nchar(fwd.primer)*0.9) & nmatch(rr) > nchar(rev.primer)*0.9)
	{ warning("Reverse primer found in reference but forward primer not found") }
	
	#remove primers in reference
	if(nmatch(fr) > nchar(fwd.primer)*0.9 & nmatch(rr) > nchar(rev.primer)*0.9) { 
		#match again not constrained to beginning/end
		fr = pairwiseAlignment(ref.seq, fwd.primer, type="overlap")
		rr = pairwiseAlignment(ref.seq, rev.primer, type="overlap")
		ref.seq = substr(ref.seq, end(pattern(fr))+1, start(pattern(rr))-1)
		if(verbose) { message("[info] Primers removed from reference") }
		
		#TESTING to see what happens when primers already removed
		#ref.seq.interest = substr(ref.seq, end(pattern(fr))+1, start(pattern(rr))-1)
		#fr2 = pairwiseAlignment(substr(ref.seq.interest,1,nchar(fwd.primer)), fwd.primer, type="overlap")
		#rr2 = pairwiseAlignment(substr(ref.seq.interest,nchar(ref.seq)-nchar(rev.primer),nchar(ref.seq)), rev.primer, type="overlap")
		#could use score but score depends on sequence length
	} else { if(verbose) { message("[info] Reference sequence used as is") } }
	
	#pairwiseAlignment only returns one match
	fpsa = suppressWarnings(pairwiseAlignment(bis.seq, fwd.primer, type="local-global")) #patternOverlap works better w/Ns
	rpsa = suppressWarnings(pairwiseAlignment(bis.seq, rev.primer, type="local-global")) #patternOverlap works better w/Ns
	if(fwd.primer=="" && rev.primer=="") { 
		bis.seq.interest = bis.seq  #primers already removed
		if(verbose) { message("[info] Processing sequence without primers") }
	} else {
		bis.seq.interest = substr(bis.seq, end(pattern(fpsa))+1, start(pattern(rpsa))-1)
	}
	
	#check: user error in primer specification
	if(all(start(pattern(rpsa)) > end(pattern(fpsa))) & all(score(fpsa) > 0) & all(score(rpsa) > 0)) {
		message("[ATTN] reverse primer found before forward primer")
		warning("reverse primer found before forward primer, swapping primers")
		bis.seq.interest = substr(bis.seq, end(pattern(fpsa))+1, start(pattern(rpsa))-1)
	}
	if(any(score(fpsa) < 0) | all(score(rpsa) < 0)) {
		message("[info] Cannot find primers in the following sequences: ", appendLF=FALSE)
		message(paste(which(score(fpsa) < 0 | score(rpsa) < 0), collapse=" "))
		message("-Trying reverse complement of ALL bisulfite sequences")
		bis.seq.rc = apply(as.matrix(bis.seq),1,rc)
		fpsarc = suppressWarnings(pairwiseAlignment(bis.seq, fwd.primer, type="local-global")) #patternOverlap works better w/Ns
		rpsarc = suppressWarnings(pairwiseAlignment(bis.seq, rev.primer, type="local-global")) #patternOverlap works better w/Ns
		if(all(score(fpsarc) > score(fpsa)) & all(score(fpsarc) > 0)) { message("[info] Forward primer matches clone reverse complement better") }
		if(all(score(rpsarc) > score(rpsa)) & all(score(rpsarc) > 0)) { message("[info] Reverse primer matches clone reverse complement better") }
		if(all(score(fpsarc) > score(fpsa)) & all(score(fpsarc) > 0) & 
			all(score(rpsarc) > score(rpsa)) & all(score(rpsarc) > 0)) {
			bis.seq.interest = substr(bis.seq.rc, end(pattern(fpsarc))+1, start(pattern(rpsarc))-1)
			message("[info] Using reverse complement of clone since it matches primers better")
			message("[ATTN] set ", paste("rev.comp", !rev.comp,sep="="), " next time")
			warning(paste("Both forward and reverse primers match reverse complement of clone better"))
		}
		if(all(score(fpsarc) > score(fpsa)) & all(score(fpsarc) > 0) & 
			all(score(rpsarc) < score(rpsa)) & all(score(rpsa) > 0)) {
			#bis.seq.interest = substr(bis.seq.rc, end(pattern(fpsarc))+1, start(pattern(rpsarc))-1)
			message("[info] Forward primer matches clone reverse complement better but reverse primer does not")
		}		
		if(all(score(fpsarc) < score(fpsa)) & all(score(fpsa) > 0) & 
			all(score(rpsarc) > score(rpsa)) & all(score(rpsarc) > 0)) {
			#bis.seq.interest = substr(bis.seq.rc, end(pattern(fpsarc))+1, start(pattern(rpsarc))-1)
			message("[info] Reverse primer matches clone reverse complement better but forward primer does not\n")
		}
	}
	
	#fallback with artificial conversion and alignment
	if(any(score(fpsa) < 0 | score(rpsa) < 0)) {
		if(verbose) { message("[info] Cannot find sequence in clones using primers\n-falling back to artificial conversion") }
		refsa = suppressWarnings(pairwiseAlignment(gsub("C","T",bis.seq), gsub("C","T",ref.seq), type="local-global"))
		refsarc = suppressWarnings(pairwiseAlignment(gsub("C","T",apply(as.matrix(bis.seq),1,rc)), gsub("C","T",ref.seq), type="local-global"))
		if(mean(score(refsa)) > mean(score(refsarc)) & all(score(refsa) > 0)) { #fallback with refsa
			bis.seq.interest = substr(bis.seq, start(pattern(refsa)), end(pattern(refsa)))
			if(verbose) { message("[info] Aligning clone sequence with reference") }
		} else if(mean(score(refsarc)) > mean(score(refsa)) & all(score(refsarc) > 0)) { #fallback with refsarc
			bis.seq.interest = substr(apply(as.matrix(bis.seq),1,rc), start(pattern(refsarc)), end(pattern(refsarc)))
			if(verbose) { 
				message("[info] Aligning reverse complemented clone sequence with reference")
				message("[ATTN] set ", paste("rev.comp", !rev.comp,sep="="), " next time")
				#cat("-set "); cat(paste("rev.comp", !rev.comp,sep="=")); cat(" next time\n") 				
			}
		} else { stop("No good matches found in sample, please check primers and reference") }
	}
	
	if(verbose) { message("[info] ", paste(length(bis.seq.interest), "clones processed with average length of", 
		round(mean(nchar(bis.seq.interest)),2), collapse=" ")) }
	if(length(sampleName) < 1) { sampleName = 1:length(bis.seq.interest) }
	if(length(sampleName) != length(bis.seq.interest)) { 
		message("[ATTN] invalid number of sample names, using numbers isntead"); 
		sampleName = 1:length(bis.seq.interest)
	}
	
	#check lengths
	difflength = nchar(bis.seq.interest) - nchar(ref.seq)
	if(length(which(difflength != 0))) { 
		warning(paste("Sequence length difference (clone-ref) detected in the following samples(diff): \n", 
		paste(paste(which(difflength !=0), difflength[which(difflength !=0)], sep="("), collapse=") "), ")", sep=""))
	} #difflength used later for plotting (not anymore!)
	#use as.matrix(pairwiseAlignment(gsub("C","T",bis.seq.interest), gsub("C","T",ref.seq), type="local-global"))
	#and look for '-' to find which column is skipped over
	
	#type="global" has end gap penalty see Biostrings manual for more details
	#pwa = pairwiseAlignment(bis.seq.interest, ref.seq, type="overlap") #used for mismatchSummary
	pwa = pairwiseAlignment(bis.seq.interest, gsub("C","T",ref.seq), type="global") 
	if(all(score(pwa) <= 0)) { stop("No good matches found in sample, please check primers and reference") }
	
	#identify CpGs in reference
	gcpos = "" #length("") is 1, length(NULL) is 0. Useful since plot(x,y) requires length(x) == length(y)
	cgsite = start(matchPattern("CG", ref.seq))
	if(NOME) { 
		if(verbose) { message("[info] NOME mode enabled: also looking for GpCs") }
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
	if(reference) { 
		spacer = 2
		drawmeCircles(c(mepos,gcpos),1, cex = size, colormatrix=c(rep(col.um, length(mepos)),rep(col.gum,length(gcpos))))
	} 
	ypos = 1-seq(0,1, length.out=(ylength+spacer))*scaling
	if(!is.na(ypos[2]) && ypos[1]-ypos[2] < 0.005) { 
		ypos = 1-seq(0,by=0.005, length.out=(ylength+spacer))
		warning("[plot] Vertical distance between samples is low. Try plotting less samples next time")
	}
	
	mesum = rep(0,length(pwa))
	metotal = 0
	nosum = rep(0,length(pwa))
	nototal = 0
	for(i in 1:length(pwa))
	{ #plot from top to bottom
		#curheight = ypos[i+spacer]
		abline(h=ypos[i+spacer], lty=3, col="grey10") #only see this line if sample failed
		
		#check if current sample failed 
		#if(i %in% difflength) { next; } #length must match (very conservative)
		if(nchar(bis.seq.interest[i]) == 0 || score(pwa)[i] <= 0) {  #primer sequence not found
			message(paste("Clone", i, "skipped due to low quality match")); next; 
		} 

		#mms = mismatchSummary(pwa[i])$subject #this is what holds the interesting data
		#mms = mms[mms$Subject == "C" & mms$Pattern == "T",] #should also check for C-N, currently assume everything else methylated
		#mestat = cgsite %in% mms$SubjectPosition #methylation status of CpG sites
		#above code does not take into account seq errors / indels
		mestat = as.matrix(pwa[i])[,cgsite] == "T"
		umestat = as.matrix(pwa[i])[,cgsite] == "C"
		mesum[i] = sum(mestat)
		metotal = metotal+sum(mestat)+sum(umestat)
		
		if(NOME) {
			#nostat = gcsite %in% mms$SubjectPosition #methylation status of GpC sites
			nostat = as.matrix(pwa[i])[,cgsite] == "T"
			unostat = as.matrix(pwa[i])[,cgsite] == "C"
			nosum[i] = sum(nostat)
			nototal = nototal+sum(nostat)+sum(unostat)
		} 
		if(NOME == 1) {
			if(verbose) { message("[plot] Plotting GpCs sites in line with CpGs sites") }
			drawmeCircles(c(mepos[mestat],mepos[umestat],gcpos[nostat],gcpos[unostat]),ypos[i+spacer], cex = size, 
				colormatrix=c(rep(col.um,length(mepos[mestat])), rep(col.me,length(mepos[umestat])),
							  rep(col.gum,length(gcpos[nostat])), rep(col.gme,length(gcpos[unostat]))))
		} else {
		drawmeCircles(c(mepos[mestat],mepos[umestat]),ypos[i+spacer], cex = size, 
			colormatrix=c(rep(col.um, length(mepos[mestat])),rep(col.me,length(mepos[umestat]))))
		}
	}
	
	if(NOME == 2) { 
		if(verbose) { message("[plot] Plotting GpCs sites under CpGs sites") }
		for(i in 1:length(pwa))
		{ #plot from top to bottom
			curheight = ypos[length(pwa)+i+spacer+1] #THIS LINE IS DIFFERENT
			abline(h=curheight, lty=3, col="grey10") #only see this line if sample failed
			
			#check if current experiment failed 
			if(nchar(bis.seq.interest[i]) == 0) { next; } #primer sequence not found
			
			#mms = mismatchSummary(pwa[i])$subject #this is what holds the interesting data
			#mms = mms[mms$Subject == "C" & mms$Pattern == "T",] #should also check for C-N
			#nostat = gcsite %in% mms$SubjectPosition #methylation status of CpG sites
			nostat = as.matrix(pwa[i])[,cgsite] == "T"
			unostat = as.matrix(pwa[i])[,cgsite] == "C"
			
			drawmeCircles(c(gcpos[nostat],gcpos[unostat]),curheight, cex = size, 
				colormatrix=c(rep(col.gum, length(gcpos[nostat])),rep(col.gme,length(gcpos[unostat]))))			
		}
	}
	
	if(!noaxis) { 		#http://www.statmethods.net/advgraphs/axes.html
		#TODO: axis rescaling with plot resize
		xaxis.loc = seq(0,nchar(ref.seq),by=50)/nchar(ref.seq)
		xaxis.val = seq(0,nchar(ref.seq),by=50)[xaxis.loc < 0.97] #remove too close
		xaxis.loc = xaxis.loc[xaxis.loc < 0.97] #remove too close
		#prettylucky = xaxis.val %in% pretty(xaxis.val, n=10)
		#can also use mtext()
		axis(1, at = c(xaxis.loc, 1), labels = c(xaxis.val, nchar(ref.seq))) 

		yaxis.val = sampleName
		#yaxis.loc = NULL
		if(reference) {
			yaxis.val = c("ref", "", yaxis.val)
			
		} 
		if(NOME == 2) { yaxis.val = c(yaxis.val, "", sampleName) }
		#yaxis.loc = seq(1,0,length.out=length(yaxis.val))*scaling
		axis(2, at=ypos, labels=yaxis.val, cex.axis=0.8)
	}
	if(NOME) { 
		title(main=paste("Result:", 
		round((1-sum(mesum)/(metotal))*100,2), "% methylated",
		"| NOME", round((1-sum(nosum)/(nototal))*100,2), "% methylated"))
	}
	else { title(main=paste("Result:", round((1-sum(mesum)/(metotal))*100,2), "% methylated")) }
	
	if(verbose) { message("[info] Done!") }
	if(getAln) { pwa }
}

#plotting old way
#points(mepos[mestat],rep(curheight,length(mepos[mestat])), 
#	cex = size, pch=21, col="black", bg=col.um) #unmethylated
#points(mepos[!mestat],rep(curheight,length(mepos[!mestat])),
#	cex = size, pch=21, col="black", bg=col.me) #methylated, could replace w/pch=19

## INDEL TEST CASES
## note that score cutoffs are much more stringent after changing to type="global"
# ref.seq = "GCGACGTTATTACGGGA"
# fwd.primer = ""
# rev.primer = ""
# bis.seq = c("GTGAGTTATTACGGGA", "GTGATGTTATTATGGGA", "GTGATGTTATTAGGGA", "GTGATGTTATTAGGGGA", "GTGACGTTATTACGGGAC", "GTGACGTTTATTACGGGA")
# methylcircleplot(ref.seq, bis.seq, fwd.primer, rev.primer, scaling=1, reference=TRUE, col.um="white", col.me= "black", col.gme="#00FA9A", col.gum ="white")
# pwa = pairwiseAlignment(bis.seq, ref.seq, type="overlap")
## see biostring reference for complete set of functions
# REF:	GCGACGTTATTACG
# CLONE1:	GTGA GTTATTACGGGA	(missing first CG)
# CLONE2:	GTGATGTTATTATGGGA	(everything is unme)
# CLONE3:	GTGATGTTATTA GGGA	(missing last CG)
# CLONE4:	GTGATGTTATTAGGGGA	(mismatch at last CG)
# CLONE5:	GTGACGTTATTACGGGAC	(extra base at end)
# CLONE6:	GTGACGTTTATTACGGGA	(extra base in middle)
