#DNA methylation figure creator
#Author: Ying Wu daiyingw@usc.edu
#Current Update: Nov 2012
#License: GPLv3
#
#Step by step examples and directions can be found here:
# https://github.com/ying-w/bioinformatics-figures/tree/master/methylcircleplot
# Must download local copy for R, since you cannot source() https
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

#colors can be set, below are OLD defaults
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
	verbose = TRUE, showNumUnconverted = FALSE, cloneName=NULL, cloneOrder = NULL, getAln = FALSE)
{
	## Unimplemented ideas
	#wanted to add reference sequence in monospace characters
	#text(0,1,"AACCCTTTTGGGGG", family="mono", adj=c(0,0))
	#	could not get font= or vfont= to work
	#	adj statement makes it left justified
	#	doesnt seem to be a way to 'stretch' font out so cannot scale with plot dim changes
	#closeness variable for plotting
	#	scaling w/cex, 0 = max distance apart, 1 = closest together
	#	under default resolution, cex=2 the bubbles are about 0.5 apart
	#	however, closeness is related to plot dim, for workaround use scaling variable
	#include export highlighted HTML option since word can read HTML
	
	if(!require("Biostrings", quietly = TRUE)) { stop("Missing Biostrings library. Please install by entering:
	source(\"http://bioconductor.org/biocLite.R\")
	biocLite(\"Biostrings\")\n\n") }
	if(!is.numeric(NOME) || !(NOME == 0 || NOME == 1 || NOME == 2)) { 
		stop("Invalid NOME flag, options are: 0 (disabled), 1 (plot together), 2 (plot seperate) ") 
	} #TRUE is 1
	
	############################################################################
	## read in bis.seq
	############################################################################
	#bis.seq can take either txt/fasta file, folder, character vector, NULL (default)
	#	if NULL then begin interactive prompt for sequence
	#TODO: if end in fasta/txt read it in as file to account for dir(folder) case
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
				if(length(cloneName) < 1) { 
					cloneName = fasta[fa_header]
				}
			} #read in file as fasta
		} else if(length(list.files(bis.seq)) > 0) { #folder
			readlist = list.files(bis.seq)
			readlist = readlist[grepl("\\.txt$", readlist) | grepl("\\.fasta$", readlist) | grepl("\\.seq$", readlist)]
			if(length(readlist) < 1) { stop("No .fasta and .txt files found in folder ", bis.seq) }
			if(length(cloneName) < 1) {
				cloneName = unlist(lapply(strsplit(readlist, "\\."), "[", 1)) #split by '.' take 1st element
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
			} #for each file in folder
		} #if folder
	} else if(!length(bis.seq)) {	#interactive
		cat("Paste in sequence, seperate different sequences with newlines, two newlines to end\n") 
		#press ctrl+z (on windows) or ctrl+d (on unix)
		bis.seq = scan(what="character")
		if(length(bis.seq) < 1) { stop("No bisulfite converted sequence entered") }
	} else if(verbose) { message("[info] Assume bis.seq contains sequence") }
		
	#cleanup: remove all \n, make capital letters (for pattern matching)
	bis.seq = toupper(gsub("[\n ]","",bis.seq)) #best way may be to remove \s+ w/PERL=TRUE
	ref.seq = toupper(gsub("[\n ]","",ref.seq))
	fwd.primer = toupper(gsub("[\n ]","",fwd.primer))
	rev.primer = toupper(gsub("[\n ]","",rev.primer)) #remember, 5'->3' orientation
	if(rev.comp) { bis.seq = apply(as.matrix(bis.seq),1,rc)
		if(verbose) { message("[info] Sequence has been reverse complemented") } 
	}
	
	#check: input
	if(length(ref.seq) != 1) { stop("Please specify reference sequence for ref.seq parameter") }
	if(length(fwd.primer) > 1 || length(rev.primer) > 1) { stop("Only one set of primers allowed, check primers") }
	allseq = grepl("[^ACTGN]", (unlist(strsplit(bis.seq, split=""))))
	if(sum(ss)/length(ss) > 0.1) { stop("Too many non-nucleotide characters found in bisulfite sequence") }
	rm(allseq)
	#check: [ACTGN] in bis.seq and ref.seq
	if(grepl("[^ACTGN]", ref.seq)) { stop("Non-nucleotide character found in reference") }
	if(any(grepl("[^ACTGN]", bis.seq))) { warning("Non-nucleotide character found in bisulfite sequence ", 
		paste(which(grepl("[^ACTGN]", bis.seq)), collapse=" ")) }
	if(grepl("[^ACTGN]", fwd.primer)) { warning("Non-nucleotide character found in foward primer") }
	if(grepl("[^ACTGN]", rev.primer)) { warning("Non-nucleotide character found in reverse primer") }
	
	#check: primers, change NULL primers into "", length is different
	if(length(fwd.primer) == 0) { fwd.primer = "" } #if NULL set to blank
	if(length(rev.primer) == 0) { rev.primer = "" } #if NULL set to blank
	if((fwd.primer == "" || rev.primer == "") && fwd.primer != rev.primer) {
		message("[ATTN] Only one primer specified, check fwd.primer and rev.primer")
	}
	
	############################################################################
	## reference processing
	############################################################################	
	#check: reference for primer sequence at beginning and end 
	#supress warning and max() to catch case of blank primers
	#artificially convert ref and primer in case primer is bisulfite sensitive
	#maybe make substr longer
	fr = suppressWarnings(pairwiseAlignment(substr(gsub("C","T",ref.seq),1,max(1,nchar(fwd.primer))), 
		gsub("C","T",fwd.primer), type="global"))
	rr = suppressWarnings(pairwiseAlignment(substr(gsub("C","T",ref.seq),nchar(ref.seq)-
		max(1,nchar(rev.primer)),nchar(ref.seq)), gsub("C","T",rev.primer), type="global"))
	
	#reverse primer match not found, try again
	if(nmatch(fr) > nchar(fwd.primer)*0.9 & !(nmatch(rr) > nchar(rev.primer)*0.9)) { 
		if(verbose) { message("[info] Forward primer found in reference but reverse primer not found\n",
		"-Trying reverse complement of reverse primer") }
		rr2 = suppressWarnings(pairwiseAlignment(substr(gsub("C","T",ref.seq),
			nchar(ref.seq)-max(1,nchar(rev.primer)),nchar(ref.seq)),
			gsub("C","T",rc(rev.primer)), type="global"))
		if(nmatch(rr2) > nmatch(rr)) { #reverse complement match better
			if(verbose) { message("[info] Better match found using reverse complement of rev.primer") }
		} else { message("[info] rev.primer still not found") }
		if(nmatch(rr2) > nchar(rev.primer)*0.9) { #reverse complement matches good enough
			rev.primer = rc(rev.primer)
			if(verbose) { message("[ATTN] Reverse complementing reverse primer, please adjust rev.primer next time") }
			rr = rr2
		} else { warning("Forward primer found in reference but reverse primer not found") }		
	} #fwd found, rev not found
	#Could only find reverse primer
	if(!(nmatch(fr) > nchar(fwd.primer)*0.9) & nmatch(rr) > nchar(rev.primer)*0.9)
	{ warning("Reverse primer found in reference but forward primer not found") }
	
	#cleanup: remove primers in reference
	if(nmatch(fr) > nchar(fwd.primer)*0.9 & nmatch(rr) > nchar(rev.primer)*0.9) { 
		#realignment not constrained to start/end of ref.seq
		fr = pairwiseAlignment(ref.seq, fwd.primer, type="local-global")
		rr = pairwiseAlignment(ref.seq, rev.primer, type="local-global")
		ref.seq = substr(ref.seq, end(pattern(fr))+1, start(pattern(rr))-1)
		if(verbose) { message("[info] Primers removed from reference") }
		
		#TESTING to see what happens when primers already removed
		#ref.seq.interest = substr(ref.seq, end(pattern(fr))+1, start(pattern(rr))-1)
		#fr2 = pairwiseAlignment(substr(ref.seq.interest,1,nchar(fwd.primer)), fwd.primer, type="overlap")
		#rr2 = pairwiseAlignment(substr(ref.seq.interest,nchar(ref.seq)-nchar(rev.primer),nchar(ref.seq)), rev.primer, type="overlap")
		#could use score but score depends on sequence length
	} else { if(verbose) { message("[info] Reference sequence used as is (no primers removed)") } }
	rm(rr, fr)
	
	############################################################################
	## bisulfite sequence processing primers
	############################################################################	
	#Extract regions of interest from bis.seq
	#remember: pairwiseAlignment only returns top match
	fpsa = suppressWarnings(pairwiseAlignment(bis.seq, fwd.primer, type="local-global")) #patternOverlap works better w/Ns
	rpsa = suppressWarnings(pairwiseAlignment(bis.seq, rev.primer, type="local-global")) #patternOverlap works better w/Ns
	if(fwd.primer=="" && rev.primer=="") { #will give score==0 pairwiseAlignment
		if(verbose) { message("[info] Processing sequence without primers") }
		if(all(nchar(bis.seq) <= nchar(ref.seq))) { #partial matches only
			bis.seq.interest = bis.seq  #primers already removed
		} else { 
			if(verbose) { message("[ATTN] Found sequence longer than reference, please specify primers next time") }
			#this will trigger artificial conversion + alignment
			fpsa@score = rep(-1, length(fpsa))
			rpsa@score = rep(-1, length(rpsa))
		} 
		#TODO: test ^
		#if(!all(nchar(bis.seq) == nchar(ref.seq))) { stop("Bisulfite sequence not the same length as reference, please specify primers") }
	} else {
		bis.seq.interest = substr(bis.seq, end(pattern(fpsa))+1, start(pattern(rpsa))-1)
	}
	
	#check: user error in primer specification
	if(all(start(pattern(rpsa)) < end(pattern(fpsa))) & all(score(fpsa) > 0) & all(score(rpsa) > 0)) {
		message("[ATTN] reverse primer found before forward primer")
		warning("reverse primer found before forward primer, swapping primers")
		bis.seq.interest = substr(bis.seq, end(pattern(rpsa))+1, start(pattern(fpsa))-1)
	}
	if(fwd.primer != "" && rev.primer != "" && (any(score(fpsa) < 0) || any(score(rpsa) < 0))) {
		#TODO: selective reverse complement
		message("[info] Cannot find primers in the following sequences: ", appendLF=FALSE)
		message(paste(which(score(fpsa) < 0 | score(rpsa) < 0), collapse=" "))
		message("-Trying reverse complement of ALL bisulfite sequences")
		bis.seq.rc = apply(as.matrix(bis.seq),1,rc)
		fpsarc = suppressWarnings(pairwiseAlignment(bis.seq.rc, fwd.primer, type="local-global")) #patternOverlap works better w/Ns
		rpsarc = suppressWarnings(pairwiseAlignment(bis.seq.rc, rev.primer, type="local-global")) #patternOverlap works better w/Ns
		if(all(score(fpsarc) > score(fpsa)) & all(score(fpsarc) > 0)) { message("[info] Forward primer matches clone reverse complement better") }
		if(all(score(rpsarc) > score(rpsa)) & all(score(rpsarc) > 0)) { message("[info] Reverse primer matches clone reverse complement better") }
		if(all(score(fpsarc) > score(fpsa)) & all(score(fpsarc) > 0) & 
			all(score(rpsarc) > score(rpsa)) & all(score(rpsarc) > 0)) { #both better
			bis.seq.interest = substr(bis.seq.rc, end(pattern(fpsarc))+1, start(pattern(rpsarc))-1)
			if(verbose) { 
				message("[info] Using reverse complement of clone since it matches primers better")
				message("[ATTN] Set ", paste("rev.comp", !rev.comp,sep="="), " next time")
			}
			warning(paste("Both forward and reverse primers match reverse complement of clone better"))
			fpsa = fpsarc
			rpsa = rpsarc
		}
		if(all(score(fpsarc) > score(fpsa)) & all(score(fpsarc) > 0) & 
			all(score(rpsarc) < score(rpsa)) & all(score(rpsa) > 0)) { #fwd better
			#bis.seq.interest = substr(bis.seq.rc, end(pattern(fpsarc))+1, start(pattern(rpsarc))-1)
			rpsarc2 = suppressWarnings(pairwiseAlignment(bis.seq.rc, rc(rev.primer), type="local-global"))
			if(all(score(rpsarc2) > 0)) { 
				if(verbose) {
					message("[info] Using reverse complement of clone with reverse complement of reverse primer")
					message("[ATTN] Set ", paste("rev.comp", !rev.comp,sep="="), " next time")
					message("[ATTN] Change reverse primer into 5'->3' orientation (reverse complement) next time")
				}
				rev.primer = rc(rev.primer)
				bis.seq.interest = substr(bis.seq.rc, end(pattern(fpsarc))+1, start(pattern(rpsarc2))-1)
				fpsa = fpsarc
				rpsa = rpsarc2
			} else { message("[info] Forward primer matches clone reverse complement better but reverse primer does not") }
		}
		if(all(score(fpsarc) < score(fpsa)) & all(score(fpsa) > 0) & 
			all(score(rpsarc) > score(rpsa)) & all(score(rpsarc) > 0)) { #rev better
			#bis.seq.interest = substr(bis.seq.rc, end(pattern(fpsarc))+1, start(pattern(rpsarc))-1)
			message("[info] Reverse primer matches clone reverse complement better but forward primer does not\n")
		}
	} #if primer not found, try to switch primer around
	
	#primer still not found, ignore primer and 
	#fallback with artificial conversion and alignment
	if(any(score(fpsa) < 0 || score(rpsa) < 0)) {
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
			warning(paste("Both forward and reverse primers match reverse complement of clone better"))
		} else { stop("No good matches found in bis.seq, please check primers and reference sequence") }
	}
	
	if(verbose) { message("[info] ", paste(length(bis.seq.interest), "clones processed with average length of", 
		round(mean(nchar(bis.seq.interest)),2), collapse=" ")) }
		
	#Assign sample names and sample order
	if(length(cloneName) < 1) { cloneName = 1:length(bis.seq.interest) }
	if(length(cloneName) != length(bis.seq.interest)) { 
		message("[ATTN] invalid number of sample names, using numbers instead"); 
		cloneName = 1:length(bis.seq.interest)
	}
	if(length(cloneOrder) == length(bis.seq.interest) && all(1:length(bis.seq.interest) %in% cloneOrder) ) {
		bis.seq.interest = bis.seq.interest[cloneOrder]
		cloneName = cloneName[cloneOrder]
		if(verbose) { message("[info] Reordering rows") }
	} else if(length(cloneOrder) > 0) { message("[ATTN] cloneOrder ignored due to incorrect input") }
	
	############################################################################
	## bisulfite sequence alignment and mC site identification
	############################################################################	
	#Align bisulfite sequence to reference
	#type="global" has end gap penalty see Biostrings manual for more details
	#pwa = pairwiseAlignment(bis.seq.interest, ref.seq, type="overlap") #used for mismatchSummary
	if(all(nchar(bis.seq.interest) < nchar(ref.seq) * 0.9)) {
		#note the switch in match type here, bis.seq.interest is expected to be shorter
		pwa = pairwiseAlignment(bis.seq.interest, gsub("C","T",ref.seq), type="global-local")
	} else {
		pwa = pairwiseAlignment(bis.seq.interest, gsub("C","T",ref.seq), type="global")
	}
	badpwa = score(pwa) <= 0
	if(all(badpwa)) { stop("No good matches found in sample, please check primers and reference") }
	
	if(any(badpwa)) {  #look at reverse complement of low scoring alignments
		pwarc = pairwiseAlignment(apply(as.matrix(bis.seq.interest[badpwa]),1,rc), gsub("C","T",ref.seq), type="global") 
		for(i in which(badpwa)[score(pwarc) > score(pwa)[badpwa]]) { 
			message("[info] Clone ", cloneName[i], " replaced with reverse complement") 
		}
		#pwa[badpwa][score(pwarc) > score(pwa)[badpwa]] = pwarc[score(pwarc) > score(pwa)[badpwa]] #not allowed
		bis.seq.interest[badpwa][score(pwarc) > score(pwa)[badpwa]] = 
			apply(as.matrix(bis.seq.interest[badpwa][score(pwarc) > score(pwa)[badpwa]]),1,rc)
		pwa = pairwiseAlignment(bis.seq.interest, gsub("C","T",ref.seq), type="global")
	}

	#check lengths
	difflength = nchar(bis.seq.interest) - nchar(ref.seq)
	if(length(which(difflength != 0))) { 
		warning(paste("Sequence length difference (clone-ref) detected in the following samples(diff): \n", 
		paste(paste(which(difflength !=0), difflength[which(difflength !=0)], sep="("), collapse=") "), ")", sep=""))
	} #difflength used later for plotting (not anymore!)
	#use as.matrix(pairwiseAlignment(gsub("C","T",bis.seq.interest), gsub("C","T",ref.seq), type="local-global"))
	#and look for '-' to find which column is skipped over
		
	#identify CpGs/GpCs in reference
	gcpos = "" #length("") is 1, length(NULL) is 0. Useful since plot(x,y) requires length(x) == length(y)
	csite = start(matchPattern("C", ref.seq)) #used to detect incomplete conversion
	cgsite = start(matchPattern("CG", ref.seq))
	csite = csite[!csite %in% cgsite]
	if(NOME) { 
		if(verbose) { message("[info] NOME mode enabled: also looking for GpCs") }
		gcsite = end(matchPattern("GC", ref.seq))
		csite = csite[!csite %in% gcsite]
		ambiguousite = intersect(cgsite, gcsite)
		cgsite = cgsite[!cgsite %in% ambiguousite]
		gcsite = gcsite[!gcsite %in% ambiguousite]
		gcpos = gcsite/nchar(ref.seq)
		if(length(gcsite) < 1) { message("[ATTN] No GCs found, NOMe plotting disabled") }
	}
	mepos = cgsite/nchar(ref.seq) #Methylation positions
	
	#Error checking on CG and C sites
	if(length(cgsite) < 1) { stop("No CpGs found in reference, please check ref.seq") }
	if(length(csite) < 1) { 
		warning("No non-CG cysteines found in reference, check that non-bisulfite converted sequence used for reference")
		message("[ATTN] No useable cysteines found to check for incomplete bisulfite conversion")
	}
	
	############################################################################
	## plotting
	############################################################################
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

	if(NOME == 1 && verbose) { message("[plot] Plotting GpCs sites in line with CpGs sites") } 
	
	mesum = rep(0,length(pwa))
	metotal = 0
	nosum = rep(0,length(pwa))
	nototal = 0
	for(i in 1:length(pwa))
	{ #plot from top to bottom
		#curheight = ypos[i+spacer]
		abline(h=ypos[i+spacer], lty=3, col="grey10") #only see this line if sample failed
		
		#check for incomplete conversion
		convertedC = as.matrix(pwa[i])[,csite] == "T"
		unconvertedC = as.matrix(pwa[i])[,csite] == "C"
		if(sum(unconvertedC)/max(sum(convertedC)+sum(unconvertedC),1) > 0.9) { #90%+ unconverted
			score(pwa)[i] = 0
			message("[ATTN] Clone ", cloneName[i], " will be skipped due to incomplete bisulfite conversion")
		} else if (sum(unconvertedC) > round(length(csite)/10)) { #10%+ unconverted
			message("[ATTN] Clone ", cloneName[i], " has ", sum(unconvertedC), "/", 
			sum(convertedC)+sum(unconvertedC), " unconverted Cs, check bisulfite conversion") 
		} else if (showNumUnconverted && sum(unconvertedC)) { #any unconverted
			message("[ATTN] Clone ", cloneName[i], " has ", sum(unconvertedC), "/", 
			sum(convertedC)+sum(unconvertedC), " unconverted Cs, check bisulfite conversion") 		
		}

		#check: if current sample failed in matching
		#if(i %in% difflength) { next; } #length must match (very conservative)
		if(nchar(bis.seq.interest[i]) == 0 || score(pwa)[i] <= 0) {  #primer sequence not found
			message("[ATTN] Clone ", cloneName[i], " skipped due to low quality match"); next; 
		} 		
		
		#mms = mismatchSummary(pwa[i])$subject #this is what holds the interesting data
		#mms = mms[mms$Subject == "C" & mms$Pattern == "T",] #should also check for C-N, currently assume everything else methylated
		#mestat = cgsite %in% mms$SubjectPosition #methylation status of CpG sites
		#above code does not take into account seq errors / indels
		umestat = as.matrix(pwa[i])[,cgsite] == "T"
		mestat = as.matrix(pwa[i])[,cgsite] == "C"
		mesum[i] = sum(mestat)
		metotal = metotal+sum(mestat)+sum(umestat)
		
		#catch non T/C csites and cgsites (1+)
		if(showNumUnconverted && length(cgsite) - (sum(umestat)+sum(mestat)) > 0) {
			message("[info] Sequencing / alignment errors: ", 
				length(cgsite)-sum(umestat)-sum(mestat), " non T/C found at CpGs in ", cloneName[i], appendLF=FALSE)
			print(table(as.matrix(pwa[i])[,csite])) #DEBUG 
		}
		if(showNumUnconverted && length(csite) - (sum(unconvertedC)+sum(convertedC)) > 0) {
			message("[info] Sequencing / alignment errors: ", 
				length(csite)-sum(unconvertedC)-sum(convertedC), " non T/C found at Cs in ", cloneName[i], appendLF=FALSE)
			print(table(as.matrix(pwa[i])[,csite])) #DEBUG
		}
		
		if(NOME) {
			#nostat = gcsite %in% mms$SubjectPosition #methylation status of GpC sites
			unostat = as.matrix(pwa[i])[,gcsite] == "T"
			nostat = as.matrix(pwa[i])[,gcsite] == "C"
			nosum[i] = sum(nostat)
			nototal = nototal+sum(nostat)+sum(unostat)
			
			if(showNumUnconverted && length(gcsite) - (sum(unostat)+sum(nostat)) > 0) {
				message("[info] Sequencing / alignment errors: ", 
					length(gcsite)-sum(nostat)-sum(unostat), " non T/C found at GpCs in ", cloneName[i], appendLF=FALSE)
				print(table(as.matrix(pwa[i])[,csite])) #DEBUG 		
			}
		} 
		if(NOME == 1) {
			drawmeCircles(c(mepos[umestat],mepos[mestat],gcpos[unostat],gcpos[nostat]),ypos[i+spacer], cex = size, 
				colormatrix=c(rep(col.um,length(mepos[umestat])), rep(col.me,length(mepos[mestat])),
							  rep(col.gum,length(gcpos[unostat])), rep(col.gme,length(gcpos[nostat]))))
		} else {
			drawmeCircles(c(mepos[umestat],mepos[mestat]),ypos[i+spacer], cex = size, 
				colormatrix=c(rep(col.um, length(mepos[umestat])),rep(col.me,length(mepos[mestat]))))
		}
	} #for each alignment
	
	if(NOME == 2) { 
		if(verbose) { message("[plot] Plotting GpCs sites under CpGs sites") }
		for(i in 1:length(pwa))
		{ #plot from top to bottom
			curheight = ypos[length(pwa)+i+spacer+1] #THIS LINE IS DIFFERENT
			abline(h=curheight, lty=3, col="grey10") #only see this line if sample failed
			
			#check if current experiment failed (see above)
			if(nchar(bis.seq.interest[i]) == 0 || score(pwa)[i] <= 0) {  #primer sequence not found
				next; 
			} 
			
			#mms = mismatchSummary(pwa[i])$subject #this is what holds the interesting data
			#mms = mms[mms$Subject == "C" & mms$Pattern == "T",] #should also check for C-N
			#nostat = gcsite %in% mms$SubjectPosition #methylation status of CpG sites
			unostat = as.matrix(pwa[i])[,gcsite] == "T"
			nostat = as.matrix(pwa[i])[,gcsite] == "C"
			
			drawmeCircles(c(gcpos[unostat],gcpos[nostat]),curheight, cex = size, 
				colormatrix=c(rep(col.gum, length(gcpos[unostat])),rep(col.gme,length(gcpos[nostat]))))			
		}
	} #if NOME=2
	
	if(!noaxis) { 		#http://www.statmethods.net/advgraphs/axes.html
		#TODO: axis rescaling with plot resize
		xaxis.loc = seq(0,nchar(ref.seq),by=50)/nchar(ref.seq)
		xaxis.val = seq(0,nchar(ref.seq),by=50)[xaxis.loc < 0.97] #remove too close
		xaxis.loc = xaxis.loc[xaxis.loc < 0.97] #remove too close
		#prettylucky = xaxis.val %in% pretty(xaxis.val, n=10)
		axis(1, at = c(xaxis.loc, 1), labels = c(xaxis.val, nchar(ref.seq))) #can also use mtext()

		yaxis.val = cloneName
		#yaxis.loc = NULL
		if(reference) {
			yaxis.val = c("ref", "", yaxis.val)
			
		} 
		if(NOME == 2) { yaxis.val = c(yaxis.val, "", cloneName) }
		#yaxis.loc = seq(1,0,length.out=length(yaxis.val))*scaling
		#TODO: dynamic scaling of cex based on nchar of yaxis.val
		axis(2, at=ypos, labels=yaxis.val, cex.axis=0.8, las=2) 
	}
	if(NOME) { 
		title(main=paste("Result:", 
		round((sum(mesum)/(metotal))*100,2), "% methylated",
		"| NOME", round((sum(nosum)/(nototal))*100,2), "% methylated"))
	}
	else { title(main=paste("Result:", round((sum(mesum)/(metotal))*100,2), "% methylated")) }
	
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
## see biostring manual for complete set of functions
#REFERENCE:	GCGACGTTATTACGGGA
# CLONE1:	GTGA GTTATTACGGGA	(missing first CG)
# CLONE2:	GTGATGTTATTATGGGA	(everything is unme)
# CLONE3:	GTGATGTTATTA GGGA	(missing last CG)
# CLONE4:	GTGATGTTATTAGGGGA	(mismatch at last CG)
# CLONE5:	GTGACGTTATTACGGGAC	(extra base at end)
# CLONE6:	GTGACGTTTATTACGGGA	(extra base in middle)
