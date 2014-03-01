# Output IGV snapshots/screenshots from R
# Author: Ying Wu
# Current version: July 2013
#
# Take genomic coordinates (from BED file) or GRanges objects and pass them to
# IGV using default port that IGV is listening on and export pic of region
#
# Possible to run this remotely by using
# ssh -Y and `cat` the shortcut, should look something like this:
# javaws "/home/<username>/.icedtea/cache/0/http/www.broadinstitute.org/igv/projects/current/igv_lm.jnlp"
#
# For a list of commands see: http://www.broadinstitute.org/igv/PortCommands
# Example: http://plindenbaum.blogspot.com/2011/07/controlling-igv-through-port-my.html
#
# Sometimes IGV will die after being open for too long for no reason, when that happens, this script will hang
#
# Note:
# New version of bedtools has igv command that will generate a script that will dump IGV figures
#
# TODO:
# need to be more robust to different delimiters
# figure out a way to first detect if IGV is running

r2igv = function(..., dist=200, igv=FALSE, prefix="tmp") {
  # Assume IGV is already open and proper session is loaded with snapshot directory set
  # This script will take bed files or vectors of whitespace seperated choordinates and 
  # make screenshots of IGV at these locations with default +/-200 window around coordinates
  # This script only works for linux probably
  # To load session, use: echo "load mysavedsession.xml" | nc 127.0.0.1 60151
  # To set ouput dir: echo "snapshotDirectory /home/usr/igvout" | nc 127.0.0.1 60151  
  
  #ret = paste(..., collapse=" ")
  #ret = strsplit(ret, split=" ")[[1]]	
  
  if(!is.null(dim(...))) {
      ret = apply(..., 1, paste, collapse=" ")
      #ret = unlist(lapply(strsplit(ret, split=" +"), paste, collapse=" ")) #remove extra spaces
  } else { ret = as.vector(...) }
  ret = t(as.data.frame(strsplit(ret, split=" +")))
  ret[,2] = as.numeric(ret[,2]) - dist
  ret[,3] = as.numeric(ret[,3]) + dist
  region = paste0(ret[,1], ":", ret[,2], "-", ret[,3])

  if(igv) {
    #http://stackoverflow.com/questions/7014081/capture-both-exit-status-and-output-from-a-system-call-in-r
    #cannot use system2 since need piping
    
    if(nrow(ret) > 500) { stop("Over 500 regions specified!!") }
    for(i in 1:nrow(ret)) {
        #possibly put all system stuff into one line seperated by ;
        Sys.sleep(2)
        message(paste("goto", region[i], ": "), appendLF=FALSE)
        system(paste("echo \"goto",  region[i], "\"| nc 127.0.0.1 60151"))
        if(i == 1) 
            Sys.sleep(10) #account for hard drive spin-up
        else
            Sys.sleep(2)
        message(paste0("Saving ", prefix, "_", ret[i,1], "-", floor(mean(as.numeric(c(ret[i,2], ret[i,3])))), ".png : "), appendLF=FALSE)
        system(paste0("echo \"snapshot ", prefix, "_", ret[i,1], "-", floor(mean(as.numeric(c(ret[i,2], ret[i,3])))), ".png\"| nc 127.0.0.1 60151"))
    }
  }
  region
}

##example using makeVenn
#roi = large[extractOverlap("large", "small", res=res, typ=typ)] #regions of interest
#names(roi) = 1:length(roi)
#r2igv(as.data.frame(roi), dist=500, igv=TRUE, prefix="largesmall")
