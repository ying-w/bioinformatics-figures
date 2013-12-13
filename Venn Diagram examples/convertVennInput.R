convertVennInput = function(venndata, type = "matrix") {
# type can be vector/matrix/list

    type = tolower(type)
    if(type != "vector" & type != "matrix" & type != "list" ) { 
        stop("Invalid type parameter, must be either vector, matrix or list") 
    }

    # convert to matrix-type called venn
    nset = 0
    if(is.list(venndata)) {
        nset = length(venndata)
        tmp = unlist(venndata)
        tmpu = unique(tmp)
        venn = matrix(FALSE, nrow = length(tmpu), ncol = nset)
        colnames(venn) = names(venndata)
        rownames(venn) = sort(tmpu)
        for(i in 1:nset) { venn[,i] = rownames(venn) %in% venndata[[i]] }
    }
    else if(is.vector(venndata)) { #vector
        while(length(venndata) > 2^nset) { nset = nset + 1 } #find number of sets
        if(length(venndata) == 2^nset -1) { 
            venndata = c(0, venndata) #catch case forget 0
        }
        else if(length(venndata) != 2^nset) { 
            stop("Incorrect number of elements in vector") 
        }
        
        venn = matrix(FALSE, nrow = sum(venndata), ncol = nset)
        curRow = venndata[1] + 1
        
        for(i in 2:length(venndata)) {
            venn[curRow:(curRow + venndata[i] -1),c(nset:1)[which(intToBits(i-1) == "01")]] = TRUE
            curRow = curRow + venndata[i]
        }
        colnames(venn) = rev(2^(0:(nset-1)))
    }
    else if(ncol(venndata) > 0) { # matrix (or data frame) might edge cases here
        nset = ncol(venndata)
        if(is.logical(venndata)) { venn = venndata }
        else if(all(venndata == 0 | venndata == 1)) { venn = as.logical(venndata) }
        else { stop("Invalid input data type, must be either vector of numbers, matrix of logicals or list of identifiers") }
    }
    else { stop("Input data type not recognized, must be either vector of numbers, matrix of logicals or list of identifiers") }
    
    #output
    if(type == "matrix") { venn } #internal representation is as matrix
    else if(type == "list") { # convert to list using rownames
        if(is.null(rownames(venn))) { rownames(venn) = 1:nrow(venn) }
        apply(venn, 2, function(i) { as.integer(rownames(venn)[which(i)]) } )
    }
    else if(type == "vector") { # convert to vector summing values
        #left columns are bigger
        tmp = venn
        for(i in ncol(venn):1) {
            tmp[venn[,i],i] = 2^(i-1)
        }
        vt = table(rowSums(tmp))
        vec = rep(0, 2^ncol(venn))
        names(vec) = 0:(2^ncol(venn)-1)
        vec[names(vt)] = vt
        vec
    }
}

