
#' Tokenize MDZ string, breaking it down into a list of number and string elemetns
#' E.g. "MD:Z:3A5^AC8" would be broken down into list(3,"A",5,"^AC",8)
parseMDZ <- function(mdz) {
  #check for problems
  if (is.na(mdz)) return(NA) 
  if (mdz=="") return(list())
  #trim off "MD:Z:" header
  if (substr(mdz,1,5)=="MD:Z:") {
    mdz <- substr(mdz,6,nchar(mdz))
  }
  #boolean vector whether each character is a digit or not
  isDigit <- as.integer(charToRaw(mdz)) < 58
  # prepare output
  out <- list()
  #iterate over charcters and split off strings when digit status changes
  curr <- substr(mdz,1,1)
  for (i in 2:nchar(mdz)) {
    if (xor(isDigit[[i-1]],isDigit[[i]])) {
      out[[length(out)+1]] <- if (isDigit[[i-1]]) {
        as.integer(paste(curr,collapse="")) 
      } else {
        paste(curr,collapse="")
      }
      curr <- substr(mdz,i,i)
    } else {
      # curr <- c(curr,substr(mdz,i,i))
      curr[[length(curr)+1]] <- substr(mdz,i,i)
    }
  }
  out[[length(out)+1]] <- if (isDigit[[i]]) {
    as.integer(paste(curr,collapse="")) 
  } else {
    paste(curr,collapse="")
  } 
  return(out)
}

#' construct new SAM streamer object
#' 
#' can be used to stream
#' chunks of SAM files to be processed iteratively. For each chunk,
#' the individual columns as well as parsed flags, MDZ and CIGAR strings
#' are made available
#' @param sam.file the input sam file or connection object
#' @param chunkSize how many lines to process at once
#' @return a sam.streamer object with methods to access the content
#' @export
new.sam.streamer <- function(sam.file, chunkSize=100) {

  library(bitops)
  library(yogitools)
  options(stringsAsFactors=FALSE)

  flagMasks <- c(
    multiSegment=0x1, allSegmentsOK=0x2, segmentUnmapped=0x4,
    nextSegmentUnmapped=0x8, revComp=0x10, nextRevComp=0x20, 
    firstSegment=0x40, lastSegment=0x80, secondary=0x100, 
    failQC=0x200, duplicate=0x400, supplementary=0x800
  )
  fieldnames <- c(
    "cname","flag","rname","pos","mapq","cigar","mrnm","mpos",
    "isize","seq","qual"
  )
  
  #custom table parser that deals with unpredictable column numbers in SAM
  lines2df <- function(lines) {
    #list of lists
    lol <- strsplit(lines,"\t")
    ncols <- max(sapply(lol,length))
    #this is inconveniently slow, but there's not really a way around it.
    tagnames <- Reduce(union,lapply(lol, function(elems) substr(elems[12:length(elems)],1,2)))
    #supplement missing columns
    df <- as.data.frame(do.call(rbind,lapply(lol, function(elems) {
        foundTags <- substr(elems[12:length(elems)],1,2)
        ord <- sapply(foundTags, function(tag) which(tagnames==tag))
        tagarray <- rep(NA_character_,length(tagnames))
        tagarray[ord] <- sapply(elems[12:length(elems)], function(el) substr(el,6,nchar(el)))
        c(elems[1:11],tagarray)
    })))
    #convert strings to numbers as appropriate
    for (i in 1:ncols) {
      if (!any(is.na(suppressWarnings(as.numeric(na.omit(df[,i])))))) {
        df[,i] <- as.numeric(df[,i])
      }
    }
    #set column names
    colnames(df) <- c(fieldnames,tagnames)
    #output
    return(df)
  }

  if (inherits(sam.file,"connection")) {
    samcon <- sam.file
  } else if (inherits(sam.file,"character")) {
    samcon <- file(sam.file,open="r")
  } else {
    stop("sam.file must be a file or connection!")
  }

  samChunk <- flagChunk <- cigarChunk <- mdzChunk <- NULL


  nextChunk <- function() {

    if (!isOpen(samcon)) {
      samChunk <<- flagChunk <<- cigarChunk <<- mdzChunk <<- NULL
      return(0)
    }

    lines <- readLines(samcon,chunkSize)

    if (length(lines)==0) {
      samChunk <<- flagChunk <<- cigarChunk <<- mdzChunk <<- NULL
      close(samcon)
      return(0)
    }

    sam <- lines2df(lines)
    if (nrow(sam) < chunkSize) {
      close(samcon)
    }

    samChunk <<- sam

    flags <- do.call(rbind,lapply(sam$flag,function(x)bitAnd(x,flagMasks)>0))
    colnames(flags) <- names(flagMasks)
    flagChunk <<- to.df(flags)

    #parse CIGAR strings
    cigarChunk <<- global.extract.groups(sam$cigar,"(\\d+)([SHNMDIP]{1})")

    #parse MD string
    mdzChunk <<- lapply(sam$MD, parseMDZ)

    return(nrow(sam))
  }

  return(list(
    nextChunk=nextChunk,
    getSamRows=function()samChunk,
    getFlags=function()flagChunk,
    getCigars=function()cigarChunk,
    getMDZs=function()mdzChunk,
    getSamElement=function(i,j)samChunk[i,j],
    getFlagRow=function(i)flagChunk[i,],
    getCigar=function(i)cigarChunk[[i]],
    getMDZ=function(i)mdzChunk[[i]]
  ))
  
}


#' Find variant call candidates from a sam file line
#' 
#' Interpret the CIGAR and MDZ strings to construct a 
#' table of variant base calls
#' 
#' @param flagV vector of flags (from samStreamer)
#' @param mdzStr MD:Z string from sam
#' @param cigarStr CIGAR string from sam
#' @param startPos alignment start position from sam
#' @param readSeq read sequence from sam
#' @param readQual quality string from sam
#' @return a table of all variant base calls
#' @export
varcallCandidates <- function(flagV, mdzStr, cigarStr, startPos, readSeq, readQual) {

  if (is.na(flagV$segmentUnmapped) || flagV$segmentUnmapped) {
    return(data.frame())
  }
  # if (grepl("^\\d+$",mdzStr) && grepl("^\\d+M$",cigarStr)) {
  if (all(sapply(mdzStr,is.numeric)) && all(cigarStr[,2]=="M")) {
    return(data.frame())
  }

  #here we derive the position indices for each MDZ element
  #first: how much does each MDZ element add to the reference cursor?
  refAdds <- sapply(mdzStr,function(x) {
    if (is.numeric(x)) x else
    if (substr(x,1,1)=="^") nchar(x)-1 else nchar(x)
  })
  readAdds <- sapply(mdzStr,function(x) {
    if (is.numeric(x)) x else
    if (substr(x,1,1)=="^") 0 else nchar(x)
  })
  #now we can use those to get the actual positions
  #note that these are not necessarily correct yet, as they 
  #don't take into account any clipping or any insertions.
  #we will get those from the CIGAR string below
  # refPos <- startPos-1+sapply(1:length(refAdds),function(j)sum(refAdds[1:j]))
  # readPos <- sapply(1:length(readAdds),function(j)sum(readAdds[1:j]))
  refPos <- startPos-1+cumsum(refAdds)
  readPos <- cumsum(readAdds)

  #Calculate the reference and read positions for CIGAR elements
  #in the same way
  cRefAdds <- as.integer(sapply(1:nrow(cigarStr), function(j) {
    if (cigarStr[j,2] %in% c("M","D")) cigarStr[j,1] else 0
  }))
  cReadAdds <- as.integer(sapply(1:nrow(cigarStr), function(j) {
    if (cigarStr[j,2] %in% c("M","I","S","H")) cigarStr[j,1] else 0
  }))
  # cRefPos <- startPos+sapply(1:length(cRefAdds),function(j)sum(cRefAdds[1:j]))
  # cReadPos <- 1+sapply(1:length(cReadAdds),function(j)sum(cReadAdds[1:j]))
  cRefPos <- startPos+cumsum(cRefAdds)
  cReadPos <- 1+cumsum(cReadAdds)
  cigarPos <- data.frame(
    refStart=c(startPos,cRefPos[1:length(cRefPos)-1]),
    refEnd=cRefPos-1,
    op=cigarStr[,2],
    readStart=c(1,cReadPos[1:length(cReadPos)-1]),
    readEnd=cReadPos-1
  )
  
  #construct a summary table from the MDZ information
  mmidx <- which(sapply(mdzStr,is.character))
  vars <- data.frame(
    refpos=refPos[mmidx],
    op=do.call(c,mdzStr[mmidx]),
    readpos=readPos[mmidx]
  )
  vars$refbase <- sub("\\^","",vars$op)
  vars$op <- sapply(vars$op,function(o) if (substr(o,1,1)=="^") "del" else "sub")
  vars$indelen <- nchar(vars$refbase)
  vars$indelen[vars$op=="sub"] <- NA
  #correct for deletion position offset
  vars[vars$op=="del","readpos"] <- vars[vars$op=="del","readpos"]+1

  #now for the brutal task of adjusting the table based on the CIGAR string.
  #we iterate over cigarPos table and account for insertions and clips
  for (j in 1:nrow(cigarPos)) {
    if (cigarPos[j,"op"] %in% c("I","S","H")) {
      #length of the insertion or clip
      ilen <- as.integer(cigarStr[j,1])
      if (cigarPos[j,"op"] == "I") {#insertion
        tableRow <- data.frame(
          refpos=cigarPos[j,"refStart"],
          op="ins",
          readpos=cigarPos[j,"readStart"],
          refbase="-",
          indelen=ilen
        )
        #find the appropriate spot in the vars table to add our insertion
        #by checking which positions are below
        below <- vars$refpos < cigarPos[j,"refStart"]
        if (all(below)) {
          #add at the end
          vars <- rbind(vars, tableRow)
        } else if (!any(below)) {
          #shift all read positions accordingly
          vars$readpos <- vars$readpos+ilen
          #then add at the start
          vars <- rbind(tableRow, vars)
        } else {
          #find correct row to split the table at
          cutPoint <- max(which(below))
          #shift positions after that cut point
          # vars$readpos[1:cutPoint] <- vars$readpos[1:cutPoint]+ilen
          vars$readpos[(cutPoint+1):nrow(vars)] <- vars$readpos[(cutPoint+1):nrow(vars)]+ilen
          #and insert row in the cut point
          vars <- rbind(vars[1:cutPoint,], tableRow, vars[(cutPoint+1):nrow(vars),])
        }
      } else {#clipping
        if (!any(vars$refpos < cigarPos[j,"refStart"])) {
          #if the clip is at the start, then we have to shift our
          #read positions
          vars$readpos <- vars$readpos+ilen
        } else {
          #otherwise the clip is at the end and it doesn't affect anything
          #so we do nothing here
        }
      }
    }
  }

  #if at this point we still have no variants, then we can quit
  if (nrow(vars) == 0) {
    return(data.frame())
  }

  #now we can add readbase and readqual information
  vars$readbase <- sapply(1:nrow(vars), function(j) {
    switch(vars[j,"op"],
      sub=substr(readSeq,vars[j,"readpos"],vars[j,"readpos"]),
      ins=substr(readSeq,vars[j,"readpos"],vars[j,"readpos"]+vars[j,"indelen"]-1),
      del="-",
    )
  })
  vars$qual <- sapply(1:nrow(vars), function(j) {
    switch(vars[j,"op"],
      sub=substr(readQual,vars[j,"readpos"],vars[j,"readpos"]),
      ins=substr(readQual,vars[j,"readpos"],vars[j,"readpos"]+vars[j,"indelen"]-1),
      del=substr(readQual,vars[j,"readpos"]-1,vars[j,"readpos"]),
    )
  })

  vars[,c("op","refbase","refpos","readbase","qual")]

}
