
#' string to character vector
#' 
#' Turns a string into a vector of single characters
#' @param str input string
#' @return a vector of single characters
to.char.array <- function (str) sapply(1:nchar(str), function(i) substr(str,i,i))

#' return a character froma string
#' 
#' return's the i'th character in a string
#' @param str input string
#' @param i the character index
#' @return the character
char.at <- function(str,i) substr(str,i,i)


#' constructor to create a new sequence object
#' 
#' creates new yogiseq object
#' @param sequence the sequence
#' @param qual the quality string
#' @param id an identifier for the sequence
#' @return the sequence object
#' @export
new.sequence <- function(sequence, qual=NULL, id=NULL) {

	.seq <- sequence
	.qual <- qual
	.id <- id

	toString <- function() {
		.seq
	}
	getQuality <- function(is) {
		if (!is.null(.qual)) {
			.qual[is]
		} else {
			# warning("This sequence has no quality track.")
			NULL
		}
	}
	getID <- function() {
		if (!is.null(.id)) {
			.id
		} else {
			# warning("This sequence has no ID field.")
			NULL
		}
	}

	structure(list(
		toString=toString,
		getQuality=getQuality,
		getID=getID
	),class="yogiseq")
}

#' print method for yogiseq
#' 
#' @param s the yogiseq object
#' @export
print.yogiseq <- function(s) print(paste("<YogiSeq:",s$getID(),">"))

#' summary method for yogiseq
#' 
#' @param s the yogiseq object
#' @export
summary.yogiseq <- function(s) c(id=s$getID(),sequence=s$toString(),phred=paste(s$getQuality(),collapse=","))

#' length method for yogiseq
#' 
#' @param s the yogiseq object
#' @export
length.yogiseq <- function(s) nchar(s$toString())

#' reverse complement of sequence
#' 
#' @param seq the yogiseq object
#' @return the reverse complement of the sequence
#' @export
reverseComplement <- function(seq) {		
	trans <- c(A='T',C='G',G='C',T='A',N='N',R='Y',Y='R',S='S',W='W',K='M',M='K')
	if (any(class(seq) == "yogiseq")) {
		revSeq <- paste(rev(sapply(to.char.array(seq$toString()), function(nc) trans[nc])),collapse="")
		revQual <- rev(seq$getQuality())
		new.sequence(revSeq,qual=revQual,id=seq$getID())
	} else {
		paste(rev(sapply(to.char.array(seq), function(nc) trans[nc])),collapse="")
	}
}

#' subsequence
#' 
#' @param seq the yogiseq object
#' @param from from index
#' @param to end index
#' @return given subsequece
#' @export
subseq <- function(s,from,to) {
	if (!any(class(s) == "yogiseq")) stop("First argument must be a YogiSeq object")
	new.sequence(
		substr(s$toString(),from,to), 
		if (!is.null(s$getQuality())) s$getQuality(from:to) else NULL,
		s$getID()
	)
}


#' write sequences to fasta file
#' 
#' @param con output connection
#' @param seqs list of yogiseq objects
#' @return given subsequece
#' @export
writeFASTA <- function(con,seqs) {
	if (inherits(con,"character")) {
		con <- file(con,open="w")
	} else if (!inherits(con,"connection")) {
		stop("con must be filename or writeable connection")
	}
	for (i in 1:length(seqs)) {
		s <- seqs[[i]]
		if (class(s) == "yogiseq") {
			writeLines(c(
				paste(">",s$getID(),sep=""),
				s$toString()
			),con)
		} else if (class(s) == "character") {
			writeLines(c(
				paste(">",names(seqs)[i],sep=""),
				s
			),con)
		} else {
			warning("Skipping unsupported data type",class(s))
		}
	}
}

#' read sequences from fasta file
#' 
#' @param con input file or connection
#' @return the sequences
#' @export
readFASTA <- function(con) {
	if (inherits(con,"character")) {
		con <- file(con,open="r")
	} else if (!inherits(con,"connection")) {
		stop("con must be filename or readable connection")
	}
	out <- list()
	id <- NULL
	seq <- NULL
	i <- 0
	while(length(line <- readLines(con, n=1)) > 0) {
		if (substr(line,1,1)==">") {
			#if old sequence exists, add it to the output
			if (!is.null(id)) {
				out[[length(out)+1]] <- new.sequence(seq,id=id)
				cat(paste("\r Read",i <- i+1,"sequences    "))
			}
			#new sequence
			id <- substr(line,2,nchar(line))
			seq <- ""
		} else {
			seq <- paste(seq,line,sep="")
		}
	}
	#add last sequence to output
	if (!is.null(id)) {
		out[[length(out)+1]] <- new.sequence(seq,id=id)
		cat(paste("\r Read",i <- i+1,"sequences    \n"))
	}
	out
}


#' write sequences to fastq file
#' 
#' @param con file or connection
#' @param seqs list of sequences
#' @export
writeFASTQ <- function(con, seqs) {

	if (inherits(con,"character")) {
		con <- file(con,open="w")
	} else if (!inherits(con,"connection")) {
		stop("con must be filename or writeable connection")
	}

	#function for decoding phred quality scores
	qualScale <- to.char.array("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~")
	qual2string <- function(qual) paste(qualScale[qual-32],collapse="")

	writeLines(unlist(lapply(seqs, function(s) {
		c(
			paste("@",s$getID(),sep=""),
			s$toString(),
			"+",
			qual2string(s$getQuality())
		)
	})),con)

}


#' creates a new fastq parser object
#' 
#' @param file or connection
#' @return new fastq parser with method 
#'     <code>parseNext(n=10,ignore.quality=FALSE)</code>
#' @export
new.fastq.parser <- function(con) {

	if (inherits(con,"character")) {
		con <- file(con,open="r")
	} else if (!inherits(con,"connection")) {
		stop("con must be filename or readable connection")
	}

	.con <- con

	#function for decoding phred quality scores
	qualScale <- to.char.array("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~")
	string2phred <- function(string) {
		out <- sapply(to.char.array(string), function(x) which(qualScale == x))
		names(out) <- NULL
		out+32
	}

	#function for parsing the next n entries from the open fastq file (or less if less than n remain)
	parse.next <- function(n=10,ignore.quality=FALSE) {

		contents <- list()
		i <- 0

		while ((i <- i+1) <= n && length(lines <- readLines(.con, n=4)) > 0) {

			if (length(lines) < 4 || substr(lines[1],1,1) != "@" || substr(lines[3],1,1) != "+") {
				stop("Corrupt read:\n",paste(lines,collapse="\n"))
			}

			id <- strsplit(substr(lines[1],2,nchar(lines[1])), " ", fixed=TRUE)[[1]][1]
			sequence <- lines[2]

			quality <- if (ignore.quality) NULL else string2phred(lines[4])

			contents[[length(contents)+1]] <- new.sequence(sequence,id=id,qual=quality)

		}

		contents
	}

	structure(list(parse.next=parse.next),class="yogi.fastq.parser")
}

##
# Needleman-Wunsch global alignment algorithm
#


#' creates a new Needleman-Wunsch global alignment
#' 
#' @param s1 sequence #1
#' @param s2 sequence #2
#' @return new alignment with methods:
#'     <code>getMatrix()</code>
#'     <code>getDistance()</code>
#'     <code>getMutations()</code>
#'     <code>getMappings()</code>
#'     <code>printAlignment()</code>
#'     <code>getMatrix()</code>
#' @export
new.alignment <- function(s1, s2) {

	#alignment algorithm requires bitwise operations
	library(bitops)

	if (any(class(s1)=="yogiseq")) {
		c1 <- c("$",to.char.array(s1$toString()))
	} else {
		c1 <- c("$",to.char.array(s1))
	}
	if (any(class(s2)=="yogiseq")) {
		c2 <- c("$",to.char.array(s2$toString()))
	} else {
		c2 <- c("$",to.char.array(s2))
	}

	#init score matrix
	mat <- matrix(nrow=length(c1), ncol=length(c2))
	mat[1,] <- 1:length(c2) - (c1[1] == c2[1])
	mat[,1] <- 1:length(c1) - (c1[1] == c2[1])

	#init trace matrix
	trace <- matrix(0, nrow=length(c1), ncol=length(c2))
	trace[1,] <- 4
	trace[,1] <- 2
	trace[1,1] <- 0

	#compute alignment matrix
	for (i in 2:length(c1)) {
		for (j in 2:length(c2)) {
			options <- c(
				rep = mat[i-1,j-1] + (c1[i] != c2[j]),
				del = mat[i-1,j] + 1,
				ins = mat[i,j-1] + 1
			)
			mat[i,j] <- min(options)

			tr.bitmasks <- 2^(which(options == min(options))-1)
			for (mask in tr.bitmasks) {
				trace[i,j] <- bitOr(trace[i,j],mask)
			}
		}
	}

	getMatrix <- function() {
		mat
	}

	getDistance <- function() {
		mat[length(c1),length(c2)]
	}

	.mutations <- NULL
	.mapping <- NULL
	run.trace <- function() {

		rep <- 1
		del <- 2
		ins <- 4

		muts <- list()
		map <- list()

		i <- length(c1)
		j <- length(c2)

		while (i > 1 || j > 1) {
			if (bitAnd(trace[i,j], rep) > 0) {
				if (c1[i] != c2[j]) {
					muts[[length(muts)+1]] <- c(c1[i], i-1, j-1, c2[j])
				}
				map[[length(map)+1]] <- c(i-1, j-1)
				i <- i-1
				j <- j-1
			} else if (bitAnd(trace[i,j], del)) {
				muts[[length(muts)+1]] <- c(c1[i], i-1, j-1, "-")
				map[[length(map)+1]] <- c(i-1, NA)
				i <- i-1
			} else if (bitAnd(trace[i,j], ins)) {
				muts[[length(muts)+1]] <- c("-", i-1, j-1, c2[j])
				map[[length(map)+1]] <- c(NA, j-1)
				j <- j-1
			} else {
				stop("uninitialized trace at ",i,j)
			}
		}
		# if (c1[1] != c2[1]) {
		# 	muts[[length(muts)+1]] <- c(c1[1],i,c2[1])
		# }
		.mapping <<- do.call(rbind,rev(map))
		.mutations <<- do.call(rbind,rev(muts))
	}

	getMutations <- function() {
		if (is.null(.mutations)) run.trace()
		.mutations
	}

	getMappings <- function() {
		if (is.null(.mapping)) run.trace()
		.mapping
	}

	printAlignment <- function() {

		if (is.null(.mapping)) run.trace()

		chars <- do.call(cbind,lapply(1:nrow(.mapping), function(k) {
			i <- .mapping[k,1]
			j <- .mapping[k,2]
			char1 <- if (is.na(c1[i+1])) '-' else c1[i+1]
			char2 <- if (is.na(c2[j+1])) '-' else c2[j+1]
			matchChar <- if (is.na(c1[i+1]) || is.na(c2[j+1])) " " 
				else if (c1[i+1] == c2[j+1]) "|" else "."
			c(char1,matchChar,char2)
		}))

		cat("\nLevenstein distance:",getDistance(),"\n")

		for (wrap in 0:(ncol(chars)/70)) {
			startcol <- wrap*70 + 1
			endcol <- if (startcol+69 > ncol(chars)) ncol(chars) else startcol+69
			cat("\n",paste(apply(chars[,startcol:endcol],1,paste,collapse=""),collapse="\n"),"\n",sep="")

		}

	}

	structure(list(
		getMatrix=getMatrix,
		getDistance=getDistance,
		getMutations=getMutations,
		getMappings=getMappings,
		printAlignment=printAlignment
	),class="yogialign")

}


###
# A class for searching strings against an index of k-mers
#
#' K-mer search
#' 
#' @param k length of k-mers
#' @return the k-mer search object with methods:
#'    <code>build.index(fasta.file)</code>
#'    <code>load.index(index.file)</code>
#'    <code>search(queries,min.hits=3,max.d=Inf)</code>
#' @export
new.kmer.search <- function(k=5) {

	library("hash")

	#extract all k-mers from given sting
	kmers <- function(s) {
		if (class(s)=="yogiseq") {
			sapply(1:(length(s)-k+1),function(i) subseq(s,i,i+k-1)$toString())
		} else {
			sapply(1:(nchar(s)-k+1),function(i) substr(s,i,i+k-1))
		}
	}

	#Fields storing the index (a hash mapping kmers to template indices)
	.kmer.index <- NULL
	#Template names 
	.kmer.index.names <- NULL
	#Template sequences
	.template.seqs <- NULL

	build.index <- function(fasta.file) {
		tryCatch({
			con <- file(fasta.file, open="r")
			seqs <- readFASTA(con)
		},
		error = function(ex) {
			# logger$fatal(ex)
			stop(ex)
		},
		finally = {
			if (exists("con") && isOpen(con)) {
				close(con)
			}
		})

		kmer.index <- hash()
		kmer.index.names <- sapply(seqs,function(s)s$getID())

		for (j in 1:length(seqs)) {
			s <- seqs[[j]]
			kms <- kmers(s)
			for (i in 1:length(kms)) {
				kmer.index[[kms[[i]]]] <- c(kmer.index[[kms[[i]]]],j)
			}
		}

		index.file <- sub(".fa","_index.rdata",fasta.file)
		save(kmer.index,kmer.index.names,seqs,file=index.file)

		.kmer.index <<- kmer.index
		.kmer.index.names <<- kmer.index.names
		.template.seqs <<- seqs
	}

	load.index <- function(index.file) {
		load(index.file)
		.kmer.index <<- kmer.index
		.kmer.index.names <<- kmer.index.names
		.template.seqs <<- seqs
	}

	search <- function(queries,min.hits=3,max.d=Inf) {
		yogitools::as.df(lapply(queries, function(s) {

			if (is.null(s) || is.na(s) || length(s) == 0 || nchar(s) == 0) {
				return(NA)
			}
			# cat(s$getID(),"\n")
			kms <- kmers(s)
			#Filter out kmers that don't occur in library
			kms <- kms[kms %in% keys(.kmer.index)]
			#No result if no kmers occur in library
			if (length(kms)==0) {
				return(NA)
			} 
			#table showing number of hits per template id
			nhits <- table(do.call(c,values(.kmer.index,kms,simplify=FALSE)))
			#filter out best match(es) if it fulfills minimum #hits requirement
			if (is.finite(max.d)) {

				#top 3 matches
				top.nhits <- head(sort(nhits[nhits >= min.hits],decreasing=TRUE),3)
				idxs <- as.integer(names(top.nhits))
				#perform alignments for top 3 matches and report distance
				d <- sapply(idxs,function(idx) {
					tseq <- .template.seqs[[idx]]
					new.alignment(tseq,s)$getDistance()
				})
				top.match <- idxs[which(d <= max.d & d==min(d))]
				if (length(top.match) == 1) {
					list(seq=.kmer.index.names[[top.match]],hits=top.nhits[d==min(d)],dist=min(d))
				} else  {
					list(seq=NA,hits=NA,dist=NA)
				}
		
			} else {

				top.nhits <- nhits[nhits >= min.hits & nhits==max(nhits)]
				if (length(top.nhits) == 1) {
					.kmer.index.names[[as.integer(names(top.nhits))]]
					list(seq=.kmer.index.names[[as.integer(names(top.nhits))]],
						hits=top.nhits,
						dist=NA
					)
				} else {
					#in case nothing gets over minimum or there are multiple choices
					list(seq=NA,hits=NA,dist=NA)
				}
			}
		
		}))
	}

	list(
		build.index=build.index,
		load.index=load.index,
		search=search
	)
}

#' Create new barcode matcher
#' 
#' @param lib library of barcode sequences to match against
#' @return a new barcode matcher objects with the method:
#'    findMatches(queries)
#' @export
new.bc.matcher <- function(lib,errCutoff=2,strictMode=FALSE) {
	library(hash)

	#safety checks
	if (!inherits(lib,"character")) {
		stop("`lib' must be of type `character'")
	}
	
	#convert library to integer matrix	
	# charMatrix <- do.call(rbind,lapply(lib,function(l) as.integer(charToRaw(l))))
	#create hash of library
	libHash <- hash(lib,1:length(lib))

	distantMatches <- function(query,errCutoff=2) {
		ncs <- c("A","C","G","T")
		dist1Hits <- do.call(c,lapply(1:nchar(query), function(i) {
			do.call(c,lapply(setdiff(ncs,substr(query,i,i)),function(nc) {
				q1 <- query
				substr(q1,i,i) <- nc
				libHash[[q1]]
			}))
		}))
		if (length(dist1Hits) > 0) {
			return(list(hits=dist1Hits,dist=1))
		} else if (errCutoff < 2) {
			return(list(hits=integer(0),dist=NA))
		}
		dist2Hits <- do.call(c,lapply(2:nchar(query), function(i) {
			do.call(c,lapply(setdiff(ncs,substr(query,i,i)),function(nc) {
				q1 <- query
				substr(q1,i,i) <- nc
				do.call(c,lapply(1:(i-1), function(j) {
					do.call(c,lapply(setdiff(ncs,substr(query,j,j)),function(nc) {
						q2 <- q1
						substr(q2,j,j) <- nc
						libHash[[q2]]
					}))
				}))
			}))
		}))
		if (length(dist2Hits) > 0) {
			return(list(hits=unique(dist2Hits),dist=2))
		} else {
			return(list(hits=unique(dist2Hits),dist=NA))
		}
	}

	findMatches <- function(queries,queries2=NA) {

		if (!inherits(queries,"character")) {
			stop("`queries' must be of type `character'")
		}
		if (!all(is.na(queries2))) {
			if (!inherits(queries2,"character")) {
				stop("`queries2' must be of type `character'")
			}
			if (length(queries) != length(queries2)) {
				stop("queries1 and queries2 must be of equal length!")
			}
			#find disagreeing reads
			disagree <- queries != queries2
			disagree[is.na(disagree)] <- TRUE
			#if there are any cases where only R2 contains a read, flip
			#it into R1 (to use preferentially)
			flips <- which(is.na(queries) & !is.na(queries2))
			if (length(flips) > 0) {
				queries[flips] <- queries2[flips]
				is.na(queries2[flips]) <- TRUE
			}
			#mark those reads with no backup
			nosecond <- is.na(queries2)
		} else {
			disagree <- rep(FALSE,length(queries))
			nosecond <- rep(TRUE,length(queries))
		}

		#iterate over queries
		# system.time({
		do.call(rbind,lapply(1:length(queries), function(i) {
			query <- queries[[i]]
			if ((is.na(query)&&nosecond[[i]]) || (strictMode && disagree[[i]]) ) {
				return(list(hits=NA,diffs=NA,nhits=0))
			}
			#look for perfect matches in hash
			if (!is.na(query)) {
				hashHit <- libHash[[query]]
			} else {
				hashHit <- NULL
			}
			#if available, try read2 as well
			if (disagree[[i]] && !nosecond[[i]]) {
				hashHit2 <- libHash[[queries2[[i]]]]
			} else {
				hashHit2 <- NULL
			}
			#return hits, if there are any
			hits <- union(hashHit,hashHit2)
			if (length(hits) > 0) {
				return(list(hits=hits[[1]],diffs=0,nhits=length(hits)))
			}
			#if no errors are tolerated, we're still done here
			if (errCutoff < 1) {
				return(list(hits=NA,diffs=NA,nhits=0))
			}
			#otherwise, run (slow) inexact search
			#convert query to integer vector
			# qChars <- as.integer(charToRaw(query))
			# #count differences in every row
			# misMatches <- apply(charMatrix,1,function(row) sum(row != qChars))
			# #find row with smallest number of differences
			# minErr <- min(misMatches)
			# hits <- which(misMatches==minErr)
			dmatch <- distantMatches(query,errCutoff)
			hits <- dmatch$hits
			minErr <- dmatch$dist

			#do same for mismatching R2
			if (disagree[[i]] && !nosecond[[i]]) {

				dmatch <- distantMatches(queries2[[i]],errCutoff)
				minErr <- min(c(minErr,dmatch$dist))
				hits <- if (minErr < dmatch$dist) hits else if (minErr > dmatch$dist) dmatch$hits else union(hits,dmatch$hits)
				# qChars2 <- as.integer(charToRaw(queries2[[i]]))
				# #count differences in every row
				# misMatches2 <- apply(charMatrix,1,function(row) sum(row != qChars2))
				#find row with smallest number of differences
				# minErr <- min(c(minErr,misMatches2))
				# hits <- union(which(misMatches==minErr),which(misMatches2==minErr))
			}

			if (is.na(minErr) || minErr > errCutoff) {
				list(hits=NA,diffs=NA,nhits=0)
			} else {
				list(hits=hits[[1]],diffs=minErr,nhits=length(hits))
			}
		}))
		# })
	}

	list(findMatches=findMatches)
}


# read.sam <- function(sam.file) {
# 	tryCatch({
# 		sam.con <- file(sam.file,open="r")
# 		lines <- readLines(sam.con)
# 		lines <- lines[substr(lines,1,1)!="@"]
# 		split <- strsplit(lines,"\t")
# 		ncol <- max(sapply(split,length))
# 		sam <- do.call(rbind,lapply(split,function(row) c(row,rep(NA,ncol-length(row)))))
# 		colnames(sam) <- c(
# 			"cname","flag","rname","pos","mapq","cigar","mrnm","mpos",
# 			"isize","seq","qual","tags",13:ncol
# 		)
# 		sam <- to.df(sam)
# 		sam$flag <- as.integer(sam$flag)
# 		sam$pos <- as.integer(sam$pos)
# 		sam$mapq <- as.integer(sam$mapq)
# 		sam$mpos <- as.integer(sam$mpos)
# 		sam$isize <- as.integer(sam$isize)
# 		sam$mapq <- as.integer(sam$mapq)
# 		sam
# 	},
# 	error=function(e) {
# 		logger$fatal(e)
# 		stop(e)
# 	},
# 	finally={
# 		if (exists("sam.con") && isOpen(sam.con)) {
# 			close(sam.con)
# 		}
# 	})
	
# }

# sam2pileup <- function(sam.file,ref.file) {
# 	tryCatch({
# 		ref.con <- file(ref.file,open="r")
# 		ref.seq <- readFASTA(ref.con)[[1]]
# 	},
# 	error=function(e) {
# 		logger$fatal(e)
# 		stop(e)
# 	},
# 	finally={
# 		if (exists("ref.con") && isOpen(ref.con)) {
# 			close(ref.con)
# 		}
# 	})

# 	# sam <- read.delim(sam.file,header=FALSE,stringsAsFactors=FALSE,skip=3)
# 	sam <- read.sam(sam.file)

# 	flagMasks <- c(
# 		multiSegment=0x1, allSegmentsOK=0x2, segmentUnmapped=0x4,
# 		nextSegmentUnmapped=0x8, revComp=0x10, nextRevComp=0x20, 
# 		firstSegment=0x40, lastSegment=0x80, secondary=0x100, 
# 		failQC=0x200, duplicate=0x400, supplementary=0x800
# 	)
# 	flags <- do.call(rbind,lapply(sam$flag,function(x)bitAnd(x,flagMasks)>0))
# 	colnames(flags) <- names(flagMasks)
# 	flags <- to.df(flags)

# 	#CIGAR: S=Soft clip, H=Hard clip, N=Intron skip, M=Match, D=Deletion, I=Insertion, P=Padded
# 	cigar <- global.extract.groups(sam$cigar,"(\\d+)([SHNMDIP]{1})")

# 	start.stop <- do.call(rbind,lapply(1:nrow(sam),function(i) {
# 		if (flags$segmentUnmapped[[i]]) {
# 			return(c(start=NA,end=NA))
# 		}
# 		l <- sum(as.integer(cigar[[i]][cigar[[i]][,2] %in% c("M","D"),1]))
# 		c(start=sam$pos[[i]],sam$pos[[i]]+l)
# 	}))
# 	pcr.dup <- apply(is.na(start.stop),1,any) | duplicated(start.stop)

# 	out.sam <- sub("/[^/]+\\.sam$","/nodup.sam",sam.file)
# 	write.table(sam[!pcr.dup,],out.sam,sep="\t",quote=FALSE,row.names=FALSE)

# 	pileup <- list(
# 		bases=replicate(length(ref.seq),character()),
# 		qual=replicate(length(ref.seq),numeric()),
# 		ins=replicate(length(ref.seq),character())
# 	)
# 	for (i in 1:nrow(sam)) {

# 		if (flags$segmentUnmapped[[i]] || pcr.dup[[i]]) {
# 			next
# 		}

# 		qtrack <- as.integer(charToRaw(sam$qual[[i]]))-33
# 		read <- to.char.array(sam$seq[[i]])
# 		tp <- sam$pos[[i]] #template position
# 		rp <- 1 #read position
# 		for (cigrow in 1:nrow(cigar[[i]])) {
# 			k <- as.integer(cigar[[i]][cigrow,1])
# 			op <- cigar[[i]][cigrow,2]
# 			if (op=="M") {
# 				mstart <- rp
# 				while (rp < mstart+k) {
# 					pileup$bases[[tp]][[length(pileup$bases[[tp]])+1]] <- read[[rp]]
# 					pileup$qual[[tp]][[length(pileup$qual[[tp]])+1]] <- qtrack[[rp]]
# 					rp <- rp+1
# 					tp <- tp+1
# 				}
# 			} else if (op=="D") {
# 				mstart <- rp
# 				for (.dummy in 1:k) {
# 					pileup$bases[[tp]][[length(pileup$bases[[tp]])+1]] <- "*"
# 					pileup$qual[[tp]][[length(pileup$qual[[tp]])+1]] <- sam$mapq[[i]]
# 					tp <- tp+1
# 				}
# 			} else if (op=="I") {
# 				ins.bases <- paste(read[rp:(rp+k-1)],collapse="")
# 				pileup$ins[[tp]][[length(pileup$ins[[tp]])+1]] <- ins.bases
# 				rp <- rp + k
# 			} else if (op %in% c("S","H")) {
# 				# tp <- tp + k
# 				rp <- rp + k
# 			} else {
# 				warning("Unsupported cigar character: ",op, sam$cigar[[i]])
# 				tp <- tp + k
# 			}
# 		}
# 	}
		

# 	pu <- mapply(function(bases, qual){
# 		data.frame(base=bases,p=10^(-qual/10))
# 	},bases=pileup$bases,qual=pileup$qual,SIMPLIFY=FALSE)
# 	names(pu) <- 1:length(ref.seq)
	
# 	list(pileup=pu,indel.track=pileup$ins)
# }

# var.call <- function(piles, ref, indel.track, ref.length, threshold=.05) {
#     bases <- c("A","C","G","T","*")
#     freqs <- do.call(rbind,lapply(piles, function(pile.i) {
#         fpile <- pile.i[pile.i$p < threshold,]
#         table(factor(fpile$base,levels=bases))
#     }))
#     d <- apply(freqs,1,sum) + sapply(indel.track,length)
#     names(d) <- names(piles)

#     #check indels
#     indel.track <- lapply(indel.track, toupper)
#     indel.idxs <- which(sapply(indel.track,length) > 0)
#     called.indels <- to.df(do.call(rbind,lapply(indel.idxs, function(i) {
#     	indel.freqs <- table(toupper(indel.track[[i]]))
#     	do.call(rbind,lapply(names(indel.freqs), function(indel) {
#     		f <- indel.freqs[[indel]]
#     		if (f > 1 && f/d[[i]] > threshold) {
#     			list(ref=ref[[i]],pos=names(piles)[[i]],alt=indel,freq=f/d[[i]])
#     		} else {
#     			NULL
#     		}
#     	}))
#     })))

#     #check SNVs
#     skimmed.freqs <- apply(freqs,c(1,2), function(x)if(x < 2) 0 else x)
#     calls <- lapply(1:nrow(freqs), function(i) {
#         f <- skimmed.freqs[i,]/d[[i]]
#         nonref <- f[setdiff(bases,ref[[i]])]
#         nonref[!is.na(nonref) & nonref > threshold]
#     })

#     idxs <- which(sapply(calls,length) > 0)
#     called.snvs <- to.df(do.call(rbind,lapply(idxs, function(i) {
#         pos <- as.numeric(names(piles)[[i]])
#         vars <- calls[[i]]
#         ref <- ref[[i]]
#         do.call(rbind,lapply(names(vars), function(base) 
#             list(ref=ref,pos=pos,alt=base,freq=vars[[base]])
#         ))
#     })))

#     #create depth vector for all positions (including those not in alignment)
#     d.all <- sapply(as.character(1:ref.length), function(pos) if (pos %in% names(d)) d[[pos]] else 0)

#     list(calls=rbind(called.snvs,called.indels),depth=d.all)
# }

# base.posteriors <- function(piles) {
#     do.call(rbind,lapply(piles, function(pile.i) {
#         #possible bases
#         qis <- c("A","C","G","T","*")
#         posteriors <- sapply(qis, function(qi) {
#             #skip impossible bases
#             if (!(qi %in% pile.i$base)) {
#                 return (0)
#             }
#             # compute the log-odds by iterating over all symbols at the pileup position
#             lo.i <- sum(sapply(1:nrow(pile.i), function(j){
#                 #the base symbol
#                 bij <- pile.i$base[[j]]
#                 #the error probability
#                 pij <- pile.i$p[[j]]
#                 #calculate the Bayes factor
#                 if (qi==bij) {
#                     log(1-pij) - log(pij/3)
#                 } else {
#                     log(pij/3) - log(1/3)
#                 }
#             })) + log(1/length(qis))
#             # then transform the log-odds to the probability
#             # being careful to avoid NaNs
#             if (lo.i > 38) 1 else exp(lo.i)/(1+exp(lo.i))
#         })
#     }))
# }




# default.error <- function(ex) {
# 	print(ex)
# 	traceback(ex)
# 	stop()
# }
# protect <- function(filename, fun, mode="r",error=default.error) {
# 	tryCatch({
# 		con <- file(filename, open=mode)
# 		fun(con)
# 	},
# 	error = error,
# 	finally = {
# 		if (exists("con") && isOpen(con)) {
# 			close(con)
# 		}
# 	})
# }







