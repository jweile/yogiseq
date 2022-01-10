
#' Determine the length of the a common prefix between
#' two strings. For example: Banana and Band have a common
#' prefix of length 3.
#' 
#' @param x first string
#' @param y second string
#' @param wildcard a wildcard character that always counts as a match
#'        or NULL to disable. Default: NULL
#' @return the length of the common prefix
#' @export
commonPrefix <- function(x,y,wildcard=NULL) {

  x <- as.integer(charToRaw(x))
  y <- as.integer(charToRaw(y))
  n <- min(length(x),length(y))
  if (n==0) {
    return(0)
  }

  if (!is.null(wildcard)) {
    wc <- as.integer(charToRaw(wildcard))
  } else {
    wc <- -1
  }

  for (i in 1:n) {
    if (x[[i]] != y[[i]]) {
      if (wc == -1 || x[[i]] != wc && y[[i]] != wc) {
        return(i-1)
      }
    }
  }
  return(n)
}

#' Construct a suffix tree for a given string.
#'
#' @param str the input string
#' @return the suffix tree (implemented as list of lists)
#' @export
suffixTree <- function(str) {
    
  stInsert <- function(suff,tree,i) {
    found <- FALSE
    if (length(tree) > 0) {
      for (j in 1:length(tree)) {
        brlabel <- names(tree)[[j]]
        p <- commonPrefix(suff,brlabel)
        if (p > 0) {
          found <- TRUE
          #bump existing content to subtree
          if (is.numeric(tree[[j]]) || p < nchar(brlabel)) {
            sublabel <- substr(brlabel,p+1,nchar(brlabel))
            tree[[j]] <- setNames(list(tree[[j]]),sublabel)
          }
          #rename node
          label <- substr(suff,1,p)
          names(tree)[[j]] <- label
          tree[[j]] <- stInsert(substr(suff,p+1,nchar(suff)),tree[[j]],i)
          #no need to search the rest of the branches
          break
        }
      }
    }
    if (!found) {
      tree[[paste0(suff,"$")]] <- i
    }
    return(tree)
  }

  tree <- list()
  n <- nchar(str)
  for (i in n:1) {
    suff <- substr(str,i,n)
    tree <- stInsert(suff,tree,i)
  }

  return(tree)

}

#' Search for a query string within a suffix tree
#' with or without wildcard characters.
#' 
#' @param tree the suffix tree
#' @param query the query string
#' @param wildcard A wildcard character that will match anything
#'        or NULL to disable. Default: NULL
#' @return An integer vector indicating the start positions of
#'         any matches
#' @export
searchSuffixTree <- function(tree,query,wildcard=NULL) {

  out <- integer()
  for (i in 1:length(tree)) {

    label <- names(tree)[[i]]
    #if this is a leaf, we can't match anything
    if (is.null(label)) {
      next
    }

    llen <- nchar(label)
    if (substr(label,llen,llen)=="$") {
      llen <- llen-1
    }

    if (llen > 0) {
      p <- commonPrefix(label,query,wildcard)
      if (p == 0 || (p < llen && nchar(query) > p)) {
          #then it's not a match
      } else {
        if (p==llen && nchar(query) > p) {
          #then we have to confirm the match further down
          out <- c(out, searchSuffixTree(
            tree[[i]], 
            substr(query,llen+1,nchar(query)),
            wildcard=wildcard
          ))
          #if no wildcard allowed, only one path can match
          if (is.null(wildcard)) {
            break
          }
        } else if (p <= llen && p == nchar(query)) {
          #full match
          out <- c(out, as.vector(unlist(tree[[i]])) )
          #if no wildcard allowed, only one path can match
          if (is.null(wildcard)) {
            break
          }
        }
      }
    }
  }
  return(out)
}

# matchSeqPos <- function(query,template,wildcard="N") {
#   searchSuffixTree(suffixTree(template),query,wildcard)
# }
