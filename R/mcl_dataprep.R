#' dataprep.mcl
#'
#' Prepares input matrices for MCL function based on otu count matrix.
#' Rows of the otu table should be subjects.
#'
#' @param otu.table Matrix of otu counts; rows are subjects, columns are OTUs.
#' @param groups Vector indicating group membership
#' @param rare.count OTU count considered "rare"; OTUs with total read count <= rare.count will be excluded. Default is 1 (exclude singletons). If <= 0, no filter is applied for rare OTUs based on counts.
#' @param rare.prop Alternative to rare.count; minimum proportion of samples in which the OTU should appear. Default is 0. If <= 0, no filter is applied for rare OTUs based on proportions.
#' @param pseudocount Count to replace zeros; default is 0.5.
#' @return Returns a list with elements Z and W for input into MCL function.
#' @export
#'
dataprep.mcl <- function(otu.table, groups, rare.count = 1, rare.prop = 0, pseudocount = 0.5) {
  if (rare.count > 0) {
    rare1 <- which(apply(otu.table, 2, sum) <= rare.count)
  } else {
    rare1 <- c()
  }
  if (rare.prop > 0) {
    rare2 <- which(apply(otu.table, 2, FUN = function(x) sum(x != 0)) <= rare.prop)
  } else {
    rare2 <- c()
  }
  if (length(rare1) > 0 | length(rare2) > 0) {
    rare <- unique(c(rare1, rare2))
    otus.common <- otu.table[, -rare]
    groups <- groups[-rare]
  } else {
    otus.common <- otu.table
  }

  if (pseudocount <= 0) {
    stop("Pseudocount must be a positive number.")
  }
  otus.common[which(otus.common == 0)] <- pseudocount
  X <- otus.common

  ## X: overall rel. abund.
  for (i in 1:nrow(X)) {
    X[i, ] <- X[i, ] / sum(X[i, ])
  }

  ## Z: group-wide proportions
  Z.tilde <- t(apply(X, 1, FUN = function(xx) aggregate(xx, by = list(groups), FUN = function(zz) sum(zz))$x))
  colnames(Z.tilde) <- aggregate(X[1,], by = list(groups), FUN = function(zz) sum(zz))$Group.1
  Z <- log(Z.tilde)

  ## W: within-group proportions
  pj.vec <- unname(table(groups))
  if (length(unique(pj.vec)) == 1) {
    W.tilde <- t(apply(X, 1, FUN = function(yy) {
      c(t(aggregate(yy, by = list(groups), FUN = function(zz) zz/sum(zz))$x))
    }))
  } else {
    W.tilde <- t(apply(X, 1, FUN = function(yy) {
      unname(unlist(aggregate(yy, by = list(groups), FUN = function(zz) zz/sum(zz))$x))
    }))
  }
  W <- log(W.tilde)

  colnames(W) <- colnames(X)
  return(list(Z = Z, W = W))
}





#' group.sums
#'
#' Adds values within groups
#'
#' @param xx Vector of values to sum within groups
#' @param groups Vector indicating group membership
#' @export
#'
group.sums <- function(xx, groups) {
    return( aggregate(xx, by = list(groups), sum)$x )
}
