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
#' @return Returns a list with elements Z, W, and groups for input into MCL function, and a vector with the names of groups with only one element (excluded from W matrix).
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

  ## W: within-group proportions
  pj.vec <- table(groups)
  if (length(unique(pj.vec)) == 1) {
    W.tilde <- t(apply(X, 1, FUN = function(yy) {
      c(t(aggregate(yy, by = list(groups), FUN = function(zz) zz/sum(zz))$x))
    }))
  } else {
    W.tilde <- t(apply(X, 1, FUN = function(yy) {
      unname(unlist(aggregate(yy, by = list(groups), FUN = function(zz) zz/sum(zz))$x))
    }))
  }
  colnames(W.tilde) <- colnames(X)[order(groups)]

  groups0 <- c()
  if (length(unique(pj.vec)) == 1) {
    groups0 <- rep(aggregate(X[1,], by = list(groups), FUN = function(zz) zz/sum(zz))$Group.1, each = unique(pj.vec))
  } else {
    for (i in 1:length(unique(groups))) {
      groups0 <- c(groups0,
                   rep(aggregate(X[1,], by = list(groups), FUN = function(zz) zz/sum(zz))$Group.1[i],
                       length(unlist(aggregate(X[1,], by = list(groups), FUN = function(xx) {
                         xx/sum(xx) } )$x[i]))))
    }
  }
  singletons <- names(pj.vec)[which(pj.vec == 1)]
  if (length(singletons) > 0) {
    singleton.features.idx <- which(groups0 %in% singletons)
    singleton.feature.names <- colnames(W.tilde)[singleton.features.idx]
    W.tilde <- W.tilde[ , -singleton.features.idx]
    groups0 <- groups0[-singleton.features.idx]
  } else {
    singleton.feature.names = NA
  }

  W <- log(W.tilde)

  ## Z: group-wide proportions
  Z.tilde <- t(apply(X, 1, FUN = function(xx) aggregate(xx, by = list(groups), FUN = function(zz) sum(zz))$x))
  colnames(Z.tilde) <- aggregate(X[1,], by = list(groups), FUN = function(zz) sum(zz))$Group.1
  Z <- log(Z.tilde)

  # groups of size > 1 are the first q1, then the q0 singleton groups are after that
  Z <- Z[ , c(which(table(groups) != 1), which(table(groups) == 1))]

  return(list(Z = Z, W = W, groups = groups0, singletons = singleton.feature.names))
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
