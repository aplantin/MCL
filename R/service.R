#' dataprep.mcl
#'
#' Prepares input matrices for MCL function based on otu count matrix.
#' Rows of the otu table should be subjects.
#'
#' @param otu.table Matrix of otu counts; rows are subjects, columns are OTUs.
#' @param groups Vector indicating group membership
#' @return Returns a list with elements Z and W for input into MCL function 
#' @export
#'
dataprep.mcl <- function(otu.table, groups) {
    ## X: overall rel. abund.
    X <- otu.table
    for (i in 1:nrow(X)) {
        X[i, ] <- X[i, ] / sum(X[i, ])
    }
    
    ## Z: group-wide proportions
    Z.tilde <- t(apply(X, 1, FUN = function(xx) aggregate(xx, by = list(groups), FUN = function(zz) sum(zz))$x))
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
