#' mcl
#'
#' Calculates MCL coefficient estimates over a grid of (lam1, lam2) pairs
#'
#'
#' @useDynLib MCL
#' @importFrom Rcpp sourceCpp
#' @importFrom stats aggregate sd
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @param Z Otu matrix at group level
#' @param W Otu matrix within groups
#' @param y Outcome vector
#' @param mu Lagrange scaling parameter; default 1
#' @param groups Vector of group membership per feature
#' @param lam.seq Sequence of lambda values; default NULL
#' @param lam.max Maximum lambda value; default NULL
#' @param nlam Vector of length 2 indicating number of lambda values to test for lam1 and lam2; default (10, 10)
#' @param min.frac Fraction of max.lam that the minimum lambda value should be; default 0.01
#' @param maxit Maximum number of iterations; default 10000
#' @param thresh Threshold for convergence; default 0.0001
#' @param std Logical indicating whether data should be standardized; default TRUE
#' @param verbose Print lambdas as indicator of progress?; default FALSE
#' @export
#'
mcl <- function(Z, W, y, mu = 1, groups, lam.seq = NULL, lam.max = NULL,
                nlam = c(10,10), min.frac = 0.1, maxit = 1e4, thresh = 1e-4,
                std = TRUE, verbose = FALSE) {
  if (std == TRUE) {
    y <- scale(y, scale = FALSE)
    Z <- scale(Z, scale = apply(Z, 2, sd)*sqrt(nrow(Z)-1))
    facZ <- 1/attr(Z, "scaled:scale")
    W <- scale(W, scale = apply(W, 2, sd)*sqrt(nrow(W)-1))
    facW <- 1/attr(W, "scaled:scale")
  } else {
    facZ <- rep(1, ncol(Z))
    facW <- rep(1, ncol(W))
  }

  groupIdx <- list()
  for (j in 1:length(unique(groups))) {
    groupIdx[[j]] <- which(groups == names(table(groups))[j]) - 1 ## -1 for 0-indexing
  }
  groups = as.numeric(as.factor(groups)) - 1


  if (is.null(lam.seq)) {
    if (is.null(lam.max)) {
      lam.max <- findMaxLams(Z, W, y, facZ, facW, groups, groupIdx, mu, maxit, thresh)
    }
    lam.seq <- list()
    lam.seq[[1]] <- exp(seq(from = log(lam.max[1]),
                            to = log(min.frac*lam.max[1]),
                            length.out = nlam[1]))
    lam.seq[[2]] <- exp(seq(from = log(lam.max[2]),
                            to = log(min.frac*lam.max[2]),
                            length.out = nlam[2]))
  }

  betas <- matrix(nrow = ncol(Z), ncol = (nlam[1]*nlam[2]))
  gammas <- matrix(nrow = ncol(W), ncol = (nlam[1]*nlam[2]))
  if (verbose) {
    pbcount <- 0
    pb <- txtProgressBar(min = 0, max = (nlam[1]*nlam[2]), style = 3)
  }
  for (i in 1:nlam[1]) {
    for (j in 1:nlam[2]) {
      if (verbose) {
        pbcount <- pbcount + 1; setTxtProgressBar(pb, pbcount)
      }
      res <- mclC(Z, W, y, facZ, facW, groups, groupIdx, mu,
                  lam1 = lam.seq[[1]][i], lam2 = lam.seq[[2]][j], thresh, maxit)
      betas[,((i-1)*nlam[2] + j)] <- res$beta
      gammas[,((i-1)*nlam[2] + j)] <- res$gamma
    }
  }
  if (verbose) {
    close(pb)
  }
  if (std == TRUE) {
    int <- attr(y, "scaled:center") -
      drop(crossprod(betas, attr(Z, "scaled:center"))) -
      drop(crossprod(gammas, attr(W, "scaled:center")))
  } else {
    int = 0
  }
  res <- list(int = int, betas = betas, gammas = gammas, lam.seq = lam.seq)
  return(res)
}

