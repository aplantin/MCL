#' mcl.grid
#'
#' Calculates MCL coefficient estimates over a grid of (lam1, lam2) pairs
#'
#'
#' @useDynLib MCL
#' @importFrom Rcpp sourceCpp
#' @importFrom stats aggregate sd
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
#' @export
#'
mcl.grid <- function(Z, W, y, mu = 1, groups, lam.seq = NULL, lam.max = NULL, nlam = c(10,10), min.frac = 0.01, maxit = 1e4, thresh = 1e-4, std = TRUE) {
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
    lam.seq[[1]] <- exp(seq(from = log(lam.max[1]), to = log(min.frac*lam.max[1]), length.out = nlam[1]))
    lam.seq[[2]] <- exp(seq(from = log(lam.max[2]), to = log(min.frac*lam.max[2]), length.out = nlam[2]))
  }

  betas <- matrix(nrow = ncol(Z), ncol = (nlam[1]*nlam[2]))
  gammas <- matrix(nrow = ncol(W), ncol = (nlam[1]*nlam[2]))
##  pbcount <- 0 ; pb <- txtProgressBar(min = 0, max = (nlam[1]*nlam[2]), style = 3)
  for (i in 1:nlam[1]) {
    for (j in 1:nlam[2]) {
##      pbcount <- pbcount + 1; setTxtProgressBar(pb, pbcount)
      res <- mclC(Z, W, y, facZ, facW, groups, groupIdx, mu, lam1 = lam.seq[[1]][i], lam2 = lam.seq[[2]][j], thresh, maxit)
      betas[,((i-1)*nlam[2] + j)] <- res$beta
      gammas[,((i-1)*nlam[2] + j)] <- res$gamma
    }
  }
##  close(pb)
  if (std == TRUE) {
    int <- attr(y, "scaled:center") - drop(crossprod(betas, attr(Z, "scaled:center"))) - drop(crossprod(gammas, attr(W, "scaled:center")))
  } else {
    int = 0
  }
  res <- list(int = int, betas = betas, gammas = gammas, lam.seq = lam.seq)
  return(res)
}

#' get.folds
#'
#' Returns folds of data for cross validation
#'
#' @param n Number of subjects
#' @param nfolds Number of folds for cross validation
#' @export
#'
get.folds <- function(n, nfolds) {
  perFold <- floor(n/nfolds)
  leftOver <- n %% nfolds
  reOrd <- sample(n)
  folds <- list()
  last <- 1
  for (i in 1:nfolds) {
    if (i <= leftOver) {
      folds[[i]] <- reOrd[last:(last+perFold)]
      last <- last + perFold + 1
    } else {
      folds[[i]] <- reOrd[last:(last + perFold - 1)]
      last <- last + perFold
    }
  }
  return(folds)
}

#' mcl.cv
#'
#' Cross-validated choice of lam1 and lam2 (over a grid of lambdas)
#'
#' @param Z Otu matrix at group level
#' @param W Otu matrix within groups
#' @param y Outcome vector
#' @param mu Lagrange scaling parameter, default 1
#' @param groups Vector indicating group membership for each feature
#' @param lam.seq Sequence of lambda values, default NULL
#' @param lam.max Maximum lambda value, default NULL
#' @param nlam 2-vector of number of lambdas to test (lam1 and lam2), default (10, 10)
#' @param min.frac Fraction of max.lam that the minimum lambda value should be, default 0.01
#' @param maxit Maximum number of iterations, default 10000
#' @param thresh Threshold for convergence, default 0.0001
#' @param std Logical indicating whether data should be standardized, default TRUE
#' @param nfolds Number of folds for cross validation, default 5
#' @export
#'
mcl.cv <- function(Z, W, y, mu = 1, groups, lam.seq = NULL, lam.max = NULL, nlam = c(10,10), min.frac = 0.01, maxit = 1e4, thresh = 1e-4, std = TRUE, nfolds = 5) {

  groupIdx <- list()
  for (j in 1:length(unique(groups))) {
    groupIdx[[j]] <- which(groups == names(table(groups))[j]) - 1 ## -1 for 0-indexing
  }
  groups = as.numeric(as.factor(groups)) - 1

  if (std == TRUE) {
    y2 <- scale(y, scale = FALSE)
    Z2 <- scale(Z, scale = apply(Z, 2, sd)*sqrt(nrow(Z)-1))
    facZ <- 1/attr(Z2, "scaled:scale")
    W2 <- scale(W, scale = apply(W, 2, sd)*sqrt(nrow(W)-1))
    facW <- 1/attr(W2, "scaled:scale")
  } else {
    y2 <- y
    Z2 <- Z
    W2 <- W
    facZ <- rep(1, ncol(Z))
    facW <- rep(1, ncol(W))
  }

  if (is.null(lam.seq)) {
    if (is.null(lam.max)) {
      lam.max <- findMaxLams(Z, W, y, mu, groups, groupIdx, facZ, facW, maxit, thresh)
    }
    lam.seq <- list()
    lam.seq[[1]] <- exp(seq(from = log(lam.max[1]), to = log(min.frac*lam.max[1]), length.out = nlam[1]))
    lam.seq[[2]] <- exp(seq(from = log(lam.max[2]), to = log(min.frac*lam.max[2]), length.out = nlam[2]))
  }

  folds <- get.folds(n = length(y), nfolds)

  ## cross-validated
  dev <- matrix(nrow = length(y), ncol = (nlam[1]*nlam[2]))
  for (k in 1:nfolds) {
    print(paste("This is fold", k))
    betas <- matrix(nrow = ncol(Z), ncol = (nlam[1]*nlam[2]))
    gammas <- matrix(nrow = ncol(W), ncol = (nlam[1]*nlam[2]))

    ## training data
    train.Z = Z[-folds[[k]], ]
    train.W = W[-folds[[k]], ]
    train.y = y[-folds[[k]]]

    ## standardize training data
    if (std == TRUE) {
      train.y <- scale(train.y, scale = FALSE)
      train.Z <- scale(train.Z, scale = apply(train.Z, 2, sd)*sqrt(nrow(train.Z)-1))
      train.facZ <- 1/attr(train.Z, "scaled:scale")
      train.W <- scale(train.W, scale = apply(train.W, 2, sd)*sqrt(nrow(train.W)-1))
      train.facW <- 1/attr(train.W, "scaled:scale")
    } else {
      train.facZ <- rep(1, ncol(train.Z))
      train.facW <- rep(1, ncol(train.W))
    }

    ## run for each (lam1, lam2) on training data
    for (i in 1:nlam[1]) {
      for (j in 1:nlam[2]) {
        res <- mclC(train.Z, train.W, train.y, train.facZ, train.facW,
                    groups, groupIdx, mu,
                    lam1 = lam.seq[[1]][i], lam2 = lam.seq[[2]][j],
                    thresh, maxit)
        betas[,((i-1)*nlam[2] + j)] <- res$beta
        gammas[,((i-1)*nlam[2] + j)] <- res$gamma
      }
    }
    if (std == TRUE) {
      int <- attr(train.y, "scaled:center") -
        drop(crossprod(betas, attr(train.Z, "scaled:center"))) -
        drop(crossprod(gammas, attr(train.W, "scaled:center")))
    } else {
      int = 0
    }

    ## calculate prediction error on test set
    test.W <- W[folds[[k]], ]
    test.Z <- Z[folds[[k]], ]
    test.y <- y[folds[[k]]]
    for (i in 1:(nlam[1]*nlam[2])) {
      fit.y <- (int[i] + test.Z %*% betas[,i] + test.W %*% gammas[,i])
      dev[folds[[k]], i] <- test.y - fit.y
    }
  }

  pe <- apply(dev, 2, FUN = function(x) sum(x^2)/length(y))
  pe.best <- which.min(pe)

  pe.lam1 <- lam.seq[[1]][ceiling(pe.best/nlam[2])]
  if (pe.best %% nlam[2] == 0) {
    pe.lam2 <- lam.seq[[2]][nlam[2]]
  } else {
    pe.lam2 <- lam.seq[[2]][pe.best %% nlam[2]]
  }

  ## fit to full data with selected lambdas
  fullfit <- mclC(Z2, W2, y2, facZ, facW, groups, groupIdx, mu, lam1 = pe.lam1, lam2 = pe.lam2, thresh, maxit)
  if (std == TRUE) {
    res.int <- attr(y2, "scaled:center") - drop(crossprod(fullfit$beta, attr(Z2, "scaled:center"))) - drop(crossprod(fullfit$gamma, attr(W2, "scaled:center")))
  } else {
    res.int = 0
  }

  res <- list(int = res.int, beta = fullfit$beta, gamma = fullfit$gamma, lam.seq = lam.seq, lam1 = pe.lam1, lam2 = pe.lam2, pe = pe)
  return(res)
}


