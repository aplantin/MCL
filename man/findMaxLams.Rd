% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{findMaxLams}
\alias{findMaxLams}
\title{findMaxLams}
\usage{
findMaxLams(Z, W, y, facZ, facW, groups, groupIdx, mu, maxit, thresh)
}
\arguments{
\item{Z}{OTU matrix (group level)}

\item{W}{OTU matrix (within-group level)}

\item{y}{Vector of outcomes}

\item{facZ}{Scaling factors for Z}

\item{facW}{Scaling factors for W}

\item{groups}{Vector indicating group membership for each feature}

\item{groupIdx}{List indicating feature indices in each group}

\item{mu}{Lagrange scaling parameter}

\item{maxit}{Maximum number of iterations allowed}

\item{thresh}{Parameter defining convergence}
}
\description{
Finds maximum lambda value (smallest lambda such
that all coefficients are zero) for both lam1 and lam2
}
