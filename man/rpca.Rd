% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rpca.R
\name{rpca}
\alias{rpca}
\title{Estimation of Approximate Factor Models}
\usage{
rpca(X, kmax, standardize = FALSE, tau = 0)
}
\arguments{
\item{X}{a matrix of size T by N.}

\item{kmax}{integer, indicating the maximum number of factors.}

\item{standardize}{logical, indicating Whether or not X should be centered and scaled.}

\item{tau}{numeric, specifying the parameter in the rank-regularized estimation.
If \code{tau = 0}, then rank regularization is not used.}
}
\value{
a list of elements:
\item{X}{}
\item{kmax}{}
\item{standardize}{}
\item{tau}{}
\item{ic2}{}
\item{pc2k}{}
\item{pc20}{}
\item{Fhat}{}
\item{Lamhat}{}
\item{Chat}{}
\item{Sigma}{}
\item{IC2}{}
\item{PC2k}{}
\item{PC20}{}
\item{fhat}{}
\item{lamhat}{}
\item{d}{}
\item{d0}{}
}
\description{
\code{rpca} estimates the approximate factor models of the given matrix.
}
\references{
Jushan Bai and Serena Ng (2002), \emph{Determining the number of factors in approximate factor models}.
\url{https://doi.org/10.1111/1468-0262.00273}

Jushan Bai and Serena Ng (2019), \emph{Rank regularized estimation of approximate factor models}.
\url{https://doi.org/10.1016/j.jeconom.2019.04.021}
}
\author{
Yankang (Bennie) Chen <yankang.chen@yale.edu>

Serena Ng <serena.ng@columbia.edu>

Jushan Bai <jushan.bai@columbia.edu>
}
