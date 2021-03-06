% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tp_apc.R
\name{tp_apc}
\alias{tp_apc}
\title{Tall-Project Imputation of Missing Value in Panel Data}
\usage{
tp_apc(X, kmax, center = FALSE, standardize = FALSE, re_estimate = TRUE)
}
\arguments{
\item{X}{a matrix of size T by N with missing values.}

\item{kmax}{integer, indicating the maximum number of factors.}

\item{center}{logical, indicating whether or not X should be demeaned}

\item{standardize}{logical, indicating whether or not X should be scaled.}

\item{re_estimate}{logical, indicating whether or not output factors,
`Fhat`, `Lamhat`, and `Chat`, should be re-estimated from the imputed data.}
}
\value{
a list of elements:
\item{Fhat}{estimated F}
\item{Lamhat}{estimated Lambda}
\item{Chat}{euqals Fhat x Lamhat'}
\item{data}{X with missing data imputed}
\item{X}{the original data with missing values}
\item{kmax}{the maximum number of factors}
\item{center}{logical, indicating whether or not X was demeaned in the algorithm}
\item{standardize}{logical, indicating whether or not X was scaled in the algorithm}
\item{re_estimate}{logical, indicating whether or not output factors,
`Fhat`, `Lamhat`, and `Chat`, were re-estimated from the imputed data}
}
\description{
\code{tp_apc} imputates the missing values in a given panel data using the
method of "Tall-Project".
}
\references{
Cahan, E., Bai, J. and Ng, S. 2019, Factor Based Imputation of Missing Data
and Covariance Matrix Estimation. unpublished manuscript, Columbia University
}
\author{
Yankang (Bennie) Chen <yankang.chen@yale.edu>

Serena Ng <serena.ng@columbia.edu>

Jushan Bai <jushan.bai@columbia.edu>
}
