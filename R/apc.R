#' @title Factor Model of Balanced Panel Data
#'
#' @description
#' \code{apc} estimates the factor model of a given balanced panel data.
#'
#' @export
#'
#' @param X a matrix of size T by N.
#' @param kmax integer, indicating the maximum number of factors.
#'
#' @return a list of elements:
#' \item{X}{the original data}
#' \item{kmax}{the maximum number of factors}
#' \item{Fhat}{estimated F}
#' \item{Lamhat}{estimated Lambda}
#' \item{Chat}{equals Fhat x Lamhat'}
#' \item{Dhat}{estimated diagonal matrix D, of dim kmax by kmax}
#' \item{d}{first kmax elements of Dhat}
#' \item{d0}{diagonal elements of Dhat}
#' \item{ehat}{equals X - Chat}
#'
#' @author Yankang (Bennie) Chen <yankang.chen@@yale.edu>
#' @author Serena Ng <serena.ng@@columbia.edu>
#' @author Jushan Bai <jushan.bai@@columbia.edu>
#'
#' @references
#' Jushan Bai and Serena Ng (2002), \emph{Determining the number of factors in approximate factor models}.
#' \url{https://doi.org/10.1111/1468-0262.00273}



apc <- function(X, kmax){
  # Error checking
  if (! is.matrix(X))
    try(X <- as.matrix(X))
  if (! is.numeric(kmax))
    stop("'kmax' must be an integer.")
  if (kmax > min(nrow(X), ncol(X)))
    stop("'kmax' must be smaller than the size of X.")

  # Create output object
  out <- list()
  out$X <- X
  out$kmax <- kmax

  T <- nrow(X)
  N <- ncol(X)
  svd_model <- svd(X)
  U <- svd_model$u
  d <- svd_model$d
  V <- svd_model$v
  D <- diag(d)
  D <- D / (sqrt(N*T))
  Dr <- D[1:kmax, 1:kmax, drop = FALSE]

  out$Fhat <- sqrt(T) * U[, 1:kmax, drop = FALSE]
  out$Lamhat <- sqrt(N) * V[, 1:kmax, drop = FALSE] %*% Dr
  out$Chat <- out$Fhat %*% t(out$Lamhat)
  out$d <- d[1:kmax]
  out$d0 <- d
  out$Dhat <- diag(out$d)
  out$ehat <- out$X - out$Chat

  return(out)
}
