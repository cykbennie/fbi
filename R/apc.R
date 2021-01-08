#' @title Factor Model of Balanced Panel Data
#'
#' @description
#' \code{apc} estiamtes the factor model of a given balanced panel data.
#'
#' @export
#'
#' @param X a matrix of size T by N.
#' @param kmax integer, indicating the maximum number of factors.
#'
#' @return a list of elements:
#' \item{Fhat}{}
#' \item{Lamhat}{}
#' \item{Chat}{}
#' \item{d}{}
#' \item{d0}{}
#' \item{ehat}{}
#' \item{Chat}{euqals Fhat x Lamhat'}
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

  T <- nrow(X)
  N <- ncol(X)
  svd_model <- svd(X)
  U <- svd_model$u
  d <- svd_model$d
  V <- svd_model$v
  D <- diag(d)
  D <- D / (sqrt(N*T))
  Dr <- D[1:kmax, 1:kmax]

  out$Fhat <- sqrt(T) * U[, 1:kmax]
  out$Lamhat <- sqrt(N) * V[, 1:kmax] %*% Dr
  out$Chat <- out$Fhat %*% t(out$Lamhat)
  out$d <- d[1:kmax]
  out$d0 <- d
  out$ehat <- X - out$Chat
  out$kmax <- kmax

  return(out)
}
