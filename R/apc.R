#' @title Factor Model of Balanced Panel Data
#'
#' @description
#' \code{apc} estiamtes the factor model of a given balanced panel data.
#'
#' @export
#'
#' @param X a matrix of size T by N.
#' @param r integer, indicating the maximum number of factors.
#'
#' @return a list of elements:
#' \item{Fhat}{}
#' \item{Lamhat}{}
#' \item{d}{}
#' \item{d0}{}
#' \item{ehat}{}
#'
#' @author Yankang (Bennie) Chen <yankang.chen@@columbia.edu>
#' @author Serena Ng <serena.ng@@columbia.edu>
#' @author Jushan Bai <jushan.bai@@columbia.edu>
#'
#' @references
#' Jushan Bai and Serena Ng (2019), \emph{Matrix Completion, Counterfactuals, and Factor Analysis of Missing Data}.
#' \url{https://arxiv.org/abs/1910.06677}
#'
#' @examples
#' results <- apc(X, r)



apc <- function(X, r){
  # Error checking
  if (! is.matrix(X))
    try(X <- as.matrix(X))
  if (! is.numeric(r))
    stop("'r' must be an integer.")
  if (r > min(nrow(X), ncol(X)))
    stop("'r' must be smaller than the size of X.")

  # Clear memory and create output object
  gc()
  out <- list()

  T <- nrow(X)
  N <- ncol(X)
  svd_model <- svd(X)
  U <- svd_model$u
  d <- svd_model$d
  V <- svd_model$v
  D <- diag(d)
  D <- D / (sqrt(N*T))
  Dr <- D[1:r, 1:r]

  out$Fhat <- sqrt(T) * U[, 1:r]
  out$Lamhat <- sqrt(N) * V[, 1:r] %*% Dr
  out$d <- d[1:r]
  out$d0 <- d
  out$ehat <- X - out$Fhat %*% t(out$Lamhat)

  return(out)
}
