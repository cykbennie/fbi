#' @title Tall-Project Imputation of Missing Value in Panel Data
#'
#' @description
#' \code{tp_apc} imputates the missing values in a given panel data using the
#' method of "Tall-Project".
#'
#' @import stats
#' @export
#'
#' @param X a matrix of size T by N with missing values.
#' @param kmax integer, indicating the maximum number of factors.
#' @param center logical, indicating whether or not X should be demeaned
#' @param standardize logical, indicating whether or not X should be scaled.
#' @param re_estimate logical, indicating whether or not output factors,
#' `Fhat`, `Lamhat`, and `Chat`, should be re-estimated from the imputed data.
#'
#' @return a list of elements:
#' \item{Fhat}{estimated F}
#' \item{Lamhat}{estimated Lambda}
#' \item{Chat}{euqals Fhat x Lamhat'}
#' \item{data}{X with missing data imputed}
#' \item{X}{the original data with missing values}
#' \item{kmax}{the maximum number of factors}
#' \item{center}{logical, indicating whether or not X was demeaned in the algorithm}
#' \item{standardize}{logical, indicating whether or not X was scaled in the algorithm}
#' \item{re_estimate}{logical, indicating whether or not output factors,
#' `Fhat`, `Lamhat`, and `Chat`, were re-estimated from the imputed data}
#'
#' @author Yankang (Bennie) Chen <yankang.chen@@yale.edu>
#' @author Serena Ng <serena.ng@@columbia.edu>
#' @author Jushan Bai <jushan.bai@@columbia.edu>
#'
#' @references
#' Cahan, E., Bai, J. and Ng, S. 2019, Factor Based Imputation of Missing Data
#' and Covariance Matrix Estimation. unpublished manuscript, Columbia University



tp_apc <- function(X, kmax, center = FALSE, standardize = FALSE,
                   re_estimate = TRUE) {
  # Error checking
  if (! is.logical(center))
    stop("'center' must be logical.")
  if (! is.logical(standardize))
    stop("'standardize' must be logical.")
  if (! is.logical(re_estimate))
    stop("'re_estimate' must be logical.")
  if ((! center) & (standardize))
    stop("The option 'center = FALSE, standardize = TRUE' is not available.")
  if (! is.matrix(X))
    try(X <- as.matrix(X))
  if (! is.numeric(kmax))
    stop("'kmax' must be an integer.")
  if (kmax > min(nrow(X), ncol(X)))
    stop("'kmax' must be smaller than the size of X.")

  # Clear memory and create output object
  gc()
  out <- list()

  out$X <- X
  out$kmax <- kmax
  out$center <- center
  out$standardize <- standardize
  out$re_estimate <- re_estimate


  T <- nrow(X)
  N <- ncol(X)
  rownames(X) <- 1:T
  colnames(X) <- 1:N

  missing <- is.na(X)
  goodT <- rowSums(is.na(X)) == 0
  goodN <- colSums(is.na(X)) == 0
  T1 <- sum(goodT)
  N1 <- sum(goodN)
  mu1 <- matrix(rep(colMeans(X, na.rm = TRUE), T), nrow = T, ncol = N, byrow = TRUE)
  sd1 <- matrix(rep(apply(X, 2, stats::sd, na.rm = TRUE), T), nrow = T, ncol = N, byrow = TRUE)


  if (center & standardize){
    # demean and standardize

    XT <- (X[,goodN] - mu1[,goodN]) / sd1[,goodN]
    XN <- (X - mu1) / sd1

    bnXT <- fbi::apc(XT, kmax)
    Fhat <- bnXT$Fhat
    Lamhat <- matrix(rep(0, N*kmax), nrow = N, ncol = kmax)

    for (i in 1:N) {
      goodTi <- missing[, i] == FALSE
      Fn1 <- Fhat[goodTi, ]
      Reg <- Fn1
      P <- solve(t(Reg) %*% Reg) %*% t(Reg)
      Lamhat[i, ] <- P %*% XN[goodTi, i]
    }

    Chat <- Fhat %*% t(Lamhat)
    Xhat <- X   # estimated data
    Xhat[missing] <- Chat[missing] * sd1[missing] + mu1[missing]

    if (re_estimate){
      reest <- fbi::apc(Xhat, kmax)
      out$data <- Xhat
      out$Fhat <- reest$Fhat
      out$Lamhat <- reest$Lamhat
      out$Chat <- (out$Fhat %*% t(out$Lamhat))*sd1 + mu1
    } else {
      out$data <- Xhat
      out$Fhat <- Fhat
      out$Lamhat <- Lamhat
      out$Chat <- Chat * sd1 + mu1
    }


  } else if (center & (!standardize)){
    # only demean, do not standardize

    XT <- X[,goodN] - mu1[,goodN]
    XN <- X - mu1

    bnXT <- fbi::apc(XT, kmax)
    Fhat <- bnXT$Fhat
    Lamhat <- matrix(rep(0, N*kmax), nrow = N, ncol = kmax)

    for (i in 1:N) {
      goodTi <- missing[, i] == FALSE
      Fn1 <- Fhat[goodTi, ]
      Reg <- Fn1
      P <- solve(t(Reg) %*% Reg) %*% t(Reg)
      Lamhat[i, ] <- P %*% XN[goodTi, i]
    }

    Chat <- Fhat %*% t(Lamhat)
    Xhat <- X   # estimated data
    Xhat[missing] <- Chat[missing] + mu1[missing]

    if (re_estimate){
      reest <- fbi::apc(Xhat, kmax)
      out$data <- Xhat
      out$Fhat <- reest$Fhat
      out$Lamhat <- reest$Lamhat
      out$Chat <- (out$Fhat %*% t(out$Lamhat)) + mu1
    } else {
      out$data <- Xhat
      out$Fhat <- Fhat
      out$Lamhat <- Lamhat
      out$Chat <- Chat + mu1
    }


  } else {
    # no demeaning or standardizing
    XT <- X[,goodN]
    XN <- X

    bnXT <- fbi::apc(XT, kmax)
    Fhat <- bnXT$Fhat
    Lamhat <- matrix(rep(0, N*(kmax+1)), nrow = N, ncol = kmax+1)

    for (i in 1:N) {
      goodTi <- missing[, i] == FALSE
      Fn1 <- Fhat[goodTi, ]
      lenTi <- sum(goodTi, na.rm = TRUE)
      Reg <- cbind(matrix(rep(1, lenTi), nrow = lenTi, ncol = 1), Fn1)
      P <- solve(t(Reg) %*% Reg) %*% t(Reg)
      Lamhat[i, ] <- P %*% XN[goodTi, i]
    }

    Chat <- cbind(matrix(rep(1, T), nrow = T, ncol = 1), Fhat) %*% t(Lamhat)
    Xhat <- X   # estimated data
    Xhat[missing] <- Chat[missing]

    if (re_estimate){
      reest <- fbi::apc(Xhat, kmax)
      out$data <- Xhat
      out$Fhat <- reest$Fhat
      out$Lamhat <- reest$Lamhat
      out$Chat <- out$Fhat %*% t(out$Lamhat)
    } else {
      out$data <- Xhat
      out$Fhat <- Fhat
      out$Lamhat <- Lamhat[,2:(kmax+1)]
      out$Chat <- Chat
    }

  }

  class(out) <- c("list", "tptw")
  return(out)

}
