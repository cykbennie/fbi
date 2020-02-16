#' @title Tall-Project Imputation of Missing Value in Panel Data
#'
#' @description
#' \code{tp_apc} imputates the missing values in a given panel data using the
#' method of "Tall-Project".
#'
#' @import stats
#' @export
#'
#' @param X1 a matrix of size T by N.
#' @param r1 integer, indicating the maximum number of factors.
#' @param center logical, indicating Whether or not X1 should be demeaned
#' @param standardize logical, indicating Whether or not X1 should be scaled.
#' @param re_estimate logical, indicating Whether or not output factors,
#' `Fhat`, `Lamhat`, and `Chat`, should be re-estimated from the imputed data.
#'
#' @return a list of elements:
#' \item{Fhat}{}
#' \item{Lamhat}{}
#' \item{Chat}{}
#' \item{data}{}
#'
#' @author Yankang (Bennie) Chen <yankang.chen@@columbia.edu>
#' @author Serena Ng <serena.ng@@columbia.edu>
#' @author Jushan Bai <jushan.bai@@columbia.edu>
#'
#' @references
#' Cahan, E., Bai, J. and Ng, S. 2019, Factor Based Imputation of Missing Data
#' and Covariance Matrix Estimation. unpublished manuscript, Columbia University



tp_apc <- function(X1, r1, center = FALSE, standardize = FALSE,
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
  if (! is.matrix(X1))
    try(X1 <- as.matrix(X1))
  if (! is.numeric(r1))
    stop("'r1' must be an integer.")
  if (r1 > min(nrow(X1), ncol(X1)))
    stop("'r1' must be smaller than the size of X1.")

  # Clear memory and create output object
  gc()
  out <- list()


  T <- nrow(X1)
  N <- ncol(X1)
  rownames(X1) <- 1:T
  colnames(X1) <- 1:N

  missing <- is.na(X1)
  goodT <- rowSums(is.na(X1)) == 0
  goodN <- colSums(is.na(X1)) == 0
  T1 <- sum(goodT)
  N1 <- sum(goodN)
  mu1 <- matrix(rep(colMeans(X1, na.rm = TRUE), T), nrow = T, ncol = N, byrow = TRUE)
  sd1 <- matrix(rep(apply(X1, 2, stats::sd, na.rm = TRUE), T), nrow = T, ncol = N, byrow = TRUE)


  if (center & standardize){
    # demean and standardize

    XT <- (X1[,goodN] - mu1[,goodN]) / sd1[,goodN]
    XN <- (X1 - mu1) / sd1

    bnXT <- fbi::apc(XT, r1)
    Fhat <- bnXT$Fhat
    Lamhat <- matrix(rep(0, N*r1), nrow = N, ncol = r1)

    for (i in 1:N) {
      goodTi <- missing[, i] == FALSE
      Fn1 <- Fhat[goodTi, ]
      Reg <- Fn1
      P <- solve(t(Reg) %*% Reg) %*% t(Reg)
      Lamhat[i, ] <- P %*% XN[goodTi, i]
    }

    Chat <- Fhat %*% t(Lamhat)
    Xhat <- X1   # estimated data
    Xhat[missing] <- Chat[missing] * sd1[missing] + mu1[missing]

    if (re_estimate){
      reest <- fbi::apc(Xhat, r1)
      out$data <- Xhat
      out$Fhat <- reest$Fhat
      out$Lamhat <- reest$Lamhat
      out$Chat <- (out$Fhat %*% t(out$Lamhat))*sd1 + mu1
    } else {
      out$data <- Xhat
      out$Fhat <- Fhat
      out$Lamhat <- Lamhat
      out$Chat <- (out$Fhat %*% t(out$Lamhat))*sd1 + mu1
    }


  } else if (center & (!standardize)){
    # only demean, do not standardize

    XT <- X1[,goodN] - mu1[,goodN]
    XN <- X1 - mu1

    bnXT <- fbi::apc(XT, r1)
    Fhat <- bnXT$Fhat
    Lamhat <- matrix(rep(0, N*r1), nrow = N, ncol = r1)

    for (i in 1:N) {
      goodTi <- missing[, i] == FALSE
      Fn1 <- Fhat[goodTi, ]
      Reg <- Fn1
      P <- solve(t(Reg) %*% Reg) %*% t(Reg)
      Lamhat[i, ] <- P %*% XN[goodTi, i]
    }

    Chat <- Fhat %*% t(Lamhat)
    Xhat <- X1   # estimated data
    Xhat[missing] <- Chat[missing] + mu1[missing]

    if (re_estimate){
      reest <- fbi::apc(Xhat, r1)
      out$data <- Xhat
      out$Fhat <- reest$Fhat
      out$Lamhat <- reest$Lamhat
      out$Chat <- (out$Fhat %*% t(out$Lamhat)) + mu1
    } else {
      out$data <- Xhat
      out$Fhat <- Fhat
      out$Lamhat <- Lamhat
      out$Chat <- (out$Fhat %*% t(out$Lamhat)) + mu1
    }


  } else {
    # no demeaning or standardizing
    XT <- X1[,goodN]
    XN <- X1

    bnXT <- fbi::apc(XT, r1)
    Fhat <- bnXT$Fhat
    Lamhat <- matrix(rep(0, N*(r1+1)), nrow = N, ncol = r1+1)

    for (i in 1:N) {
      goodTi <- missing[, i] == FALSE
      Fn1 <- Fhat[goodTi, ]
      lenTi <- sum(goodTi, na.rm = TRUE)
      Reg <- cbind(matrix(rep(1, lenTi), nrow = lenTi, ncol = 1), Fn1)
      P <- solve(t(Reg) %*% Reg) %*% t(Reg)
      Lamhat[i, ] <- P %*% XN[goodTi, i]
    }

    Chat <- cbind(matrix(rep(1, T), nrow = T, ncol = 1), Fhat) %*% t(Lamhat)
    Xhat <- X1   # estimated data
    Xhat[missing] <- Chat[missing]

    if (re_estimate){
      reest <- fbi::apc(Xhat, r1)
      out$data <- Xhat
      out$Fhat <- reest$Fhat
      out$Lamhat <- reest$Lamhat
      out$Chat <- out$Fhat %*% t(out$Lamhat)
    } else {
      out$data <- Xhat
      out$Fhat <- Fhat
      out$Lamhat <- Lamhat
      out$Chat <- out$Fhat %*% t(out$Lamhat)
    }

  }

  return(out)

}
