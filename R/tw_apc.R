#' @title Tall-Wide Imputation of Missing Value in Panel Data
#'
#' @description
#' \code{tw_apc} imputates the missing values in a given panel data using the
#' method of "Tall-Wide".
#'
#' @import stats
#' @import pracma
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
#'
#' @author Yankang (Bennie) Chen <yankang.chen@@yale.edu>
#' @author Serena Ng <serena.ng@@columbia.edu>
#' @author Jushan Bai <jushan.bai@@columbia.edu>
#'
#' @references
#' Jushan Bai and Serena Ng (2019), \emph{Matrix Completion, Counterfactuals, and Factor Analysis of Missing Data}.
#' \url{https://arxiv.org/abs/1910.06677}



tw_apc <- function(X, kmax, center = FALSE, standardize = FALSE,
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
    XN <- (X[goodT,] - mu1[goodT,]) / sd1[goodT,]

    bnXT <- fbi::apc(XT, kmax)
    bnXN <- fbi::apc(XN, kmax)

    HH <- pracma::mldivide(bnXN$Lamhat[1:N1, 1:kmax], bnXT$Lamhat[1:N1, 1:kmax])
    Lamhat <- bnXN$Lamhat %*% HH
    Fhat <- bnXT$Fhat
    Chat <- bnXT$Fhat %*% t(Lamhat)
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
      out$Chat <- Chat*sd1 + mu1
    }


  } else if (center & (!standardize)){
    # only demean, do not standardize

    XT <- X[,goodN] - mu1[,goodN]
    XN <- X[goodT,] - mu1[goodT,]

    bnXT <- fbi::apc(XT, kmax)
    bnXN <- fbi::apc(XN, kmax)

    HH <- pracma::mldivide(bnXN$Lamhat[1:N1, 1:kmax], bnXT$Lamhat[1:N1, 1:kmax])
    Lamhat <- bnXN$Lamhat %*% HH
    Fhat <- bnXT$Fhat
    Chat <- bnXT$Fhat %*% t(Lamhat)
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
    XN <- X[goodT,]

    bnXT <- fbi::apc(XT, kmax)
    bnXN <- fbi::apc(XN, kmax)

    HH <- pracma::mldivide(bnXN$Lamhat[1:N1, 1:kmax], bnXT$Lamhat[1:N1, 1:kmax])
    Lamhat <- bnXN$Lamhat %*% HH
    Fhat <- bnXT$Fhat
    Chat <- bnXT$Fhat %*% t(Lamhat)
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
      out$Lamhat <- Lamhat
      out$Chat <- Chat
    }

  }

  class(out) <- c("list", "tptw")
  return(out)

}
