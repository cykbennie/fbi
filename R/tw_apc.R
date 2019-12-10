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
#' @param X1 a matrix of size T by N.
#' @param r1 integer, indicating the maximum number of factors.
#' @param center logical, indicating Whether or not X1 should be demeaned
#' @param standardize logical, indicating Whether or not X1 should be scaled.
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
#' Jushan Bai and Serena Ng (2019), \emph{Matrix Completion, Counterfactuals, and Factor Analysis of Missing Data}.
#' \url{https://arxiv.org/abs/1910.06677}



tw_apc <- function(X1, r1, center = FALSE, standardize = FALSE) {
  # Error checking
  if (! is.logical(center))
    stop("'center' must be logical.")
  if (! is.logical(standardize))
    stop("'standardize' must be logical.")
  if ((! center) & (standardize))
    stop("The option 'center = FALSE, standardize = TRUE' is not available.")
  if (! is.matrix(X1))
    try(X1 <- as.matrix(X1))
  if (! is.numeric(r1))
    stop("'r1' must be an integer.")
  if (r1 > min(nrow(X1), ncol(X1)))
    stop("'r1' must be smaller than the size of X1.")

  # Create output object
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
    XN <- (X1[goodT,] - mu1[goodT,]) / sd1[goodT,]

    bnXT <- fbi::apc(XT, r1)
    bnXN <- fbi::apc(XN, r1)

    HH <- pracma::mldivide(bnXN$Lamhat[1:N1, 1:r1], bnXT$Lamhat[1:N1, 1:r1])
    Lamhat <- bnXN$Lamhat %*% HH
    Fhat <- bnXT$Fhat
    Chat <- bnXT$Fhat %*% t(Lamhat)
    Xhat <- X1   # estimated data
    Xhat[missing] <- Chat[missing] * sd1[missing] + mu1[missing]

    reest <- fbi::apc(Xhat, r1)
    out$data <- Xhat
    out$Fhat <- reest$Fhat
    out$Lamhat <- reest$Lamhat
    out$Chat <- (out$Fhat %*% t(out$Lamhat))*sd1 + mu1


  } else if (center & (!standardize)){
    # only demean, do not standardize

    XT <- X1[,goodN] - mu1[,goodN]
    XN <- X1[goodT,] - mu1[goodT,]

    bnXT <- fbi::apc(XT, r1)
    bnXN <- fbi::apc(XN, r1)

    HH <- pracma::mldivide(bnXN$Lamhat[1:N1, 1:r1], bnXT$Lamhat[1:N1, 1:r1])
    Lamhat <- bnXN$Lamhat %*% HH
    Fhat <- bnXT$Fhat
    Chat <- bnXT$Fhat %*% t(Lamhat)
    Xhat <- X1   # estimated data
    Xhat[missing] <- Chat[missing] + mu1[missing]

    reest <- fbi::apc(Xhat, r1)
    out$data <- Xhat
    out$Fhat <- reest$Fhat
    out$Lamhat <- reest$Lamhat
    out$Chat <- (out$Fhat %*% t(out$Lamhat)) + mu1


  } else {
    # no demeaning or standardizing

    XT <- X1[,goodN]
    XN <- X1[goodT,]

    bnXT <- fbi::apc(XT, r1)
    bnXN <- fbi::apc(XN, r1)

    HH <- pracma::mldivide(bnXN$Lamhat[1:N1, 1:r1], bnXT$Lamhat[1:N1, 1:r1])
    Lamhat <- bnXN$Lamhat %*% HH
    Fhat <- bnXT$Fhat
    Chat <- bnXT$Fhat %*% t(Lamhat)
    Xhat <- X1   # estimated data
    Xhat[missing] <- Chat[missing]

    reest <- fbi::apc(Xhat, r1)
    out$data <- Xhat
    out$Fhat <- reest$Fhat
    out$Lamhat <- reest$Lamhat
    out$Chat <- out$Fhat %*% t(out$Lamhat)

  }

  return(out)

}
