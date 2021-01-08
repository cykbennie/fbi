#' @title Residual Overlay
#'
#' @description
#' \code{res_overlay.twtp} estimates the covariance and correlation matrix of the
#' unbalanced panel data using the method of residual overlay.
#'
#' @import stats
#' @import pracma
#' @export
#'
#' @param object an object of class `tptw`, i.e. the output of
#' \code{\link{tp_apc}} or \code{\link{tw_apc}}.
#' @param method integer 1 to 4, indicating which residual overlay
#' method to use. They correspond to the four methods described in the paper.
#' @param S the number of iterations.
#'
#' @return a list of elements:
#' \item{method}{the method of residual overlay}
#' \item{S}{the number of iterations}
#' \item{cov}{estimated covariance matrix}
#' \item{cor}{estimated correlation matrix}
#'
#' @author Yankang (Bennie) Chen <yankang.chen@@yale.edu>
#' @author Serena Ng <serena.ng@@columbia.edu>
#' @author Jushan Bai <jushan.bai@@columbia.edu>
#'
#' @references
#' Cahan, E., Bai, J. and Ng, S. 2019, Factor Based Imputation of Missing Data
#' and Covariance Matrix Estimation. unpublished manuscript, Columbia University



res_overlay.twtp <- function(object, method = 1, S = 500) {
  # Error checking
  if (!inherits(object, "tptw"))
    stop("Object must be of class 'tptw', i.e. the output of tp_apc or
         tw_apc.")
  if (!method %in% c(1,2,3,4))
    stop("'method' must be of integer 1, 2, 3, or 4.")
  if (! is.numeric(S))
    stop("'S' must be an integer.")



  # Preparation
  TW0 <- object
  missingX <- is.na(TW0$X)
  notmissingX <- !missingX
  N <- ncol(TW0$data)
  goodT <- which(rowSums(missingX)==0)
  goodN <- which(colSums(missingX)==0)
  badN <- which(colSums(missingX)!=0)
  badT <- which(rowSums(missingX)!=0)

  TW0$ehat <- TW0$data - TW0$Chat
  TW0$cov <- cov(TW0$data)
  TW0$cor <- cor(TW0$data)

  dum1 <- pracma::zeros(N,N)
  dum2 <- pracma::zeros(N,N)
  RW0 <- list()


  # Residual overlay
  if (method == 1){
    # residual overlay 1
    # iid sampling
    uhat <- list()
    uhat$obs <- TW0$ehat[notmissingX]
    No <- length(uhat$obs)
    Nm <- length(TW0$ehat[missingX])

    for (s in 1:S){
      RW0$uhat <- TW0$ehat
      trial <- pracma::randsample(No, Nm)
      RW0$uhat[missingX] <- uhat$obs[trial]
      RW0$data <- TW0$data
      RW0$data[missingX] <- TW0$data[missingX] + RW0$uhat[missingX]
      dum1 <- dum1 + cov(RW0$data)
      dum2 <- dum2 + cor(RW0$data)
    } # end s

    RW0$cov <- dum1/S
    RW0$cor <- dum2/S


  } else if (method == 2){
    # residual overlay 2
    # sampling by columns
    for (s in 1:S){
      RW0$data <- TW0$data
      RW0$uhat <- TW0$ehat
      for (j in 1:length(badN)){
        jj <- badN[j]
        UU1 <- TW0$ehat[,jj]
        UUgood <- which(!missingX[,jj])
        UUbad <- which(missingX[,jj])
        UU2 <- UU1[UUgood]
        trial <- pracma::randsample(length(UUgood), length(UUbad),
                                    replacement = TRUE)
        RW0$uhat[UUbad,jj] <- UU2[trial]
        RW0$data[UUbad,jj] <- RW0$data[UUbad,jj] + RW0$uhat[UUbad,jj]
      } # end j
      dum1 <- dum1 + cov(RW0$data)
      dum2 <- dum2 + cor(RW0$data)
    } # end s

    RW0$cov <- dum1/S
    RW0$cor <- dum2/S


  } else if (method == 3){
    # residual overlay 3
    # iid sampling
    uhat <- list()
    uhat$obs <- TW0$ehat[notmissingX]
    uhat$mu <- mean(uhat$obs)
    uhat$sd <- sd(uhat$obs)
    No <- length(uhat$obs)
    Nm <- length(TW0$ehat[missingX])

    for (s in 1:S) {
      RW0$uhat <- TW0$ehat
      RW0$uhat[missingX] <- pracma::randn(Nm,1) * uhat$sd + uhat$mu
      RW0$data <- TW0$data
      RW0$data[missingX] <- TW0$data[missingX] + RW0$uhat[missingX]
      dum1 <- dum1 + cov(RW0$data)
      dum2 <- dum2 + cor(RW0$data)
    } # end s

    RW0$cov <- dum1/S
    RW0$cor <- dum2/S


  } else {
    # residual overlay 4
    # sampling by columns
    for (s in 1:S){
      RW0$data <- TW0$data
      RW0$uhat <- TW0$ehat
      for (j in 1:length(badN)){
        jj <- badN[j]
        UU1 <- TW0$ehat[,jj]
        UUgood <- which(!missingX[,jj])
        UUbad <- which(missingX[,jj])
        UU2 <- UU1[UUgood]
        uu2.mu <- mean(UU2)
        uu2.sd <- std(UU2)
        nbad <- length(UUbad)
        RW0$uhat[UUbad,jj] <- randn(nbad,1)*uu2.sd + uu2.mu
        RW0$data[UUbad,jj] <- RW0$data[UUbad,jj] + RW0$uhat[UUbad,jj]
      } # end j
      dum1 <- dum1 + cov(RW0$data)
      dum2 <- dum2 + cor(RW0$data)
    } # end s

    RW0$cov <- dum1/S
    RW0$cor <- dum2/S

  } # end residual overlay


  out <- list()
  out$method <- method
  out$S <- S
  out$cov <- RW0$cov
  out$cor <- RW0$cor

  return(out)

}
