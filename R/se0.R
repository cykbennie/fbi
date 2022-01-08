#' @title Standard Error of Selected Points (Baseline)
#'
#' @description
#' \code{se0} produces the estimated standard error of C^hat
#' produced by the \code{\link{tw_apc}} or \code{\link{tp_apc}} function.
#'
#' @import pracma
#' @export
#'
#' @param object an object of class `tw` or `tp`.
#' @param tpoints integer or vector of integers, indicating t of the
#' (i,t) pair(s) of interest.
#' @param npoints integer or vector of integers, indicating i of the
#' (i,t) pair(s) of interest.
#' @param qq placeholder.
#'
#' @return a list of elements:
#' \item{tpoints}{t's of the (i,t) pair(s) of interest}
#' \item{npoints}{i's of the (i,t) pair(s) of interest}
#' \item{Fhat}{estimated F}
#' \item{Lamhat}{estimated Lambda}
#' \item{Chat}{euqals Fhat x Lamhat'}
#' \item{SigmaC}{estimated variance of C}
#' \item{SigmaF}{estimated variance of F}
#' \item{SigmaL}{estimated variance of L}
#'
#' @author Yankang (Bennie) Chen <yankang.chen@@yale.edu>
#' @author Serena Ng <serena.ng@@columbia.edu>
#' @author Jushan Bai <jushan.bai@@columbia.edu>
#'
#' @references
#' Jushan Bai and Serena Ng (2002), \emph{Determining the number of factors in approximate factor models}.
#' \url{https://doi.org/10.1111/1468-0262.00273}
#'
#' Jushan Bai and Serena Ng (2019), \emph{Rank regularized estimation of approximate factor models}.
#' \url{https://doi.org/10.1016/j.jeconom.2019.04.021}
#'
#' Jushan Bai and Serena Ng (2021), \emph{Matrix Completion, Counterfactuals, and Factor Analysis of Missing Data}.
#' \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2021.1967163}
#'
#' Ercument Cahan, Jushan Bai, and Serena Ng (2021),
#' \emph{Factor-Based Imputation of Missing Values and Covariances in Panel Data of Large Dimensions}.
#' \url{https://arxiv.org/abs/2103.03045}



se0 <- function(object, npoints, tpoints, qq){
  # Error checking
  if (!(inherits(object, "tw") | inherits(object, "tp")))
    stop("Object must be of class 'tw' or 'tp', i.e. the output of tw_apc or
         tp_apc.")

  Fhat <- object$Fhat
  Lhat <- object$Lamhat
  Chat <- object$Chat
  ehat <- object$ehat
  Dhat <- object$Dhat
  T <- nrow(Fhat)
  r <- ncol(Fhat)
  N <- nrow(Lhat)

  out <- list()

  out$tpoints <- tpoints
  out$npoints <- npoints
  out$Fhat <- c()
  out$Lhat <- c()
  out$Chat <- c()
  out$SigmaC <- c()
  out$SigmaF <- c()
  out$SigmaL <- c()

  for (ipoints in 1:length(tpoints)) {
    ii <- npoints[ipoints]
    tt <- tpoints[ipoints]

    Lhat_ii <- Lhat[ii,]
    Fhat_tt <- Fhat[tt,]

    LEhat_t <- pracma::repmat(matrix(ehat[tt,], byrow = FALSE),1,r) * Lhat
    var.Lehat_t <- t(LEhat_t) %*% LEhat_t / N
    var.Lhat <- t(Lhat) %*% Lhat / N

    FEhat_i <- pracma::repmat(matrix(ehat[,ii], byrow = FALSE),1,r) * Fhat
    var.Fehat_i <- t(FEhat_i) %*% FEhat_i / T

    for (k in 1:qq) {
      var.Fehat_i <- var.Fehat_i + t(FEhat_i[(k+1):nrow(FEhat_i),]) %*% FEhat_i[1:(nrow(FEhat_i)-k),] / T
    }

    # variance of common component
    V0 <- Lhat_ii %*% solve(var.Lhat) %*% var.Lehat_t %*% solve(var.Lhat) %*% t(Lhat_ii)
    W0 <- Fhat_tt %*% var.Fehat_i %*% t(Fhat_tt)

    # Store output
    out$Fhat <- c(out$Fhat, Fhat_tt[1])   # only save the first factor
    out$Lhat <- c(out$Lhat, Lhat_ii[1])   # only save the first loading
    out$Chat <- c(out$Chat, Chat[tt,ii])
    out$SigmaC <- c(out$SigmaC, V0/N + W0/T)
    temp.F <- diag(solve(Dhat) %*% var.Lehat_t %*% solve(Dhat))
    temp.L <- diag(var.Fehat_i)
    out$SigmaF <- c(out$SigmaF, temp.F[1,1])
    out$SigmaL <- c(out$SigmaL, temp.L[1,1])
  }

  return(out)
}
