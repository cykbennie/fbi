#' @title Standard Error of Selected Points (TP)
#'
#' @description
#' \code{se.tp} produces the estimated standard error of C^hat
#' produced by the \code{\link{tp_apc}} function.
#'
#' @import pracma
#' @export
#'
#' @param object an object of class `tp`.
#' @param tpoints integer or vector of integers, indicating t of the
#' (i,t) pair(s) of interest.
#' @param npoints integer or vector of integers, indicating i of the
#' (i,t) pair(s) of interest.
#' @param qq placeholder.
#' @param re_estimate logical. If `FALSE`, use first pass estimation (Lemma 2 of
#' Cahan, Bai, and Ng (2021)). If `TRUE`, use re-estimation (Proposition 1).
#'
#' @return a list of elements:
#' \item{tpoints}{t's of the (i,t) pair(s) of interest}
#' \item{npoints}{i's of the (i,t) pair(s) of interest}
#' \item{re_estimate}{logical. If `FALSE`, use first pass estimation; if `TRUE`,
#' use re-estimation.}
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



se.tp <- function(object, npoints, tpoints, qq, re_estimate = TRUE){
  # Error checking
  if (!inherits(object, "tp"))
    stop("Object must be of class 'tp', i.e. the output of tp_apc.")

  Fhat <- object$Fhat
  Lhat <- object$Lamhat
  Chat <- object$Chat
  ehat <- object$ehat
  Dhat <- object$Dhat
  T <- nrow(Fhat)
  r <- ncol(Fhat)
  N <- nrow(Lhat)
  missingX <- is.na(object$X)   # 1 not observed; 0 observed
  goodN <- which(colSums(missingX)==0)   # goodN: no missing for any period
  No <- length(goodN)

  out <- list()

  out$tpoints <- tpoints
  out$npoints <- npoints
  out$Fhat <- c()
  out$Lhat <- c()
  out$Chat <- c()
  out$SigmaC <- c()
  out$SigmaF <- c()
  out$SigmaL <- c()

  if (re_estimate) {
    # re-estimation using X tilde
    for (ipoints in 1:length(tpoints)) {
      ii <- npoints[ipoints]
      tt <- tpoints[ipoints]
      N_tt <- N - sum(missingX[tt,])
      T_ii <- T - sum(missingX[,ii])
      obsi.tt <- which(missingX[tt,]==0)
      obst.ii <- which(missingX[,ii]==0)
      notobsi.tt <- which(missingX[tt,]==1)
      notobst.ii <- which(missingX[,ii]==1)

      var_Lhato <- t(Lhat[obsi.tt,]) %*% Lhat[obsi.tt,] / No
      var_Lhatm <- t(Lhat[notobsi.tt,]) %*% Lhat[notobsi.tt,] / No
      A_Lam <- var_Lhatm %*% solve(var_Lhato)

      # this assumes that the data are ordered. Otherwise
      # we need to use for i in JT (see check_goodNJ)
      if (i<No) {
        B_Lam <- (diag(r) + A_Lam) %*% (N_tt/N)
      }
      if (i>No & i<N_tt) {
        B_Lam <- diag(r) %*% (N_tt/N)
      }

      Lhat_ii <- matrix(Lhat[ii,], nrow = 1)
      LEhat_t <- pracma::repmat(matrix(ehat[tt,obsi.tt], byrow = FALSE),1,r) * Lhat[obsi.tt,]
      var.Lehat_t <- B_Lam %*% t(LEhat_t) %*% LEhat_t %*% t(B_Lam) / N_tt
      var.Lhat <- t(Lhat) %*% Lhat / N

      Fhat_tt <- matrix(Fhat[tt,], nrow = 1)
      FEhat_i <- pracma::repmat(matrix(ehat[obst.ii,ii], byrow = FALSE),1,r) * Fhat[obst.ii,]
      var.Fehat_i <- t(FEhat_i) %*% FEhat_i / T_ii

      for (k in 1:qq) {
        var.Fehat_i <- var.Fehat_i + t(FEhat_i[(k+1):nrow(FEhat_i),]) %*% FEhat_i[1:(nrow(FEhat_i)-k),] / T_ii
      }

      V0 <- (N_tt/N)^2 * Lhat_ii %*% solve(var.Lhat) %*% var.Lehat_t %*% solve(var.Lhat) %*% t(Lhat_ii)
      W0 <- Fhat_tt %*% var.Fehat_i %*% t(Fhat_tt)

      # Store output
      out$Fhat <- c(out$Fhat, Fhat_tt[1])   # only save the first factor
      out$Lhat <- c(out$Lhat, Lhat_ii[1])   # only save the first loading
      out$Chat <- c(out$Chat, Chat[tt,ii])
      out$SigmaC <- c(out$SigmaC, V0/N_tt + W0/T_ii)
      temp.F <- diag(solve(Dhat) %*% var.Lehat_t %*% solve(Dhat))
      temp.L <- diag(var.Fehat_i)
      out$SigmaF <- c(out$SigmaF, temp.F[1])
      out$SigmaL <- c(out$SigmaL, temp.L[1])
    }

    } else {
    # first pass estimation
    for (ipoints in 1:length(tpoints)) {
      ii <- npoints[ipoints]
      tt <- tpoints[ipoints]
      N_tt <- N - sum(missingX[tt,])
      T_ii <- T - sum(missingX[,ii])
      obsi.tt <- which(missingX[tt,]==0)
      obst.ii <- which(missingX[,ii]==0)

      Lhat_ii <- matrix(Lhat[ii,], nrow = 1)
      LEhat_t <- pracma::repmat(matrix(ehat[tt,obsi.tt], byrow = FALSE),1,r) * Lhat[obsi.tt,]
      var.Lehat_t <- t(LEhat_t) %*% LEhat_t / No
      var.Lhat <- t(Lhat) %*% Lhat / N

      Fhat_tt <- matrix(Fhat[tt,], nrow = 1)
      FEhat_i <- pracma::repmat(matrix(ehat[obst.ii,ii], byrow = FALSE),1,r) * Fhat[obst.ii,]
      var.Fehat_i <- t(FEhat_i) %*% FEhat_i / T_ii

      for (k in 1:qq) {
        var.Fehat_i <- var.Fehat_i + t(FEhat_i[(k+1):nrow(FEhat_i),]) %*% FEhat_i[1:(nrow(FEhat_i)-k),] / T_ii
      }

      # variance of common component
      V0 <- Lhat_ii %*% solve(var.Lhat) %*% var.Lehat_t %*% solve(var.Lhat) %*% t(Lhat_ii)
      W0 <- Fhat_tt %*% var.Fehat_i %*% t(Fhat_tt)

      # Store output
      out$Fhat <- c(out$Fhat, Fhat_tt[1])   # only save the first factor
      out$Lhat <- c(out$Lhat, Lhat_ii[1])   # only save the first loading
      out$Chat <- c(out$Chat, Chat[tt,ii])
      out$SigmaC <- c(out$SigmaC, V0/No + W0/T_ii)
      temp.F <- diag(solve(Dhat) %*% var.Lehat_t %*% solve(Dhat))
      temp.L <- diag(var.Fehat_i)
      out$SigmaF <- c(out$SigmaF, temp.F[1])
      out$SigmaL <- c(out$SigmaL, temp.L[1])
      }
    }
  return(out)
}
