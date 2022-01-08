#' @title Standard Error of Selected Points (TW)
#'
#' @description
#' \code{se.tw} produces the estimated standard error of C^hat
#' produced by the \code{\link{tw_apc}} function.
#'
#' @import pracma
#' @export
#'
#' @param object an object of class `tw`.
#' @param tpoints integer or vector of integers, indicating t of the
#' (i,t) pair(s) of interest.
#' @param npoints integer or vector of integers, indicating i of the
#' (i,t) pair(s) of interest.
#' @param qq placeholder.
#' @param re_estimate logical. If `FALSE`, use first pass estimation. If `TRUE`,
#' use re-estimation.
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



se.tw <- function(object, npoints, tpoints, qq, re_estimate){
  # Error checking
  if (!inherits(object, "tw"))
    stop("Object must be of class 'tw', i.e. the output of tw_apc.")

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
  goodT <- which(rowSums(missingX)==0)   # goodT: no missing individuals
  To <- length(goodT)
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

      if (ii<=No) {# tall block
        TT <- T
        NN <- No
      }

      if (tt<=To) {# wide block
        TT <- To
        NN <- N
      }

      if (ii<=No & tt<=To) {
        TT <- T
        NN <- N
      }

      if (ii>No & tt>To) {# miss
        TT <- T - To
        NN <- N - No
      }

      Lhat_ii <- matrix(Lhat[ii,], nrow = 1)
      Fhat_tt <- matrix(Fhat[tt,], nrow = 1)
      LEhat_t <- pracma::repmat(matrix(ehat[tt,], byrow = FALSE),1,r) * Lhat
      var.Lehat_t <- t(LEhat_t) %*% LEhat_t / N
      var.Lhat <- t(Lhat) %*% Lhat / N

      FEhat_i <- pracma::repmat(matrix(ehat[,ii], byrow = FALSE),1,r) * Fhat
      var.Fehat_i <- t(FEhat_i) %*% FEhat_i / T

      for (k in 1:qq) {
        var.Fehat_i <- var.Fehat_i + t(FEhat_i[(k+1):nrow(FEhat_i),]) %*% FEhat_i[1:(nrow(FEhat_i)-k),] / T
      }

      B.Lam <- diag(r)
      B.F <- diag(r)

      var_Lhato <- t(Lhat[1:No,]) %*% Lhat[1:No,] / No
      var_Lhatm <- t(Lhat[(No+1):N,]) %*% Lhat[(No+1):N,] / (N-No)
      var_Fhato <- t(Fhat[1:To,]) %*% Fhat[1:To,] / To
      var_Fhatm <- t(Fhat[(To+1):T,]) %*% Fhat[(To+1):T,] / (T-To)

      if (ii<No) {
        B.Lam <- (No/N) * diag(r) + (1-No/N) * var_Lhatm %*% solve(var_Lhato)
      }
      if (tt<To) {
        B.F <- (To/T) * diag(r) + (1-To/T) * var_Fhatm %*% solve(var_Fhato)
      }

      # variance of common component
      V0 <- Lhat_ii %*% solve(var.Lhat) %*% B.Lam %*% var.Lehat_t %*% t(B.Lam) %*% solve(var.Lhat) %*% t(Lhat_ii)
      W0 <- Fhat_tt %*% B.F %*% var.Fehat_i %*% t(B.F) %*% t(Fhat_tt)

      # Store output
      out$Fhat <- c(out$Fhat, Fhat_tt[1])   # only save the first factor
      out$Lhat <- c(out$Lhat, Lhat_ii[1])   # only save the first loading
      out$Chat <- c(out$Chat, Chat[tt,ii])
      out$SigmaC <- c(out$SigmaC, V0/NN + W0/TT)
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

      if (ii<=No) {# tall not in bal block
        TT <- T
        NN <- No
      }

      if (tt<=To) {# wide not in bal block
        TT <- To
        NN <- N
      }

      if (ii<=No & tt<=To) {# bal block
        if (min(No,T) > min(N,To)) {
          TT <- T
          NN <- No
        }
        if (min(No,T) <= min(N,To)) {
          TT <- To
          NN <- N
        }
      }

      if (ii>No & tt>To) {# miss
        TT <- T - To
        NN <- N - No
      }

      Lhat_ii <- matrix(Lhat[ii,], nrow = 1)
      Fhat_tt <- matrix(Fhat[tt,], nrow = 1)
      LEhat_t <- pracma::repmat(matrix(ehat[tt,obsi.tt], byrow = FALSE),1,r) * Lhat[obsi.tt,]
      var.Lehat_t <- t(LEhat_t) %*% LEhat_t / NN
      var.Lhat <- t(Lhat) %*% Lhat / N

      FEhat_i <- pracma::repmat(matrix(ehat[obst.ii,ii], byrow = FALSE),1,r) * Fhat[obst.ii,]
      var.Fehat_i <- t(FEhat_i) %*% FEhat_i / TT

      for (k in 1:qq) {
        var.Fehat_i <- var.Fehat_i + t(FEhat_i[(k+1):nrow(FEhat_i),]) %*% FEhat_i[1:(nrow(FEhat_i)-k),] / TT
      }

      # variance of common component
      V0 <- Lhat_ii %*% solve(var.Lhat) %*% var.Lehat_t %*% solve(var.Lhat) %*% t(Lhat_ii)
      W0 <- Fhat_tt %*% var.Fehat_i %*% t(Fhat_tt)

      # Store output
      out$Fhat <- c(out$Fhat, Fhat_tt[1])   # only save the first factor
      out$Lhat <- c(out$Lhat, Lhat_ii[1])   # only save the first loading
      out$Chat <- c(out$Chat, Chat[tt,ii])
      out$SigmaC <- c(out$SigmaC, V0/NN + W0/TT)
      temp.F <- diag(solve(Dhat) %*% var.Lehat_t %*% solve(Dhat))
      temp.L <- diag(var.Fehat_i)
      out$SigmaF <- c(out$SigmaF, temp.F[1])
      out$SigmaL <- c(out$SigmaL, temp.L[1])
      }
    }
  return(out)
}
