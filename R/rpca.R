#' @title Estimation of Approximate Factor Models
#'
#' @description
#' \code{rpca} estimates the approximate factor models of the given matrix.
#'
#' @export
#'
#' @param X a matrix of size T by N.
#' @param kmax integer, indicating the maximum number of factors.
#' @param standardize logical, indicating Whether or not X should be centered and scaled.
#' @param tau numeric, specifying the parameter in the rank-regularized estimation.
#' If \code{tau = 0}, then rank regularization is not used.
#'
#' @return a list of elements:
#' \item{X}{}
#' \item{kmax}{}
#' \item{standardize}
#' \item{tau}{}
#' \item{ic2}{}
#' \item{pc2k}{}
#' \item{pc20}{}
#' \item{Fhat}{}
#' \item{Lamhat}{}
#' \item{Chat}{}
#' \item{Sigma}{}
#' \item{IC2}{}
#' \item{PC2k}{}
#' \item{PC20}{}
#' \item{fhat}{}
#' \item{lamhat}{}
#' \item{d}{}
#' \item{d0}{}
#'
#' @author Yankang (Bennie) Chen <yankang.chen@@columbia.edu>
#' @author Serena Ng <serena.ng@@columbia.edu>
#' @author Jushan Bai <jushan.bai@@columbia.edu>
#'
#' @references
#' Jushan Bai and Serena Ng (2002), \emph{Determining the number of factors in approximate factor models}.
#' \url{https://onlinelibrary.wiley.com/doi/pdf/10.1111/1468-0262.00273}
#'
#' Jushan Bai and Serena Ng (2017), \emph{Rank regularized estimation of approximate factor models}.
#' \url{https://www.sciencedirect.com/science/article/pii/S0304407619300764}
#'
#' @examples
#' results <- rpca(X, kmax, standardize = FALSE, tau = None)
#' summary(results)
#' results$d



rpca <- function(X, kmax, standardize = FALSE, tau = 0) {
  # Error checking
  if (! is.numeric(kmax))
    stop("'kmax' must be an integer.")
  if (! is.numeric(tau))
    stop("'tau' must be numeric.")
  if (! is.matrix(X))
    try(X <- as.matrix(X))
  if (kmax > min(nrow(X), ncol(X)))
    stop("'kmax' must be smaller than the size of X.")

  # Clear memory and create output object
  gc()
  out <- list()
  out$X <- X
  out$kmax <- kmax
  out$standardize <- standardize
  out$tau <- tau

  if (standardize)
    X <- scale(X, center = TRUE, scale = TRUE)

  T <- nrow(X)
  N <- ncol(X)
  ii <- 1:kmax
  ii[kmax+1] <- 0   # 0 factors is put last
  NT <- N*T
  NT1 <- N+T
  CT <- (NT1/NT) * log(min(N,T)) * ii
  CT[kmax+1] <- 0


  if (tau != 0) {   # baing19 (rank regularized)
    print("Bai Ng 2019 (rank regularized)", quote = FALSE)

    svdx <- svd(X)
    U <- svdx$u
    d0 <- svdx$d
    D <- diag(d0)
    V <- svdx$v
    d0 <- diag(D)
    total <- sum(d0^2)

    Du <- D[1:kmax,1:kmax]
    Dr <- Du - tau * diag(kmax)
    du <- diag(Du)
    dr <- diag(Dr)
    dr <- dr * (dr>0)

    if (tau>0){
      d <- dr
      D <- diag(dr)
    }else{
      d <- du
      D <- Du
    }

    explained <- cumsum(d^2)
    explained[kmax+1] <- 0.0
    Sigma <- total * rep(1, kmax+1) - explained
    IC2 <- log(Sigma) + CT
    ic2 <- which.min(IC2)
    PC2k <- Sigma + CT*Sigma[kmax]
    PC20 <- Sigma + CT*Sigma[kmax+1]
    pc2k <- which.min(PC2k)
    pc20 <- which.min(PC20)

    out$ic2 <- ic2 * (ic2 <= kmax)
    out$pc2k <- pc2k * (pc2k <= kmax)
    out$pc20 <- pc20 * (pc20 <= kmax)
    out$Fhat <- U[,1:kmax] %*% D[1:kmax,1:kmax]^(1/2)
    out$Lamhat <- V[,1:kmax] %*% D[1:kmax,1:kmax]^(1/2)
    out$Chat <- out$Fhat %*% t(out$Lamhat)
    out$Sigma <- Sigma
    out$IC2 <- IC2
    out$PC2k <- PC2k
    out$PC20 <- PC20
    out$fhat <- out$Fhat * sqrt(T)
    out$lamhat <- out$Lamhat * sqrt(N)
    out$d <- d
    out$d0 <- d0


  } else {   # baing02 (not rank regularized)
    print("Bai Ng 2002 (not rank regularized)", quote = FALSE)

    IC2 <- rep(0, kmax+1)
    PC20 <- rep(0, kmax+1)
    PC2k <- rep(0, kmax+1)
    Sigma <- rep(0, kmax+1)

    svdxx <- svd(X %*% t(X))
    ev <- svdxx$u
    eigval <- svdxx$d
    ev1 <- svdxx$v

    sumeigval <- cumsum(eigval)/sum(eigval)
    Fhat0 <- sqrt(T) * ev
    Lambda0 <- t(X) %*% Fhat0 / T

    Sigma[kmax+1] <- mean(colSums(X^2/T))
    IC2[kmax+1] <- log(Sigma[kmax+1])
    PC2k[kmax+1] <- Sigma[kmax+1]
    PC20[kmax+1] <- Sigma[kmax+1]

    for (i in kmax:1) {
      Fhat <- Fhat0[,1:i]
      lambda <- Lambda0[,1:i]
      chat <- Fhat %*% t(lambda)
      ehat <- X - chat
      Sigma[i] <- mean(colSums(ehat^2/T))
      IC2[i] <- log(Sigma[i]) + CT[i]
      PC2k[i] <- Sigma[i] + CT[i]*Sigma[kmax]
      PC20[i] <- Sigma[i] + CT[i]*Sigma[kmax+1]
    }

    ic2 <- which.min(IC2)
    pc2k <- which.min(PC2k)
    pc20 <- which.min(PC20)

    out$ic2 <- ic2 * (ic2 <= kmax)
    out$pc2k <- pc2k * (pc2k <= kmax)
    out$pc20 <- pc20 * (pc20 <= kmax)
    out$Fhat <- Fhat0[,1:kmax]
    out$Lamhat <- Lambda0[,1:kmax]
    out$Chat <- out$Fhat %*% t(out$Lamhat)
    out$d <- eigval[1:kmax]
    out$IC2 <- IC2
    out$PC2k <- PC2k
    out$PC20 <- PC20
    out$Sigma <- Sigma
  }

  class(out) <- c("rpca", "list")

  print("Available output:", quote = FALSE)
  print(names(out))
  return(out)
}
