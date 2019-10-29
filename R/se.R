#' @title Standard Error of C^hat
#'
#' @description
#' \code{se.rpca} produces the estimated standard error of C^hat
#' produced by the \code{\link{rpca}} function.
#'
#' @export
#'
#' @param object an object of class \code{\link{rpca}}.
#' @param xpoints placeholder.
#' @param qq placeholder.
#'
#' @return standard error of C^hat
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
#' results <- rpca(X, kmax, standardize = FALSE, rank_reg = FALSE, tau = 0)
#' coef(results)


se.rpca <- function(object, xpoints, qq){
  if (!inherits(object, "rpca"))
    stop("Object must be of class 'rpca'")

  Chat_se <- NA

  for (ipoints in 1:length(xpoints)) {
    #tt <- xpoints[ipoints,1]
    #ii <- xpoints[ipoints,2]
    #TT <- xpoints[ipoints,3]
    #NN <- xpoints[ipoints,4]

    temp_xpoints <- xpoints[[ipoints]]
    tt <- temp_xpoints[1]
    ii <- temp_xpoints[2]
    TT <- temp_xpoints[3]
    NN <- temp_xpoints[4]

    Lamhat <- object$Lamhat
    r <- ncol(Lamhat)
    Fhat <- object$Fhat
    Lamhat_ii <- matrix(Lamhat[ii,], nrow = 1)
    Fhat_tt <- matrix(Fhat[tt,], nrow = 1)
    Chat <- object$Chat
    X <- object$X  # clarify
    ehat <- X - Chat
    LamEhat_t <- matrix(rep(ehat[tt,], r), ncol = r,
                        nrow = ncol(ehat), byrow = FALSE) * Lamhat
    FEhat_i <- matrix(rep(ehat[,ii], r), ncol = r,
                      nrow = nrow(ehat), byrow = FALSE) * Fhat
    bread_t <- t(LamEhat_t) %*% LamEhat_t / NN
    bread_i <- t(FEhat_i) %*% FEhat_i / TT

    for (k in 1:qq) {
      bread_i <- bread_i + t(FEhat_i[(k+1):nrow(FEhat_i),]) %*%
        FEhat_i[1:(nrow(FEhat_i)-k),] / TT
    }

    Sigma_Lamhat <- t(Lamhat) %*% Lamhat / nrow(Lamhat)
    V0 <- Lamhat_ii %*% solve(Sigma_Lamhat) %*% bread_t %*%
      solve(Sigma_Lamhat) %*% t(Lamhat_ii)
    W0 <- Fhat_tt %*% bread_i %*% t(Fhat_tt)
    temp <- V0/NN + t(W0)/TT

    Chat_se <- cbind(Chat_se, sqrt(temp))

  }

  Chat_se <- Chat_se[,2:ncol(Chat_se)]  # Remove the first column of NAs
  return(Chat_se)

}
