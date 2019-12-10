#' @title Estimate Treatment Effect
#'
#' @description
#' \code{tnt} estimates the treatment effect.
#'
#' @import pracma
#' @export
#'
#' @param data list containing x1, x2, y0, y1, N0, N1, T0, and T1.
#' @param param list containing K, r, do_FE, do_IFE, and maxit1.
#'
#' @return a list of elements:
#' \item{est}{}
#' \item{SE}{}
#' \item{V}{}
#' \item{it1}{}
#'
#' @author Yankang (Bennie) Chen <yankang.chen@@columbia.edu>
#' @author Serena Ng <serena.ng@@columbia.edu>
#' @author Jushan Bai <jushan.bai@@columbia.edu>
#'
#' @references
#' Jushan Bai and Serena Ng (2019), \emph{Matrix Completion, Counterfactuals, and Factor Analysis of Missing Data}.
#' \url{https://arxiv.org/abs/1910.06677}



tnt <- function(data, param){
  # Error checking
    if (! all(c("x1", "x2", "y0", "y1", "N0", "N1", "T0", "T1") %in% names(data)))
      stop("'data' must be a list containing x1, x2, y0, y1, N0, N1, T0, and T1.")
    if (! all(c("K", "r", "do_FE", "do_IFE", "maxit1") %in% names(param)))
      stop("'param' must be a list containing K, r, do_FE, do_IFE, and maxit1.")



  # Clear memory and create output object
  gc()
  output <- list()

  x1 <- data$x1
  x2 <- data$x2
  y0 <- data$y0
  y1 <- data$y1
  N0 <- data$N0
  N1 <- data$N1
  T0 <- data$T0
  T1 <- data$T1
  K <- param$K
  N <- N0 + N1
  T <- T0 + T1
  r <- param$r

  missing <- is.na(y0)

  if (param$do_FE == 0) {
  	X1 <- x1
    X2 <- x2
    Y0 <- y0
    Y1 <- y1
  }

  if (param$do_FE == 1) {
    X1 <- removeFE(x1,N,T,N0,T0)$X1
    X2 <- removeFE(x2,N,T,N0,T0)$X1
    Y0 <- removeFE(y0,N,T,N0,T0)$X1
    FE <- removeFE(y0,N,T,N0,T0)$FE
    Y1 <- y1 - FE$i
  }

  if (param$do_FE == 2) {
    X1 <- demeanXY(x1,N,T,N0,T0)$X1
    X2 <- demeanXY(x2,N,T,N0,T0)$X1
    Y0 <- demeanXY(y0,N,T,N0,T0)$X1
    FE <- demeanXY(y0,N,T,N0,T0)$FE
    Y1 <- y1 - FE$i - FE$t
  }


  X1vecN0 <- X1[, 1:N0]
  X2vecN0 <- X2[, 1:N0]
  Y0vecN0 <- Y0[, 1:N0]
  X1vecN <- X1
  X2vecN <- X2
  oneN0 <- pracma::ones(nrow(X1vecN0), 1)
  reg0 <- cbind(X1vecN0, X2vecN0, oneN0)

  out1 <- list()
  out2 <- list()
  out1$beta <- pracma::zeros(ncol(reg0), 1)

  if (param$do_IFE == 1) {
    out1$beta <- solve(t(reg0) %*% reg0) %*% t(reg0) %*% Y0vecN0
    it1 <- 1
    done <- 0

    while ((done == 0) & (it1 < param$maxit1)) {
    	Zhat <- Y0[,1:N0] - out1$beta[1] * X1[,1:N0] - out1$beta[2] * X2[,1:N0] - out1$beta[3]
      itFE <- fbi::apc(Zhat, r)
      Y00 <- Y0[,1:N0] - itFE$Fhat %*% t(itFE$Lamhat)
      Y00vec <- Y00
      out2$beta <- solve(t(reg0) %*% reg0) %*% t(reg0) %*% Y00vec
      gap <- out1$beta - out2$beta

      if (sqrt(sum(gap^2)) < 1e-5) {
				done <- 1
				break
      } else {
       	out1$beta <- out2$beta
        it1 <- it1+1
    	}
    } # end iteractive fixed effect

  } else {
  	it1 <- 0
  }


  Y0tilde <- Y0 - out1$beta[1] * X1 - out1$beta[2] * X2 - out1$beta[3]
  Y1hat <- Y1 - out1$beta[1] * X1 - out1$beta[2] * X2 - out1$beta[3]

  Alw3 <- fbi::tw_apc(Y0tilde, r, center = FALSE, standardize = FALSE)
  Fhat <- Alw3$Fhat
  Lhat <- Alw3$Lamhat

  Y0hat <- Fhat %*% t(Lhat)
  ehat <- Y0tilde - Y0hat

  Lam <- list()
  Lam$Sigma <- (t(Lhat) %*% Lhat) / N

  if (N1 > 1) {
    Lam$treated <- matrix(colSums(Lhat[(N0+1):N, ]),
                          ncol = 1) / N1
  } else {
    Lam$treated <- matrix(Lhat[N, ],
                          ncol = 1)
  }

  denom <- T * N0 - r * (T+N0) + r^2 - K
  etilde <- ehat[, 1:N0]
  sig2e <- sum(etilde^2) / denom

  V <- list()
  V$Att_t <- pracma::zeros(T, 1)
  V$Att_it <- pracma::zeros(T, 1) # individual effect

  for (t in ((T0+1):T)) {
  	ee <- matrix(ehat[t, 1:N0], nrow = 1)
  	Lehat <- Lhat[1:N0, ] * pracma::repmat(t(ee), 1, r)
  	Gamt <- t(Lehat) %*% Lehat / N0
  	bread <- solve(Lam$Sigma) %*% Gamt %*% solve(Lam$Sigma)
  	temp1 <- t(Lam$treated) %*% bread %*% Lam$treated
  	V$Att_t[t, 1] <- temp1/N0 + sig2e/N1
  	temp2 <- matrix(Lhat[1, ], nrow = 1) %*% bread %*% matrix(Lhat[1, ], ncol = 1)
  	V$Att_it[t, 1] <- temp2 + sig2e/N1
  }

  est <- list()
  est$Att_it <- Y1hat[(T0+1):T, (N0+1):N] - Y0hat[(T0+1):T, (N0+1):N]

  if (N1 > 1) {
  	est$Att_t <- colMeans(t(est$Att_it))
  	est$Att <- mean(est$Att_it)
  } else {
	est$Att_t <- est$Att_it
	est$Att <- mean(est$Att_t)
  }

  V$Att_t <- matrix(V$Att_t[(T0+1):T, 1], nrow = 1)
  V$Att_it <- matrix(V$Att_it[(T0+1):T, 1], ncol = 1)

  SE <- list()
  SE$Att_it <- sqrt(V$Att_it)
  SE$Att_t <- sqrt(V$Att_t)

  V$Att <- mean(V$Att_t)
  SE$Att <- mean(SE$Att_t)
  est$beta <- out1$beta


  output$est <- est
  output$SE <- SE
  output$V <- V
  output$it1 <- it1

  return(output)
}
