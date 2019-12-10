#' @title Demean Panel Data
#'
#' @description
#' \code{demeanXY} demeans the panel data.
#'
#' @import pracma
#'
#' @param X detaframe or matrix of the original panel data.
#' @param N integer, total number of columns of the panel data.
#' @param T integer, total number of rows of the panel data.
#' @param N0 integer, the number of columns in the panel data with full data availability.
#' @param T0 integer, the number of rows in the panel data with full data availability.
#'
#' @return a list of elements:
#' \item{X1}{demeaned data}
#' \item{FE}{estimated fixed effects matrix}
#'
#' @author Yankang (Bennie) Chen <yankang.chen@@columbia.edu>
#' @author Serena Ng <serena.ng@@columbia.edu>
#' @author Jushan Bai <jushan.bai@@columbia.edu>



demeanXY <- function(X, N, T, N0, T0){
	aihat <- pracma::zeros(1, N)
	xihat <- pracma::zeros(T, 1)

	aihat[1:N0] <- colMeans(X[, 1:N0])
	aihat[(N0+1):N] <- colMeans(X[1:T0, (N0+1):N])

	xihat[1:T0, 1] <- rowMeans(X[1:T0, ])
	xihat[(T0+1):T, 1] <- rowMeans(X[(T0+1):T, 1:N0])

	cc <- mean(X[1:T0, 1:N0])
	FE <- list()
	FE$i <- pracma::repmat(aihat, T, 1);
	FE$t <- pracma::repmat(xihat, 1, N);
	X1 <- X - FE$i - FE$t

	output <- list()
	output$X1 <- X1
	output$FE <- FE
	return(output)
}
