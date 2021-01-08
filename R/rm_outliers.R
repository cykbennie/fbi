#' @title Remove outliers of the FRED-MD Data Set
#'
#' @description
#' \code{rm_outliers.fredmd} removes outliers of the FRED-MD data set produced by
#' the \code{\link{fredmd}} function.
#'
#' @import stats
#' @export
#'
#' @param object an object of class \code{\link{fredmd}}.
#'
#' @return FRED-MD data of class \code{fredmd} with outliers removed.
#'
#' @author Yankang (Bennie) Chen <yankang.chen@@yale.edu>
#'
#' @references
#' Michael W. McCracken and Serena Ng (2015), \emph{FRED-MD and FRED-QD: Monthly and Quarterly Databases for Macroeconomic Research}.
#' \url{https://research.stlouisfed.org/econ/mccracken/fred-databases/}
#'
#' @examples
#' library(fbi)
#' data <- fredmd(date_start = NULL, date_end = NULL, transform = TRUE)
#' newdata <- rm_outliers.fredmd(data)



rm_outliers.fredmd <- function(object) {
  # Error checking
  if (!inherits(object, "fredmd"))
    stop("Object must be of class 'fredmd'")


  data <- object
  N <- ncol(data)
  X <- data[, 2:N]

  # Calcualte median of each series
  median_X <- apply(X, 2, stats::median, na.rm = TRUE)

  # Repeat median of each series over all data points in the series
  median_X_mat <- matrix(rep(median_X, nrow(X)), nrow = nrow(X),
                         ncol = ncol(X), byrow = TRUE)

  # Calculate quartiles
  Q <- apply(X, 2, stats::quantile, probs = c(0.25, 0.75), na.rm = TRUE)

  # Calculate interquartile range (IQR) of each series
  IQR <- Q[2, ] - Q[1, ]

  # Repeat IQR of each series over all data points in the series
  IQR_mat <- matrix(rep(IQR, nrow(X)), nrow = nrow(X),
                    ncol = ncol(X), byrow = TRUE)

  # Determine outliers
  Z <- abs(X - median_X_mat)
  outlier <- (Z > (10 * IQR_mat))

  # Replace outliers with NaN
  Y <- X
  Y[outlier] <- NA

  # Cleaned data
  outdata <- data
  outdata[, 2:N] <- Y
  class(outdata) <- c("data.frame", "fredmd")
  return(outdata)

  # Print the number of outliers
  print("Number of outliers:", quote = FALSE)
  print(sum(outlier, na.rm = TRUE), quote = FALSE)

}
