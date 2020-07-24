#' @title Describe selected variables in the FRED-QD Data Set
#'
#' @description
#' \code{describe_qd} provides a description of the selected variables
#' in the FRED-QD data set.
#'
#' @export
#'
#' @param varname string or a vector strings of the format \code{"X1"} to
#' \code{"X135"}.
#' @param name.only logical. If \code{TRUE}, return a dataframe with variable
#' names and types of transformation only;
#' if \code{FALSE}, return a dataframe with more details.
#' @param verbose logical, indicating whether or not descriptions should be printed.
#'
#' @return a vector of variable names, or a data frame with detailed descriptions.
#'
#' @author Yankang (Bennie) Chen <yankang.chen@@columbia.edu>
#'
#' @references
#' Michael W. McCracken and Serena Ng (2020), \emph{FRED-QD Updated Appendix}.
#' \url{https://s3.amazonaws.com/files.fred.stlouisfed.org/fred-md/FRED-QD_appendix.pdf}
#'
#' @examples
#' library(fbi)
#' varnames <- describe_qd(c("X32", "X56"), name.only = TRUE, verbose = FALSE)



describe_qd <- function(varname, name.only = TRUE, verbose = FALSE) {
  # Error checking
  for (name in varname) {
    if (substr(name, 1, 1) != "X")
      stop("varname must be a string or a vector strings
           of the format 'X1' to 'X248'")

    tempnum <- as.integer(substr(name, 2, 4))
    if ((tempnum < 1) | (tempnum > 248))
      stop("varname must be a string or a vector strings
           of the format 'X1' to 'X248'")
  }

  fredqd_description <- fredqd_description

  index <- c()
  for (name in varname) {
    temp <- which.max(fredqd_description$varname == name)
    index <- c(index, temp)
  }

  data <- fredqd_description[index, -(ncol(fredqd_description)-1)]

  if (name.only) {
    outdata <- data[, c("description", "ttype")]
  } else {
    outdata <- data
  }

  if (verbose)
    print(t(data))

  return(outdata)

}
