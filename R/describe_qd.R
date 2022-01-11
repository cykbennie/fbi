#' @title Describe selected variables in the FRED-QD Data Set
#'
#' @description
#' \code{describe_qd} provides a description of the selected variables
#' in the FRED-QD data set.
#'
#' @export
#'
#' @param varname string or a vector strings of the FRED variable name,
#' such as \code{GDPC1}.
#' @param name.only logical. If \code{TRUE}, return a dataframe with variable
#' names and types of transformation only;
#' if \code{FALSE}, return a dataframe with more details.
#'
#' @return a vector of variable names, or a data frame with detailed descriptions.
#'
#' @author Yankang (Bennie) Chen <yankang.chen@@yale.edu>
#'
#' @references
#' Michael W. McCracken and Serena Ng (2020), \emph{FRED-QD Updated Appendix}.
#' \url{https://s3.amazonaws.com/files.fred.stlouisfed.org/fred-md/FRED-QD_appendix.pdf}
#'
#' @examples
#' library(fbi)
#' varnames <- describe_qd(c("GDPC1", "Y033RC1Q027SBEAx"), name.only = TRUE)



describe_qd <- function(varname, name.only = TRUE) {
  fredqd_description <- fredqd_description

  index <- c()
  for (name in varname) {
    temp <- which.max(fredqd_description$fred == name)
    index <- c(index, temp)
  }

  data <- fredqd_description[index,]

  if (name.only) {
    outdata <- data[, c("fred", "description", "tcode", "ttype")]
  } else {
    outdata <- data
  }

  return(outdata)

}
