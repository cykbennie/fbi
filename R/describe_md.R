#' @title Describe selected variables in the FRED-MD Data Set
#'
#' @description
#' \code{describe_md} provides a description of the selected variables
#' in the FRED-MD data set.
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
#' Michael W. McCracken and Serena Ng (2015), \emph{FRED-MD Updated Appendix}.
#' \url{https://s3.amazonaws.com/files.fred.stlouisfed.org/fred-md/Appendix_Tables_Update.pdf}
#'
#' @examples
#' library(fbi)
#' varnames <- describe_md(c("RPI", "RETAILx"), name.only = TRUE)



describe_md <- function(varname, name.only = TRUE) {
  fredmd_description <- fredmd_description

  index <- c()
  for (name in varname) {
    temp <- which.max(fredmd_description$fred == name)
    index <- c(index, temp)
  }

  data <- fredmd_description[index,]

  if (name.only) {
    outdata <- data[, c("fred", "description", "tcode", "ttype")]
  } else {
    outdata <- data
  }

  return(outdata)

}
