#' @title FRED-MD Data Set Description
#'
#' @description A description of the FRED-MD data set.
#' @docType data
#' @keywords datasets
#' @name fredmd_description
#' @usage data(fredmd_description)
#'
#' @format A data frame with 135 rows and 9 variables. The variables are as follows:
#' \describe{
#'   \item{id}{series ID number}
#'   \item{tcode}{code of transformation}
#'   \item{ttype}{type of transformation}
#'   \item{fred}{variable name used in the FRED-MD data set}
#'   \item{description}{description of the series}
#'   \item{gsi}{variable name used in the Global Insights Basic Economics Database (GSI)}
#'   \item{gsi:description}{description of the series in GSI}
#'   \item{group}{group of the series}
#'   \item{edited}{logical, indicating if the data has been editted}
#' }
#'
#' @source The \code{fredmd_description} data were obtained from
#' Michael W. McCracken and Serena Ng (2015), \emph{FRED-MD Updated Appendix}.
#' \url{https://s3.amazonaws.com/files.fred.stlouisfed.org/fred-md/Appendix_Tables_Update.pdf}
NULL
