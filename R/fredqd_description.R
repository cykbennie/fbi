#' @title FRED-QD Data Set Description
#'
#' @description A description of the FRED-QD data set.
#' @docType data
#' @keywords datasets
#' @name fredqd_description
#' @usage data(fredqd_description)
#'
#' @format A data frame with 248 rows and 10 variables. The variables are as follows:
#' \describe{
#'   \item{id}{series ID number}
#'   \item{sw_id}{series ID number in SW (2012)}
#'   \item{tcode}{code of transformation}
#'   \item{ttype}{type of transformation}
#'   \item{sw_factors}{logical, indicating whether a series was used in SW (2012) when constructing factors}
#'   \item{fred_mnemonic}{mnemonic in FRED-QD}
#'   \item{sw_mnemonic}{mnemonic used in SW (2012)}
#'   \item{description}{a brief definition of the series}
#'   \item{group}{group of the series}
#'   \item{varname}{"X" + id}
#' }
#'
#' @source The \code{fredqd_description} data were obtained from
#' Michael W. McCracken and Serena Ng (2020), \emph{FRED-QD Updated Appendix}.
#' \url{https://s3.amazonaws.com/files.fred.stlouisfed.org/fred-md/FRED-QD_appendix.pdf}
NULL
