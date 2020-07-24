#' @title Loading FRED-MD Data Set
#'
#' @description \code{fredmd} loads the official FRED-MD data set and provides
#' a few tools to manipulate the data set.
#'
#' @import readr
#' @export
#'
#' @param date_start Date or \code{NULL}, the start date (included) of the data selection.
#' If \code{NULL}, select till the latest data available.
#' @param date_end Date or \code{NULL}, the end date (included) of the data selection.
#' If \code{NULL}, select up to the earliest data available.
#' @param transform logical, indicating Whether or not the FRED-MD data set
#' should be transformed according to the transformation code.
#' @param local logical, indicating Whether or not the FRED-MD data set
#' should be loaded from the local files or downloaded online
#' @return a subset of the (transformed) FRED-MD data of class \code{fredmd}.
#'
#' @author Yankang (Bennie) Chen <yankang.chen@@columbia.edu>
#'
#' @references
#' Michael W. McCracken and Serena Ng (2015), \emph{FRED-MD and FRED-QD: Monthly and Quarterly Databases for Macroeconomic Research}.
#' \url{https://research.stlouisfed.org/econ/mccracken/fred-databases/}
#'
#' @examples
#' library(fbi)
#' data <- fredmd(date_start = NULL, date_end = NULL, transform = TRUE, local = FALSE)



fredmd <- function(date_start = NULL, date_end = NULL, transform = TRUE,
                   local = FALSE) {
  # Error checking
  if (!is.logical(transform))
    stop("'transform' must be logical.")
  if ((class(date_start) != "Date") && (!is.null(date_start)))
    stop("'date_start' must be Date or NULL.")
  if ((class(date_end) != "Date") && (!is.null(date_end)))
    stop("'date_end' must be Date or NULL.")

  if (class(date_start) == "Date") {
    if (as.numeric(format(date_start, "%d")) != 1)
      stop("'date_start' must be Date whose day is 1.")
    if (date_start < as.Date("1959-01-01"))
      stop("'date_start' must be later than 1959-01-01.")
  }

  if (class(date_end) == "Date") {
    if (as.numeric(format(date_end, "%d")) != 1)
      stop("'date_end' must be Date whose day is 1.")
  }



  # Prepare raw data
  if (local) {
    # local files
    rawdata <- fredmd_2020_04

  } else {
    # download online
    dataurl <- "https://s3.amazonaws.com/files.fred.stlouisfed.org/fred-md/monthly/current.csv"
    rawdata <- readr::read_csv(url(dataurl), col_names = FALSE, col_types = cols(X1 = col_date(format = "%m/%d/%Y")),
                               skip = 2)
    rawdata <- rawdata[1:(nrow(rawdata) - 1), ] # remove NA rows
    rawdata <- as.data.frame(rawdata)
    header <- c("date", colnames(rawdata))[1:ncol(rawdata)]
    colnames(rawdata) <- header
  }



  # Import tcode tcodes is an internal data of the R package
  tcode <- tcodes_md


  # Subfunction transxf: data transformation based on tcodes
  transxf <- function(x, tcode) {
    # Number of observations (including missing values)
    n <- length(x)

    # Value close to zero
    small <- 1e-06

    # Allocate output variable
    y <- rep(NA, n)
    y1 <- rep(NA, n)

    # TRANSFORMATION: Determine case 1-7 by transformation code
    if (tcode == 1) {
      # Case 1 Level (i.e. no transformation): x(t)
      y <- x

    } else if (tcode == 2) {
      # Case 2 First difference: x(t)-x(t-1)
      y[2:n] <- x[2:n] - x[1:(n - 1)]

    } else if (tcode == 3) {
      # case 3 Second difference: (x(t)-x(t-1))-(x(t-1)-x(t-2))
      y[3:n] <- x[3:n] - 2 * x[2:(n - 1)] + x[1:(n - 2)]

    } else if (tcode == 4) {
      # case 4 Natural log: ln(x)
      if (min(x, na.rm = TRUE) > small)
        y <- log(x)

    } else if (tcode == 5) {
      # case 5 First difference of natural log: ln(x)-ln(x-1)
      if (min(x, na.rm = TRUE) > small) {
        x <- log(x)
        y[2:n] <- x[2:n] - x[1:(n - 1)]
      }

    } else if (tcode == 6) {
      # case 6 Second difference of natural log:
      # (ln(x)-ln(x-1))-(ln(x-1)-ln(x-2))
      if (min(x, na.rm = TRUE) > small) {
        x <- log(x)
        y[3:n] <- x[3:n] - 2 * x[2:(n - 1)] + x[1:(n - 2)]
      }

    } else if (tcode == 7) {
      # case 7 First difference of percent change:
      # (x(t)/x(t-1)-1)-(x(t-1)/x(t-2)-1)
      y1[2:n] <- (x[2:n] - x[1:(n - 1)])/x[1:(n - 1)]
      y[3:n] <- y1[3:n] - y1[2:(n - 1)]
    }

    return(y)
  }


  # Transform data
  if (transform) {
    # Apply transformations
    N <- ncol(rawdata)
    data <- rawdata
    data[, 2:N] <- NA

    # Perform transformation using subfunction transxf (see below for
    # details)
    for (i in 2:N) {
      temp <- transxf(rawdata[, i], tcode[i - 1])
      data[, i] <- temp
    }

  } else {
    data <- rawdata
  }


  # Null case of date_start and date_end
  if (is.null(date_start))
    date_start <- as.Date("1959-01-01")
  if (is.null(date_end))
    date_end <- data[, 1][nrow(data)]


  # Subset data
  index_start <- which.max(data[, 1] == date_start)
  index_end <- which.max(data[, 1] == date_end)

  outdata <- data[index_start:index_end, ]
  class(outdata) <- c("data.frame", "fredmd")
  return(outdata)

}
