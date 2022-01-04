## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(fbi)

filepath <- "https://files.stlouisfed.org/files/htdocs/fred-md/monthly/2021-08.csv"
data <- fredmd(filepath ,date_start = NULL, date_end = NULL, transform = TRUE)
N <- ncol(data)

# View the head lines of data
head(data)

## -----------------------------------------------------------------------------
data_clean <- rm_outliers.fredmd(data)

## -----------------------------------------------------------------------------
col_na_prop <- apply(is.na(data_clean), 2, mean)
data_select <- data_clean[, (col_na_prop < 0.05)]
data_bal <- na.omit(data_select)
X_bal <- data_bal[,2:ncol(data_bal)]
rownames(X_bal) <- data_bal[,1]

# View balanced data
head(X_bal)

# Run rpca
out <- rpca(X_bal, kmax = 8, standardize = FALSE, tau = 0)

## -----------------------------------------------------------------------------
head(out$Fhat)
head(out$Lamhat)
head(out$Chat)

## -----------------------------------------------------------------------------
xpt <- list(c(1,2,3,4), c(2,3,4,5), c(3,4,5,6), c(4,5,6,7))

se.rpca(object = out, xpoints = xpt, qq = 50)

