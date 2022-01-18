# fbi: Factor-Based Imputation for Missing Data
[![Build Status](https://travis-ci.com/cykbennie/fbi.svg?branch=master)](https://travis-ci.com/cykbennie/fbi)

The codes in the `fbi` package are developed based on the following papers:

* [Bai, Jushan and Ng, Serena (2021), "Matrix Completion, Counterfactuals, and Factor Analysis of Missing Data"](https://doi.org/10.1080/01621459.2021.1967163).

* [Cahan, E., Bai, J. and Ng, S. (2021), "Factor-Based Imputation of Missing Values and Covariances in Panel Data of Large Dimensions"](https://arxiv.org/abs/2103.03045).

* [Bai, Jushan and Ng, Serena (2002), "Determining the number of factors in approximate factor models"](https://onlinelibrary.wiley.com/doi/pdf/10.1111/1468-0262.00273).

* [Bai, Jushan and Ng, Serena (2019), "Rank regularized estimation of approximate factor models"](https://www.sciencedirect.com/science/article/pii/S0304407619300764).

* [McCracken, Michael W. and Ng, Serena (2015), "FRED-MD and FRED-QD: Monthly and Quarterly Databases for Macroeconomic Research"](https://research.stlouisfed.org/econ/mccracken/fred-databases/).

## Requirements

The `fbi` package requires the following three R packages:

* `stats`
* `readr`
* `pracma`

They should be installed prior to the installation of the `fbi` package:

``` r
install.packages("stats")
install.packages("readr")
install.packages("pracma")
```

## Installation

The `fbi` package can be installed from github:

``` r
devtools::install_github("cykbennie/fbi")
```

## Reference Manual

Refer to the [fbi.pdf](fbi.pdf) file for details.

## Vignette

After installing the package, use the code:
``` r
browseVignettes("fbi")
```
or

``` r
vignette("factor_fred", package = "fbi")
```
to see an example using the FRED-MD dataset (https://research.stlouisfed.org/econ/mccracken/fred-databases/).

## Author

* Yankang (Bennie) Chen (<yankang.chen@yale.edu>)

## License

This project is licensed under the GNU General Public License v3.0 License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

* Serena Ng (<serena.ng@columbia.edu>)
* Jushan Bai (<jushan.bai@columbia.edu>)
