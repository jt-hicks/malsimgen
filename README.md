
# malsimgen

<!-- badges: start -->
<!-- badges: end -->

This package helps to generate simulated runs of the deterministic version of malariasimulation (https://github.com/mrc-ide/deterministic-malaria-model).

## Installation

You can install the development version of malsimgen (and mamasante, which has required functions and houses the model) from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jt-hicks/malsimgen")
devtools::install_github("jt-hicks/mamasante")

```

## Example

Basic usage:

``` r
library(malsimgen)
sim_data <- malsimgen::data_gen(volatility=1,
                                init_EIR=100,
                                min_EIR=0.01,
                                max_EIR=500,
                                plot_instance=TRUE)
```

