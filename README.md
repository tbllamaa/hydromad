
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hydromad

<!-- badges: start -->

<!-- badges: end -->

The goal of hydromad is to provide a modelling framework for
environmental hydrology: water balance accounting and flow routing in
spatially aggregated catchments.

Hydromad supports simulation, estimation, assessment and visualisation
of flow response to time series of rainfall and other drivers. A minimal
unit hydrograph framework is used, where areal rainfall is passed
through a soil moisture accounting (SMA) model to estimate effective
rainfall; this is then passed through a routing model to estimate
streamflow. Included are several implementations of common hydrological
models consistent with this framework.

## Installation

<!--- You can install the released version of hydromad from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("hydromad")
```
--->

Currently you can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JosephGuillaume/hydromad")
```

## Contributing

The maintainers of Hydromad are always keen to receive contributions.
However to help us we would like you to consider the following:

  - Please read and use our [CONTRIBUTING
    document](https://github.com/josephguillaume/hydromad/blob/master/docs/CONTRIBUTING.md),
    this will really help us integrating your solution into Hydromad.  
  - This is a list of the \[current issues\] (which you might be able to
    help with).  
  - As part of contributing we would like you to consider the ropensci
    [Code of Conduct](https://ropensci.org/code-of-conduct/), which we
    use for Hydromad.

## Usage

This is a basic example which shows you how to solve a common problem:

``` r
library(hydromad)
## basic example code
data(Cotter)
## IHACRES CWI model with exponential unit hydrograph
## an unfitted model, with ranges of possible parameter values
modx <- hydromad(Cotter[1:1000], sma = "cwi", routing = "expuh",
                 tau_s = c(2, 100), v_s = c(0, 1))
modx
#> 
#> Hydromad model with "cwi" SMA and "expuh" routing:
#> Start = 1966-05-01, End = 1969-01-24
#> 
#> SMA Parameters:
#>       lower upper     
#> tw        0   100     
#> f         0     8     
#> scale    NA    NA     
#> l         0     0 (==)
#> p         1     1 (==)
#> t_ref    20    20 (==)
#> Routing Parameters:
#>       lower upper  
#> tau_s     2   100  
#> v_s       0     1
```

## License

This project is licensed under the terms of the [MIT
license](https://github.com/josephguillaume/hydromad/blob/master/LICENSE.txt)
