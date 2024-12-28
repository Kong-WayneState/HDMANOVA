
# HDMANOVA

<!-- badges: start -->
<!-- badges: end -->

The goal of HDMANOVA is to test the equality of multi-sample high dimensional mean vectors.


## Installation

You can install the development version of HDMANOVA like so:

``` r
remotes::install_github("Kong-WayneState/HDMANOVA")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(HDMANOVA)
## basic example code

set.seed(123)
X.list <- rdatalist(g = 3)
str(X.list)
#> List of 3
#>  $ : num [1:50, 1:100] -0.5605 -0.7104 2.1988 -0.7152 -0.0736 ...
#>  $ : num [1:50, 1:100] -0.494 0.254 1.068 -0.108 -0.604 ...
#>  $ : num [1:50, 1:100] 2.371 -0.12 1.008 0.508 0.953 ...

SKtest(X.list)
#> $statistic
#> [1] 0.07726602
#> 
#> $p.value
#> [1] 0.469206

CXtest(X.list)
#> $statistic
#> [1] 15.10297
#> 
#> $p.value
#> [1] 0.4785073

CXtest(X.list, A = "Omega")
#> $statistic
#> [1] 19.25311
#> 
#> $p.value
#> [1] 0.1506127

HBtest(X.list)
#> $Statistic
#> [1] -0.001000385
#> 
#> $p.value
#> [1] 0.5003991

KVtest(X.list)
#> $Statistic
#> [1] 0.1878805
#> 
#> $p.value
#> [1] 0.4254852
#> 
#> $window
#> [1] 10

KVtest(X.list, method = "MPT", L = NA)
#> $Statistic
#> [1] 0.2170077
#> 
#> $p.value
#> [1] 0.4141012
#> 
#> $window
#> [1] 1

```

## License
This package is licensed under MIT license.
