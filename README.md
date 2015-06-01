<!-- README.md is generated from README.Rmd. Please edit that file -->
frailtyr
--------

frailtyr is a library for fitting shared frailty models with support for general frailty distributions. **This is currently a work in progress**, part of Google Summer of Code 2015.

Installation
------------

To get started, you can install the package from github using `devtools`.

``` r
library(devtools)
install_github("vmonaco/frailtyr")
```

Clone and build
---------------

Instead of installing from github, you can clone and build the project from source. To create a project in RStudio from this repository, you must have both RStudio and git installed.

-   In RStudio, go to File -\> New Project -\> Version Control -\> Git
-   Name the project (eg. frailtyr), choose location, and specify the Repository URL as `https://github.com/vmonaco/frailtyr`

More information about [RStudio and git](https://support.rstudio.com/hc/en-us/articles/200532077-Version-Control-with-Git-and-SVN)

Also see the chapter in the [R packages book](http://r-pkgs.had.co.nz/git.html)

### Build everything

To load the package (simulate installing and loading with `library(frailty)`):

``` r
devtools::load_all()
```

or shortuct: **Ctrl+Shift+L**

Use this command to build a source package (binary = FALSE by default)

``` r
devtools::build()
```

Other usefull commands and shortcuts:

-   Build and reload everything: **Ctrl+Shift+B**
-   Check the package, `devtools::check()`: **Ctrl+Shift+E**
-   Run all tests, `devtools::test()`: **Ctrl+Shift+T**

### Build vignettes

Build the vignettes with the command:

``` r
devtools::build_vignettes()
```

Compiled vignettes will reside in `inst/doc`. In particular, see frailtyr.pdf

### Run tests

frailtyr uses the `testthat` library for tests. Test files are all located in `tests/testthat`. To run all tests, use the command `devtools::test()` or shortcut **Ctrl+Shift+T**

### git pull

Use RStudio to pull the latest code if you've already cloned the repository. In the "git" pane, look for the Pull button. You can also view the commit history from there and see the diff between modified files, if any.

Example
-------

The following code shows how to generate data and fit the model. For now, the baseline hazard estimator is very slow. Estimating parameters for 60 data points takes about 90 seconds, as shown in the example.

``` r
set.seed(1234)
library(frailtyr)
#> Loading required package: survival
#> Loading required package: nleqslv
#> Welcome to frailtyr
data <- genfrail(beta=c(log(2)), theta=2, N=30, K=2, frailty="gamma", covariates="normal")
ptm <- proc.time()
fitfrail(Surv(time, status) ~ Z1 + cluster(family), data, frailty="gamma")
#> Iteration  1  :  beta <-  0.3715843 , theta <-  1 
#> Iteration  2  :  beta <-  0.3715843 , theta <-  1 
#> Iteration  3  :  beta <-  0.3715843 , theta <-  1 
#> Iteration  4  :  beta <-  0.3715843 , theta <-  1 
#> Iteration  5  :  beta <-  0.4351054 , theta <-  0.413516 
#> Iteration  6  :  beta <-  0.4310277 , theta <-  0.7082956 
#> Iteration  7  :  beta <-  0.4296597 , theta <-  0.6397397 
#> Iteration  8  :  beta <-  0.4294886 , theta <-  0.6182223 
#> Iteration  9  :  beta <-  0.4295137 , theta <-  0.6201941 
#> Iteration  10  :  beta <-  0.4295115 , theta <-  0.6201543 
#> Iteration  11  :  beta <-  0.4295117 , theta <-  0.620154 
#> Iteration  12  :  beta <-  0.4295116 , theta <-  0.620154
#> $beta
#>        Z1 
#> 0.4295116 
#> 
#> $theta
#>          
#> 0.620154 
#> 
#> $method
#> [1] "fitfrail"
#> 
#> $call
#> fitfrail(formula = Surv(time, status) ~ Z1 + cluster(family), 
#>     data = data, frailty = "gamma")
#> 
#> attr(,"class")
#> [1] "fitfrail"
proc.time() - ptm
#>    user  system elapsed 
#>  91.591   0.436  92.335
```

These results can be compared to existing shared frailty models.

``` r
library(survival)
library(frailtypack)
```

``` r
coxph(Surv(time, status) ~ Z1 + frailty.gamma(family), data=data)
#> Call:
#> coxph(formula = Surv(time, status) ~ Z1 + frailty.gamma(family), 
#>     data = data)
#> 
#>                       coef se(coef) se2   Chisq DF   p    
#> Z1                    0.44 0.202    0.178  4.72  1.0 0.030
#> frailty.gamma(family)                     21.48 13.4 0.075
#> 
#> Iterations: 7 outer, 43 Newton-Raphson
#>      Variance of random effect= 0.658   I-likelihood = -122.5 
#> Degrees of freedom for terms=  0.8 13.4 
#> Likelihood ratio test=38.3  on 14.2 df, p=0.000536  n= 60
frailtyPenal(Surv(time, status) ~ Z1 + cluster(family), data=data, n.knots=10, kappa=2)
#> 
#> Be patient. The program is computing ... 
#> The program took 0.06 seconds
#> Call:
#> frailtyPenal(formula = Surv(time, status) ~ Z1 + cluster(family), 
#>     data = data, n.knots = 10, kappa = 2)
#> 
#> 
#>   Shared Gamma Frailty model parameter estimates  
#>   using a Penalized Likelihood on the hazard function 
#> 
#>        coef exp(coef) SE coef (H) SE coef (HIH)       z       p
#> Z1 0.464124   1.59062    0.213063      0.213063 2.17834 0.02938
#> 
#>     Frailty parameter, Theta: 0.798327 (SE (H): 0.514862 ) p = 0.060503 
#>  
#>       penalized marginal log-likelihood = -180.75
#>       Convergence criteria: 
#>       parameters = 4.68e-05 likelihood = 0.000103 gradient = 2.31e-09 
#> 
#>       LCV = the approximate likelihood cross-validation criterion
#>             in the semi parametrical case     = 3.23753 
#> 
#>       n= 60
#>       n events= 36  n groups= 30
#>       number of iterations:  11 
#> 
#>       Exact number of knots used:  10 
#>       Value of the smoothing parameter:  2, DoF:  8.97
```
