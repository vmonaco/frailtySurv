frailtySurv
--------

# Version 1.3 (June 2016)

## Bug fixes
  * Made both loglik and score fit methods use the same convergence control parameters
  * Improved loglik parameter estimation and convergence
  * Fixed the labels of residuals boxplot x-axis
  * Added check for Hmisc package in plotting
  * Fixed clang warning with static casts
  
## New features
  * Added init.beta and init.theta control parameters for parameter initialization
  * Added censor.time parameter to genfrail for user-defined censorship times
  * Added option to specify Lambda.times when plotting residuals
  
# Version 1.2.2 (September 2015)

## Bug fixes
  * Fixed ambiguous C++ calls to pow

# Version 1.2.1 (September 2015)

## Bug fixes
  * Resolved undefined global vars
  * Added `Lambda.all` to the fitted model object for estimated cumulative baseline hazard at all observed times
  * Minor changes to man pages

# Version 1.2.0 (September 2015)

## New features
  * Renamed package to `frailtySurv`
  * Added uniform censorship to `genfrail`

# Version 1.1.1 (September 2015)

## Bug fixes
  * Fixed bug in `summary.fitfrail`
  * Fixed `censor.param` parameter in `genfrail`
  
# Version 1.1.0 (September 2015)

## New features
  * Control parameters for numerical integration `(int.iter, int.reltol, int.abstol)`
  * Added `summary.fitfrail` for summarizing the survival curve 
  * Added parameter to `fitfrail` to compute the SE of parameter estimates

## Bug fixes
  * Renamed `Lambda.time` to `Lambda.times` (not backwards compatible)
  * Improved print.fitfrail
  * Remove missing observations
  * Improved print.summary.genfrail
  * Added sanity checks to frailty distribution params
  * Fixed parallel support on Windows (runs in serial with a warning for now)
  
# Version 1.0.0 (August 2015)

  **Initial release**