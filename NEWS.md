frailtySurv
--------

# Version 1.3.6 (April 2019)
  * Added a call to set.seed before each example in genfrail

# Version 1.3.5 (August 2018)
  * Added citation to JSS article

# Version 1.3.4 (June 2018)

## Bug fixes
  * Changed default genfrail control params to constants instead derived from .Machine
  * Added tryCatch to simfrail to avoid failing completely and generate a warning

# Version 1.3.3 (June 2018)

## Bug fixes
  * Fixed rifcond error
  * Fixed missing figures in online docs

# Version 1.3.2 (February 2017)

## Bug fixes
  * Renamed "hazard" "cumhaz" in plot/summary functions to avoid confusion with the baseline hazard rate.
  * Fixed typos in fitfrail.control docs.
  * Fixed a bug in lognormal 2nd derivative, resulting in biased estimates using fitmethod="score".
  * Added checks for object and param matching in plot/summary functions

# Version 1.3.1 (December 2016)
  
## New features
  * Added numeric integration control parameters to genfrail (see genfrail.control)

## Bug fixes
  * Removed unnecssary imports

# Version 1.3 (June 2016)
  
## New features
  * Added init.beta and init.theta control parameters for parameter initialization
  * Added censor.time parameter to genfrail for user-defined censorship times
  * Added option to specify Lambda.times when plotting residuals

## Bug fixes
  * Made both loglik and score fit methods use the same convergence control parameters
  * Improved loglik parameter estimation and convergence
  * Fixed the labels of residuals boxplot x-axis
  * Fixed clang warning with static casts
  
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