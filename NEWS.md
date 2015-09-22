frailtyr
--------

# Version 1.2.0 (September 2015)

## Bug fixes
  * Renamed package to `frailtySurv`

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