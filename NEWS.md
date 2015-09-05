frailtyr
--------

# Version 1.0.0 (August 2015)

  **Initial release**
  
# Version 1.1.0 (September 2015)

## New features
  * Control parameters for numerical integration `(int.iter, int.reltol, int.abstol)`
  * Added `summary.fitfrail` for summarizing the survival curve 
  * Added parameter to `fitfrail` to compute the SE estimates

## Bug fixes
  * Renamed `Lambda.time` to `Lambda.times` (not backwards compatible)
  * Improved print.fitfrail
  * Remove missing observations
  * Improved print.summary.genfrail
  * Added sanity checks to frailty distribution params
  * Fixed parallel support on Windows (run in serial with a warning for now)