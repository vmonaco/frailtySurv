#' Run functions in parallel with reproducible results
#' 
#' 
plapply <- function(reps, FUN, cores=0) {
  # To make the result reproducable when parallelized, seed each run 
  seeds <- sample(1:1e7, reps, replace=FALSE)
  
  use.parallel <- cores != 1
  have.parallel <- requireNamespace("parallel", quietly = TRUE)
  if (use.parallel && have.parallel) {
    # Run in parallel
    if (cores <= 0) {
      cores <- parallel::detectCores() + cores
    }
    result <- parallel::mclapply(seeds, FUN, mc.preschedule=FALSE, 
                               mc.set.seed=TRUE, mc.cores=cores, mc.silent=FALSE)
  } else if (use.parallel && !have.parallel) {
    # Run in serial with a warning
    warning('Cannot run plapply in parallel: requires the "parallel" package.')
    result <- lapply(seeds, FUN)
  } else {
    # Run in serial
    result <- lapply(seeds, FUN)
  }
  
  result
}