genfrail.control <- function(censor.reltol = .Machine$double.eps^0.25,
                             censor.subdivisions = 1000L,
                             crowther.reltol = .Machine$double.eps^0.25,
                             crowther.subdivisions = 1000L
                             ) {
  
  if (censor.reltol < 0) 
    stop ("Invalid censor relative tolerance")
  if (censor.subdivisions < 0) 
    stop ("Invalid censor subdivisions")
  if (crowther.reltol < 0) 
    stop("Invalid Crowther relative tolerance")
  if (crowther.subdivisions < 0) 
    stop ("Invalid Crowther integration subdivisions")
  
  list(censor.reltol=censor.reltol,
       censor.subdivisions=censor.subdivisions,
       crowther.reltol=crowther.reltol,
       crowther.subdivisions=crowther.subdivisions
       )
}
