.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to frailtySurv v1.3.3")
}

.onLoad <- function(libname, pkgname) {
  invisible()
}

.onUnload <- function (libpath) {
  library.dynam.unload("frailtySurv", libpath)
}