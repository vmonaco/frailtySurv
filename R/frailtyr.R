.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to frailtyr")
}

.onLoad <- function(libname, pkgname) {
  invisible()
}

.onUnload <- function (libpath) {
  library.dynam.unload("frailtyr", libpath)
}