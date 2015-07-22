summary.genfrail <- function(dat, ...) {
  
  s <- append(attributes(dat), list(
      n.obs=length(dat$time),
      n.clusters=length(unique(dat$cluster)),
      mean.cluster=mean(table(dat$cluster)),
      censor.rate=1 - sum(dat$status)/length(dat$time)
  ))
  
  class(s) <- "summary.genfrail"
  
  s
}