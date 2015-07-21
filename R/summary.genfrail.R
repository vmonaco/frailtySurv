summary.genfrail <- function(gen, ...) {
  
  s <- append(attributes(gen), list(
      n.obs=length(gen$time),
      n.clusters=length(unique(gen$cluster)),
      mean.cluster=mean(table(gen$cluster)),
      censor.rate=1 - sum(gen$status)/length(gen$time)
  ))
  
  class(s) <- "summary.genfrail"
  
  s
}