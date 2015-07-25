\name{genfrail}
\alias{genfrail}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate survival data with a hidden shared frailty among clusters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
genfrail(beta = c(log(2)), frailty = c("gamma", "lognormal", "posstab"), censor.mu = 130, censor.sigma = 15, theta = 2, N = 300, K = 2, tau = 4.6, C = 0.01, covariates = c("normal", "uniform"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{frailty}{
%%     ~~Describe \code{frailty} here~~
}
  \item{censor.mu}{
%%     ~~Describe \code{censor.mu} here~~
}
  \item{censor.sigma}{
%%     ~~Describe \code{censor.sigma} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
}
  \item{C}{
%%     ~~Describe \code{C} here~~
}
  \item{covariates}{
%%     ~~Describe \code{covariates} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ kwd1 }
\keyword{ kwd2 }% __ONLY ONE__ keyword per line
