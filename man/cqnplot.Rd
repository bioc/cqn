\name{cqnplot}
\alias{cqnplot}
\title{
  Plot the systematic effect estimated as part of a CQN normalization.~
}
\description{
  This function plots the estimated systematic effect which are removed
  suring CQN normalization.
}
\usage{
cqnplot(x, n = 1, col = "grey60", ylab = "QR fit", xlab = "", type = "l", lty = 1, ...)
}
\arguments{
  \item{x}{ The result of a call to \code{cqn}; an object of class \code{cqn}.}
  \item{n}{ Which systematic effect is plotted.}
  \item{col}{ A vector of colors, as in \code{plot}.}
  \item{ylab}{ y-label as in \code{plot}.}
  \item{xlab}{ x-label as in \code{plot}.}
  \item{type}{ type, as in \code{plot}.}
  \item{lty}{ line type, as in \code{plot}.}
  \item{\dots}{ These arguments are passed to \code{matplot}}
}
\value{
  This function is invoked for its side effect.
}
\author{
  Kasper Daniel Hansen
}
\examples{
data(montgomery.subset)
data(sizeFactors.subset)
data(uCovar)
cqn.subset <- cqn(montgomery.subset, lengths = uCovar$length, 
                  x = uCovar$gccontent, sizeFactors = sizeFactors.subset,
                  verbose = TRUE)
cqnplot(cqn.subset, n = 1)
}
\keyword{hplot}
