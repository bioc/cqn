\name{cqn}
\alias{cqn}
\alias{cqn.fixedlength}

\title{
CQN (conditional quantile normalization) for RNA-Seq data
}
\description{
This function implements CQN (conditional quantile normalization) for RNA-Seq data.
}
\usage{
cqn(counts, x, lengths, sizeFactors = NULL, subindex = NULL, tau = 0.5, sqn = TRUE, verbose = FALSE)
cqn.fixedlength(counts, x, lengths, sizeFactors = NULL, subindex = NULL, tau = 0.5, sqn = TRUE, verbose = FALSE)
}
\arguments{
  \item{counts}{
    An object that can be coerced to a \code{matrix} of region by sample
    counts.  Ought to have integer values.
}
  \item{x}{
    This is a covariate whose systematic influence on the counts will be
    removed.  Typically the GC content.  Has to have the same length as
    the number of rows of counts.
}
  \item{lengths}{
    The lengths (in bp) of the regions in counts.  Has to have the same
    length as the number of rows of counts.
}
  \item{sizeFactors}{
    An optional vector of sizeFactors, ie. the sequencing effort of the
    various samples.  If \code{NULL} this is calculated as the column
    sums of \code{counts}.
}
  \item{subindex}{
    An optional vector of indices into the rows of \code{counts}.  If
    not given, this becomes the indices of genes with row means of
    \code{counts} greater then 50.
}
  \item{tau}{
    This argument is passed to \code{rq}, it indicates what quantile is
    being fit.  The default should only be changed by expert users..
}
  \item{sqn}{
    This argument indicates whether the residuals from the systematic
    fit are (subset) quantile normalized.  The default should only be
    changed by expert users.
}
  \item{verbose}{
    Is the function verbose?
}
}

\details{
  These functions implement the CQN (conditional quantile normalization)
  for RNA-Seq data.  The functions remove a single systematic effect,
  contained in the argument \code{x}, which will typicall be GC
  content.  The effect of \code{lengths} will either be modelled as a
  smooth function (which we recommend), if you are using \code{cqn} or
  as an offset (equivalent to modelling using RPKMs), if you are using
  \code{cqn.fixedlength}.

  Final corrected values are equal to \code{value$y + value$offset}.
}

\value{
  A \code{list} with the following components
  \item{counts}{The value of argument \code{counts}.}
  \item{x}{The value of argument \code{x}.}
  \item{lengths}{The value of argument \code{lengths}.}
  \item{sizeFactors}{The value of argument \code{sizeFactors}.  In case
    the argument was \code{NULL}, this is the value used internally.}
  \item{subindex}{The value of argument \code{subindex}.   In case
    the argument was \code{NULL}, this is the value used internally.}
  \item{y}{The dependent value used in the systematic effect fit.  Equal
  to log2 tranformed reads per millions.}
  \item{offset}{The estimated offset.}
  \item{offset0}{A single number used internally for identifiability.}
  \item{func1}{The estimated effect of function 1 (argument \code{x}).  This is a matrix of
    function values on a grid.  Columns are samples and rows are grid points.}
  \item{grid1}{The grid points on which function 1 (argument \code{x}) was evaluated.}
  \item{knots1}{The knots used for function 1 (argument \code{x}).}
  \item{func2}{The estimated effect of function 2 (lengths).  This is a matrix of
    function values on a grid.  Columns are samples and rows are grid points.}
  \item{grid2}{The grid points on which function 2 (lengths) was evaluated.}
  \item{knots2}{The knots used for function 2 (lengths).}
}

\references{ Hansen, K.D., Irizarry, R.A. and Wu Z., Removing technical
  variability in RNA-seq data using conditional quantile normalization,
  Johns Hopkins, Dept of Biostatistics Working Papers. Working Paper
  227, \url{http://www.bepress.com/jhubiostat/paper227} }

\author{ Kasper Daniel Hansen, Zhijin Wu}
\note{
  Internally, the function uses a custom implementation of subset
  quantile normalization, contained in the (not exported) \code{SQN2} function.
}

\seealso{
  The package vignette.
}
\examples{
data(montgomery.subset)
data(sizeFactors.subset)
data(uCovar)
cqn.subset <- cqn(montgomery.subset, lengths = uCovar$length, 
                  x = uCovar$gccontent, sizeFactors = sizeFactors.subset,
                  verbose = TRUE)
}
\keyword{models}
