\name{montgomery.subset}
\alias{montgomery.subset}
\alias{sizeFactors.subset}
\alias{uCovar}
\docType{data}
\title{
  Mongtomery RNA-seq data.
}
\description{
  A gene by sample count matrix for 10 samples from from Montgomery et al.
  Also included is information about these genes (length and gc content)
  as well as sequencing depth for each of the samples.
}
\usage{
data(montgomery.subset)
data(sizeFactors.subset)
data(uCovar)
}
\format{
  \code{montgomery.subset} is a data frame with 23552 observations on 10
  different samples, the column names are the sample ids.
  \code{sizeFactors.subset} a a named vector of length 10 containing the
  number of mapped reads for each of the 10 samples.  \code{uCovar} is a
  data frame with 23552 observations on 2 different covariates: gc
  content and genic length in bp.
}
\details{
  Gene models are union models based on Ensembl 61.  These gene models
  were constructed using Genominator.  Genes that have zero counts in
  all 10 samples were excluded.  
}
\references{
  SB Montgomery, M Sammeth, M Gutierrez-Arcelus, RP Lach, C Ingle, J
  Nisbett, R Guigo, ET Dermitzakis, (2010) 
  ``Transcriptome genetics using second generation sequencing in a Caucasian
  population''. Nature 464(7289), 773-777.
}
\keyword{datasets}
