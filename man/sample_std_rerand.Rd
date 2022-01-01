\name{sample_std_rerand}
\alias{sample_std_rerand}
\title{Draw from Rerandomization Distribution}
\description{
\code{sample_std_rerand} is an internal function that draws from the non-Normal distribution under rerandomization; this is a mixture of a Normal and a truncated Normal.
}
\usage{
sample_std_rerand(n, K, pa, R2)
}
\arguments{
  \item{n}{The number of draws.}
  \item{K}{The number of covariates.}
  \item{pa}{The acceptance probability used for rerandomization.}
  \item{R2}{The R-squared between covariates and potential outcomes.}
}
\value{
Random draws from the non-Normal distribution under rerandomization.
}

\author{Xinran Li}

\references{
Branson, Z., Li, X., and Ding, P. (2022). Power and Sample Size Calculations for Rerandomized Experiments.
}


\keyword{internal}
