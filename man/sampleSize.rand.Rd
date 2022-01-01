\name{sampleSize.rand}
\alias{sampleSize.rand}
\title{Compute Sample Size under Complete Randomization}
\description{
\code{sampleSize.rand} computes the sample size needed to obtain a desired level of power when we use the mean-difference estimator for a completely randomized experiment with two groups. The sample size depends on the proportion of subjects assigned to each group, the potential outcome variation in each group, the variation of individual treatment effects (i.e., treatment effect heterogeneity), and the magnitude of the average treatment effect.
}
\usage{
sampleSize.rand(power = 0.8,
  p1 = 0.5, p0 = 0.5,
  s1, s0, s.tau = 0,
  tau, alpha = 0.05)
}
\arguments{
  \item{power}{The desired level of power, ranging between 0 and 1. The default is 0.8.}
  \item{p1}{The proportion of subjects in the treatment group. The default is 0.5.}
  \item{p0}{The proportion of subjects in the control group. The default is 0.5.}
  \item{s1}{Standard deviation of the potential outcomes in the treatment group.}
  \item{s0}{Standard deviation of the potential outcomes in the control group.}
  \item{s.tau}{Standard deviation of individual treatment effects. The default is 0, denoting an additive treatment effect.}
  \item{tau}{Magnitude of the average treatment effect.}
  \item{alpha}{Level at which we reject the null. The default is 0.05.}
}

\value{
The sample size needed to obtain a desired level of power when we use the mean-difference estimator for a completely randomized experiment with two groups. In other words, given a particular \code{power}, this yields the sample size such that the probability we reject the null hypothesis of no treatment effect (when there truly is a treatment effect) is equal to \code{power}.
}

\examples{
#The sample size needed for 0.8 power when
#the standard deviation is 4 in each group,
#and the average treatment effect is 2.
sampleSize.rand(s1 = 4, s0 = 4, tau = 2)
#same as before, but when
#the average treatment effect is 0.8.
sampleSize.rand(s1 = 4, s0 = 4, tau = 0.8)
  
#The same examples as above,
#but now with treatment effect heterogeneity.
#We set the standard deviation of treatment effects
#to be that of potential outcomes.
#Note that, compared to the previous examples,
#sample size always decreases.
#This will always happen when power > 0.5;
#this is discussed in Branson, Li, and Ding (2022).
sampleSize.rand(s1 = 4, s0 = 4, s.tau = 4, tau = 2)
#same as before, but when
#the average treatment effect is 0.8.
sampleSize.rand(s1 = 4, s0 = 4, s.tau = 4, tau = 0.8)
}

\author{Zach Branson
}

\references{
Branson, Z., Li, X., and Ding, P. (2022). Power and Sample Size Calculations for Rerandomized Experiments.
}


\keyword{sample size}
\keyword{complete randomization}
