\name{sampleSize.rerand}
\alias{sampleSize.rerand}
\title{Compute Sample Size under Rerandomization}
\description{
\code{sampleSize.rerand} computes the sample size needed to obtain a desired level of power when we use the mean-difference estimator for a rerandomized experiment with two groups, where the Mahalanobis distance is constrained to be below a prespecified threshold. The sample size depends on the proportion of subjects assigned to each group, the potential outcome variation in each group, the variation of individual treatment effects (i.e., treatment effect heterogeneity), the magnitude of the average treatment effect, the number of covariates, the linear association between covariates and potential outcomes, and the acceptance probability used for rerandomization. 
}
\usage{
sampleSize.rerand(power = 0.8,
  p1 = 0.5, p0 = 0.5,
  s1, s0, s.tau = 0,
  tau, alpha = 0.05,
  K, pa, R2)
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
  \item{K}{The number of covariates.}
  \item{pa}{The acceptance probability used for rerandomization.}
  \item{R2}{The R-squared between covariates and potential outcomes.}
}

\value{
The sample size needed to obtain a desired level of power when we use the mean-difference estimator for a rerandomized experiment with two groups, where the Mahalanobis distance is constrained to be below a prespecified threshold. In other words, given a particular \code{power}, this yields the sample size such that the probability we reject the null hypothesis of no treatment effect (when there truly is a treatment effect) is equal to \code{power}.
}

\examples{
#The sample size needed for 0.8 power when
#the standard deviation is 4 in each group,
#the average treatment effect is 2,
#there are 50 covariates, and R2 = 0.3.
sampleSize.rerand(s1 = 4, s0 = 4, tau = 2,
  K = 50, R2 = 0.3, pa = 0.01)
#same as before, but when
#the average treatment effect is 0.8.
sampleSize.rerand(s1 = 4, s0 = 4, tau = 0.8,
  K = 50, R2 = 0.3, pa = 0.01)
  
#The same examples as above,
#but now with treatment effect heterogeneity.
#We set the standard deviation of treatment effects
#to be that of potential outcomes.
#Note that, compared to the previous examples,
#sample size always decreases.
#This will always happen when power > 0.5;
#this is discussed in Branson, Li, and Ding (2022).
sampleSize.rerand(s1 = 4, s0 = 4, s.tau = 4, tau = 2,
  K = 50, R2 = 0.3, pa = 0.01)
#same as before, but when
#the average treatment effect is 0.8.
sampleSize.rerand(s1 = 4, s0 = 4, s.tau = 4, tau = 0.8,
  K = 50, R2 = 0.3, pa = 0.01)
}

\author{Zach Branson
}

\references{
Branson, Z., Li, X., and Ding, P. (2022). Power and Sample Size Calculations for Rerandomized Experiments.
}


\keyword{sample size}
\keyword{rerandomization}
