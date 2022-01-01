\name{power.rand}
\alias{power.rand}
\title{Compute Power under Complete Randomization}
\description{
\code{power.rand} computes the power of the mean-difference estimator for a completely randomized experiment with two treatment groups. The power depends on the sample size in each group, the potential outcome variation in each group, the variation of individual treatment effects (i.e., treatment effect heterogeneity), and the magnitude of the average treatment effect.
}
\usage{
power.rand(N1, N0,
  s1, s0, s.tau = 0,
  tau, alpha = 0.05,
  exact = FALSE)
}
\arguments{
  \item{N1}{The number of subjects in the treatment group.}
  \item{N0}{The number of subjects in the control group.}
  \item{s1}{Standard deviation of the potential outcomes in the treatment group.}
  \item{s0}{Standard deviation of the potential outcomes in the control group.}
  \item{s.tau}{Standard deviation of individual treatment effects. The default is 0, denoting an additive treatment effect.}
  \item{tau}{Magnitude of the average treatment effect.}
  \item{alpha}{Level at which we reject the null. The default is 0.05.}
  \item{exact}{Whether power is computed exactly or approximately. As shown in Theorem 1 of Branson, Li, and Ding (2022), power is the sum of the right-tail of a Normal distribution (which is usually large) and the left-tail (which is usually small). When exact = FALSE, power equals just the right-tail, which corresponds to the right-hand bound in Theorem 1 of Branson, Li, and Ding (2022). When exact = TRUE, power equals the sum, which corresponds to the left-hand side of Theorem 1. The default is FALSE.}
}
\value{
The power of the mean-difference estimator for a completely randomized experiment with two treatment groups. In other words, the probability that we reject the null hypothesis of no treatment effect (when there truly is a treatment effect), given particular settings for a randomized experiment.
}

\examples{
#Power when 100 subjects are in each group,
#the standard deviation is 4 in each group,
#and the average treatment effect is 2.
power.rand(N1 = 100, N0 = 100,
s1 = 4, s0 = 4, tau = 2)
#same as before, but when
#the average treatment effect is 0.8.
power.rand(N1 = 100, N0 = 100,
s1 = 4, s0 = 4, tau = 0.8)
  
#The same examples as above,
#but now with treatment effect heterogeneity.
#We set the standard deviation of treatment effects
#to be that of potential outcomes.
#Note that, compared to the previous examples,
#power increases for large treatment effects but
#decreases for small treatment effects.
#This phenomenon is discussed in Branson, Li, and Ding (2022).

power.rand(N1 = 100, N0 = 100,
s1 = 4, s0 = 4, s.tau = 4, tau = 2)
#same as before, but when
#the average treatment effect is 0.8.
power.rand(N1 = 100, N0 = 100,
s1 = 4, s0 = 4, s.tau = 4, tau = 0.8)
}

\author{Zach Branson
}

\references{
Branson, Z., Li, X., and Ding, P. (2022). Power and Sample Size Calculations for Rerandomized Experiments.
}


\keyword{power}
\keyword{complete randomization}
