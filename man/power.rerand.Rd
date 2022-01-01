\name{power.rerand}
\alias{power.rerand}
\title{Compute Power under Rerandomization}
\description{
\code{power.rerand} computes the power of the mean-difference estimator for a rerandomized experiment with two treatment groups, where the Mahalanobis distance is constrained to be below a prespecified threshold. The power depends on the sample size in each group, the potential outcome variation in each group, the variation of individual treatment effects (i.e., treatment effect heterogeneity), the magnitude of the average treatment effect, the number of covariates, the linear association between covariates and potential outcomes, and the acceptance probability used for rerandomization. 
}
\usage{
power.rerand(N1, N0,
  s1, s0, s.tau = 0,
  tau, alpha = 0.05,
  K, pa, R2,
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
  \item{K}{The number of covariates.}
  \item{pa}{The acceptance probability used for rerandomization.}
  \item{R2}{The R-squared between covariates and potential outcomes.}
  \item{exact}{Whether power is computed exactly or approximately. As shown in Theorem 3 of Branson, Li, and Ding (2022), power is the sum of the right-tail of a non-Normal distribution (which is usually large) and the left-tail (which is usually small). When exact = FALSE, power equals just the right-tail, which corresponds to the right-hand bound in Theorem 3 of Branson, Li, and Ding (2022). When exact = TRUE, power equals the sum, which corresponds to the left-hand side of Theorem 3. The default is FALSE.}
}
\value{
The power of the mean-difference estimator for a rerandomized experiment with two treatment groups, where the Mahalanobis distance is constrained to be below a prespecified threshold. In other words, the probability that we reject the null hypothesis of no treatment effect (when there truly is a treatment effect), given particular settings for a rerandomized experiment.
}

\examples{
#Power when 50 subjects are in each group,
#the standard deviation is 4 in each group,
#the average treatment effect is 2,
#there are 10 covariates,
#covariates are moderately related with outcomes,
#and the acceptance probability is 0.01.
power.rerand(N1 = 50, N0 = 50,
s1 = 4, s0 = 4, tau = 2,
K = 10, pa = 0.01, R2 = 0.3)
#same as before, but when
#the average treatment effect is 0.8.
power.rerand(N1 = 50, N0 = 50,
s1 = 4, s0 = 4, tau = 0.8,
K = 10, pa = 0.01, R2 = 0.3)
  
#The same examples as above,
#but now with treatment effect heterogeneity.
#We set the standard deviation of treatment effects
#to be that of potential outcomes.
#Note that, compared to the previous examples,
#power increases for large treatment effects but
#decreases for small treatment effects.
#This phenomenon is discussed in Branson, Li, and Ding (2022).

power.rerand(N1 = 50, N0 = 50,
s1 = 4, s0 = 4, s.tau = 4, tau = 2,
K = 10, pa = 0.01, R2 = 0.3)
#same as before, but when
#the average treatment effect is 0.8.
power.rerand(N1 = 50, N0 = 50,
s1 = 4, s0 = 4, s.tau = 4, tau = 0.8,
K = 10, pa = 0.01, R2 = 0.3)
}

\author{Zach Branson
}

\references{
Branson, Z., Li, X., and Ding, P. (2022). Power and Sample Size Calculations for Rerandomized Experiments.
}


\keyword{power}
\keyword{rerandomization}
