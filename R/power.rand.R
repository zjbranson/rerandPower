#Computes the power of the mean-difference estimator
#for a completely randomized experiment, where
#N1 subjects are assigned to one group and
#N0 subjects assigned to another group.
#The power is provided by Theorem 1 of Branson, Li, and Ding (2022);
#thus, this function simply computes what's in Theorem 1.

#Under the null of no effect,
# \hat{\tau} \sim N(0, S1^2/N1 + S0^2/N0 - Stau^2/N)
#where S1^2 and S0^2 are the variances in each group,
#and N = N1 + N0.
#Under the alternative,
# \hat{\tau} \sim N(\tau, S1^2/N1 + S0^2/N0 - Stau^2/N)

#Using these, we can compute power, i.e.,
#the probability we reject the null
#under the alternative for a given tau.

#This function takes the following inputs:
# N1: number of subjects in group 1 (e.g. treatment)
# N0: number of subjects in group 0 (e.g. control)
# s1: standard deviation in group 1
# s0: standard deviation in group 0
# s.tau: standard deviation of treatment effects (default 0)
# tau: average treatment effect
# alpha: level at which we reject the null (default 0.05)
# exact: whether power is computed exactly or approximately.
#        When exact = FALSE, power equals the right-hand bound
#        in Theorem 1 of Branson, Li, and Ding (2022).
#        When exact = TRUE, power equals the left-hand side.
#        (default is FALSE)
#It outputs the power of the mean-difference estimator
#under complete randomization.
power.rand = function(N1, N0,
  s1, s0, s.tau = 0,
  tau, alpha = 0.05,
  exact = FALSE){
  #the group sample size proportions are
  N = N1 + N0
  p1 = N1/N; p0 = N0/N
  #Compute the true variance V
  V = (1/p1)*s1^2 + (1/p0)*s0^2 - s.tau^2
  #Compute the limit of the variance estimator, V.tilde
  V.tilde = (1/p1)*s1^2 + (1/p0)*s0^2
  #The power then is (right-hand side of Theorem 1)
  power = 1 - pnorm( (qnorm(1-alpha/2)*sqrt(V.tilde/N) - tau)/sqrt(V/N) )
  #If instead we set exact = TRUE, then
  #power is the left-hand side of Theorem 1.
  if(exact){
    power = 1 - pnorm( (qnorm(1-alpha/2)*sqrt(V.tilde/N) - tau)/sqrt(V/N) )  +
  pnorm( (qnorm(alpha/2)*sqrt(V.tilde/N) - tau)/sqrt(V/N) )
  }
  return(power)
}