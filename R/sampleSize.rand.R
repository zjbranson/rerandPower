#Computes the sample size N needed
#to obtain a desired level of power
#when we use the mean-difference estimator
#for a completely randomized experiment.
#The sample size is provided by Theorem 2 of Branson, Li, and Ding (2022);
#thus, this function simply computes what's in Theorem 2.

#This function takes the following inputs:
# power: the desired level of power (default 0.8)
# p1: the proportion of subjects in group 1 (e.g. treatment)
#     (default is 0.5)
# p0: the proportion of subjects in group 0 (e.g. control)
#     (default is 0.5)
# s1: standard deviation in group 1
# s0: standard deviation in group 0
# s.tau: standard deviation of treatment effects (default 0)
# tau: additive treatment effect
# alpha: level at which we reject the null (default 0.05)
#It outputs the sample size needed to obtain power
#using the mean-difference estimator under complete randomization.
sampleSize.rand = function(power = 0.8, p1 = 0.5, p0 = 0.5,
  s1, s0, s.tau = 0,
  tau, alpha = 0.05){
  #Compute the true variance V
  V = (1/p1)*s1^2 + (1/p0)*s0^2 - s.tau^2
  #Compute the limit of the variance estimator, V.tilde
  V.tilde = (1/p1)*s1^2 + (1/p0)*s0^2

  #the normal-distribution quantiles are
  z.alpha = qnorm(1 - alpha/2)
  z.gamma = qnorm(1-power)

  #then, the desired sample size is
  N = ((z.alpha*sqrt(V.tilde) - z.gamma*sqrt(V) )/tau)^2
  return(N)
}