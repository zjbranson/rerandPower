#Computes the sample size N needed
#to obtain a desired level of power
#when we use the mean-difference estimator
#for a rerandomized experiment.
#The sample size is provided by Theorem 5 of Branson, Li, and Ding (2022);
#thus, this function simply computes what's in Theorem 5.

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
# K: the number of covariates
# pa: the acceptance probability
# R2: the R-squared between covariates and outcome
#It outputs the sample size needed to obtain power
#using the mean-difference estimator under rerandomization.
sampleSize.rerand = function(power = 0.8, p1 = 0.5, p0 = 0.5,
  s1, s0, s.tau = 0,
  tau, alpha = 0.05,
  K, pa, R2){
  #Compute the true variance V
  V = (1/p1)*s1^2 + (1/p0)*s0^2 - s.tau^2
  #Compute the limit of the variance estimator, V.tilde
  V.tilde = (1/p1)*s1^2 + (1/p0)*s0^2
  #Compute the limit of the R2 estimator, R2.tilde
  R2.tilde = V*R2/V.tilde

  #Now we need to approximate quantiles of
  #the non-Normal distribution under rerandomization.
  #Draws are obtained using internal function sample_std_rerand();
  #by default, this uses one million draws.
  #Note that the sample size depends on two different quantiles:
  #one where R2 = R2, and one where R2 = R2.tilde
  nuDraws = sample_std_rerand(n=10^6, R2 = R2, K = K, pa = pa)
  nuDraws.tilde = sample_std_rerand(n=10^6, R2 = R2.tilde, K = K, pa = pa)
  #Then, the (approximate) quantiles are
  nu.alpha = as.numeric(quantile(nuDraws.tilde, prob = 1 - alpha/2 ))
  nu.gamma = as.numeric(quantile(nuDraws, prob = 1 - power))
  #Note that, indeed, the alpha quantile depends on R2.tilde,
  #but the gamma quantile depends on R2.

  #then, the desired sample size is
  N = (( nu.alpha*sqrt(V.tilde) - nu.gamma*sqrt(V) )/tau)^2
  return(N)
}