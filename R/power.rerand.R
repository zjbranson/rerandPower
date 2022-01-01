#Computes the power of the mean-difference estimator
#for a rerandomized experiment, where
#N1 subjects are assigned to one group and
#N0 subjects assigned to another group,
#where the Mahalanobis distance is below
#some threshold a.
#The power is provided by Theorem 3 of Branson, Li, and Ding (2022);
#thus, this function simply computes what's in Theorem 3.

#Note that this function requires the library Runuran
#This function takes the following inputs:
# N1: number of subjects in group 1 (e.g. treatment)
# N0: number of subjects in group 0 (e.g. control)
# s1: standard deviation in group 1
# s0: standard deviation in group 0
# s.tau: standard deviation of treatment effects (default 0)
# tau: additive treatment effect
# alpha: level at which we reject the null (default 0.05)
# K: the number of covariates
# pa: the acceptance probability
# R2: the R-squared between covariates and outcome
# exact: whether power is computed exactly or approximately.
#        When exact = FALSE, power equals the right-hand bound
#        in Theorem 3 of Branson, Li, and Ding (2022).
#        When exact = TRUE, power equals the left-hand side.
#        (default is FALSE)
#It outputs the power of the mean-difference estimator
#under rerandomization.
power.rerand = function(N1, N0,
  s1, s0, s.tau = 0,
  tau, alpha = 0.05,
  K, pa, R2,
  exact = FALSE){
  #the group sample size proportions are
  N = N1 + N0
  p1 = N1/N; p0 = N0/N
  #Compute the true variance V
  V = (1/p1)*s1^2 + (1/p0)*s0^2 - s.tau^2
  #Compute the limit of the variance estimator, V.tilde
  V.tilde = (1/p1)*s1^2 + (1/p0)*s0^2
  #Compute the limit of the R2 estimator, R2.tilde
  R2.tilde = V*R2/V.tilde

  #Now we need to approximate the 1 - alpha/2-quantile
  #of the non-Normal distribution under rerandomization.
  #Draws are obtained using internal function sample_std_rerand();
  #by default, this uses one million draws.
  nuDraws = sample_std_rerand(n=10^6, K = K, pa = pa, R2 = R2.tilde)
  #Then, the (approximate) quantile is
  nu.quantile = as.numeric(quantile( nuDraws, prob = 1 - alpha/2 ))

  #Then, the power is defined with the empirical survival function
  #of the draws from the non-Normal distribution.
  #(This is the right-hand side of Theorem 3.)
  power = mean(nuDraws > (nu.quantile * sqrt(V.tilde/N) - tau)/sqrt(V/N))
  #If instead we set exact = TRUE, then
  #power is the left-hand side of Theorem 3.
  if(exact){
    #First we need to approximate the alpha/2 quantile
    nu.quantile.alpha2 = as.numeric(quantile( nuDraws, prob = alpha/2 ))
    #Then, the power is defined with the empirical survival function
    #and the empirical cumulative distribution function
    #of the draws from the non-Normal distribution.
    #(This is the left-hand side of Theorem 3.)
    power = mean(nuDraws > (nu.quantile * sqrt(V.tilde/N) - tau)/sqrt(V/N)) +
      mean(nuDraws <= (nu.quantile.alpha2 * sqrt(V.tilde/N) - tau)/sqrt(V/N))
  }
  return(power)
}