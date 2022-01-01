
#Internal function to draw from the distribution
#of the mean-difference estimator under rerandomization.
#This distribution is a mixture of a standard normal
#and a truncated normal distribution.
#This function takes as inputs:
# n: the number of draws
# K: the number of covariates
# pa: the acceptance probability
# R2: the R-squared between covariates and outcome
sample_std_rerand <- function(n, K, pa, R2){
  #standard normal
  epsilon0 = rnorm(n)
  #rerandomization threshold
  a = qchisq(pa, df=K)
  #draw from the square-root of a truncated chi-squared
  chi_aK = sqrt(Runuran::urchisq(n, df=K, lb=0, ub=a))
  #draw a random sign
  S = 2*rbinom(n, 1, prob=0.5) - 1
  #draw the \beta_K random variable
  if(K>=2){
    beta_K = rbeta(n, shape1=1/2, shape2=(K-1)/2)
  }else{
    beta_K = 1
  }
  #finally, draw from the distribution of the mean-difference
  #estimator under rerandomization.
  #This distribution is a mixture of a standard normal
  #and a truncated normal distribution.
  draw = sqrt(1-R2) * epsilon0 + sqrt(R2)*chi_aK*S*sqrt(beta_K)
  return(draw)
}
