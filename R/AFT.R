
# ----------------------------------------------------------------- #
#  AFT  Estimation  #
#
#  Input:
#     Y0: observations for survival time
#     Delta0: censoring or not. 1 represents non-censoring.
#     Z0: covariates
#     dist: accelerate failure time model, 'lognormal', 'weibull'
#  Output:
#     beta_est: the estimation of coefficients
# ----------------------------------------------------------------- #

AFT <- function(Y0, Delta0, Z0, dist='lognormal')
{
  Y = sort(Y0);
  n = length(Y);
  Delta = Delta0[order(Y0)];
  Z = matrix(Z0[order(Y0),], nrow=n);

  # fitAFT = survreg(Surv(Y, Delta) ~ Z, dist='weibull');
  fitAFT = survreg(Surv(Y, Delta) ~ Z, dist=dist);
  beta_new = fitAFT$coefficients;

  F_est = ecdf( exp(cbind( rep(1,n), matrix(Z0[order(Y0),], nrow=n) ) %*% beta_new) );

  return( list( beta_est = beta_new, F_est=F_est) );
}


