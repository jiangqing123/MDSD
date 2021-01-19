
# -------------------------------------------------------------------------- #
#    OLS & MQE     #
#
#  Applying the MQE method by Sgouropoulos et al. (2015) for the estimation
# -------------------------------------------------------------------------- #

Naive1.MQE <- function(Y0, Delta0, Z0, maxiter=1000, ep=1e-4)
{

  Y = sort(Y0);
  n = length(Y);
  Delta = Delta0[order(Y0)];
  Z = matrix(Z0[order(Y0),], nrow=n);

  beta_0 = lm(log(Y)~Z)$coefficients[-1];
  Converge=0;  iter=1;  error_0 = 0;

  while( iter <= maxiter )
  {
    bx = c( Z %*% beta_0 );
    Z.perm = matrix(Z[order(bx),], nrow=n);
    beta_new = lm(log(Y) ~ Z.perm)$coefficients;

    error_new = mean(( log(Y) - Z.perm %*% beta_new[-1] )^2);
    if( abs(error_0 - error_new) < ep )
    {
      Converge = 1;
      break;
    }else{
      beta_0 = beta_new[-1];
      error_0 = error_new;
      iter = iter + 1;
    }
  }

  # F_est = ecdf( c(Z %*% beta_New) );
  F_est = ecdf( exp(cbind( rep(1,n), matrix(Z0[order(Y0),], nrow=n) ) %*% beta_new) );

  return( list( beta_est = beta_new, F_est=F_est) );

}
