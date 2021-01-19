
# ----------------------------------------------------------------------------- #
#  Naive OLS & MQE
#
#  Applying the same MQE method as in Naive1, but it abandons all the censored
#  observations and only adopts those uncensored ones.
# ----------------------------------------------------------------------------- #

Naive2.MQE <- function(Y0, Delta0, Z0, maxiter=1000, ep=1e-4)
{
   Y = sort(Y0);
   n = length(Y);
   Delta = Delta0[order(Y0)];
   Z = matrix(Z0[order(Y0),], nrow=n);

   Y.e = Y[Delta == 1];
   n.e = length(Y.e);
   Z.e = matrix(Z[Delta == 1, ], nrow=n.e );

   beta_0 = lm(log(Y.e)~Z.e)$coefficients[-1];
   Converge=0;  iter=1;  error_0 = 0;

   while( iter <= maxiter )
   {
      bx = c( Z.e %*% beta_0 );
      Z.e.perm = matrix(Z.e[order(bx),], nrow=n.e);
      beta_new = lm(log(Y.e)~Z.e.perm)$coefficients;

      error_new = mean( ( log(Y.e) - Z.e.perm %*% beta_new[-1] )^2 );
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

   F_est = ecdf( exp(cbind( rep(1,n), matrix(Z0[order(Y0),], nrow=n) ) %*% beta_new) );

   return(list( beta_est = beta_new, F_est=F_est))
}




