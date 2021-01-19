

KM.MQE.unpair <- function(Y0, Delta0, Z0,  maxiter=1000, ep=1e-4)
{

  Y = sort(Y0);
  n = length(Y);
  m = nrow(Z0);
  Z = cbind(rep(1, m), Z0);
  Delta = Delta0[order(Y0)];

  beta_ols = rep( 0.1, ncol(Z) );
  # lm(Y~Z-1)$coefficients;

 if(sum(1-Delta) == 0)
 {
    beta_0 = beta_ols;
    Converge=0;  iter=1;  error_0 = 0;

    S.T = survfit( Surv(Y, Delta)~ 1 )$surv;
    F.T = unique(1-S.T)[ 1: min(m, sum(Delta))];
    Y_tau = quantile(x = Y, probs = F.T);

    #while( iter <= maxiter )
    #{
      bx = c( Z %*% beta_0 );
      Z.perm = matrix(Z[order(bx),], nrow=m);

      Obj_fun <- function(x)
                 {
                   B = c(Z.perm %*% x);
                   Zx_tau = quantile(x = B, probs = F.T);
                   return( sum( (log(Y_tau) - Zx_tau)^2 ) );
                 }

      beta_new = optim( par = beta_0, fn = Obj_fun, method = "Nelder-Mead" )$par;

      #Zx_tau = quantile(x = exp(Z.perm %*% beta_new), probs = F.T);
      #error = mean( (Y_tau - Zx_tau )^2);

      #if( abs(error_0 - error) < ep )
      #{
      #  Converge = 1;
      #  break;
      #}else{
      #  beta_0 = beta_new;
      #  error_0 = error;
      #  iter = iter + 1;
      #}
    #}

 } else
 {
    S.T = survfit( Surv(Y, Delta)~ 1 )$surv;
    Y.e = Y[Delta==1];
    F.T = unique(1-S.T)[ 1: min(m, length(unique(Y.e)) ) ];

    d.tau = abs( F.T - c(0, F.T[-length(F.T)]) );

    beta_0 = beta_ols;  Converge=0;  iter=1;  error_0 = 0;
    Y_tau = unique(Y.e)[ 1: min(m, length(unique(Y.e)) ) ];

    #while( iter <= maxiter )
    #{
      bx = c( Z %*% beta_0 );
      Z.perm = matrix(Z[order(bx),], nrow=m);

      Obj_fun <- function(x)
                 {
                   B = c(Z.perm %*% x);
                   Zx_tau = quantile(x = B, probs = F.T);
                   return( sum( d.tau * ( log(Y_tau) - Zx_tau)^2 ) );
                 }

      beta_new = optim( par = beta_0, fn = Obj_fun, method = "Nelder-Mead" )$par;
      #Zx_tau = quantile(x = exp(Z.perm %*% beta_new), probs = F.T);

      #error = mean(d.tau * ( Y_tau - Zx_tau )^2);
      #if( abs(error_0 - error) < ep )
      #{
      #  Converge = 1;
      #  break;
      #}else{
      #  beta_0 = beta_new;
      #  error_0 = error;
      #  iter = iter + 1;
      #}
    #}

 }

 F_est = ecdf( exp(Z %*% beta_new) );

 return( list(beta_est=beta_new, F_est=F_est) );
}


