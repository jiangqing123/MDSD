
KM.MQE <- function(Y0, Delta0, Z0, maxiter=1000, ep=1e-4)   # k.d,
{

  Y = sort(Y0);
  n = length(Y);
  Z = cbind( rep(1,n), matrix(Z0[order(Y0),], nrow=n) );
  Delta = Delta0[order(Y0)];
  beta_ols = lm(log(Y)~Z-1)$coefficients;

 if(sum(1-Delta) == 0)
 {
    beta_0 = beta_ols;
    Converge=0;  iter=1;  error_0 = 0;

    while( iter <= maxiter )
    {
      bx = c( Z %*% beta_0 );
      Z.perm = matrix(Z[order(bx),], nrow=n);
      beta_new = lm(log(Y) ~ Z.perm-1)$coefficients;
      error = mean((log(Y) - Z.perm %*% beta_new )^2);
      if( abs(error_0 - error) < ep )
      {
        Converge = 1;
        break;
      }else{
        beta_0 = beta_new;
        error_0 = error;
        iter = iter + 1;
      }
    }

 } else
 {
    S.T = survfit( Surv(log(Y), Delta)~ 1 )$surv; # survfit( Surv(Y, Delta)~ 1 )$surv;
    logY.e = log(Y[Delta==1]); # Y.e = Y[Delta==1];
    F.T = (1-S.T)[Delta == 1];
    d.tau = abs( F.T - c(0, F.T[-length(F.T)]) );

    beta_0 = beta_ols;  Converge=0;  iter=1;  error_0 = 0;

    while( iter <= maxiter )
    {
      bx = c( Z %*% beta_0 );
      Z.perm = matrix(Z[order(bx),], nrow=n);

      Obj_fun <- function(x)
                 {
                   # B = exp(c(Z.perm %*% x));
                   B = c(Z.perm %*% x);
                   Zx_tau = quantile(x = B, probs = F.T);
                   # return( sum( (Y.e - Zx_tau)^2 ) );
                   return( sum( d.tau * (logY.e - Zx_tau)^2 ) );
                 }

      # beta_new = optimize( f=Obj_fun, interval=c(-1,4) )$minimum;
      beta_new = optim( par = beta_0, fn = Obj_fun, method = "Nelder-Mead" )$par;

      error = mean(( log(Y)-Z.perm %*% beta_new )[Delta == 1]^2);  #  max(abs(beta_0 - beta_new));
      if( abs(error_0 - error) < ep )
      {
        Converge = 1;
        break;
      }else{
        beta_0 = beta_new;
        error_0 = error;
        iter = iter + 1;
      }
    }

 }

 # bx = Z %*% beta_new;
 # error_new = mean((log(Y) - bx[order(bx)] )[Delta ==1]^2);
 # F_n = ecdf( log(Y[Delta ==1]) );
 # Ui = F_n(Z %*% beta_new);
 # rho = 1;
 # k = n * k.d;
 # for( j in 1:floor(n/k) )
 # {
 #   Cj = sum( ( Ui > (j-1)*k.d ) & ( Ui <= j*k.d ) )/n;
 #   rho = rho - 0.5 * abs( Cj - k.d);
 # }
 # return( list(beta_est=beta_new, rho_est=rho, MSE=error_new) );

  # F_est = ecdf( c(Z %*% beta_New) );
  F_est = ecdf( exp(Z %*% beta_new) );

 return( list(beta_est=beta_new, F_est=F_est) );

}


