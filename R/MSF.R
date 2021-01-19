
# -------------------------------------- #
#    MSF: Matching Survival Function     #
# -------------------------------------- #

MSF <- function(Y0, Delta0, Z0, maxiter=1000, ep=1e-4)
{
  Y = sort(Y0);
  n = length(Y);
  Delta = Delta0[order(Y0)];
  Z = cbind( rep(1, n),  matrix(Z0[order(Y0),], nrow=n) );

  if( sum(1-Delta) < 1 )
  {
     #print("Censoring Rates equal to 0.");

     beta_0 = lm(log(Y)~Z-1)$coefficients;
     S = ecdf(Y);
     S.T = S(Y);

     Obj.f <- function(x){
                 S.Bf = ecdf( exp(c( Z %*% x )) );
                 S.B = S.Bf(Y);
                 return( sum( (S.T - S.B)^2 ) )
              }

     grad <- function(x){
               S.Bf = ecdf( exp(c( Z %*% x )) );
               S.B = S.Bf(Y);
               g.sum = 0;

               for(i in 1:n)
               {
                 g.sum = g.sum + 2*(S.T[i]- S.B[i])*( ( (exp(c(Z%*%x))<=Y[i])/n * exp(c(Z%*%x)) ) %*% Z );
               }

               return(g.sum);
             }

     # beta_new = optimize(f=Obj.f, interval = c(-1, 4) )$minimum;
      beta_new = optim( par = beta_0, fn = Obj.f, method = "Nelder-Mead" )$par;

  } else {

     beta_0 = lm(log(Y)~Z-1)$coefficients;
     S.T = 1-survfit( Surv(Y, Delta)~ 1 )$surv;

     Obj.f <- function(x)
              {
                S.Bf = ecdf( exp(c( Z %*% x )) );
                S.B = S.Bf(Y);
                return( sum( (S.T - S.B)^2 ) )
              }

     beta_new = optim( par = beta_0, fn = Obj.f, method = "Nelder-Mead" )$par;
 }

  F_est = ecdf( exp(Z %*% beta_new) );

  return( list(beta_est=beta_new, F_est=F_est) );

}


