

MSF.unpair <- function(y0, Delta0, Z0,  maxiter=1000, ep=1e-4)
{
  n = length(y0);
  Y0 = sort(y0);
  m = nrow(Z0);
  Z = cbind(rep(1, m), Z0);
  Delta = Delta0[order(y0)];

  Lam = 0; # seq(0, 3, 0.5);
  ERR = rep(0, length(Lam) );
  BETA_est = matrix(0, nrow=length(Lam), ncol = ncol(Z));

 for(l in 1:length(Lam) )
 {
    # print(Lam[l]);

   if( Lam[l] == 0 )
   {
      Y = log(Y0);
   } else{
      Y = ( (Y0)^{Lam[l]} - 1 )/(Lam[l]);
   }

   # print( sum( is.na(Y) ) );

   if( sum( is.na(Y) ) >= 1 )
   {
      ERR[l] = 1e+30;  next;
   } else{

    beta_ols = rep( 0.1, ncol(Z) ); # lm(Y~Z-1)$coefficients;   #  rep(0, ncol(Z));   #

    if(sum(1-Delta) == 0)
    {

      F_Y = ecdf(Y);
      F.T = F_Y(Y);

      R_fun = function(x, k.d=0.05) {
                Ui = F_Y(Z %*% x);

                rho = 1;
                k = n * k.d;

                for( j in 1:floor(n/k) )
                {
                  Cj = sum( ( Ui > (j-1)*k.d ) & ( Ui <= j*k.d ) )/n;
                  rho = rho - 0.5 * abs(Cj - k.d);
                }

                return(rho);
             }

        Obj_fun <- function(x){
                   F.Bf = ecdf( c(Z %*% x) );
                   F.B = F.Bf(Y);
                   return( sum( (F.T - F.B)^2 ) )
                 }

      beta_new = optim( par = beta_ols, fn = Obj_fun, method = "Nelder-Mead" )$par;

      BETA_est[l, ] = beta_new;
      ERR[l] = R_fun(x=beta_new, k.d=0.05);

    } else
    {
      KM.fit = survfit( Surv(Y, Delta) ~ 1 );
      S.T = KM.fit$surv;
      Y.e = Y[Delta == 1];

      F.T = 1 - S.T[Delta == 1];
      F_Y <- function(x)
             {
               if(sum(Y.e <= x) >= 1)
               {
                 return( max(F.T[Y.e <= x]) )
               } else{
                 return(0)
               }
             }

      R_fun <- function(x, k.d=0.05)
               {
                 Ui = apply( matrix(Z %*% x, ncol=1), 1, F_Y );

                 rho = 1;
                 k = n * k.d;

                 for( j in 1:floor(n/k) )
                 {
                   Cj = sum( ( Ui > (j-1)*k.d ) & ( Ui <= j*k.d ) )/n;
                   rho = rho - 0.5 * abs(Cj - k.d);
                 }

                 return(rho);
               }

      Obj.f <- function(x){
                 F.Bf = ecdf( c( Z %*% x ) );
                 F.B = F.Bf(Y);
                 return( sum( (1 - S.T - F.B)^2 ) )
               }

      beta_new = optim( par = beta_ols, fn = Obj.f, method = "Nelder-Mead" )$par;

      ERR[l] = R_fun(x=beta_new, k.d=0.05);
      BETA_est[l,] = beta_new;
    }

  }

 }

 # cbind(Lam, ERR, BETA_est);

 beta_New = BETA_est[which.max(ERR),];
 # F_est = ecdf( c(Z %*% beta_New) );
 F_est = ecdf( exp(Z %*% beta_new) );

 return( list(beta_est=beta_New, F_est=F_est) ); #, opt = Lam[which.max(ERR)]

}


