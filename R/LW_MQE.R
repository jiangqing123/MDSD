
# ------------------------------- #
#  Locally Weighted & MQE         #
# ------------------------------- #

LW.MQE = function(Y0, Delta0, Z0, maxiter=1000, ep=1e-4)
{

  Y = sort(Y0);
  n = length(Y);
  Delta = Delta0[order(Y0)];
  Z = cbind( rep(1, n),  matrix(Z0[order(Y0),], nrow=n) );
  p = ncol(Z);

  beta_ols = lm(log(Y)~Z-1)$coefficients;

  if( sum(1-Delta) < 1 )
  {
    #print("Censoring Rates equal to 0.");
    beta_0 =  beta_ols;
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

  } else {

     n.c = sum(1 - Delta);

     # Y.e = Y[Delta == 1];
     # Y.c = Y[Delta == 0];
     # Y.r = Y;  Y.r[Delta == 0] = max(Y)^{100};

     # -- Biquadratic Kernal Fun. -- #

     h = 80* n^{-1/3}; # sd( log(Y) ) ; 50; 100;

   B_nk.f = function(x0, x, h, kernel.type="4th")
   {
     # the kernel weight function Bnk(x0, x), where x0 is a scalar, and x is a vector
     # returns a vector
     # h is the bandwidth

     if(!is.vector(x0)) x0 = as.vector(x0);

     xx = (x-x0)/h;
     if(kernel.type == "4th"){
       xx[abs(xx) >= 1] = 1;
       w = 15/16 * (1 - xx^{2})^{2};  #biquadratic kernel
     }
     w = w/sum(w);
     return(w);
   }

   Kernel.f = function(U0, U, h, kernel.type="4th")
   {
     # U: n*k matrix
     # U0: 1*k matrix
     # return: K((U-U0)/h)
     if(!is.vector(U0)) U0 = as.vector(U0);
     n = nrow(U);
     if(kernel.type=="4th"){
        tt = rbind(U, U0);
        tmp = apply(tt, 2, function(x) {
                B_nk.f(x0 = x[n+1], x = x[1:n], h = h, kernel.type = kernel.type)
            });
        tmp = apply(tmp, 1, prod);
        tmp = tmp/sum(tmp);
     }
     return(tmp);
   }


   F.T = function(y0, x0, y, x, delta, h, kernel.type = "4th")
   {
     # tau0(y0, x0) = F(T<y0|x0); so y0 is the C_i, and x0 is the xi in the paper
     # x0: k-dimensional covariate vector
     # y: n-vector of observed survival time = T^C
     # x: k-dimensional covariate matrix
     # delta: the censoring indicator function
     # h: bandwidth parameter

     n = length(y);

     # --   kernel weights   -- ##
     p = qr(x)$rank;

     if(p >1) Bn = Kernel.f(U0=x0, U=x, h=h, kernel.type)
     else Bn = B_nk.f(x0=x0, x=x, h=h, kernel.type);

     if (y0 < max(y)) {
        # sort the data y, and the delta, Bn correspondingly to the order of sorted y
        y2 = sort(y);
        Order = order(y); # so z[Order] = z2
        Bn2 = Bn[Order];
        delta2 = delta[Order];
        eta = which(delta2==1 & y2<=y0); # the index of those observations satisfying delta2==1 & z2<=y0
        Bn3 = Bn2[n:1];  # change the order of Bn2, make the first obs of Bn2 to be the last of Bn3
        tmp = 1 - Bn2 /cumsum(Bn3)[n:1];
        out = 1 - prod(tmp[eta], na.rm=T); # na.rm=T, as some of those tmp=NA as the denom =0
     }
     else out = 1;
     return(out)
   }

 # ---------------------------------- #
   F_logY = rep(0, n);

   for(i in 1:n)
   {
     F_logY[i] = F.T(y0=log(Y[i]), x0=Z[i,, drop=FALSE], y=log(Y), x=Z, delta=Delta, h=h, kernel.type = "4th");
   }

   LCRQ <- function(y, delta, tau, kernel.type = "4th"){
      # Locally weighted censored quantile regression method
      # x is a design matrix
      # y is the observed survival time = min(T, C)
      # delta is the censoring indicator function with 1 standing for uncensored, and 0 censored
      # tau is the quantile level of interest

      n = length(y);
      ind = which(delta == 0);
      w = rep(1, n); # the weight vector

      if(length(ind) >= 1){
        for(i in 1:length(ind)){
            tau.star = F_logY[ind[i]];
            if (tau > tau.star) w[ind[i]] = (tau - tau.star) / (1-tau.star)
            else w[ind[i]] = 1;
        }
        # pseudo observations
        ind2 = which(w != 1);
        y.pse = rep(max(y)+100, length(ind2));

        yy = c(y, y.pse);
        ww = c(w, 1-w[ind2]);
     }
     else{
        yy = y;
        ww = w;
    }

    rq1 = rq(yy~1, weights=ww, tau=tau);
    return(rq1$coeff);
  }

   # ---------------------------------- #
   #logY.pse = rep( max(log(Y))+100, n.c );
   #logYY = c( log(Y), logY.pse );

   tau_i = F_logY[Delta==1];
   Q_logY = rep(0, length(tau_i) );

   for( ti in 1:length(tau_i) )
   {
      Q_logY[ti] = LCRQ(y=log(Y), delta=Delta, tau=tau_i[ti], kernel.type = "4th");
   }

   # ---------------------------------- #
    beta_0 = beta_ols;

    Obj_fun <- function(x)
    {
      zx = Z %*% x;
      Q_zx = as.numeric( quantile(x=zx, probs=tau_i) );
      Obj = sum( (Q_logY - Q_zx)^2 );
      return(Obj)
    }

    beta_new = optim( par=beta_0, fn=Obj_fun, method="Nelder-Mead")$par;    ##

  }

  # F_est = ecdf( c(Z %*% beta_New) );
  F_est = ecdf( exp(Z %*% beta_new) );

 return( list( beta_est = beta_new, F_est=F_est) );
}


