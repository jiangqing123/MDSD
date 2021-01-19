\name{Naive1.MQE}
\alias{Naive1.MQE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matching a Distribution by Matching Quantiles Estimation with All Data
}
\description{
Implementing a semiparametric estimation procedure, proposed by Sgouropoulos et al.,
to match a target distribution by that of a linear combination of a set of random variables, by matching quantiles estimation.
}
\usage{
Naive1.MQE(Y0, Delta0, Z0, maxiter=1000, ep=1e-4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y0}{observations for survival time}
  \item{Delta0}{indicator for censoring or not. 1 represents non-censoring, while 0 represents censoring.}
  \item{Z0}{covariates for matching survival function}
  \item{maxiter}{the maximum number of iterations. Default to 1000.}
  \item{ep}{the absolute convergence tolerance. Default to 1e-4.}
}
\details{
No matter the observations are censoring or not, all data are regarded as non-censoring. Then the method proposed by Sgouropoulos et al. is applied to estimate the coefficients of linear combination.
}
\value{
A list of outputs.
\item{beta_est}{the estimation of coefficients}
\item{F_est}{the estimation of distribution function}
}
\references{
Sgouropoulos, N., Yao, Q. and Yastremiz, C. (2015), Matching a Distribution by Matching Quantiles Estimation, Journal of the American Statistical Association, 110(510), 742-759.
}
\author{
Baosheng Liang, Qing Jiang
}
\note{
%%  ~~further notes~~
}

\seealso{

}

\examples{
##---- Generate right censored data ----
n <- 100
beta <- c(sqrt(2), 1);
Z1 <- rnorm(n);
Z2 <- rnorm(n);
Z0 <- cbind(Z1, Z2);  # 2-dimension covariates
Epsil <- rnorm(n);
Ti <- exp( c(Z0 \%*\% beta) + Epsil );  # accelerate failure time model
Ci <- runif(n=n, min=0.5, max=10);
Y0 <- pmin(Ci, Ti);     # observations
Delta0 <- 1 * (Ti <= Ci);   # indicator for censoring

##---- Estimation procedure ----
Naive1.MQE(Y0, Delta0, Z0)
}

