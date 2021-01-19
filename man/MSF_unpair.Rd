\name{MSF.unpair}
\alias{MSF.unpair}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict the Survival Rates by Matching Survival Functions for Unpaired Data
}
\description{
Implementing a semiparametric estimation procedure, based on section 2.1 of Jiang et al.,
to predict the survival rates based on matching survival functions and non-informative censoring mechanism for unpaired data.
}
\usage{
MSF.unpair(y0, Delta0, Z0,  maxiter=1000, ep=1e-4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y0}{observations for survival time}
  \item{Delta0}{indicator for censoring or not. 1 represents non-censoring, while 0 represents censoring.}
  \item{Z0}{covariates for matching survival function}
  \item{maxiter}{the maximum number of iterations. Default to 1000.}
  \item{ep}{the absolute convergence tolerance. Default to 1e-4.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of outputs.
\item{beta_est}{the estimation of coefficients}
\item{F_est}{the estimation of distribution function}
}
\references{
Jiang, Q., Xia, Y. and Liang, B., Matching distributions for survival data.
}
\author{
Baosheng Liang, Qing Jiang
}
\note{
%%  ~~further notes~~
}

\seealso{
MSF
}
\examples{
##---- Generate right censored data ----
n <- 100
beta <- c(sqrt(2), 1);
Z1 <- rnorm(n);
Z2 <- rnorm(n);
Z00 <- cbind(Z1, Z2);   # 2-dimension covariates
Epsil <- rnorm(n);
Ti <- exp( c(Z00 \%*\% beta) + Epsil );  # accelerate failure time model
Ci <- runif(n=n, min=0.5, max=10);
Y0 <- pmin(Ci, Ti);     # observations
Delta0 <- 1 * (Ti <= Ci);   # indicator for censoring
Ind <- sample(x = 1:n, size= floor(n * 0.6), replace = FALSE);  # 1:n; #
Z0 <- Z00[Ind, ];  # missing data

##---- Estimation procedure ----
MSF.unpair(Y0, Delta0, Z0)
}



