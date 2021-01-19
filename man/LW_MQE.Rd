\name{LW.MQE}
\alias{LW.MQE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict the Survival Rates by Mathcing Censored Quantiles for Conditional independent censoring
}
\description{
Implementing a semiparametric estimation procedure, based on section 2.2 of Jiang et al.,
to predict the survival rates based on matching censored quantiles for conditional independent censoring.
}
\usage{
LW.MQE(Y0, Delta0, Z0, maxiter=1000, ep=1e-4)
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
KM.MQE
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
Ci <- exp( c(Z0 \%*\% c(0.2, 0.3)) + runif(n=n, min=0, max=3.5) );  # conditional independent
Y0 <- pmin(Ci, Ti);     # observations
Delta0 <- 1 * (Ti <= Ci);   # indicator for censoring

##---- Estimation procedure ----
LW.MQE(Y0, Delta0, Z0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line



