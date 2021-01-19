\name{AFT}
\alias{AFT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Regression for Accelerate Failure Time Model}
\description{
Fit an accelerate failure time model.
}
\usage{
AFT(Y0, Delta0, Z0, dist='lognormal')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y0}{observations for survival time}
  \item{Delta0}{indicator for censoring or not. 1 represents non-censoring, while 0 represents censoring.}
  \item{Z0}{covariates for accelerate failure time model}
  \item{dist}{assumed distribution for Y0 variable.  If the argument is a character string, then it is assumed to name an element from survreg.distributions. These include "weibull", "exponential", "gaussian", "logistic","lognormal" and "loglogistic". Otherwise, it is assumed to be a user defined list conforming to the format described in survreg.distributions.}
}
\details{
Implementing a parametric estimation procedure, based on chapter 2.2 of Kalbfleisch and Prentice,
to estimate the coefficients of an accelerate failure model.
Then the distribution of Y0 could be predicted.
}
\value{
A list of outputs.
\item{beta_est}{the estimation of coefficients}
\item{F_est}{the estimation of distribution function}
}
\references{
Kalbfleisch, J. D. and Prentice, R. L., The statistical analysis of failure time data, Wiley, 2002.
}
\author{
Baosheng Liang, Qing Jiang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
AFT(Y0, Delta0, Z0)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~AFT }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~AFT }% __ONLY ONE__ keyword per line
