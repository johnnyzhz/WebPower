\name{wp.mc.t}
\alias{wp.mc.t}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power analysis for t-test based on Monte Carlo simulation
}
\description{
Power analysis for t-test based on Monte Carlo simulation
}
\usage{
wp.mc.t(n = NULL, R0 = 1e+05, R1 = 1000, mu0 = 0, mu1 = 0, 
sd = 1, skewness = 0, kurtosis = 3, alpha = 0.05, 
type = c("two.sample", "one.sample", "paired"), 
alternative = c("two.sided", "less", "greater"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Sample size
}
  \item{R0}{
Number of replications under the null
}
  \item{R1}{
Number of replications
}
  \item{mu0}{
Population mean under the null
}
  \item{mu1}{
Population mean under the alternative
}
  \item{sd}{
Standard deviation
}
  \item{skewness}{
Skewness
}
  \item{kurtosis}{
kurtosis
}
  \item{alpha}{
Significance level
}
  \item{type}{
Type of anlaysis
}
  \item{alternative}{
alternative hypothesis
}
}
\examples{
\donttest{
set.seed(20220722)
########## Chapter 16. Monte Carlo t-test #############
wp.mc.t(n=20 , mu0=0, mu1=0.5, sd=1, skewness=0, 
kurtosis=3, type = c("one.sample"), alternative = c("two.sided"))

wp.mc.t(n=40 , mu0=0, mu1=0.3, sd=1, skewness=1, 
kurtosis=6, type = c("paired"), alternative = c("greater"))

wp.mc.t(n=c(15, 15), mu1=c(0.2, 0.5), sd=c(0.2, 0.5), 
skewness=c(1, 2), kurtosis=c(4, 6), type = c("two.sample"), alternative = c("less"))
}

}
\references{
Zhang, Z., & Yuan, K.-H. (2018). Practical Statistical Power Analysis Using Webpower and R (Eds). Granger, IN: ISDSA Press.
}