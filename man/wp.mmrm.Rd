\name{wp.mmrm}
\alias{wp.mmrm}
\alias{wp.mmrm.ar1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power analysis for longitudinal data analysis
}
\description{
The two functions are adapted from the R package longpower by Michael C. Donohue. More will be added later.
}
\usage{
wp.mmrm(N = NULL, Ra = NULL, ra = NULL, sigmaa = NULL, 
Rb = NULL, rb = NULL, sigmab = NULL, lambda = 1, 
delta = NULL, alpha = 0.05, power = NULL, 
alternative = c("two.sided", "one.sided"))
wp.mmrm.ar1(N = NULL, rho = NULL, ra = NULL, sigmaa = NULL, 
rb = NULL, sigmab = NULL, lambda = 1, times = 1:length(ra), 
delta = NULL, alpha = 0.05, power = NULL, 
alternative = c("two.sided", "one.sided"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
N
}
  \item{Ra}{
Ra
}
  \item{ra}{
ra
}
  \item{sigmaa}{
sigmaa
}
  \item{Rb}{
Rb
}
  \item{rb}{
rb
}
  \item{sigmab}{
sigmab
}
  \item{lambda}{
lambda
}
  \item{delta}{
delta
}
  \item{alpha}{
Significance level
}
  \item{power}{
Statistical power
}
  \item{alternative}{
alternative hypothesis
}
\item{rho}{rho}
\item{times}{times}
}

\references{
Lu, K., Luo, X., Chen, P.-Y. (2008) Sample size estimation for repeated measures analysis in randomized clinical trials with missing data. International Journal of Biostatistics, 4, (1)
}

