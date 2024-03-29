% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/webpower.R
\name{sem.effect.size}
\alias{sem.effect.size}
\title{Calculate the Effect Size for SEM}
\usage{
sem.effect.size(full.model.pop, reduced.model)
}
\arguments{
\item{full.model.pop}{Full model (under the alternative hypothesis) with population parameters.}
\item{reduced.model}{Reduced model (under the null hypothesis) lavaan specification.}
}
\value{
An object of the power analysis.
\item{delta}{Effect size.}
\item{df}{Degrees of freedom}
\item{RMSEA}{RMSEA}
}
\description{
This function is for calculating SEM effect size.
}
\examples{
full.model.pop <-'
y1 ~ 0.4*x
y2 ~ 0.5*x + 0.2*y1
y3 ~ 0.4*x
y4 ~ 0.4*y1 + 0.4*y2 + 0.4*y3
y1 ~~ 0.84*y1
y2 ~~ 0.61*y2
y3 ~~ 0.84*y3
y4 ~~ 0.27*y4
'

reduced.model <-'
y1 ~ x
y2 ~ x 
y3 ~ x
y4 ~ y1 + y3
'

sem.effect.size(full.model.pop, reduced.model)
}
\references{
Demidenko, E. (2007). Sample size determination for logistic regression revisited. Statistics in medicine, 26(18), 3385-3397.

Zhang, Z., & Yuan, K.-H. (2018). Practical Statistical Power Analysis Using Webpower and R (Eds). Granger, IN: ISDSA Press.
}
