\name{wp.kanova}
\alias{wp.kanova}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power analysis for two-way, three-way and k-way ANOVA
}
\description{
Power analysis for two-way, three-way and k-way ANOVA
}
\usage{
wp.kanova(n = NULL, ndf = NULL, f = NULL, ng = NULL, alpha = 0.05, power = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Sample size
}
  \item{ndf}{
Numerator degrees of freedom
}
  \item{f}{
Effect size
}
  \item{ng}{
Number of groups
}
  \item{alpha}{
Significance level
}
  \item{power}{
Statistical power
}
}

\examples{
########## Chapter 6. two-way and three-way ANOVA  #############
## Main effect of two-way ANOVA
wp.kanova(n=120, ndf=2, f=0.2, alph=0.05, ng=6)

## Interaction effect of two-way ANOVA
wp.kanova(n=120, ndf=2, f=0.4, alph=0.05, ng=6)

## Interaction effect of three-way ANOVA
wp.kanova(n=360, ndf=4, f=0.3, alph=0.05, ng=18)

}

\references{
Zhang, Z., & Yuan, K.-H. (2018). Practical Statistical Power Analysis Using Webpower and R (Eds). Granger, IN: ISDSA Press.
}
