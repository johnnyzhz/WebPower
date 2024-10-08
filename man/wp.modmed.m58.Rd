% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model58.R
\name{wp.modmed.m58}
\alias{wp.modmed.m58}
\title{model58}
\usage{
wp.modmed.m58(
  c1,
  a1,
  c2,
  d1,
  b1,
  b2,
  cp,
  sige12,
  sige22,
  sigx_w,
  n,
  sigx2 = 1,
  sigw2 = 1,
  nrep_power = 1000,
  alpha = 0.05,
  b = 1000,
  nb = n,
  w_value = 0,
  power_method = "product",
  ncore = 1,
  simulation_method = "percentile",
  pop.cov = NULL,
  mu = NULL,
  varnames = c("x", "w", "m", "xw", "mw", "y")
)
}
\arguments{
\item{c1}{regression coefficient of outcome (m) on moderator (w)}

\item{a1}{regression coefficient of mediator (m) on predictor (x)}

\item{c2}{regression coefficient of outcome (m) on the product (xw)}

\item{d1}{regression coefficient of outcome (y) on moderator (w)}

\item{b1}{regression coefficient of outcome (y) on mediator (m)}

\item{b2}{regression coefficient of outcome (y) on the product (mw)}

\item{cp}{regression coefficient of outcome (y) on predictor (x)}

\item{sige12}{variance of error in the first regression equation}

\item{sige22}{variance of error in the second regression equation}

\item{sigx_w}{covariance between predictor (x) and moderator (w)}

\item{n}{sample size}

\item{sigx2}{variance of predictor (x)}

\item{sigw2}{variance of moderator (w)}

\item{nrep_power}{number of replications for finding power}

\item{alpha}{type 1 error rate}

\item{b}{number of bootstrap iterations}

\item{nb}{bootstrap sample size, default to n, used when simulation method is "percentile"}

\item{w_value}{moderator level, value of w}

\item{power_method}{"product" for using the indirect effect value in power calculation, or "joint" for using joint significance in power calculation}

\item{ncore}{number of cores to use for the percentile bootstrap method, default is 1, when ncore > 1, parallel is used}

\item{simulation_method}{"percentile" for using percentile bootstrap CI in finding significance of mediation, or "MC" for using Monte Carlo CI in finding significance of mediation}

\item{pop.cov}{covariance matrix, default to NULL if using the regression coefficient approach}

\item{mu}{mean vector, default to NULL if using the regression coefficient approach}

\item{varnames}{name of variables for the covariance matrix}
}
\value{
power of indirect effect, direct effect, and moderation
}
\description{
power analysis of model 58 in Introduction to Mediation, Moderation, and Conditional Process Analysis. Powers are obtained through either the percentile bootstrap method or the Monte Carlo method. The conditional indirect effect value is (a1 + c2w)(b1 + b2w); the index of moderated mediation cannot be applied here since the conditional indirect effect is not linear.
}
\examples{
test = wp.modmed.m58(c1 = 0.2, a1 = 0.2, c2 = 0.1, b2 = 0.1,
     b1 = 0.2, cp = 0.2, d1 = 0.2, w_value = 0.3, simulation_method = "MC",
     sigx2 = 1, sigw2 = 1, sige12 = 1, sige22 = 1, sigx_w = 0.5,
     n = 50, nrep_power = 1000, b = 1000, alpha = 0.05, ncore = 1)
print(test)
}
