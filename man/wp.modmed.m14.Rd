% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model14.R
\name{wp.modmed.m14}
\alias{wp.modmed.m14}
\title{model14}
\usage{
wp.modmed.m14(
  a1,
  cp,
  b1,
  d1,
  b2,
  sige12,
  sige22,
  n,
  sigx_w,
  sigx2 = 1,
  sigw2 = 1,
  nrep_power = 1000,
  alpha = 0.05,
  b = 1000,
  nb = n,
  w_value = 0,
  power_method = "product",
  simulation_method = "percentile",
  ncore = 1,
  pop.cov = NULL,
  mu = NULL,
  varnames = c("y", "x", "w", "m", "mw")
)
}
\arguments{
\item{a1}{regression coefficient of mediator (m) on predictor (x)}

\item{cp}{regression coefficient of outcome (y) on predictor (x)}

\item{b1}{regression coefficient of outcome (y) on mediator (m)}

\item{d1}{regression coefficient of outcome (y) on moderator (w)}

\item{b2}{regression coefficient of outcome (y) on the product (mw)}

\item{sige12}{variance of error in the first regression equation}

\item{sige22}{variance of error in the second regression equation}

\item{n}{sample size}

\item{sigx_w}{covariance between predictor (x) and moderator (w)}

\item{sigx2}{variance of predictor (x)}

\item{sigw2}{variance of moderator (w)}

\item{nrep_power}{number of replications for finding power}

\item{alpha}{type 1 error rate}

\item{b}{number of bootstrap iterations}

\item{nb}{bootstrap sample size, default to n, used when simulation method is "percentile"}

\item{w_value}{moderator level, value of w}

\item{power_method}{"product" for using the indirect effect value in power calculation, or "joint" for using joint significance in power calculation}

\item{simulation_method}{"percentile" for using percentile bootstrap CI in finding significance of mediation, or "MC" for using Monte Carlo CI in finding significance of mediation}

\item{ncore}{number of cores to use for the percentile bootstrap method, default is 1, when ncore > 1, parallel is used}

\item{pop.cov}{covariance matrix, default to NULL if using the regression coefficient approach}

\item{mu}{mean vector, default to NULL if using the regression coefficient approach}

\item{varnames}{name of variables for the covariance matrix}
}
\value{
power of indirect effect, direct effect, moderation, and the index of moderated mediation
}
\description{
power analysis of model 14 in Introduction to Mediation, Moderation, and Conditional Process Analysis. Powers are obtained through either the percentile bootstrap method or the Monte Carlo method. The conditional indirect effect value is a1(b1 + b2w); the index of moderated mediation is a1b2.
}
\examples{
test = wp.modmed.m14(a1 = 0.2, cp = 0.2, b1 = 0.5, d1 = 0.5, b2 = 0.2, sigx2 = 1,
                    sigw2 = 1, sige12 = 1, sige22 = 1, sigx_w = 0.5, n = 50,
                    w_value = 0.5, simulation_method = "MC",
                    nrep_power = 1000, alpha = 0.05, b = 1000, ncore = 1)
print(test)
}
