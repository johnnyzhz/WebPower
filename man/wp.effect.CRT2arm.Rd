% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/webpower.R
\name{wp.effect.CRT2arm}
\alias{wp.effect.CRT2arm}
\title{Effect size calculatator based on raw data for Cluster Randomized Trials with 2 Arms}
\usage{
wp.effect.CRT2arm(file)
}
\arguments{
\item{file}{The input data set.
The first column of the data is the ID variable, the second column represents cluster,
the third column is the outcome variable, and the fourth column is the condition variable (0 for control, 1 for condition).}
}
\value{
A list including effect size f and ICC.
}
\description{
This function is for effect size and ICC calculation for CRT with 2 arms based on empirical data.
Cluster randomized trials (CRT) are a type of multilevel design for the situation when the entire cluster is randomly assigned to either a treatment arm or a contral arm \cite{(Liu, 2013)}.
The data from CRT can be analyzed in a two-level hierachical linear model, where the indicator variable for treatment assignment is included in second level. If a study contains multiple treatments, then mutiple indicators will be used.
This function is for designs with 3 arms (i.e., two treatments and a control).
Details leading to power calculation can be found in \cite{Raudenbush (1997)} and \cite{Liu (2013)}.
The Effect size f specifies the main effect of treatment, the mean difference between the treatment clusters and the control clusters.
This function is used to calculate the effect size with a input data set.
}
\examples{
#Empirical data set CRT2:
CRT2
#ID cluster score group
#1 1 6 0
#2 1 2 0
#3 2 6 1
#4 2 5 1
#5 3 1 0
#6 3 4 0
#7 4 6 1
#8 4 4 1

#To calculate the effect size and ICC based on empirical data
wp.effect.CRT2arm (CRT2)
#  Effect size for CRT2arm
#  
#         f  ICC
#  1.264911 -0.5
#  
#  NOTE: f is the effect size.
#  URL: http://psychstat.org/crt2arm

}
\references{
Liu, X. S. (2013). Statistical power analysis for the social and behavioral sciences: basic and advanced techniques. Routledge.

Raudenbush, S. W. (1997). Statistical analysis and optimal design for cluster randomized trials. Psychological Methods, 2(2), 173.
}
