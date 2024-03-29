% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic_functions.R
\name{get_survival_splines_cloglog}
\alias{get_survival_splines_cloglog}
\title{Estimate risk and survival function splines at initial site, cloglog link}
\usage{
get_survival_splines_cloglog(
  data,
  time = "time",
  event = "delta",
  degree = 2,
  knots = 5,
  max_time = max(data[, time]),
  bknots = c(0, max_time),
  weights = rep(1, nrow(data)),
  confint_times = NULL,
  confint_restriction = 0
)
}
\arguments{
\item{data}{dataframe containing censoring and survival time data}

\item{time}{name of survival time column (default is "time")}

\item{event}{name of censoring indicator column (default is "delta")}

\item{degree}{spline degrees of freedom}

\item{knots}{spline knot locations}

\item{max_time}{time to estimate spline to}

\item{bknots}{boundary knots (two values)}

\item{weights}{vector of weights}

\item{confint_times}{times to store squared sum of influence functions for variance estimation in weighted case}

\item{confint_restriction}{pivot for hybrid integration (should be >0 if steep hazard observed near 0) in weighted case where confidence intervals desired}
}
\value{
spline parameters for risk and survival functions
}
\description{
Estimate risk and survival function splines at initial site, cloglog link
}
