% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UKF.R
\name{propagate_model}
\alias{propagate_model}
\title{propagate_model
Used inside UKF_dT, previously named kura_ode_model
Take the ode model and current augmented state
and propagate ind vars to next dT using runge-kutta}
\usage{
propagate_model(t, ode_model, dt, dT, N_p, x)
}
\arguments{
\item{t}{dummy time variable, because ode models don't have explcit time in our examples}

\item{ode_model}{model with ind vars y and N_p parameters}

\item{dt}{smaller time step size within dT}

\item{dT}{time step size that comes from time series data step}

\item{N_p}{number of unknown model params}

\item{x}{augmented state vector}
}
\value{
propagated augmented column vector
}
\description{
propagate_model
Used inside UKF_dT, previously named kura_ode_model
Take the ode model and current augmented state
and propagate ind vars to next dT using runge-kutta
}
\examples{
Example
}
