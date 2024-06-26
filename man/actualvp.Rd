% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_weather.R
\name{actualvp}
\alias{actualvp}
\title{Calculate Actual Vapor Pressure}
\usage{
actualvp(Ta, h, ...)
}
\arguments{
\item{Ta}{vector of air temperature (◦C).}

\item{h}{vector of relative humidity (\%).}

\item{...}{additional parameters accepted by function \code{satvp}.}
}
\value{
vector of actual vapor pressure (kPa).
}
\description{
calculates actual vapor pressure (kPa). Only valid for temperature above
0◦C.
}
\examples{
data("p1_tah")
temperature <- p1_tah$Ta
humidity <- p1_tah$h
actualvp(temperature, humidity)
}
