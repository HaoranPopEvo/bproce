% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_weather.R
\name{sunRS}
\alias{sunRS}
\title{Calculate Sun Rise/Set Time}
\usage{
sunRS(doy, long, lat)
}
\arguments{
\item{doy}{day of the year. A scalar.}

\item{long}{longtitude of the site (-180 to 180 degree). A scalar.}

\item{lat}{latitude of the site (-90 to 90 degree). A scalar.}
}
\value{
vector of two values. The first one is sunrise time and the second one
is sunset time.
}
\description{
calculates sunrise and sunset time (h).
}
\examples{
doy <- DOY(2016, 9, 22)
sunRS(doy, 85.44675833, 42.47670833)
}
\references{
Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
Higher Education Press. 33pp
}
