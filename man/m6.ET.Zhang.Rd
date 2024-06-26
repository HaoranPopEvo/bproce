% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_evapotranspiration.R
\name{m6.ET.Zhang}
\alias{m6.ET.Zhang}
\title{Calculate ET (Method 6)}
\usage{
m6.ET.Zhang(PET, P, w = 2)
}
\arguments{
\item{PET}{potential evapotranspiration (mm month-1) estimated by Priestley-Taylor method (function \code{"m6.PET.PriestleyTaylor"}).}

\item{P}{predicitation (mm month-1)}

\item{w}{plant-available water coefficient that represents the relative difference in plant water use for transpiration. Default 2 for forests while 0.5 for grassland when PET is estimated using the Priestley–Taylor model (function "\code{m6.PET.PriestleyTaylor}").}
}
\value{
vector of monthly evapotranspiration (mm month-1)
}
\description{
Calculates evapotranspiration (mm month−1) using the model proposed by Zhang et al. (2001)
}
\examples{
PET <- c(14.26, 24.65, 64.79, 74.40, 137.02, 153.30, 164.61, 164.30, 111.30, 60.14, 40.50, 10.54)
P <- c(38.0, 42.0, 83.0, 73.0, 96.0, 46.0, 115.0, 173.0, 70.0, 85.0, 58.0, 72.0)
m6.ET.Zhang(PET, P)
}
\references{
Zhang, L., Dawes, W. R., and Walker, G. R. (2001). Response of mean annual evapotranspiration to vegetation changes at catchment scale. Water Resources Research, 37(3), 701–708.
Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
Higher Education Press. 104pp
}
