% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_evapotranspiration.R
\name{m8.PET.HargreavesSamani}
\alias{m8.PET.HargreavesSamani}
\title{Calculate PET using Hargreaves-Samani Model}
\usage{
m8.PET.HargreavesSamani(ta, tmax, tmin, Ra)
}
\arguments{
\item{ta}{daily mean air temperature (◦C)}

\item{tmax}{maximum air temperature (◦C)}

\item{tmin}{minimum air temperature (◦C)}

\item{Ra}{extraterrestrial solar radiation (MJ m−2 d −1)}
}
\value{
vector of daily potential evapotranspiration (mm d-1)
}
\description{
Calculates daily potential evapotranspiration (mm d-1) based on Hargreaves-Samani model.
}
\examples{
ta <- c(-3.05, -1.53, 5.31, 8.01, 15.99, 21.17, 23.17, 23.65,
19.66, 13.25, 7.61, -2.39)
tmax <- ta + 10
tmin <- ta - 5
Ra <- c(-0.87, 27.20, 78.73, 95.00, 149.09, 163.70, 174.69,
185.77, 125.35, 66.43, 33.21, 0.22) # W m-2
Ra <- Ra * (1e-6*24*3600) #convert to MJ m-2 d-1
m8.PET.HargreavesSamani(ta, tmax, tmin, Ra)
}
\references{
Hargreaves, G. H., and Samani, Z. A. (1982). Estimating potential evapotranspiration. Journal of the Irrigation and Drainage Division, 108(3), 225–230\cr
Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
Higher Education Press. 103pp
}
