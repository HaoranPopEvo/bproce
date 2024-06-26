% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_evapotranspiration.R
\name{m1.ET.Sun}
\alias{m1.ET.Sun}
\title{Calculate ET (Method 1)}
\usage{
m1.ET.Sun(PET, LAI, P)
}
\arguments{
\item{PET}{monthly potential evapotranspiration (mm month-1) estimated by FAO method (function \code{"m1.PET.FAO"}).}

\item{LAI}{leaf area index}

\item{P}{precipitation (mm month-1)}
}
\value{
vector of monthly evapotranspiration (mm month-1)
}
\description{
Calculates evapotranspiration (mm month−1) using the model proposed by Sun et al. (2011)
}
\examples{
LAI <- c(0.0, 0.0, 0.2, 0.2, 0.2, 0.3, 1.3, 1.9, 1.2, 0.5, 0.3, 0.0)
PET <- c(14.26, 24.65, 64.79, 74.40, 137.02, 153.30, 164.61, 164.30, 111.30, 60.14, 40.50, 10.54)
P <- c(38.0, 42.0, 83.0, 73.0, 96.0, 46.0, 115.0, 173.0, 70.0, 85.0, 58.0, 72.0)
m1.ET.Sun(PET, LAI, P)
}
\references{
Sun, G., Alstad, K., Chen, J., Chen, S., Ford, C. R., Lin, G., Liu, C., Lu, N., McNulty, S. G., Miao, H., Noormets, A., Vose, J. M., Wilske, B., Zeppel, M., Zhang, Y., and Zhang, Z. (2011). A general predictive model for estimating monthly ecosystem evapotranspiration. Ecohydrology, 4(2), 245–255.\cr
Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
Higher Education Press. 103pp
}
