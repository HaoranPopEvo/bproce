##Script info----------------------
# Created on 2022-06-04 11:57:51 CST
#
# author: Hao-Ran Wu
# College of Life Sciences
# Zhejiang University
# Hangzhou 310012, China
# E-mail: haoranwu@zju.edu.cn
#
# This script provides function to calculate evapotranspiration (ET) based on
#  different models. In most of cases, potential evapotranspiration (PET) is
#  calculated first before calculating ET.
#
#  Step1: Calculate PET
#       function                   model
#    -- 'm1.PET.FAO'               FAO Penman-Monteith
#    -- 'm2.PET.Thornthwaite'      Thornthwaite, 1948
#    -- 'm3.PET.Hamon'             Hamon, 1963
#    -- 'm4.PET.BlaneyCriddle'     Blaney-Criddle, 1950s
#    -- 'm5.PET.Turc'              Turc, 1961
#    -- 'm6.PET.PriestleyTaylor'   Priestley-Taylor, 1972
#    -- 'm7.PET.Makkink'           Makkink, 1957
#    -- 'm8.PET.HargreavesSamani'  Samani, 1982
#
#  Step2: Calculate ET
#    (For some methods, PET is not required)
#    -- 'm1.ET.Sun'                Sun et al., 2011a
#    -- 'm2.ET.Sun'                Sun et al., 2011b
#    -- 'm3.ET.Fang'               Fang et al., 2016
#    -- 'm4.ET.type1'              Fang et al., 2016
#    -- 'm5.ET.type2'              Fang et al., 2016
#    -- 'm6.ET.Zhang'              Zhang et al., 2001

##m1.PET.FAO--------------------------
#' Calculate PET using Penman-Monteith Model
#'
#' Calculates potential evapotranspiration (to be more exactly,
#' grass reference ET, mm). Only valid for temperature above 0◦C. In
#' additional, this model assumes a stand that has a 0.12 m canopy height,
#' a leaf area index (LAI) of 4.8, a bulk surface resistence of 70 s m-1,
#' and an albedo of 0.23.
#' @param ta vector of air temperature (◦C).
#' @param vpd vapor pressure deficient (kPa).
#' @param Rn vector of net radiation (MJ m-2)
#' @param u2 vector of wind speed (m s-1) at 2 m height.
#' @param G vector of soil heat flux (MJ m−2).
#' @param EL elevation (m).
#' @param Cn numerator constant thta changes with reference surface and
#' calculation time step (900 ◦C mm s−3 Mg−1 d−1 for 24 h time steps, and
#' 37 ◦C mm s−3 Mg−1 d−1 for hourly time steps).
#' @param Cd denominator constant that changes with reference surface and
#' calculation steps (0.34 s m−1 for 24 h time steps, 0.24 s m-1 for hourly
#' time steps during daytime, and 0.96 s m-1 for hourly nighttime for grass
#' reference surface).
#'
#' @return vector of potential evapotranspiration (mm d-1).
#' @references
#' Allen, R. G., Smith, M., Perrier, A., and Pereira, L. S. (1994). An update for the definition of reference evapotranspiration. International Commission on Irrigation and Drainage Bulletin, 43(2), 1–34.\cr
#' Djaman, K., Koudahe, K., Lombard, K., and O’Neill, M. (2018). Sum of hourly vs daily Penman–Monteith grass-reference evapotranspiration under semiarid and arid climate. Irrigation and Drainage Engineering, 7(1), doi: 10.4172/2168-9768.1000202.\cr
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 98pp
#' @export
#' @examples
#' ta <- c(21.17, 23.17, 23.65, 19.66)
#' vpd <- c(1.15, 0.95, 0.84, 0.67)
#' Rn <- c(163.70, 174.69, 185.77, 125.35)*24*3600/1e6 #convert from W m-2 to MJ m-2 d-1
#' u2 <- c(2.04, 1.79, 1.52, 1.75)
#' G <- c(-5.52, -6.64, -8.13, -1.17)*24*3600/1e6 #convert from W m-2 to MJ m-2 d-1
#' EL <- 0
#' m1.PET.FAO(ta, vpd, Rn, u2, G, EL)
m1.PET.FAO <- function(ta, vpd, Rn, u2, G, EL, Cn=900, Cd=0.34){
  delta <- (4098 * 0.6108) * exp(17.27*ta/(ta+237.3)) / (ta + 237.3)^2
  #slope of the saturation water vapor pressure at air temperature (T, kPa ◦C−1)

  #λ is the amount of energy (enthalpy) that must be added to
  #   transform a given quantity of water into a gas (i.e., vapor)
  p <- 101.3 - 0.01055*EL
  ga <- (0.001013 * p)/(0.622 * (2.501 - 0.002361 * ta))
  #thermodynamic psychrometric constant (0.000666 ◦C-1)

  (0.408 * delta * (Rn-G) + ga * Cn/(ta+273) * u2 * vpd) /
    (delta + ga * (1+Cd * u2))
}

##m2.PET.Thornthwaite--------------------------
#' Calculate PET using Thornthwaite Model
#'
#' Calculates monthly potential evapotranspiration (cm) based on method
#' proposed by Thornthwaite (1948), the most widely used temperature-based monthly
#' scale PET model.
#' @param month_ta a vector of 12 elements, each showing the mean air temperature for a month within a year.
#' @param Ld a vector of 12 elements, each showing the mean daytime length (12·h), it is time from sunrise to sunset in multiples of 12 hours.
#' @return vector of monthly potential evapotranspiration (cm).
#' @references
#' Thornthwaite, C. W. (1948). An approach toward a Rational classification of climate. Geographical Review, 38(1), 55–94.\cr
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 99pp
#' @export
#' @examples
#' ta <- c(1.05, 4.53, 9.31, 12.01, 19.99, 25.17, 27.17, 27.65, 23.66, 17.25, 11.61, 2.39)
#' doy <- DOY(2016, 9, 22)
#' daylength <- diff(sunRS(doy, 85.44675833, 42.47670833))[[1]]
#' daylength <- daylength/12 #hours to multiples of 12 hours
#' m2.PET.Thornthwaite(ta, daylength)
m2.PET.Thornthwaite <- function(month_ta, Ld){
  if(length(month_ta)!=12) stop("'month_ta' should be a vector of 12 elements.")
  if(length(Ld)!=12) stop("'Ld' should be a vector of 12 elements.")
  if(any(length(month_ta)<0)) warning("temperature of some months is lower than 0, which is not allowed.")
  heatI <- sum((month_ta/5)^1.514) #annual heat index
  a <- 6.75 * 10^(-7) * heatI^3 - 7.71 * 10^(-5) * heatI^2 + 0.01792 * heatI + 0.49239
  1.6 * Ld * (10 * month_ta/heatI)^a
}

##m3.PET.Hamon--------------------------
#' Calculate PET using Hamon's model
#'
#' Calculates daily potential evapotranspiration (mm) based on Hamon's PET
#' model.
#' @param ta mean air temperature (◦C)
#' @param Ld mean daytime length (12·h), it is time from sunrise to sunset in multiples of 12 hours.
#' @return vector of daily potential evapotranspiration (mm).
#' @references
#' Hamon, W. R. (1963). Computation of direct runoff amounts from storm rainfall. International Association of Scientific Hydrology Publication, 63, 52–62.\cr
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 100pp
#' @export
#' @examples
#' ta <- c(-3.05, -1.53, 5.31, 8.01, 15.99, 21.17, 23.17, 23.65, 19.66, 13.25, 7.61, -2.39)
#' doy <- DOY(2016, 9, 22)
#' daylength <- diff(sunRS(doy, 85.44675833, 42.47670833))[[1]]
#' daylength <- daylength/12 #hours to multiples of 12 hours
#' m3.PET.Hamon(ta, daylength)
m3.PET.Hamon <- function(ta, Ld){
  es <- satvp(ta)*10 #from kPa to mb
  0.1651 * Ld * (216.7 * es/(ta + 273.3))
}

##m4.PET.BlaneyCriddle--------------------------
#' Calculate PET using Blaney-Criddle model
#'
#' Calculates daily potential evapotranspiration (mm) based on Blaney-Criddle model.
#' It is used to estimate ET for the western USA. Calculated PET is the potential
#' water use for a reference crop.
#'
#' @param ta daily mean temperature (◦C). A vector.
#' @param dhour the total daytime hours for the period used. A vector.
#' @return vector of potential water use (PET) for a reference crop (mm).
#' @references
#' Doorenbos, J., and Pruitt, W. O. (1977). Guidelines for Predicting Crop Water Requirements. FAO Irrigation and Drainage Paper, 24 (revised). Rome.\cr
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 100pp
#' @export
#' @examples
#' m4.PET.BlaneyCriddle(21.17, 12)
m4.PET.BlaneyCriddle <- function(ta, dhour){
  dhour/(365*12) * (0.46 * ta + 8.13)*100
}

##m5.PET.Turc--------------------------
#' Calculate PET using Turc model
#'
#' Calculates daily evapotranspiration (mm d-1) based on Turc model.
#'
#' @param ta daily mean air temperature (◦C)
#' @param h relative humidity (%). Value from 0-100. A vector.
#' @param Rs daily solar radiation (MJ m-2 d-1). A vector.
#' @return vector of potential water use (PET) for a reference crop (mm).
#' @references
#' Doorenbos, J., and Pruitt, W. O. (1977). Guidelines for Predicting Crop Water Requirements. FAO Irrigation and Drainage Paper, 24 (revised). Rome.\cr
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 100pp
#' @export
#' @examples
#' ta <- c(-3.05, -1.53, 5.31, 8.01, 15.99, 21.17, 23.17, 23.65,
#' 19.66, 13.25, 7.61, -2.39)
#' h <- c(73, 74, 68, 63, 58, 60, 71, 75, 74, 77, 74, 81)
#' Rs <- c(-0.87, 27.20, 78.73, 95.00, 149.09, 163.70, 174.69,
#' 185.77, 125.35, 66.43, 33.21, 0.22) # W m-2
#' Rs <- Rs * (1e-6*24*3600) #convert to MJ m-2 d-1
#' m5.PET.Turc(ta, h, Rs)
m5.PET.Turc <- function(ta, h, Rs){
  Rs <- 100*Rs/4.1868
  ifelse(
    h<50,
    0.013*(ta/(ta+15))*(Rs+50)*(1+(50-h)/70),
    0.013*(ta/(ta+15))*(Rs+50)
  )
}

##m6.PET.PriestleyTaylor--------------------------
#' Calculate PET using Priestley-Taylor model
#'
#' Calculates daily potential evapotranspiration (mm d-1) based on Priestley-Taylor model.
#'
#' @param ta daily mean air temperature (◦C). A vector.
#' @param EL elevation (m). A vector.
#' @param G heat flux density to the ground (MJ m−2 d −1). A vector.
#' @param Rn vector of net radiation (MJ m-2). A vector.
#' @param alpha calibration constant, \code{alpha} = 1.26 for wet or humid conditions.
#' @return vector of daily potential evapotranspiration (mm d-1)
#' @references
#' Priestley, C. H. B., and Taylor, R. J. (1972). On the assessment of surface heat flux and evaporation using large-scale parameters. Monthly Weather Review, 100(2), 81–92.\cr
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 101pp
#' @export
#' @examples
#' ta <- c(-3.05, -1.53, 5.31, 8.01, 15.99, 21.17, 23.17, 23.65,
#' 19.66, 13.25, 7.61, -2.39)
#' Rn <- c(-0.87, 27.20, 78.73, 95.00, 149.09, 163.70, 174.69,
#' 185.77, 125.35, 66.43, 33.21, 0.22) # W m-2
#' Rn <- Rn * (1e-6*24*3600) #convert to MJ m-2 d-1
#' G <- -5.52 * (1e-6*24*3600) #convert to MJ m-2 d-1
#' m6.PET.PriestleyTaylor(ta, 0, G, Rn)
m6.PET.PriestleyTaylor <- function(ta, EL, G, Rn, alpha=1.26){
  lambda <- 2.501 - 0.002361 *ta #latent heat of vaporization (MJ kg−1)
  delta <- 0.2*(0.00738*ta+0.8072)^7-0.000116 #slope of the saturation vapor pressure-temperature curve (kPa ◦C−1)
  cp <- 0.001013 #specific heat of moist air at constant pressure (MJ kg−1 ◦C −1)
  P <- 101.3 - 0.01055 * EL # atmospheric pressure (kPa)
  ga <- cp*P/(0.622*lambda) # psychrometric constant modified by the ratio of canopy resistance to atmospheric resistance (kPa ◦C−1)
  alpha * delta/(delta + ga)*(Rn-G)/lambda
}

##m7.PET.Makkink--------------------------
#' Calculate PET using Makkink model
#'
#' Calculates daily potential evapotranspiration (mm d-1) based on Makkink model.
#'
#' @param ta daily mean air temperature (◦C)
#' @param EL elevation (m)
#' @param Rs solar radiation (MJ m-2 d-1)
#' @return vector of daily potential evapotranspiration (mm d-1)
#' @references
#' Makkink, G. F. (1957). Testing the Penman formula by means of lysimeters. Journal of the Institution of Water Engineerrs, 11, 277–288.
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 102pp
#' @export
#' @examples
#' ta <- c(-3.05, -1.53, 5.31, 8.01, 15.99, 21.17, 23.17, 23.65,
#' 19.66, 13.25, 7.61, -2.39)
#' Rs <- c(-0.87, 27.20, 78.73, 95.00, 149.09, 163.70, 174.69,
#' 185.77, 125.35, 66.43, 33.21, 0.22) # W m-2
#' Rs <- Rs * (1e-6*24*3600) #convert to MJ m-2 d-1
#' m7.PET.Makkink(ta, 0, Rs)
m7.PET.Makkink <- function(ta, EL, Rs){
  delta <- 0.2*(0.00738*ta+0.8072)^7-0.000116
  lambda <- 2.501 - 0.002361 *ta
  P <- 101.3 - 0.01055 * EL
  cp <- 0.001013
  ga <- cp*P/(0.622*lambda)

  0.61 *(delta/(delta+ga))*Rs/lambda - 0.12
}

##m8.PET.HargreavesSamani--------------------------
#' Calculate PET using Hargreaves-Samani Model
#'
#' Calculates daily potential evapotranspiration (mm d-1) based on Hargreaves-Samani model.
#'
#' @param ta daily mean air temperature (◦C)
#' @param tmax maximum air temperature (◦C)
#' @param tmin minimum air temperature (◦C)
#' @param Ra extraterrestrial solar radiation (MJ m−2 d −1)
#' @return vector of daily potential evapotranspiration (mm d-1)
#' @references
#' Hargreaves, G. H., and Samani, Z. A. (1982). Estimating potential evapotranspiration. Journal of the Irrigation and Drainage Division, 108(3), 225–230\cr
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 103pp
#' @export
#' @examples
#' ta <- c(-3.05, -1.53, 5.31, 8.01, 15.99, 21.17, 23.17, 23.65,
#' 19.66, 13.25, 7.61, -2.39)
#' tmax <- ta + 10
#' tmin <- ta - 5
#' Ra <- c(-0.87, 27.20, 78.73, 95.00, 149.09, 163.70, 174.69,
#' 185.77, 125.35, 66.43, 33.21, 0.22) # W m-2
#' Ra <- Ra * (1e-6*24*3600) #convert to MJ m-2 d-1
#' m8.PET.HargreavesSamani(ta, tmax, tmin, Ra)
m8.PET.HargreavesSamani <- function(ta, tmax, tmin, Ra){
  lambda <- 2.501 - 0.002361 *ta
  0.0023 * Ra *(tmax-tmin)^0.5 *(ta + 17.8)/lambda
}

##m1.ET.Sun--------------------------
#' Calculate ET (Method 1)
#'
#' Calculates evapotranspiration (mm month−1) using the model proposed by Sun et al. (2011)
#'
#' @param PET monthly potential evapotranspiration (mm month-1) estimated by FAO method (function \code{"m1.PET.FAO"}).
#' @param LAI leaf area index
#' @param P precipitation (mm month-1)
#' @return vector of monthly evapotranspiration (mm month-1)
#' @references
#' Sun, G., Alstad, K., Chen, J., Chen, S., Ford, C. R., Lin, G., Liu, C., Lu, N., McNulty, S. G., Miao, H., Noormets, A., Vose, J. M., Wilske, B., Zeppel, M., Zhang, Y., and Zhang, Z. (2011). A general predictive model for estimating monthly ecosystem evapotranspiration. Ecohydrology, 4(2), 245–255.\cr
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 103pp
#' @export
#' @examples
#' LAI <- c(0.0, 0.0, 0.2, 0.2, 0.2, 0.3, 1.3, 1.9, 1.2, 0.5, 0.3, 0.0)
#' PET <- c(14.26, 24.65, 64.79, 74.40, 137.02, 153.30, 164.61, 164.30, 111.30, 60.14, 40.50, 10.54)
#' P <- c(38.0, 42.0, 83.0, 73.0, 96.0, 46.0, 115.0, 173.0, 70.0, 85.0, 58.0, 72.0)
#' m1.ET.Sun(PET, LAI, P)
m1.ET.Sun <- function(PET, LAI, P){
  11.94 + 4.76 * LAI + PET * (0.032 * LAI + 0.0026 * P + 0.15)
}

##m2.ET.Sun--------------------------
#' Calculate ET (Method 2)
#'
#' Calculates evapotranspiration (mm month−1) using the model proposed by Sun et al. (2011)
#'
#' @param PET monthly potential evapotranspiration (mm month-1) estimated by Hamon's method (function \code{"m3.PET.Hamon"}).
#' @param LAI leaf area index
#' @param P precipitation (mm month-1)
#' @return vector of monthly evapotranspiration (mm month-1)
#' @references
#' Sun, G., Caldwell, P., Noormets, A., McNulty, S. G., Cohen, E., Moore Myers, J. M., Domec, J-C., Treasure, E., Mu, Q., Xiao, J., John, R., and Chen, J. (2011). Upscaling key ecosystem functions across the conterminous United States by a watercentric ecosystem model. Journal of Geophysical Research: Biogeosciences, 116(G3), https://doi.org/10.1029/2010JG001573.
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 103pp
#' @export
#' @examples
#' LAI <- c(0.0, 0.0, 0.2, 0.2, 0.2, 0.3, 1.3, 1.9, 1.2, 0.5, 0.3, 0.0)
#' PET <- c(14.26, 24.65, 64.79, 74.40, 137.02, 153.30, 164.61, 164.30, 111.30, 60.14, 40.50, 10.54)
#' P <- c(38.0, 42.0, 83.0, 73.0, 96.0, 46.0, 115.0, 173.0, 70.0, 85.0, 58.0, 72.0)
#' m2.ET.Sun(PET, LAI, P)
m2.ET.Sun <- function(PET, LAI, P){
  0.174 * P + 0.502 * PET + 5.31 * LAI + 0.0222 * PET * LAI
}

##m3.ET.Fang--------------------------
#' Calculate ET (Method 3)
#'
#' Calculates evapotranspiration (mm month−1) using the model proposed by Fang et al. (2016)
#'
#' @param PET monthly potential evapotranspiration (mm month-1) estimated by Hamon's method (function \code{"m3.PET.Hamon"}).
#' @param VPD vapor pressure deficit (kPa)
#' @param Rn solar radiation (MJ m-2 month-1)
#' @return vector of monthly evapotranspiration (mm month-1)
#' @references
#' Fang, Y., Sun, G., Caldwell, P., McNulty, S. G., Noormets, A., Domec, J. C., King, J., Zhang, Z., Zhang, X., Lin, G., Zhou, G., Xiao, J., Chen, J., and Zhou, G. (2016). Monthly land cover-specific evapotranspiration models derived from global eddy flux measurements and remote sensing data. Ecohydrology, 9(2), 248–266.
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 104pp
#' @export
#' @examples
#' PET <- c(-2.26, 70.5, 204.07, 246.24, 386.44, 424.31, 452.8, 481.52, 324.91, 172.19, 86.08, 0.57)
#' VPD <- c(0.15, 0.17, 0.34, 0.55, 0.91, 1.15, 0.95, 0.84, 0.67, 0.40, 0.32, 0.11)
#' Rn <- c(-0.87, 27.20, 78.73, 95.00, 149.09, 163.70, 174.69, 185.77, 125.35, 66.43, 33.21, 0.22)
#' m3.ET.Fang(PET, VPD, Rn)
m3.ET.Fang <- function(PET, VPD, Rn){
  0.42 + 0.74 * PET - 2.73 * VPD + 0.10 * Rn
}

##m4.ET.type1--------------------------
#' Calculate ET (Method 4)
#'
#' Calculates evapotranspiration (mm month−1) using the model proposed by Fang et al. (2016)
#'
#' @param Rn monthly net solar radiation (MJ month-1)
#' @param PET potential evapotranspiration (mm month-1) estimated by Hamon's method (function \code{"m3.PET.Hamon"}).
#' @param LAI leaf area index
#' @param P predication (mm month-1)
#' @param SWC soil water content (%)
#' @param Ta mean air temperature (◦C)
#' @param land.type type of land. "\code{shrub}" for shrubland,
#' "\code{crop}" for cropland, "\code{grass}" for grassland,
#' "\code{deciduous}" for deciduous forest, "\code{needle}" for evergreen needle leaf forest,
#' "\code{broad}" for evengreen broad leaf forest, "\code{mixed}" for mixed forest, "\code{savannas}" for savannas.
#' @return vector of monthly evapotranspiration (mm month-1)
#' @references
#' Fang, Y., Sun, G., Caldwell, P., McNulty, S. G., Noormets, A., Domec, J. C., King, J., Zhang, Z., Zhang, X., Lin, G., Zhou, G., Xiao, J., Chen, J., and Zhou, G. (2016). Monthly land cover-specific evapotranspiration models derived from global eddy flux measurements and remote sensing data. Ecohydrology, 9(2), 248–266.
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 104pp
#' @export
#' @examples
#' Rn <- c(-0.87, 27.20, 78.73, 95.00, 149.09, 163.70, 174.69, 185.77, 125.35, 66.43, 33.21, 0.22)
#' LAI <- c(0.0, 0.0, 0.2, 0.2, 0.2, 0.3, 1.3, 1.9, 1.2, 0.5, 0.3, 0.0)
#' P <- c(38.0, 42.0, 83.0, 73.0, 96.0, 46.0, 115.0, 173.0, 70.0, 85.0, 58.0, 72.0)
#' m4.ET.type1(Rn, LAI = LAI, P = P, land.type = "shrub")
m4.ET.type1 <- function(Rn, PET=NULL, LAI=NULL, P=NULL, SWC=NULL, Ta=NULL,
                        land.type = c("shrub","crop","grass","deciduous","needle","broad","mixed","savanna")){
  land.type <- land.type[1]
  if(land.type=="shrub"){
    -4.59 + 13.02 * LAI + 0.10 * Rn + 0.11 * P
  } else if(land.type=="crop"){
    0.87 + 0.19 * Rn + 13.99 * LAI + 0.06 * P
  } else if(land.type=="grass"){
    -0.87 + 0.2 * Rn + 0.1 * P + 0.24 * SWC
  } else if(land.type=="deciduous"){
    -14.22 + 0.74 * PET + 0.1 * Rn
  } else if(land.type=="needle"){
    13.47 + 0.1 * Rn + 1.35 * Ta
  } else if(land.type=="broad"){
    0.01 + 0.63 * Ta + 0.46 * SWC + 0.14 * Rn
  } else if(land.type=="mixed"){
    -8.76 + 0.95 * PET
  } else if(land.type=="savannas"){
    -8.87 + 33.46 * LAI + 0.07 * Rn
  } else{
    stop("invalid option 'land.type'")
  }
}

##m5.ET.type2--------------------------
#' Calculate ET (Method 5)
#'
#' Calculates evapotranspiration (mm month−1) using the model proposed by Fang et al. (2016)
#'
#' @param PET potential evapotranspiration (mm month-1) estimated by Hamon's method (function \code{"m3.PET.Hamon"}).
#' @param LAI leaf area index
#' @param P predicitation (mm month-1)
#' @param land.type type of land. "\code{shrub}" for shrubland,
#' "\code{crop}" for cropland, "\code{grass}" for grassland,
#' "\code{deciduous}" for deciduous forest, "\code{needle}" for evergreen needle leaf forest,
#' "\code{broad}" for evengreen broad leaf forest, "\code{mixed}" for mixed forest, "\code{savannas}" for savannas.
#' @return vector of monthly evapotranspiration (mm month-1)
#' @references
#' Fang, Y., Sun, G., Caldwell, P., McNulty, S. G., Noormets, A., Domec, J. C., King, J., Zhang, Z., Zhang, X., Lin, G., Zhou, G., Xiao, J., Chen, J., and Zhou, G. (2016). Monthly land cover-specific evapotranspiration models derived from global eddy flux measurements and remote sensing data. Ecohydrology, 9(2), 248–266.
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 104pp
#' @export
#' @examples
#' LAI <- c(0.0, 0.0, 0.2, 0.2, 0.2, 0.3, 1.3, 1.9, 1.2, 0.5, 0.3, 0.0)
#' PET <- c(14.26, 24.65, 64.79, 74.40, 137.02, 153.30, 164.61, 164.30, 111.30, 60.14, 40.50, 10.54)
#' P <- c(38.0, 42.0, 83.0, 73.0, 96.0, 46.0, 115.0, 173.0, 70.0, 85.0, 58.0, 72.0)
#' m5.ET.type2(PET, LAI = LAI, P = P, land.type="shrub")
m5.ET.type2 <- function(PET, LAI=NULL, P=NULL,
                        land.type = c("shrub","crop","grass","deciduous","needle","broad","mixed","savanna")){
  land.type <- land.type[1]
  if(land.type=="shrub"){
    -3.11 + 0.39 * PET + 0.09 * P + 11.127 * LAI
  } else if(land.type=="crop"){
    -8.15 + 0.86 * PET + 0.01 * P + 9.54 * LAI
  } else if(land.type=="grass"){
    -1.36 + 0.70 * PET + 0.04 * P + 6.56 * LAI
  } else if(land.type=="deciduous"){
    -14.82+ 0.98 * PET + 2.72 * LAI
  } else if(land.type=="needle"){
    0.10  + 0.64 * PET + 0.04 * P + 3.53 * LAI
  } else if(land.type=="broad"){
    7.71  + 0.74 * PET + 1.85 * LAI
  } else if(land.type=="mixed"){
    -8.763 + 0.95 * PET
  } else if(land.type=="savannas"){
    -5.56 + 0.18 * PET + 0.1 * P + 44.63 * LAI
  } else{
    stop("invalid option 'land.type'")
  }
}

##m6.ET.Zhang--------------------------
#' Calculate ET (Method 6)
#'
#' Calculates evapotranspiration (mm month−1) using the model proposed by Zhang et al. (2001)
#'
#' @param PET potential evapotranspiration (mm month-1) estimated by Priestley-Taylor method (function \code{"m6.PET.PriestleyTaylor"}).
#' @param P predicitation (mm month-1)
#' @param w plant-available water coefficient that represents the relative difference in plant water use for transpiration. Default 2 for forests while 0.5 for grassland when PET is estimated using the Priestley–Taylor model (function "\code{m6.PET.PriestleyTaylor}").
#' @return vector of monthly evapotranspiration (mm month-1)
#' @references
#' Zhang, L., Dawes, W. R., and Walker, G. R. (2001). Response of mean annual evapotranspiration to vegetation changes at catchment scale. Water Resources Research, 37(3), 701–708.
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 104pp
#' @export
#' @examples
#' PET <- c(14.26, 24.65, 64.79, 74.40, 137.02, 153.30, 164.61, 164.30, 111.30, 60.14, 40.50, 10.54)
#' P <- c(38.0, 42.0, 83.0, 73.0, 96.0, 46.0, 115.0, 173.0, 70.0, 85.0, 58.0, 72.0)
#' m6.ET.Zhang(PET, P)
m6.ET.Zhang <- function(PET, P, w=2){
  P * (1 + w * PET/P)/(1 + w * PET/P + P/PET)
}
