##Script info----------------------
# Created on 2022-06-04 09:17:46 CST
#
# author: Hao-Ran Wu
# College of Life Sciences
# Zhejiang University
# Hangzhou 310012, China
# E-mail: haoranwu@zju.edu.cn
#
# This script provides functions to calculate different sets of input variables
#   or parameters that are required in biophysical models. These variables are:
#   * date:
#       function 'DOY': day of the year.
#       function 'JDo': Julian Day.
#
#   * diurnal changes of air/soil temperature:
#       function 'diurnalT': diurnal temperature
#
#   * VPD and related variables:
#       function 'satvp': saturation vapor pressure
#       function 'satvd': saturation vapor density
#       function 'actualvd': actual vapor density
#       function 'actualvp': actual vapor pressure
#       function 'vpdeficit': vapor pressure deficit
#       function 'dewT': dew point temperature
#       function 'wetbulbT': wet-bulb temperature
#
#   * solar radiation:
#       function 'rn.m1', 'rn.m2': daily net solar radiation
#       function 'ra': daily extraterrestrial radiation
#       function 'rs': daily (incoming) solar radiation
#       function 'rso': daily short wave radiation
#       function 'rns': daily net shortwave radiation
#       function 'rnl': daily net longwave radiation
#
#   * daylight hours and zenith angle of the Sun:
#       function 'dh': maximal daylight hours
#       function 'sunRS': sunrise and sunset time
#       function 'zangle': solar zenith, azimuth, and elevation angle
#
#   * heat flux and storage
#       function 'heatflux': heat flux density of soil
#       function 'heatstore': heat storage in soil
#
#   * wind speed:
#       function 'windspeed': wind speed

##DOY------------------------
#' Calculate DOY
#'
#' calculates day of the year.
#'
#' @param year vector of years.
#' @param month vector of months.
#' @param day vector of days.
#'
#' @return vector of DOY (day of the years).
#' @export
#' @examples
#' ## return a day number of September 22nd in 2016
#' DOY(2016,9,22)
#'
#' ## return multiple day numbers
#' years <- c(2016,2017,2018)
#' months <- c(9,10,2)
#' days <- c(22,13,14)
#' DOY(years, months, days)
DOY <- function(year, month, day){
  floor(275*month/9)-floor((month+9)/12)*(1+floor((year-4*floor(year/4)+2)/3))+day-30
}

##JDo------------------------
#' Calculate Julian Day
#'
#' calculates Julian days of the year.
#'
#' @param year vector of years.
#' @param month vector of months.
#' @param day vector of days.
#'
#' @return vector of Julian days of the year.
#' @export
#' @examples
#' ## return a day number of September 22nd in 2016
#' JDo(2016,9,22)
#'
#' ## return multiple day numbers
#' years <- c(2016,2017,2018)
#' months <- c(9,10,2)
#' days <- c(22,13,14)
#' JDo(years, months, days)
JDo <- function(year, month, day){
  367.0*year-floor(7*(year+floor((month+9)/12))/4)+floor(275*month/9)+day+1721013.5
}

#1.2 Diurnal Changes of Air Temperature and Humidity------------------------
##hourlyT--------------------------
#' Calculate Hourly Air Temperature
#'
#' calculates hourly air temperature based on the algorithm proposed by
#' Parton and Logan (1981).
#'
#' @param hours time in hours (from 1 to 24) used for calculating air temperature
#' @param t_max daily maximum air temperature (Celsius or Fahrenheit)
#' @param t_min daily minimum temperature (Celsius or Fahrenheit)
#' @param time_sr sunrise time (hour)
#' @param time_ss sunset time (hour)
#' @param alpha time lag for maximum air temperature (hour). Default 1.5.
#' @param beta nighttime temperature coefficient. Default 4.
#' @param gamma time lag for minimum air temperature (hour). Default 1.
#' @return a vector of predicted air temperature
#' @references
#' Parton, W. J., and Logan, J. A. (1981). A model for diurnal variation in
#' soil and air temperature. Agricultural Meteorology, 23, 205–216.\cr
#' @export
#' @examples
#' lng <- 85.446758
#' lat <- 42.476708
#' t_max <- 30.2
#' t_min <- 14.3
#' doy <- DOY(2016, 9, 22)
#' hours <- seq(0, 23.5, 0.5)            # half-hour time resolution
#' time_sr <- sunRS(doy, lng, lat)[[1]]  # calculate sunrise time
#' time_ss <- sunRS(doy, lng, lat)[[2]]  # calculate sunset time
#' hourlyT(hours, t_max, t_min, time_sr, time_ss,
#'         alpha = 1.5, beta = 4, gamma = 1)
#'
hourlyT <- function(hours, t_max, t_min, time_sr, time_ss,
                    alpha = 1.5, beta = 4, gamma = 1){
  m <- hours - time_sr - gamma
  dl <- time_ss-time_sr
  Z <- 24-dl

  t_sr <- (t_max-t_min)*sin(pi*(dl-gamma)/(dl+2*alpha))+t_min
  n <- ifelse(hours>time_ss,hours-time_ss,24-time_ss+hours)
  ifelse(
    hours<=time_ss&hours>(time_sr+gamma),
    (t_max-t_min)*sin(pi*m/(dl+2*alpha))+t_min,
    (t_sr-t_min)*exp(-beta*n/Z)+t_min
  )
}

#1.3 Atmosphere Water Vapor Pressure and VPD-----------------
##satvp---------------------------
#' Calculate Saturation Vapor Pressure
#'
#' calculates saturation vapor pressure (kPa). Only valid for temperature above
#' 0◦C.
#' @param Ta vector of air temperature (◦C).
#' @param method type of estimator, either \code{"Tetens"} for model developed by
#' Monteith and Unsworth (2013) or \code{"Lowry"} for model developed by
#' Fritschen and Gay (1979).
#'
#' @return vector of saturation vapor pressure (kPa).
#' @references
#' Monteith, J., and Unsworth, M. (2013). Principles of Environmental Physics:
#' Plants, Animals, and the Atmosphere. Academic Press. 423pp.\cr
#' Lowry, W. P. (2013). Weather and Life: An Introduction to Biometeorology.
#' Elsevier. 305pp.Monteith, J., and Unsworth, M. (2013). Principles of
#' Environmental Physics: Plants, Animals, and the Atmosphere. Academic Press. 423pp\cr
#' Fritschen, L. J., and Gay, L. W. (1979). Environmental Instrumentation.
#' Springer Science & Business Media. 216pp.
#' @examples
#' data("p1_tah")
#' temperature <- p1_tah$Ta
#'
#' satvp(temperature) ##1 Tetens method
#' satvp(temperature, method="Lowry") ##2 Lowry method
#' @export
satvp <- function(Ta, method=c("Tetens","Lowry")){
  if(length(method)==2){
    if(!all(c("Tetens","Lowry")%in%method)) stop("invalid input `method`")
    method<-"Tetens"
  } else if(length(method)==1){
    if(!any(c("Tetens","Lowry")%in%method)) stop("invalid input `method`")
  } else{
    stop("invalid input `method`")
  }

  if(method=="Tetens") return(0.6118*exp(17.502*Ta/(Ta+240.97)))
  if(method=="Lowry") return(
    (6.1078+Ta*(0.44365185+Ta*(0.01428945+Ta*(0.00026506485+Ta*(3.0312404*1e-6+Ta*(2.034809*1e-8+Ta*6.1368209*1e-11))))))/10
  )
}

##satvd---------------------------
#' Calculate Saturation Vapor Density
#'
#' calculates saturation vapor density (kg m-3). Only valid for temperature above
#' 0◦C.
#' @param Ta vector of air temperature (◦C).
#' @param ... additional parameters accepted by function \code{satvp}.
#'
#' @return vector of saturation vapor density (kg m-3).
#' @export
#' @examples
#' data("p1_tah")
#' temperature <- p1_tah$Ta
#' satvd(temperature)
satvd <- function(Ta, ...){
  es <- satvp(Ta, ...)
  es/(4.62*10e-5*(Ta+273.15))
}

##actualvd---------------------------
#' Calculate Actual Vapor Density
#'
#' calculates actual vapor density (kg m-3). Only valid for temperature above
#' 0◦C.
#' @param Ta vector of air temperature (◦C).
#' @param h vector of relative humidity (%).
#' @param ... additional parameters accepted by function \code{satvp}.
#'
#' @return vector of actual vapor density (kg m-3).
#' @export
#' @examples
#' data("p1_tah")
#' temperature <- p1_tah$Ta
#' humidity <- p1_tah$h
#' actualvd(temperature, humidity)
actualvd <- function(Ta, h, ...){
  Es <- satvd(Ta, ...)
  h*Es/100
}

##actualvp--------------------------
#' Calculate Actual Vapor Pressure
#'
#' calculates actual vapor pressure (kPa). Only valid for temperature above
#' 0◦C.
#' @param Ta vector of air temperature (◦C).
#' @param h vector of relative humidity (%).
#' @param ... additional parameters accepted by function \code{satvp}.
#'
#' @return vector of actual vapor pressure (kPa).
#' @export
#' @examples
#' data("p1_tah")
#' temperature <- p1_tah$Ta
#' humidity <- p1_tah$h
#' actualvp(temperature, humidity)
actualvp <- function(Ta, h, ...){
  Ea <- actualvd(Ta, h, ...)
  Ea*(Ta+273.15)/2170
}

##vpdeficit--------------------------
#' Calculate Vapor Pressure Deficit
#'
#' calculates vapor pressure deficit (VPD, kPa). Only valid for temperature above
#' 0◦C.
#' @param Ta vector of air temperature (◦C).
#' @param h vector of relative humidity (%).
#' @param ... additional parameters accepted by function \code{satvp}.
#'
#' @return vector of actual vapor pressure (kPa).
#' @export
#' @examples
#' data("p1_tah")
#' temperature <- p1_tah$Ta
#' humidity <- p1_tah$h
#' vpdeficit(temperature, humidity)
vpdeficit <- function(Ta, h, ...){
  satvp(Ta, ...)-actualvp(Ta, h, ...)
}

##dewT--------------------------
#' calculates Dew Point Temperature
#'
#' calculates dew point temperature (◦C). Only valid for temperature above
#' 0◦C.
#' @param Ta vector of air temperature (◦C).
#' @param h vector of relative humidity (%).
#' @param ... additional parameters accepted by function \code{satvp}.
#'
#' @return vector of dew point temperature (◦C).
#' @references
#' Fritschen, L. J., and Gay, L. W. (1979). Environmental Instrumentation.
#' Springer Science & Business Media. 216pp
#' @export
#' @examples
#' data("p1_tah")
#' temperature <- p1_tah$Ta
#' humidity <- p1_tah$h
#' dewT(temperature, humidity)
dewT <- function(Ta, h, ...){
  ea <- actualvp(Ta, h, ...)
  (237.3*log10(ea/0.61078))/(17.269-log10(ea/0.61078))
}

##wetbulbT--------------------------
#' Calculate Wet-bulb Temperature
#'
#' calculates Wet-bulb temperature (◦C). Only valid for temperature above
#' 0◦C.
#' @param Ta vector of air temperature (◦C).
#' @param h vector of relative humidity (%).
#' @param ... additional parameters accepted by function \code{satvp}.
#' @param Cp specific heat of air (29.3 J mol-1 ◦C-1).
#' @param lambda latent heat vaporization of water (40.660 kJ mol-1)
#'
#' @return vector of dew point temperature (◦C).
#' @references
#' Fritschen, L. J., and Gay, L. W. (1979). Environmental Instrumentation.
#' Springer Science & Business Media. 216pp
#' @export
#' @examples
#' data("p1_tah")
#' temperature <- p1_tah$Ta
#' humidity <- p1_tah$h
#' wetbulbT(temperature, humidity)
wetbulbT <- function(Ta, h, ..., Cp = 29.3, lambda = 40660){
  #λ is the amount of energy (enthalpy) that must be added to transform a given quantity of water into a gas (i.e., vapor)
  G <- Cp/lambda  #0.000666 without calculation (error)
  #thermodynamic psychrometric constant (0.000666 ◦C-1)
  ea <- actualvp(Ta, h, ...)
  Ea <- actualvd(Ta, h, ...)
  es <- satvp(Ta, ...)
  (ea+G*Ea*(Ta+273.15))/(es+G*Ea)
}

#A word of caution is that these equations should only be applied when air
#temperature (Ta, ◦C) is greater than zero.

#The second caution is that these calculations should not be applied
#for longer time scales (e.g., more than 2 – 3 hours) because
#atmospheric vapor density or pressure can be very different


#1.4 Solar Radiation---------------------------------
##rn.m1-----------------
#' Calculate Net Solar Radiation (Method 1)
#'
#' calculates net solar radiation values (MJ m-2 d-1) based on (incoming) and extraterrestrial solar radiation.
#' @param ta mean air temperature (◦C). A vector.
#' @param tmax maximum air temperature (◦C). A vector.
#' @param tmin minimum air temperature (◦C). A vector.
#' @param Rs solar radiation (MJ m-2 d-1).
#' @param Ra extraterrestrial solar radiation (MJ m-2 d-1).
#'
#' @return vector of net solar radiation (MJ m-2 d-1).
#' @export
#' @references
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 102pp
#' @examples
#' lati <- -22.9
#' doy <- DOY(2022, 5, 15)
#' ah <- 7.1
#' Rs <- rs(lati, doy, ah)
#' Ra <- ra(lati, doy)
#' rn.m1(25, 30, 20, Rs, 120)
rn.m1 <- function(ta, tmax, tmin, Rs, Ra){
  f <- 1.2*Rs/Ra+0.1
  0.77*Rs-2.45*10^(-9)*f*(0.261*exp(-7.7710*10^(-4)*(ta+273.15)^2)-0.02)*((tmax+273.15)^4+(tmin+273.15)^4)+0.83
}

##ra-----------------
#' Calculate Extraterrestrial Radiation
#'
#' calculates daily extraterrestrial radiation (MJ m-2 d-1).
#' @param lati latitude of the site in degree.
#' @param doy day of the year. It can be calculated by the function "\code{DOY}".
#' @return vector of net solar radiation (MJ m-2 d-1).
#' @export
#' @references
#' Allen, R. G., Pereira, L. S., Raes, D. and Smith, M. (1998). Crop Evapotranspiration – Guidelines for Computing Crop Water Requirements, FAO Irrigation and Drainge Paper 56, FAO, 1998, ISBN 92-5-104219-5.
#' @examples
#' doy <- DOY(2022,9,3)
#' lati <- -20
#' ra(lati, doy)
ra <- function(lati, doy){
  lati <- lati*pi/180
  gsc <- 0.082 #solar constant (MJ m-2 min-1)
  dr <- 1 + 0.033*cos(2*pi*doy/365)
  delta <- 0.409*sin(2*pi*doy/365-1.39)
  omega <- acos(-tan(lati) * tan(delta))
  24*60/pi * gsc * dr *(omega*sin(lati)*sin(delta)+
                          cos(lati)*cos(delta)*sin(omega))
}

##dh------------------------
#' Calculate Daylight Hours
#'
#' calculates daylight hours. It is the maximum possible during sunshine or daylight hours (h).
#' @param lati latitude of the site in degree.
#' @param doy day of the year. It can be calculated by the function "\code{DOY}".
#' @return value(s) of daylight hours.
#' @export
#' @references
#' Allen, R. G., Pereira, L. S., Raes, D. and Smith, M. (1998). Crop Evapotranspiration – Guidelines for Computing Crop Water Requirements, FAO Irrigation and Drainge Paper 56, FAO, 1998, ISBN 92-5-104219-5.
#' @examples
#' doy <- DOY(2022,9,3)
#' lati <- -20
#' dh(lati, doy)
dh <- function(lati, doy){
  lati <- lati*pi/180
  delta <- 0.409*sin(2*pi*doy/365-1.39)
  omega <- acos(-tan(lati) * tan(delta))
  24/pi*omega
}

##rs---------------------
#' Calculate Solar Radiation
#'
#' calculates (incoming) solar radiation when it is not measured.
#' @param lati latitude of the site in degree.
#' @param doy day of the year. It can be calculated by the function "\code{DOY}".
#' @param ah actual duration of sunshine in a day (h).
#' @param as regression constant, expressing the fraction of extraterrestrial radiation reaching the earth on overcast days. Default 0.25.
#' @param bs fraction of extraterrestrial radiation reaching the earth on clear days. Default 0.5.
#' @return value(s) of incoming solar radiation (or shortwave radiation, MJ m-2 d-1).
#' @references
#' Allen, R. G., Pereira, L. S., Raes, D. and Smith, M. (1998). Crop Evapotranspiration – Guidelines for Computing Crop Water Requirements, FAO Irrigation and Drainge Paper 56, FAO, 1998, ISBN 92-5-104219-5.
#' @export
#' @examples
#' rs(-22.9, DOY(2022, 5, 15), 7.1)
rs <- function(lati, doy, ah, as = 0.25, bs = 0.5){
  Ra <- ra(lati, doy)
  N <- dh(lati, doy)
  (as+bs*ah/N)*Ra
}

##rso----------------
#' Calculate Short Wave Radiation
#'
#' calculates short wave radiation on a clear-sky day (MJ m-2 d-1).
#' @param lati latitude of the site in degree.
#' @param doy day of the year. It can be calculated by the function "\code{DOY}".
#' @param tau clear sky transmissivity. Default 0.75.
#' @param EL station elevation (m). Default 0.
#' @param method method to calculate short wave radiation. Valid options are:\cr
#' "\code{tau}": use clear sky transmissivity.\cr
#' "\code{EL}": use elevation.\cr
#'
#' @return value(s) of short wave radiation (MJ m-2 d-1).
#' @export
#' @references
#' Allen, R. G., Pereira, L. S., Raes, D. and Smith, M. (1998). Crop Evapotranspiration – Guidelines for Computing Crop Water Requirements, FAO Irrigation and Drainge Paper 56, FAO, 1998, ISBN 92-5-104219-5.
#' @examples
#' rso(-22.9, DOY(2022, 5, 15))
rso <- function(lati, doy, tau=0.75, EL=0,
                method = "tau"){
  Ra <- ra(lati, doy)
  if(method=="tau"){
    tau * Ra
  } else if(method=="EL"){
    (0.75 + 2e-5*EL)*Ra
  } else{
    stop("invalid option 'method'")
  }
}

##rns------------------
#' Calculate Net Shortwave Radiation
#'
#' calculates net shortwave radiation (MJ m-2 day-1)
#' @param Rs (incoming) solar radiation (MJ m-2 d-1). It can be calculated by the function "\code{rs}" if it is not measured.
#' @param albedo albedo or canopy reflection coefficient, which is 0.23 for the hypothetical grass reference crop.
#'
#' @return values of net shortwave radiatin (MJ m-2 day-1).
#' @export
#' @references
#' Allen, R. G., Pereira, L. S., Raes, D. and Smith, M. (1998). Crop Evapotranspiration – Guidelines for Computing Crop Water Requirements, FAO Irrigation and Drainge Paper 56, FAO, 1998, ISBN 92-5-104219-5.
#' @examples
#' Rs <- rs(-22.9, DOY(2022, 5, 15), 7.1)
#' rns(Rs)
rns <- function(Rs, albedo = 0.23){
  (1-albedo) * Rs
}

##rnl---------------------------------
#' Calculate Net Longwave Radiation
#'
#' calculates net longwave radiation (MJ m-2 d-1).
#' @param Rs (incoming) solar radiation (MJ m-2 d-1). It can be calculated by the function "\code{rs}" if it is not measured.
#' @param Rso short wave radiation on a clear-sky day (MJ m-2 d-1). It can be calculated by the function "\code{rso}".
#' @param ea actual vapor pressure (kPa). It can be calculated by the function "\code{actualvp}".
#' @param tmax maximum absolute temperature during a day (◦C).
#' @param tmin minimum absolute temperature during a day (◦C).
#'
#' @return value(s) of net longwave radiation (MJ m-2 d-1).
#' @export
#' @references
#' Allen, R. G., Pereira, L. S., Raes, D. and Smith, M. (1998). Crop Evapotranspiration – Guidelines for Computing Crop Water Requirements, FAO Irrigation and Drainge Paper 56, FAO, 1998, ISBN 92-5-104219-5.
#' @examples
#' lati <- -22.7
#' doy <- DOY(2022, 5, 15)
#' ah <- 7.1
#' Rs <- rs(lati, doy, ah)
#' Rso <- rso(lati, doy)
#' rnl(Rs, Rso, 2.1, 25.1, 19.1)
rnl <- function(Rs, Rso, ea, tmax, tmin){
  4.903e-09 * (((tmax + 273.16)^4 + (tmin + 273.16)^4)/2) *
    (0.34 - 0.14 * sqrt(ea)) * (1.35 * Rs/Rso - 0.35)
}

##rn.m2------------------
#' Calculate Net Radiation (Method 2)
#'
#' calculates net solar radiation (MJ m-2 day-1) based on net shortwave and longwave radiation.
#' @param Rns net shortwave radiation (MJ m-2 day-1). It can be calculated by the function "\code{rns}".
#' @param Rnl net longwave radiation (MJ m-2 day-1). It can be calculated by the function "\code{rnl}".
#'
#' @return values of net radiation (MJ m-2 day-1)
#' @export
#' @references
#' Allen, R. G., Pereira, L. S., Raes, D. and Smith, M. (1998). Crop Evapotranspiration – Guidelines for Computing Crop Water Requirements, FAO Irrigation and Drainge Paper 56, FAO, 1998, ISBN 92-5-104219-5.
#' @examples
#' lati <- -22.7
#' doy <- DOY(2022, 5, 15)
#' ah <- 7.1
#' Rs <- rs(lati, doy, ah)
#' Rso <- rso(lati, doy)
#' Rns <- rns(Rs)
#' Rnl <- rnl(Rs, Rso, 2.1, 25.1, 19.1)
#' rn.m2(Rns, Rnl)
rn.m2 <- function(Rns, Rnl){
  Rns-Rnl
}

#https://www.fao.org/3/X0490E/x0490e07.htm#chapter%203%20%20%20meteorological%20data

##sunRS--------------------------
#' Calculate Sun Rise/Set Time
#'
#' calculates sunrise and sunset time (h).
#'
#' @param doy day of the year. A scalar.
#' @param long longtitude of the site (-180 to 180 degree). A scalar.
#' @param lat latitude of the site (-90 to 90 degree). A scalar.
#'
#' @return vector of two values. The first one is sunrise time and the second one
#' is sunset time.
#' @references
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 33pp
#' @export
#' @examples
#' doy <- DOY(2016, 9, 22)
#' sunRS(doy, 85.44675833, 42.47670833)
sunRS <- function(doy, long, lat){
  if(length(doy)!=1 || length(long)!=1 || length(lat)!=1) stop("All parameters should be scalars.")
  N<-doy
  #N -- day of the year
  t1<-N+(6-long/15)/24 #for sun rise
  t2<-N+(18-long/15)/24#/*for sun set*/
  #/*compute sunrise time*/
  rise <- function(t){
    result <- tryCatch({
      K<-pi/180.0
      M<-0.9856*t-3.289
      L<-M+1.916*sin(M*K)+0.02*sin(2*K*M)+282.634
      if(L >= 360.0 || L >= 360) L<-L-360
      tanRA<-0.91746*tan(L*pi/180)
      RA<-atan(tanRA)/K
      if((L > 90.0) && (L <= 270.0)) RA<-180.0+RA
      if ((L > 270.0) && (L <= 360.0)) RA<-360.0+RA
      RA<-RA/15.0
      sind<-0.39782*sin(L*K)
      cosd<-cos(asin(sind))
      x<-(-0.01454-sind*sin(lat*K))/(cosd*cos(lat*K))
      H<-acos(x)/K

      H<-360.0-H
      Tsr<-(H/15.0)+RA-(0.06571*t)-6.622

      if(Tsr < 0.0) Tsr=Tsr+24.0
      if(Tsr > 24.0) Tsr=Tsr-24.0
      return(Tsr)
    },error = function(e) print("You could probably try to force your returned arccosine angle not within [-1,1]"))
  }
  #/*compute sunset time*/
  sset <- function(t){
    result <- tryCatch({
      K<-pi/180.0
      M<-0.9856*t-3.289
      L<-M+1.916*sin(M*K)+0.02*sin(2*M*K)+282.634

      if(L>=360.0) L<-L-360
      tanRA<-0.91746*tan(L*K)
      RA<-atan(tanRA)/K
      if((L > 90.0) && (L <= 270.0)) RA<-180.0+RA
      if((L > 270.0) && (L <= 360.0)) RA<-360.0+RA
      RA<-RA/15.0

      sind<-0.39782*sin(L*K)
      cosd<-cos(asin(sind))
      x<-(-0.01454-sind*sin(lat*K))/(cosd*cos(lat*K))
      H<-acos(x)/K

      Tst<-H/15.0+RA-0.06571*t-6.622
      if(Tst < 0.0) Tst=Tst+24.0
      return(Tst)
    },error = function(e) print("You could probably try to force your returned arccosine angle not within [-1,1]")
    )
  }
  c(rise=rise(t1),set=sset(t2))
}

##zangle--------------------------
#' Calculate Solar Zenith/Azimuth Angle
#'
#' calculates solar zenith, azimuth, and elevation angle.
#'
#' @param year a scalar.
#' @param month month of year (1-12). A scalar.
#' @param day day of the month (1-30 or 31). A scalar.
#' @param long longtitude of the site (-180 to 180 degree). A scalar.
#' @param lat latitude of the site (-90 to 90 degree). A scalar.
#' @param time_series vector of time in a day. Values within 0-24 (hour).
#'
#' @return a \code{"data.frame"} of four columns. The first column is time of a
#' day. The last three columns show solar zenith, azimuth, and solar elevation angles, respectively. Unit in degrees.
#' @references
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 33pp
#' @export
#' @examples
#' zangle(2016, 9, 22, 85.44675833, 42.47670833)
zangle <- function(year, month, day, long, lat, time_series=0:23){
  if(length(year)!=1) stop("'year' should be a scalar.")
  if(length(month)!=1) stop("'month' should be a scalar.")
  if(length(day)!=1) stop("'day' should be a scalar.")
  if(length(long)!=1) stop("'long' should be a scalar.")
  if(length(lat)!=1) stop("'lat' should be a scalar.")
  if(any(time_series>24) || any(time_series<0)) stop("'time_series' out of 0-24.")

  K<-pi/180.0
  RStime <- sunRS(DOY(year, month, day), long, lat)
  Avect <- c()
  Zvect <- c()
  Tvect <- c()
  #0:(ceiling(24/time_interval)-1)
  for(time in time_series){
    #time<-i/(1/time_interval)
    UT<-time+long/15.0 #	/*the Universal Time*/
    UT<-UT-floor(UT/24.0)*24
    JDo<-367.0*year-floor(7*(year+floor((month+9)/12))/4)+floor(275*month/9)+day+1721013.5
    #/* Julian day of the year*/
    JD<-JDo+UT/24
    Time<-(JD-2451545.0)/36525.0
    M2<-357.528+35999.050*Time
    L<-280.46+36000.772*Time+(1.915-0.0048*Time)*sin(M2*K)+0.02*sin(2*M2*K)
    #/*the true geocentric longitude of the sun*/
    L<-L-floor(L/360.0)*360
    if(L<0.0) L<-360.0+L
    Sigma<-23.439-0.013*Time
    tanRA<-cos(Sigma*K)*(sin(K*L)/cos(K*L))
    RA<-atan(tanRA)/K

    if((L >90) && (L <= 270)) RA<-180+RA
    if((L > 270) && (L <= 360)) RA<-360+RA
    RA<-RA/15.0
    sind<-0.39782*sin(L*K) #/*d -- declination of the sun*/
    cosd<-cos(asin(sind))

    N<-floor(275*month/9)-floor((month+9)/12)*(1+floor((year-4*floor(year/4)+2)/3))+day-30
    GMST<-6.6106172+0.0657098242*N+1.00273791*UT #	/*the Greenwich mean sidereal time in hours*/
    Omega<-125.04452-1934.13626*Time+0.002071*Time*Time
    E<- -0.00029*sin(Omega*K)
    GAST<-GMST+E
    GHA<-15*(GAST-RA)
    LHA<-GHA-long
    LHA<-LHA-floor(LHA/360.0)*360.0

    cosz<-sin(lat*K)*sind+cos(K*lat)*cosd*cos(K*LHA)
    Z<-(acos(cosz)/K)
    tana<-sin(LHA*K)/(cos(LHA*K)*sin(lat*K)-sind/(cosd*cos(lat*K)))

    if((LHA >= 0.0) && (LHA<=180.0)){
      if(tana>=0.0){
        A<-180.0+atan(tana)/K
      } else{
        A<-360.0+atan(tana)/K
      }
    } else{
      if(tana>=0.0){
        A<-atan(tana)/K
      } else{
        A<-180.0+atan(tana)/K
      }
    }

    if((time >= as.numeric(RStime[1])) && (time <= as.numeric(RStime[2]))){
      Tvect <- c(Tvect,time)
      Avect <- c(Avect,A)
      Zvect <- c(Zvect,Z)
    }
  }

  if(length(Tvect)==0) warning("time series not within the sunrise and sunset time, so no results left.")

  data.frame(time=Tvect,zenith=Zvect,azimuth=Avect,elevation=90-Zvect)
}

#1.5  Heat Storages in Soil, Air and Vegetation-------------------------
##heatflux--------------------------
#' Calculate Heat Flux Density
#'
#' calculates heat flux density (W m-2) of soil.
#'
#' @param detT temperature difference (◦C/K) between the top and the bottom of the soil layer.
#' @param d thickness of the soil layer (m).
#' @param kappa the thermal conductivity of the soil (W m−1 K−1). It is ~ 2.5 (W m−1 K−1)
#' for soil minerals, ~1.92 (W m−1 K−1) for organic matter, and 4.18 (W m−1 K−1) for water.
#' @return a vector of heat flux density.
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 37pp
#' @export
#' @examples
#' ## soil, delta T = 2, d = 1m
#' heatflux(2, 1)
#'
#' ## soil minerals, delta T = 2, d = 1m
#' heatflux(2, 1, kappa = 2.5)
#'
#' ## organic matter, delta T = 2, d = 1m
#' heatflux(2, 1, kappa = 1.92)
#'
#' ## water, delta T = 2, d = 1m
#' heatflux(2, 1, kappa = 4.18)
heatflux <- function(detT, d, kappa=2.5){
  kappa*(detT)/d
}

##heatstore--------------------------
#' Calculate Heat Storage
#'
#' calculates the heat storage (W m-3) in soil over a period of time.
#'
#' @param detT temperature difference (◦C/K) over a period of time.
#' @param dettime period of time (s).
#' @param theta the volumetric soil water content (%).
#' @param roub the soil bulk density (kg m−1).
#' @return a vector of heat storage.
#' @references
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 37pp
#' @export
#' @examples
#' ## delta T = 3 Celsius/Kelvin, 6 hours, and soil water content of 50%
#' heatstore(3, 3600*6, 50)
heatstore <- function(detT, dettime, theta, roub=2.7){
  cd<-890     #specific heat capacities of the dry mineral soil (J kg−1 K−1)
  cw<-4190    #specific heat capacities of the soil water (J kg−1 K−1)
  rouw<-1000  #the density of water (kg m-3)
  (roub*cd+theta/100*rouw*cw)*(detT/dettime)
}

#1.6  Vertical Profile of Wind Speed-----------------------------
##ws--------------------------
#' Calculate Vertical Profile of Wind Speed
#'
#' calculates vertical profile of wind speed, or fit the model of wind profile.
#' Note that this algorithm is valid only for the neutral atmosphere and should
#' not be applied for tall canopies or vegetation with large variation across
#' the horizontal space.
#'
#' @param z specified height (above the ground level) of the wind speed measurement (m). A vector.
#' @param uz measured wind speed at height \code{z}. A vector of the same length as parameter \code{z}.
#' @param z0 the roughness length (m) at which wind speed is near zero.
#' @param d zero plain displacement (m).
#' @param kappa the von Karmon constant with an average value of 0.35–0.43.
#' A value of 0.40 is often used in the literature.
#' @param ustar the friction velocity (m s-1) which depends on the shear stress
#' (kg m-1 s-2) at the boundary of the flow and air density (kg m-3).
#' @param windu.dmax maximal value of parameter \code{d} to be considered when fiting the
#' wind profile model.
#' @param mode three options.\cr
#' "\code{wind2}" convert wind speed measured at a certain
#' height to the wind speed at 2meters, the same as function "\code{wind2}" in package
#' \code{sirad}. Parameter \code{kappa}, \code{ustar}, and \code{windu.dmax} will
#' not be used.\cr
#' "\code{windz}" calculates wind speed at any height based on wind profile model.
#' Parameter \code{windu.dmax} will not be used.\cr
#' "\code{windu}" fits the wind profile model using measured wind speed data at
#' certain heights. Parameter \code{z0}, \code{d}, and \code{ustar} will not be used
#' and they will be estimated by regression. Models with zero plain displacement over
#' than \code{windu.dmax} will not be considered.
#'
#' @return When \code{mode = "wind2"}, it returns a vector of wind speed (m s-1) at
#' 2 meters. When \code{mode = "windz"}, it returns a vector of vertical profile of
#' wind speed (m s-1) for different height (indicated by the name of the vector).
#' When \code{mode = "windu"}, it returns a list with four components: formula,
#' estimated parameters, R square, and adjusted R square value of the regression.
#' @references
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 39pp
#' @export
#' @examples
#' ## calculate wind profile
#' z <- seq(1,10,0.1)
#' wp <- windspeed(z, mode = "windz")
#' plot(z,wp)
#
#' ## calculate wind speed at 2m based on measured speed at 9m
#' windspeed(z=9, uz=1.601216)
#'
#' ## estimate parameters based on measured wind speeds
#' md <- windspeed(z, wp, mode = "windu")
#' md$coef
ws <- function(z, uz = NULL, z0=0.01474926, d=0.07994099,
                  kappa = 0.4, ustar=0.1, windu.dmax = 2,
                  mode = "wind2"){
  if(mode=="wind2"){
    if(is.null(uz)) stop("parameter 'uz' not found")
    if(length(z)!=length(uz)) stop("unequal length between 'z' and 'uz'")

    uz * log((2-d)/z0) / log(1/z0 * z - d/z0)
  } else if(mode=="windz"){
    zad <- ifelse((z-d)/z0>=1,(z-d)/z0,1)
    prof <- ustar/kappa*log(zad)
    names(prof) <- z
    prof
  } else if(mode=="windu"){
    if(is.null(uz)) stop("parameter 'uz' not found")
    d_seq <- seq(0,windu.dmax,0.01)
    r2_seq <- do.call(c,lapply(d_seq,function(cand_d){
      uz.po <- uz[z-cand_d>0]
      z.po <- z[z-cand_d>0]
      m <- lm(uz.po~I(log(z.po-cand_d)))
      summary(m)$r.squared
    }))
    d <- d_seq[which(max(r2_seq)==r2_seq)]
    d_seq <- seq(d-0.01,d+0.01,0.0001)
    r2_seq <- do.call(c,lapply(d_seq,function(cand_d){
      uz.po <- uz[z-cand_d>0]
      z.po <- z[z-cand_d>0]
      m <- lm(uz.po~I(log(z.po-cand_d)))
      summary(m)$r.squared
    }))
    d <- d_seq[which(max(r2_seq)==r2_seq)]
    uz.po <- uz[z-d>0]
    z.po <- z[z-d>0]
    md <- lm(uz.po~I(log(z.po-d)))
    md <- lm(uz~I(log(z-d)))
    u2k <- md$coefficients[[2]]
    z0 <- exp(-md$coefficients[[1]]/u2k)
    ustar <- u2k * kappa
    list(
      formula = "uz = (ustar/kappa) * log((z-d)/z0)",
      coef = c(ustar = ustar, kappa = kappa, d = d, z0 = z0),
      r.square = summary(md)$r.square,
      adj.r.square = summary(md)$adj.r.square
    )
  }
}
