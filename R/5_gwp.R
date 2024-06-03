##Script info----------------------
# Created on 2022-06-04 12:19:52 CST
#
# author: Hao-Ran Wu
# College of Life Sciences
# Zhejiang University
# Hangzhou 310012, China
# E-mail: haoranwu@zju.edu.cn
#
# This scripts provides function to model global warming potentials (GWP) caused
#  by either greenhouse gases (GHG) or albedo change.
#   -- 'AGWP.CO2': absolute GWP caused by carbon dioxide
#   -- 'AGWP.CH4': absolute GWP caused by methane
#   -- 'AGWP.N2O': absolute GWP caused by nitrogen oxide
#   -- 'transmittance': upwelling transmittance value
#   -- 'GWP.albedo': realtive GWP caused by albedo change

#AGWP.CO2--------------------
#' AGWP of CO2
#'
#' Calculates absolute global warming potential (W m-2 kg-1 yr) of carbon dioxide for a given period.
#' @param TH chosen time horizon (i.e. 20, 100, 500 years)
#' @param mRF_CO2 marginal radiative forcing of carbon dioxide (W kg-1). It may change
#' between different years, so please keep up the pace of IPCC report.
#' @importFrom stats integrate
#'
#' @return a scalar
#' @export
#' @examples
#' ##1 AGWP CO2 for 100 years
#' AGWP.CO2()
#'
#' ##2 AGWP CO2 for 20 years
#' AGWP.CO2(TH = 20)
AGWP.CO2 <- function(TH=100, mRF_CO2 = 0.908){
  #AGWP_CO2 <- 9.351e-14 # W m-2 kg-1 yr
  #GWP_CO2 <- 1
  S_Earth <- 5.1e+14 #m2
  integrate(
    function(ti) 0.2173+0.224*(exp(-ti/394.4))+0.2824*(exp(-ti/36.54))+0.2763*(exp(-ti/4.304)),
    0,TH
  )$value * (mRF_CO2/S_Earth)
}

#AGWP.CH4------------------
#' AGWP of CH4
#'
#' Calculates absolute global warming potential (W m-2 kg-1 yr) of methane for a given period.
#' @param TH chosen time horizon (i.e. 20, 100, 500 years)
#' @param mRF_CH4 marginal radiative forcing of methane (W kg-1). It may change
#' between different years, so please keep up the pace of IPCC report.
#' @importFrom stats integrate
#'
#' @return a scalar
#' @export
#' @examples
#' ##1 AGWP CH4 for 100 years
#' AGWP.CH4()
#'
#' ##2 AGWP CH4 for 20 years
#' AGWP.CH4(TH = 20)
#'
#' ##3 relative GWP CH4 for 100 years
#' AGWP.CH4()/AGWP.CO2()
AGWP.CH4 <- function(TH=100, mRF_CH4 = 107.3805){
  #AGWP_CH4 <- 2.61e-12  # W m-2 kg-1 yr
  #GWP_CH4 <- 25
  S_Earth <- 5.1e+14 #m2
  integrate(
    function(ti) exp(-ti/12.4),
    0,TH
  )$value * (mRF_CH4/S_Earth)
}

#AGWP.N2O-------------------
#' AGWP of N2O
#'
#' Calculates absolute global warming potential (W m-2 kg-1 yr) of nitrogen oxide for a given period.
#' @param TH chosen time horizon (i.e. 20, 100, 500 years)
#' @param mRF_N2O marginal radiative forcing of nitrogen oxide (W kg-1). It may change
#' between different years, so please keep up the pace of IPCC report.
#' @importFrom stats integrate
#'
#' @return a scalar
#' @export
#' @examples
#' ##1 AGWP N2O for 100 years
#' AGWP.N2O()
#'
#' ##2 AGWP N2O for 20 years
#' AGWP.N2O(TH = 20)
#'
#' ##3 relative GWP N2O for 100 years
#' AGWP.N2O()/AGWP.CO2()
AGWP.N2O <- function(TH=100, mRF_N2O = 185.1244){
  #AGWP_N2O <- 2.43e-11  # W m-2 kg-1 yr
  #GWP_N2O <- 298
  S_Earth <- 5.1e+14 #m2
  integrate(
    function(ti) exp(-ti/121),
    0,TH
  )$value * (mRF_N2O/S_Earth)
}

#transmittance------------------
#' Calculate Upwelling Transmittance Value
#'
#' calculates upwelling transmittance value.
#' @param Sw_down the total incident sunlight (W m-2). A vector or a scalar.
#' @param doy day of the year. It can be calculated by function "\code{DOY}". A vector or a scalar.
#' @param zenith zenith angle. It can be calculated by function "\code{zangle}". A vector or a scalar.
#'
#' @return upwelling transmittance value(s).
#' @export
#' @examples
#' doy <- DOY(2018, 7, 27)
#' zenith <- zangle(2018, 7, 27, 85.44675833, 42.47670833, 9.4)[[2]]
#' transmittance(691.55, doy, zenith)
transmittance <- function(Sw_down, doy, zenith){
  Isc <- 1367 #W m-2 solar constant
  Sw_TOA <- Isc*(1+0.0334*cos(doy*2*pi/365.25))*cos(zenith*pi/180)
  Sw_down/Sw_TOA
}

#GWP.albedo------------------------
#' GWP From Albedo Change
#'
#' Calculates relative global warming potential due to albedo change.
#' @param Sw_up reflected light (W m-2).
#' @param Sw_down total incident sunlight (W m-2).
#' @param alpha_s_REF the albedo value of the reference site.
#' @param Ta upwelling transmittance value. It can be calculated by function \code{transmittance}.
#' @param N the number of days. Default 1.
#' @param TH chosen time horizon (i.e. 20, 100, 500 years). Default 100.
#' @param mRF_CO2 marginal radiative forcing of carbon dioxide (W kg-1). It may change
#' between different years, so please keep up the pace of IPCC report.
#' @param method either "\code{IRF}” or \code{"AGWP"} method.
#'
#' @return a value of relative global warming potential
#' @references
#' Chen, J. (2021). Biophysical Models and Applications in Ecosystem Analysis.
#' Higher Education Press. 138pp
#' @export
#' @examples
#' alpha_REF <- 83/482
#' GWP.albedo(134, 478, alpha_REF, 0.854)
GWP.albedo <- function(Sw_up, Sw_down, alpha_s_REF, Ta,
                       N = 1, TH = 100, mRF_CO2 = 0.908, method="IRF"){
  S_Earth <- 5.1e+14 #m2

  ti <- 1:TH #calculate AF_CO2
  AF_CO2 <- mean(0.2173+0.224*(exp(-ti/394.4))+0.2824*(exp(-ti/36.54))+0.2763*(exp(-ti/4.304)))

  alpha_s <- Sw_up/Sw_down
  delta_alpha <- alpha_s-alpha_s_REF
  mRF_alpha <- -Sw_down*delta_alpha

  if(method=="IRF"){
    mRF_TOA <- -1/N*Ta*mRF_alpha
    mRF_TOA/(AF_CO2*mRF_CO2*TH)
  } else if(method=="AGWP"){
    AGWP_TOA <- -1/N*(Ta*mRF_alpha/S_Earth)
    AGWP_TOA/AGWP.CO2(TH=TH,mRF_CO2=mRF_CO2) #method2
  }
}
