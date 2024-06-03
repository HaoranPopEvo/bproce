##Script info----------------------
# Created on 2022-06-04 11:14:04 CST
#
# author: Hao-Ran Wu
# College of Life Sciences
# Zhejiang University
# Hangzhou 310012, China
# E-mail: haoranwu@zju.edu.cn
#
# This script provides functions to model ecosystem production.
#   -- 'photo.michment': calculate ecosystem production based on Michaelis-Menten model
#   -- 'photo.michment.fit': fit Michaelis-Menten model.
#   -- 'photo.landsberg': calculate ecosystem production based on Landsberg model
#   -- 'photo.landsberg.fit': fit Landsberg model
#   -- 'photo.farquhar': calculate ecosystem production based on Farquhar model
#   -- 'photo.ballberry': model stomatal conductance based on Ball-Berry model
#   -- 'photo.ballberry.fit': fit Ball-Berry model

#2.2 Core Biophysical Models for Ecosystem Production
##photo.michment--------------------------
#' Calculate Photosynthesis Rate Based On Michaelis-Menten Model
#'
#' calculates the photosynthesis rate (Pn, μmol m-2 s-1) based on Michaelis-Menten model.
#'
#' @param PAR photosynthetically active radiation (μmol m-2 s-1).
#' @param alpha the photochemical efficiency of photosynthesis at low light.
#' @param Pm the maximum photosynthetic capacity of a leaf or an ecosystem (μmol m-2 s-1).
#' @param Rd respiration rate (μmol m-2 s-1).
#' @param beta additional shape factor.
#' @param method type of calculating methods. Option "\code{mechaelis}" adopted the
#' original Michaelis-Menten Model without using parameter \code{beta}. Option "\code{landsand}
#' adopted the non-rectangular hyperbolic model proposed by Landsberg and Sands (2011). Option "\code{peat}"
#' used the model applied by Peat (1970), another form of non-rectangular hyperbolic model.
#' @return a vector of photosynthesis rate (Pn, μmol m-2 s-1).
#' @references
#' Michaelis, L., and Menten, M. L. (1913). Die Kinetik der Invertinwirkung. Biochem Z, 49: 333–369.\cr
#' Landsberg, J. J., and Sands, P. (2011). Physiological Ecology of Forest Production:
#' Principles, Processes and Models (Vol. 4). London: Elsevier/Academic Press. 352pp.\cr
#' @export
#' @examples
#' ## Michaelis-Menten model
#' photo.michment(seq(10,100,10))
#'
#' ## Non-rectangular hyperbolic model
#' photo.michment(seq(10,100,10), method = "landsand", beta = -2)
#' photo.michment(seq(10,100,10), method = "peat", beta = -2)
photo.michment <- function(PAR, alpha=0.05, Pm=10, Rd=0, beta=NULL,
                           method=c("michaelis","landsand","peat")){
  if(length(method)==3&&all(c("michaelis","landsand","peat")%in%method)){
    method<-"michaelis"
  }
  if(length(method)!=1) stop("invalid `method`")
  if(!any(c("michaelis","landsand","peat")%in%method)) stop("invalid `method`")

  if(method=="michaelis") return(alpha*PAR*Pm/(alpha*PAR+Pm)-Rd)
  if(method=="landsand"){
    if(is.null(beta)) stop("parameter `beta` not specified.")
    return(2*alpha*PAR/(1+alpha*PAR/Pm+sqrt((1+alpha*PAR/Pm)^2-4*alpha*beta*PAR/Pm))-Rd)
  }
  if(method=="peat"){
    if(is.null(beta)) stop("parameter `beta` not specified.")
    return((1/(2*beta))*(alpha*PAR+Pm-sqrt((alpha*PAR+Pm)^2-4*alpha*PAR*Pm*beta))-Rd)
  }
}

##photo.michment.fit-----------------------
#' Fit Michaelis-Menten Models For Photosynthesis
#'
#' \code{photo.michment.fit} is used to fit Michaelis-Menten photosynthesis models. It
#' uses linear regression after linearization of the Michaelis-Menten equation.
#' @param PAR photosynthetically active radiation (μmol m-2 s-1).
#' @param Pn measured photosynthesis rate (μmol m-2 s-1).
#' @param Rd potential range of respiration rate (μmol m-2 s-1). This could be a vector, and each
#' of its elements will be used to fit Michaelis-Menten model, although only the model with maximal
#' R-square will be chosen.
#'
#' @return a list of four objects: \cr\cr
#' \code{$formula}: formula of the model.\cr
#' \code{$coef}: fitted coefficients of the model.\cr
#' \code{$r.squared}: R square of the model.\cr
#' \code{$adj.r.squared}: adjusted R square of the model.\cr
#' @export
#' @examples
#' PAR <- seq(10,100,10)
#' Pn <- photo.michment(PAR, alpha = 0.05, Pm=10, Rd=0) + jitter(0.01) # add small noises
#' photo.michment.fit(PAR, Pn)
photo.michment.fit <- function(PAR, Pn, Rd=seq(0,max(Pn),0.01)){
  r.square <- do.call(c,lapply(Rd, function(rd){
    y <- PAR/(Pn+rd)
    md <- lm(y~PAR)
    summary(md)$r.squared
  }))
  index <- which(max(r.square)==r.square)
  if(length(index)==1){
    if(index==1 && Rd[index][1]!=0) warning("it is recommended to reduce the lower boundary of vector 'Rd'")
    if(index==length(Rd)) warning("it is recommended to increase the upper boundary of vector 'Rd'")
  }
  Rd <- Rd[index][1]
  y <- PAR/(Pn+Rd)
  md <- lm(y~PAR)
  list(
    formula = "Pn = alpha * PAR * Pm / (alpha * PAR + Pm)",
    coef = c(alpha=1/md$coefficients[[1]],Pm=1/md$coefficients[[2]],Rd=Rd),
    r.squared = summary(md)$r.squared,
    adj.r.squared = summary(md)$adj.r.squared
  )
}

##photo.landsberg--------------------------
#' Calculate Photosynthesis Rate Based On Landsberg Model
#'
#' calculates the photosynthesis rate (Pn, μmol m-2 s-1) based on Landsberg model.
#'
#' @param PAR photosynthetically active radiation (μmol m-2 s-1).
#' @param alpha the slope of the change in photosynthesis rate with PAR.
#' @param Pm the maximum photosynthetic capacity of a leaf or an ecosystem (μmol m-2 s-1).
#' @param Icomp the light compensation point (μmol m-2 s-1) at which the photosynthesis rate is zero.
#' @return a vector of photosynthesis rate (Pn, μmol m-2 s-1).
#' @references
#' Landsberg, J. J. (1977). Some useful equations for biological studies. Experimental Agriculture,
#' 13(3), 273–286.
#' @export
#' @examples
#' photo.landsberg(seq(200,300,10))
photo.landsberg <- function(PAR, alpha=0.008, Pm=10, Icomp=200){
  Pm*(1-exp(-alpha*(PAR-Icomp)))
}

##photo.landsberg.fit-----------------------
#' Fit Landsberg Model For Photosynthesis
#'
#' \code{photo.landsberg.fit} is used to fit Landsberg photosynthesis model It
#' uses linear regression after linearization of the Landsberg equation.
#' @param PAR photosynthetically active radiation (μmol m-2 s-1).
#' @param Pn measured photosynthesis rate (μmol m-2 s-1).
#' @param Pm potential range of the maximum photosynthetic capacity of a leaf or an
#' ecosystem (μmol m-2 s-1). This could be a vector, and each of its elements will
#' be used to fit Landsberg model, although only the model with maximal
#' R-square will be chosen.
#'
#' @return a list of four objects: \cr\cr
#' \code{$formula}: formula of the model.\cr
#' \code{$coef}: fitted coefficients of the model.\cr
#' \code{$r.squared}: R square of the model.\cr
#' \code{$adj.r.squared}: adjusted R square of the model.\cr
#' @export
#' @examples
#' PAR <- seq(200,300,10)
#' Pn <- photo.landsberg(PAR, alpha = 0.008, Pm = 10, Icomp = 200) + jitter(0.01) # add small noises
#' photo.landsberg.fit(PAR, Pn)
photo.landsberg.fit <- function(PAR, Pn, Pm=seq(max(Pn)/2,max(Pn)*2,0.01)){
  r.square <- do.call(c,lapply(Pm, function(pm){
    valid_index <- 1-Pn/pm>0
    Pn <- Pn[valid_index]
    PAR <- PAR[valid_index]
    y <- log(1-Pn/pm)
    md <- lm(y~PAR)
    summary(md)$r.squared
  }))
  index <- which(max(r.square)==r.square)
  if(length(index)==1){
    if(index==1) warning("it is recommended to reduce the lower boundary of vector 'Pm'")
    if(index==length(Pm)) warning("it is recommended to increase the upper boundary of vector 'Pm'")
  }
  Pm <- Pm[index]
  valid_index <- 1-Pn/Pm>0
  Pn <- Pn[valid_index]
  PAR <- PAR[valid_index]
  y <- log(1-Pn/Pm)
  md <- lm(y~PAR)

  SSE <- sum((Pn-Pm*(1-exp(predict(md))))^2)
  SST <- sum((Pn-mean(Pn))^2)

  list(
    formula = "Pn = Pm * (1-exp(alpha*(PAR-Icomp)))",
    coef = c(Pm=Pm,Icomp=-md$coefficients[[1]]/md$coefficients[[2]],alpha=-md$coefficients[[2]]),
    r.square = 1-SSE/SST,
    adj.r.square = 1-(SSE/SST)*(length(PAR)-1)/(length(PAR)-2)
  )
}

##photo.farquhar--------------------------
#' Calculate Net Leaf CO2 Assimilation Based On Farquhar Model
#'
#' calculates net leaf CO2 assimilation rate (An, μmol m-2 s-1) based on Farquhar's model.
#'
#' @param ci the intercellular CO2 concentration (µmol mol−1).
#' @param Vmax the maximum activity of Rubisco (µmol m−2 s−1).
#' @param Jmax the electron transport rate (µmol m−2 s−1) which can vary with absorbed
#' photosynthetically active radiation (aPAR).
#' @param Q10 Q10 value of leaf respiration modeled with leaf temperature.
#' @param Tp the triose phosphate utilization rate (µmol m−2).
#' @param oi the oxygen (O2) concentration in the atmosphere (209 mol mol−1).
#' @param Rd respiration rate (μmol m-2 s-1).
#' @return a vector of net leaf CO2 assimilation rate (An, μmol m-2 s-1).
#' @references
#' Farquhar, G. D., von Caemmerer, S. V., and Berry, J. A. (1980). A biochemical model
#' of photosynthetic CO2 assimilation in leaves of C3 species. Planta, 149(1), 78–90.
#' @export
#' @examples
#' photo.farquhar(seq(10000,100000,1000), 20, 200, 3.18) ##NOTE: not validated by datasets
photo.farquhar <- function(ci, Vmax, Jmax, Q10, Tp=Inf, oi=209e+6, Rd=0){
  ### Q10 <- exp(10*0.1159) #NOTE!! soil respiration, not leaf   ##NOT VALIDATED!!
  tau <- 0.5*oi/(2600*0.57*Q10)
  Kc <- 30*2.1*Q10
  Ko <- 30000*1.2^Q10
  Ac <- Vmax*(ci-tau)/(ci+Kc*(1+oi/Ko))
  Aj <- Jmax*(ci-tau)/(4*ci+8*tau)
  Ap <- rep(3*Tp,length(Ac))

  do.call(c,lapply(1:length(Ac),function(ii) min(Ac[ii],Aj[ii],Ap[ii])))-Rd
}


##photo.ballberry--------------------
#' Calculate Stomatal Conductance Based On Ball-Berry Model
#'
#' \code{photo.ballberry} uses Ball-Berry model to calculate stomatal conductance
#' (µmol m−2 s-1) based on a set of variables.
#' @param data a \code{data.frame} with several columns. It contains required data
#' to compute stomatal conductance for different models. The column names must be:\cr
#' "\code{An}": for photosynthesis rate (µmol m−2 s-1)\cr
#' "\code{cs}": for CO2 concentration of leaf surface (µmol mol-1)\cr
#' "\code{hs}": for relative humidity at the leaf surface (a value between 0-1,
#' rather than percentage)\cr
#' "\code{VPD}": vapor pressure deficit (kPa).
#' @param g0 the minimum conductance value (mol m−2 s−1). Parameter.
#' @param g1 parameter.
#' @param tau parameter.
#' @param VPD_0 the value of VPD at which stomatal conductance becomes zero (kPa). Parameter.
#' @param model type of model used.\cr
#' "\code{bb}": the original Ball-Berry model.
#' "\code{launing}": model proposed by Leuning in 1990 by adding parameter \code{tau}
#' into the original Ball-Berry model.
#' "\code{c3}": model proposed by Leuning in 1995. A VPD term is used for C3 plants.
#' "\code{vpd}": model proposed by Leuning in 1995. VPD, instead of relative humidity,
#' is used in the model.
#' @return vector of stomatal conductance (µmol m−2 s-1)
#' @export
#' @references
#' Ball, J. T., Woodrow, I. E., and Berry, J. A. (1987). A model predicting stomatal conductance and its contribution to the control of photosynthesis under different environmental conditions. In Progress in Photosynthesis Research (pp. 221–224). Springer, Dordrecht. http://doi.org/10.1007/978-94-017-0519-6-48.\cr
#' Leuning, R. (1990). Modelling stomatal behaviour and photosynthesis of Eucalyptus grandis. Functional Plant Biology, 17(2), 159–175.\cr
#' Leuning, R. (1995). A critical appraisal of a combined stomatal‐photosynthesis model for C3 plants. Plant, Cell & Environment, 18(4), 339–355.
#' @examples
#' cs <- seq(100, 200, 5)
#' data <- data.frame(cs = cs, An = rep(15, length(cs)), hs = rep(0.6, length(cs)))
#' photo.ballberry(data, 0.6, 3)
photo.ballberry <- function(data, g0, g1, tau=NULL, VPD_0=NULL, model = "bb"){
  if(model=="bb"){
    if(is.null(data$An)) stop("couldn't find 'An' in 'data'")
    if(is.null(data$cs)) stop("couldn't find 'cs' in 'data'")
    if(is.null(data$hs)) stop("couldn't find 'hs' in 'data'")

    g0 + g1 * data$hs * data$An / data$cs
  } else if(model=="leuning"){
    if(is.null(data$An)) stop("couldn't find 'An' in 'data'")
    if(is.null(data$cs)) stop("couldn't find 'cs' in 'data'")
    if(is.null(data$hs)) stop("couldn't find 'hs' in 'data'")
    if(is.null(tau)) stop("parameter 'tau' not specified")

    g0 + g1 * data$An * data$hs / (data$cs - tau)
  } else if(model=="c3"){
    if(is.null(data$An)) stop("couldn't find 'An' in 'data'")
    if(is.null(data$cs)) stop("couldn't find 'cs' in 'data'")
    if(is.null(data$hs)) stop("couldn't find 'hs' in 'data'")
    if(is.null(data$VPD)) stop("couldn't find 'VPD' in 'data'")
    if(is.null(tau)) stop("parameter 'tau' not specified")
    if(is.null(VPD_0)) stop("parameter 'VPD_0' not specified")

    g0 + g1 * data$An * data$hs / ((data$cs - tau) * (1 + data$VPD/VPD_0))
  } else if(model=="vpd"){
    if(is.null(data$An)) stop("couldn't find 'An' in 'data'")
    if(is.null(data$cs)) stop("couldn't find 'cs' in 'data'")
    if(is.null(data$VPD)) stop("couldn't find 'VPD' in 'data'")

    g0 + 1.6 * (1 + g1/sqrt(data$VPD)) * data$An/data$cs
  } else{
    stop("invalid input 'model'")
  }
}


##photo.ballberry.fit------------------------
#' Fit Ball-Berry Model
#'
#' \code{photo.ballberry.fit} is used to fit Ball-Berry Model for calculating
#' stomatal conductance (µmol m−2 s-1). It uses linear regression after
#' linearization of the Ball-Berry model.
#' @param data a \code{data.frame} with several columns. It contains required data
#' to compute stomatal conductance for different models. The column names must be:\cr
#' "\code{gs}": for stomatal conductance (µmol m−2 s-1).\cr
#' "\code{An}": for photosynthesis rate (µmol m−2 s-1)\cr
#' "\code{cs}": for CO2 concentration of leaf surface (µmol mol-1)\cr
#' "\code{hs}": for relative humidity at the leaf surface (a value between 0-1,
#' rather than percentage)\cr
#' "\code{VPD}": vapor pressure deficit (kPa).
#' @param model type of model used.\cr
#' "\code{bb}": the original Ball-Berry model.
#' "\code{launing}": model proposed by Leuning in 1990 by adding parameter \code{tau}
#' into the original Ball-Berry model.
#' "\code{vpd}": model proposed by Leuning in 1995. VPD, instead of relative humidity,
#' is used in the model.
#' @param g0 potential g0 for linear regression. A vector. Only used when \code{model=="launing"}.
#'
#' @return a list of four objects: \cr\cr
#' \code{$formula}: formula of the model.\cr
#' \code{$coef}: fitted coefficients of the model.\cr
#' \code{$r.squared}: R square of the model.\cr
#' \code{$adj.r.squared}: adjusted R square of the model.\cr
#' @export
#' @references
#' Ball, J. T., Woodrow, I. E., and Berry, J. A. (1987). A model predicting stomatal conductance and its contribution to the control of photosynthesis under different environmental conditions. In Progress in Photosynthesis Research (pp. 221–224). Springer, Dordrecht. http://doi.org/10.1007/978-94-017-0519-6-48.\cr
#' Leuning, R. (1990). Modelling stomatal behaviour and photosynthesis of Eucalyptus grandis. Functional Plant Biology, 17(2), 159–175.\cr
#' Leuning, R. (1995). A critical appraisal of a combined stomatal‐photosynthesis model for C3 plants. Plant, Cell & Environment, 18(4), 339–355.
#' @examples
#' cs <- seq(100, 200, 5)
#' data <- data.frame(cs = cs, An = rep(15, length(cs)), hs = rep(0.6, length(cs)))
#' data$gs <- photo.ballberry(data, g0 = 0.6, g1 = 3, tau = 20, model = "leuning") + jitter(0.01)
#' photo.ballberry.fit(data, model = "leuning")
photo.ballberry.fit <- function(data, model = "bb",
                                g0 = seq(min(data$gs)/2,mean(data$gs),0.001)){
  if(model=="bb"){
    if(is.null(data$gs)) stop("couldn't find 'gs' in 'data'")
    if(is.null(data$An)) stop("couldn't find 'An' in 'data'")
    if(is.null(data$cs)) stop("couldn't find 'cs' in 'data'")
    if(is.null(data$hs)) stop("couldn't find 'hs' in 'data'")

    x <- with(data, hs*An/cs)
    md <- lm(data$gs~x)
    list(
      formula = "gs = g0 + g1 * h2s * An / cs",
      coef = c(g0 = md$coefficients[[1]], g1 = md$coefficients[[2]]),
      r.square = summary(md)$r.squared,
      adj.r.square = summary(md)$adj.r.squared
    )
  } else if(model=="leuning"){
    if(is.null(data$gs)) stop("couldn't find 'gs' in 'data'")
    if(is.null(data$An)) stop("couldn't find 'An' in 'data'")
    if(is.null(data$cs)) stop("couldn't find 'cs' in 'data'")
    if(is.null(data$hs)) stop("couldn't find 'hs' in 'data'")

    r.square <- do.call(c,lapply(g0, function(gg0){
      y <- with(data, cs*(gs-gg0))
      x1 <- with(data,gs-gg0)
      x2 <- with(data,An*hs)
      summary(lm(y~x1+x2-1))$r.squared
    }))

    index <- which(max(r.square)==r.square)
    if(length(index)==1){
      if(index==1) warning("it is recommended to reduce the lower boundary of vector 'g0'")
      if(index==length(g0)) warning("it is recommended to increase the upper boundary of vector 'g0'")
    }

    g0 <- g0[index[1]]
    y <- with(data, cs*(gs-g0))
    x1 <- with(data,gs-g0)
    x2 <- with(data,An*hs)
    md <- lm(y~x1+x2-1)

    list(
      formula = "gs = g0 + g1 * An * hs / (cs - tau)",
      coef = c(g0=g0,g1=md$coefficients[[2]],tau=md$coefficients[[1]]),
      r.square = summary(md)$r.squared,
      adj.r.square = summary(md)$adj.r.squared
    )
  } else if(model=="vpd"){
    y <- with(data, gs-1.6*An/cs)
    x <- with(data, 1.6*An/(cs*sqrt(VPD)))
    md <- lm(y~x)

    list(
      formula = "gs = g0 + 1.6 * (1 + g1/(sqrt(VPD))) * An / cs",
      coef = c(g0 = md$coefficients[[1]], g1 = md$coefficients[[2]]/1.6),
      r.square = summary(md)$r.squared,
      adj.r.square = summary(md)$adj.r.squared
    )
  } else{
    stop("invalid input 'model'")
  }
}
