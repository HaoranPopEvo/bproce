##Script info----------------------
# Created on 2022-06-04 11:14:04 CST
#
# author: Hao-Ran Wu
# College of Life Sciences
# Zhejiang University
# Hangzhou 310012, China
# E-mail: haoranwu@zju.edu.cn
#
# This script provides function to fit ecosystem respiration model and calculate
#  Q10 value. 'respir.fit' is the only function in this script. Use different
#  options in argument 'method' to apply different models.
#      option               model                     input variable(s)
#   -- 't1.log-linear':     Log-linear Model          temperature
#   -- 't2.linear':         Linear Model              temperature
#   -- 't3.polynomial':     Polynomial Model          temperature
#   -- 't4.arrhenius':      Arrhenius Model           temperature
#   -- 't5.logistic':       Logistic Model            temperature
#   -- 't6.gamma':          Gamma Model               temperature
#   -- 'tm1.DeForest':      DeForest et al., 2006     temperature, soil moisture
#   -- 'tm2.Xu':            Xu et al., 2011           temperature, soil moisture
#   -- 'tm3.Concilio':      Concilio et al., 2005     temperature, soil moisture
#   -- 'tmd.Xu':            Xu et al., 2011           temperature, soil moisture, DOY
#   -- 'd.fourier':         Davidson et al., 2006     DOY

#respir.fit---------------------------
#' Models for Ecosystem Respiration
#'
#' \code{respir.fit} is used to fit ecosystem respiration models.\cr
#'
#' @param R vector of respiration (μmol CO2 m-2 s-1).
#' @param ta vector of temperature (°C)
#' @param m vector of soil moisture (%)
#' @param doy vector of day of the year
#' @param method type of models used. Valid options are: \code{"t1.log-linear"},
#'  \code{"t2.linear"},  \code{"t3.polynomial"},  \code{"t4.arrhenius"},
#'  \code{"t5.logistic"},  \code{"t6.gamma"}, \code{"tm1.DeForest"},
#'  \code{"tm2.Xu"}, \code{"tm3.Concilio"}, \code{"tmd.Xu"}, \code{"d.fourier"}.
#' @param n order of polynomial model. Only used when \code{method == "t3.polynomial"}.
#' @param Q10st start temperature to calculate Q10 value.
#' @importFrom stats lm
#' @importFrom stats nls
#' @importFrom stats predict
#' @importFrom stats coefficients
#'
#' @return a list of four or five objects: \cr\cr
#' \code{$formula}: formula of the model.\cr
#' \code{$coef}: fitted coefficients of the model.\cr
#' \code{$r.squared}: R square of the model.\cr
#' \code{$adj.r.squared}: adjusted R square of the model.\cr
#' \code{$Q10}: Q10 value calculated from the model.
#'
#' @references
#' Wofsy, S. C., Goulden, M. L., Munger, J. W., Fan, S. M., Bakwin, P. S., Daube, B. C., Bassow, S. L., and Bazzaz, F. A. (1993). Net exchange of CO2 in a mid-latitude forest. Science, 260(5112), 1314–1317.\cr
#' Yu, X., Zha, T., Pang, Z., Wu, B., Wang, X., Chen, G., Li, C., Cao, J., Jia, G., Li, X., and Wu, H. (2011). Response of soil respiration to soil temperature and moisture in a 50-year-old oriental arborvitae plantation in China. PloS One, 6(12), e28397.\cr
#' Görres, C. M., Kammann, C., and Ceulemans, R. (2016). Automation of soil flux chamber measurements: Potentials and pitfalls. Biogeosciences, 13(6), 1949–1966.\cr
#' Barr, A. G., Griffis, T. J., Black, T. A., Lee, X., Staebler, R. M., Fuentes, J. D., Chen, Z., and Morgenstern, K. (2002). Comparing the carbon budgets of boreal and temperate deciduous forest stands. Canadian Journal of Forest Research, 32(5), 813–822.\cr
#' Khomik, M., Arain, M. A., Liaw, K. L., and McCaughey, J. H. (2009). Debut of a flexible model for simulating soil respiration–soil temperature relationship: Gamma model. Journal of Geophysical Research: Biogeosciences, 114(G3), doi: 10.1029/2008JG000851.\cr
#' DeForest, J. L., Noormets, A., McNulty, S. G., Sun, G., Tenney, G., and Chen, J. (2006) Phenophases alter the soil respiration-temperature relationship in an oak-dominated forest. International Journal of Biometeorology, 51(2), 135–144.\cr
#' Martin, J. G., Bolstad, P. V., Ryu, S. R., and Chen, J. Q.(2009). Modeling soil respiration based on carbon, nitrogen, and root mass across diverse Great Lake forests. Agricultural and Forest Meteorology, 149(10), 1722–172.\cr
#' Xu, M., and Qi, Y. (2001). Soil-surface CO2 efflux and its spatial and temporal variations in a young ponderosa pine plantation in northern California. Global Change Biology, 7(6), 667–677.\cr
#' Concilio, A., Ma, S., Li, Q., LeMoine, J., Chen, J., North, M., Moorhead, D., and Jensen, R. (2005). Soil respiration response to prescribed burning and thinning in mixed-conifer and hardwood forests. Canadian Journal of Forest Research, 35(7), 1581–1591.\cr
#' Davidson, E. A., Janssens, I. A., and Luo, Y. (2006). On the variability of respiration in terrestrial ecosystems: Moving beyond Q10. Global Change Biology, 12(2), 154–164.\cr
#' @export
#' @examples
#' ##1 log-linear model
#' ta <- c(-1.64,-1.59,-1.57,-1.57,-1.51,-1.50,-1.50,-1.50,
#' -1.50,-1.48,-1.47,-1.42,-1.41,-1.41,-1.41,-1.39,-1.39)
#' R <- c(0.491,0.493,0.495,0.495,0.498,0.499,0.499,0.499,
#' 0.499,0.500,0.500,0.503,0.504,0.504,0.504,0.505,0.505)
#' respir.fit(R, ta)
#'
#' ##2 logistic model
#' ta <- 0:40
#' R <- 5.8164/(1+exp(3.8513-0.3105*ta))
#' respir.fit(R, ta, method = "t5.logistic")
#'
#' ##3 lon then linear
#' tm <- expand.grid(20:40, seq(50,80,5))
#' ta <- tm[,1]
#' m <- tm[,2]
#' R <- 4.0896 * exp(0.0438 * ta) + 0.2178 * m -8.3419
#' respir.fit(R, ta, m, method = "tm1.DeForest")
#'
#' ##4 log then parabola
#' R <- 66.8272 * exp(0.1004 * ta) * 8.95373781e-6 * (m + 17.1191)^2
#' respir.fit(R, ta, m, method = "tm2.Xu")
#'
#' ##5 Concilio
#' R <- 319.9639 * exp(0.0310 * ta) * exp(-0.0027 * m) * 0.00003 * ta * m
#' respir.fit(R, ta, m, method = "tm3.Concilio")
#'
#' ##6 Xu
#' tmd <- expand.grid(5:15, seq(50,80,5), seq(100,300,25))
#' ta <- tmd[,1]
#' m <- tmd[,2]
#' doy <- tmd[,3]
#' R <- 99.29090 * exp(0.04391 * ta) * 7.37849575e-6 *
#' (m + 56.15651)^2 - 4.96519958e-7 * (doy + 2721.70967)^2
#' respir.fit(R, ta, m, doy, method = "tmd.Xu")
respir.fit <- function(R, ta=NULL, m=NULL, doy=NULL,
                       method="t1.log-linear",
                       n=2, Q10st = 5){
  use_nls <- function(formula, start_values){
    tryCatch({
      nls(formula, start = start_values)
    }, error = function(e){
      tryCatch({
        nls(formula, start = start_values, algorithm = "plinear")
      }, error = function(e){
        tryCatch({
          nls(formula, start = start_values, algorithm = "port")
        }, error = function(e){
          stop("cannot solve equation using 'nls' model.")
        })
      })
    })
  }

  if(method=="t1.log-linear"){
    md <- lm(log(R) ~ ta)
    list(
      formula = "R = a0 * exp(a1 * ta)",
      coef = c(a0=exp(md$coefficients[[1]]),a1=md$coefficients[[2]]),
      r.squared = summary(md)$r.squared,
      adj.r.squared = summary(md)$adj.r.squared,
      Q10 = as.numeric(exp(predict(md,list(ta=Q10st+10)))/exp(predict(md,list(ta=Q10st))))
    )
  } else if(method=="t2.linear"){
    md <- lm(R ~ ta)
    coef <- md$coefficients
    names(coef) <- c("a0","a1")
    list(
      formula = "R = a0 + a1 * ta",
      coef = coef,
      r.squared = summary(md)$r.squared,
      adj.r.squared = summary(md)$adj.r.squared,
      Q10 = as.numeric((coef[1] + coef[2] * (Q10st+10))/(coef[1] + coef[2] * Q10st))
    )
  } else if(method=="t3.polynomial"){
    md <- eval(parse(text = paste0("lm(R ~ ",paste(paste0("I(ta^",1:n,")"),collapse = " + ")," )")))
    coef <- md$coefficients
    names(coef) <- paste0("a",0:(length(names(coef))-1))
    list(
      formula = paste0("a0 + a1 * ta + ",paste(paste0("a",2:n," * ta^",2:n),collapse = " + ")),
      coef = coef,
      r.squared = summary(md)$r.squared,
      adj.r.squared = summary(md)$adj.r.squared,
      Q10 = as.numeric(predict(md,list(ta=Q10st+10))/predict(md,list(ta=Q10st)))
    )
  } else if(method=="t4.arrhenius"){
    md <- lm(log(R) ~ I(1/(ta+273.15-227.13)))
    E0 <- -as.numeric(md$coefficients[2])
    R10 <- exp(as.numeric(md$coefficients[1])-E0/56.03)
    list(
      formula = "R = R10 * exp(E0 * (1/56.02 - 1/(T + 273.15 - 227.13)))",
      coef = c(R10=R10, E0=E0),
      r.squared = summary(md)$r.squared,
      adj.r.squared = summary(md)$adj.r.squared,
      Q10 = as.numeric(exp(predict(md,list(ta=Q10st+10)))/exp(predict(md,list(ta=Q10st))))
    )
  } else if(method=="t5.logistic"){
    Krange <- seq(max(R)/2,max(R)*2,(max(R)*2-max(R)/2)/100)
    adj.r <- do.call(c,lapply(Krange, function(K){
      R.t <- R[(K-R)/R>0]
      ta.t <- ta[(K-R)/R>0]
      md <- lm(I(log((K-R.t)/R.t)) ~ ta.t)
      summary(md)$adj.r.squared
    }))
    index <- which(max(adj.r)==adj.r)
    K <- Krange[index]
    R.t <- R[(K-R)/R>0]
    ta.t <- ta[(K-R)/R>0]
    md <- lm(I(log((K-R.t)/R.t)) ~ ta.t)

    formula <- R ~ alpha / (1 + exp(beta0 - beta1 * ta))
    start_values <- list(
      alpha = K, beta0 = md$coefficients[[1]], beta1 = -md$coefficients[[2]]
    )

    md <- use_nls(formula, start_values)
    r.squared <- 1-sum((R-predict(md))^2)/sum((R-mean(R))^2)
    list(
      formula = "R = alpha / (1 + exp(beta0 - beta1 * ta))",
      coef = summary(md)$coefficients[,1],
      r.squared = r.squared,
      adj.r.squared = 1-(1-r.squared)*(length(R)-1)/(length(R)-1-1),
      Q10 = predict(md,list(ta=Q10st+10))/predict(md,list(ta=Q10st))
    )
  } else if(method=="t6.gamma"){
    ta.t <- ta[ta>0]
    R.t <- R[ta>0]
    md <- lm(log(R.t) ~ log(ta.t) + ta.t)
    coef <- md$coefficients
    names(coef) <- c("beta0","alpha","beta1")
    list(
      formula = "R = ta ^ alpha * exp(beta0 + beta1 * ta)",
      coef = coef,
      r.squared = summary(md)$r.squared,
      adj.r.squared = summary(md)$adj.r.squared,
      Q10 = as.numeric(exp(predict(md,list(ta.t=Q10st+10)))/exp(predict(md,list(ta.t=Q10st))))
    )
  } else if(method=="tm1.DeForest"){
    md1 <- lm(log(R) ~ ta)
    resid <- R - exp(predict(md1))
    md2 <- lm(resid ~ m)
    coef <- c(
      R10 = exp(md1$coefficients[[1]]),
      beta = md1$coefficients[[2]],
      a = md2$coefficients[[2]],
      b = md2$coefficients[[1]]
    )

    formula <- R ~ R10 * exp(beta * ta) + a * m + b
    start_values <- list(
      R10 = coef[[1]], beta = coef[[2]], a = coef[[3]], b = coef[[4]]
    )

    md <- use_nls(formula, start_values)

    r.squared <- 1-sum((R-predict(md))^2)/sum((R-mean(R))^2)

    list(
      formula = "R = R10 * exp(beta * ta) + a * m + b",
      coef = summary(md)$coefficients[,1],
      r.squared = r.squared,
      adj.r.squared = 1-(1-r.squared)*(length(R)-1)/(length(R)-2-1)
    )


  } else if(method=="tm2.Xu"){
    md1 <- lm(log(R) ~ ta)
    resid <- R/exp(predict(md1))
    md2 <- lm(resid ~ m + I(m^2))

    formula <- R ~ alpha_plus_beta1 * exp(beta0 * ta) * (m - beta2)^2
    start_values <- list(
      alpha_plus_beta1 = exp(md1$coefficients[[1]]) *  md2$coefficients[[1]],
      beta0 = md1$coefficients[[2]],
      beta2 = md2$coefficients[[2]]/(2*md2$coefficients[[1]])
    )

    md <- use_nls(formula, start_values)
    coef <- c(
      R0 = exp(md1$coefficients[[1]]),
      beta = md1$coefficients[[2]],
      a = md2$coefficients[[3]],
      b = md2$coefficients[[2]],
      c = md2$coefficients[[1]]
    )

    r.squared <- 1-sum((R-predict(md))^2)/sum((R-mean(R))^2)
    list(
      formula = "R = alpha * exp(beta0 * T) * beta1 * (m - beta2)^2",
      coef = summary(md)$coefficients[,1],
      r.squared = r.squared,
      adj.r.squared = 1-(1-r.squared)*(length(R)-1)/(length(R)-3-1)
    )
  } else if(method=="tm3.Concilio"){
    ta.t <- ta[ta>0]
    m.t <- m[ta>0]
    R.t <- R[ta>0]
    y.t <- log(R.t) - log(ta.t) - log(m.t)
    md <- lm(y.t ~ ta.t + m.t)
    coef <- md$coefficients
    names(coef) <- c("log(R0)+log(beta2)", "beta0", "beta1")
    list(
      formula = "R = R0 * exp(beta0 * ta) * exp(beta1 * m) * ta * m",
      coef = coef,
      r.squared = summary(md)$r.squared,
      adj.r.squared = summary(md)$adj.r.squared
    )
  } else if(method=="tmd.Xu"){
    md1 <- lm(R~I(doy^2)+I(doy))
    md2 <- lm(log(R) ~ ta)
    resid <- R/exp(predict(md2))
    md3 <- lm(resid ~ I(m^2)+I(m))

    formula <- R ~ alpha_plus_beta1 * exp(beta0 * ta) * (m - beta2)^2 + beta3 * (doy - beta4)^2
    start_values <- list(
      alpha_plus_beta1 = exp(md2$coefficients[[1]]) * md3$coefficients[[2]],
      beta0 = md2$coefficients[[2]],
      beta2 = md3$coefficients[[3]]/(-2*md3$coefficients[[2]]),
      beta3 = md1$coefficients[[2]],
      beta4 = md1$coefficients[[3]]/(-2*md1$coefficients[[2]])
    )
    md <- use_nls(formula, start_values)

    r.squared <- 1-sum((R-predict(md))^2)/sum((R-mean(R))^2)
    list(
      formula = "R = alpha * exp(beta0 * ta) * beta1 * (m - beta2)^2 + beta3 * (doy - beta4)^2",
      coef = summary(md)$coefficients[,1],
      r.squared = r.squared,
      adj.r.squared = 1-(1-r.squared)*(length(R)-1)/(length(R)-5-1)
    )
  } else if(method=="d.fourier"){
    F1range <- F2range <- seq(0,2*pi,0.1)
    paraSpace <- expand.grid(fai1=F1range,fai2=F2range)

    cat("optimize parameter settings...")
    adj.r <- apply(paraSpace, MARGIN=1,function(xx){
      fai1 <- as.numeric(xx[1])
      fai2 <- as.numeric(xx[2])
      md <- lm(R ~ I(sin(doy*2*pi/365+fai1)) + I(sin(2*doy*2*pi/365+fai2)))
      summary(md)$adj.r.squared
    })
    cat("done.\n")
    fai1 <- paraSpace[which(max(adj.r)==adj.r)[1],1]
    fai2 <- paraSpace[which(max(adj.r)==adj.r)[1],2]
    md <- lm(R ~ I(sin(doy*2*pi/365+fai1)) + I(sin(2*doy*2*pi/365+fai2)))

    formula <- R ~ kappa0 + kappa1 * sin(doy * 2*pi/365 + fai1) + kappa2 * sin(2 * doy * 2*pi/365 + fai2)
    start_values <- list(
      kappa0 = md$coefficients[[1]], kappa1 = md$coefficients[[2]],
      fai1 = fai1, kappa2 = md$coefficients[[3]], fai2 = fai2
    )
    md <- use_nls(formula, start_values)
    coef <- summary(md)$coefficients[,1]

    r.squared <- 1-sum((R-predict(md))^2)/sum((R-mean(R))^2)
    list(
      formula = "kappa0 + kappa1 * sin(doy * 2*pi/365 + fai1) + kappa2 * sin(2 * doy * 2*pi/365 + fai2)",
      coef = coef,
      r.squared = r.squared,
      adj.r.squared = 1-(1-r.squared)*(length(R)-1)/(length(R)-2-1)
    )
  } else{
    stop("invalid option 'method'")
  }
}
