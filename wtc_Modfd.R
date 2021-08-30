# This is a simple modification to the wtc function from the package biwavelet. With some data the original 
# wtc functions throws an error: non-stationary AR part from CSS
# I have added the method = "ML" (maximulm likelihood) as one of the arguments for the function arima.
# Dated: 30th August 2021

wtc_Modfd <- function (d1, d2, pad = TRUE, dj = 1/12, s0 = 2 * dt, J1 = NULL, 
          max.scale = NULL, mother = "morlet", param = -1, lag1 = NULL, 
          sig.level = 0.95, sig.test = 0, nrands = 300, quiet = FALSE) 
{
  mother <- match.arg(tolower(mother), MOTHERS)
  checked <- check.data(y = d1, x1 = d2)
  xaxis <- d1[, 1]
  dt <- checked$y$dt
  t <- checked$y$t
  n <- checked$y$n.obs
  if (is.null(J1)) {
    if (is.null(max.scale)) {
      max.scale <- (n * 0.17) * 2 * dt
    }
    J1 <- round(log2(max.scale/s0)/dj)
  }
  if (is.null(lag1)) {
    d1.ar1 <- arima(d1[, 2], order = c(1, 0, 0),method="ML")$coef[1] # Added method to be ML to avoid error: non-stationary AR part from CSS
    d2.ar1 <- arima(d2[, 2], order = c(1, 0, 0),method="ML")$coef[1] # Added method to be ML to avoid error: non-stationary AR part from CSS
    lag1 <- c(d1.ar1, d2.ar1)
  }
  wt1 <- wt(d = d1, pad = pad, dj = dj, s0 = s0, J1 = J1, max.scale = max.scale, 
            mother = mother, param = param, sig.level = sig.level, 
            sig.test = sig.test, lag1 = lag1[1])
  wt2 <- wt(d = d2, pad = pad, dj = dj, s0 = s0, J1 = J1, max.scale = max.scale, 
            mother = mother, param = param, sig.level = sig.level, 
            sig.test = sig.test, lag1 = lag1[2])
  d1.sigma <- sd(d1[, 2], na.rm = T)
  d2.sigma <- sd(d2[, 2], na.rm = T)
  s.inv <- 1/t(wt1$scale)
  s.inv <- matrix(rep(s.inv, n), nrow = NROW(wt1$wave))
  smooth.wt1 <- smooth.wavelet(s.inv * (abs(wt1$wave)^2), dt, 
                               dj, wt1$scale)
  smooth.wt2 <- smooth.wavelet(s.inv * (abs(wt2$wave)^2), dt, 
                               dj, wt2$scale)
  coi <- pmin(wt1$coi, wt2$coi, na.rm = T)
  CW <- wt1$wave * Conj(wt2$wave)
  CW.corr <- (wt1$wave * Conj(wt2$wave) * max(wt1$period))/matrix(rep(wt1$period, 
                                                                      length(t)), nrow = NROW(wt1$period))
  power <- abs(CW)^2
  power.corr <- (abs(CW)^2 * max.scale)/matrix(rep(wt1$period, 
                                                   length(t)), nrow = NROW(wt1$period))
  smooth.CW <- smooth.wavelet(s.inv * (CW), dt, dj, wt1$scale)
  rsq <- abs(smooth.CW)^2/(smooth.wt1 * smooth.wt2)
  phase <- atan2(Im(CW), Re(CW))
  if (nrands > 0) {
    signif <- wtc.sig(nrands = nrands, lag1 = lag1, dt = dt, 
                      ntimesteps = n, pad = pad, dj = dj, J1 = J1, s0 = s0, 
                      max.scale = max.scale, mother = mother, sig.level = sig.level, 
                      quiet = quiet)
  }
  else {
    signif <- NA
  }
  results <- list(coi = coi, wave = CW, wave.corr = CW.corr, 
                  power = power, power.corr = power.corr, rsq = rsq, phase = phase, 
                  period = wt1$period, scale = wt1$scale, dt = dt, t = t, 
                  xaxis = xaxis, s0 = s0, dj = dj, d1.sigma = d1.sigma, 
                  d2.sigma = d2.sigma, mother = mother, type = "wtc", 
                  signif = signif)
  class(results) <- "biwavelet"
  return(results)
}
