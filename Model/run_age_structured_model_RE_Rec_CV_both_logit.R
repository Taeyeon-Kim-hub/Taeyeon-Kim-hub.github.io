# Run TMB (Region setting)
run_tmb_region <- function(region_name, dll_name = "age_structured_model_RE_Rec_CV_both_logit") {
  target_yield <- data_yield_mackerel[, region_name]
  target_cpue  <- data_cpue_mackerel[, region_name]
  if (region_name == "Total") {
    target_waa_vec <- mean_WAA_Total
  } else if (region_name == "East") {
    target_waa_vec <- mean_WAA_East
  } else if (region_name == "West") {
    target_waa_vec <- mean_WAA_West
  } else if (region_name == "South") {
    target_waa_vec <- mean_WAA_South
  }
  target_WAA_matrix <- matrix(rep(target_waa_vec, nyrs), nrow = nyrs, ncol = nages, byrow = TRUE)
  data <- list(
    nyears = nyrs, nages = nages, 
    WAA_ini = target_WAA_matrix,
    CAA_obs = CAA_obs_matrix,    
    yield = as.numeric(target_yield), CPUE = as.numeric(target_cpue),
    lambda = lambda_vec, Spawn_month = Spawn_month_val, range_q = bounds_q,
    M = M_val, maturate = mat_at_age, Fref = Fref_vec, SS_age = SS_age_vec, musig_logit_h = musig_logit_h_val,
    ratio_female = ratio_female_val, 
    CV_yield = CV_yield_val,
    CV_cpue = CV_cpue_val
  )
  
  parameters <- region_settings[[region_name]]$init
  Low_bound  <- region_settings[[region_name]]$low
  Up_bound   <- region_settings[[region_name]]$up
  
  random_vars <- c("logRec","logFt")
  map_par <- list() 
  
  obj <- MakeADFun(data = data, parameters = parameters, DLL = dll_name, random = random_vars, map = map_par, silent = TRUE)
  
  lower_vec <- rep(-Inf, length(obj$par))
  upper_vec <- rep(Inf, length(obj$par))
  names(lower_vec) <- names(obj$par)
  names(upper_vec) <- names(obj$par)
  
  for(p in names(Low_bound)) if(p %in% names(lower_vec)) lower_vec[names(lower_vec) == p] <- Low_bound[[p]]
  for(p in names(Up_bound))  if(p %in% names(upper_vec)) upper_vec[names(upper_vec) == p] <- Up_bound[[p]]
  
  fit <- nlminb(
    start = obj$par, objective = obj$fn, gradient = obj$gr,
    lower = lower_vec, upper = upper_vec,
    control = list(eval.max = 10000, iter.max = 10000) 
  )
  
  sdrep <- sdreport(obj)
  return(list(region = region_name, fit = fit, obj = obj, sdrep = sdrep, reports = obj$report()))
}

# Base setting
# Age selection, maturity at age, years, sample size
nages <- 6
nyrs <- nrow(data_yield_mackerel)
mat_at_age <- c(0.6, 0.85, 1.0, 1.0, 1.0, 1.0)
CAA_obs_matrix <- as.matrix(data_age_comp_annual[, -1])
SS_age_vec <- rowSums(CAA_obs_matrix)
ratio_female_val <- 0.6
M_val <- 0.53
Fref_vec <- seq(0, 2.0, length.out = 100) 
lambda_vec <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0) 
bounds_q <- c(exp(-40.0), exp(-2.0))

# Other values (Regional specific)
Spawn_month_val <- 4
musig_logit_h_val <- c(0.6, 1.5)  # non-informative (0.64,0.58)
CV_yield_val <- 0.15
CV_cpue_val <- 0.05
#-----------------------------------------------------------------------------#
# Total (NOT regional)
region_settings <- list()
region_settings[["Total"]] <- list(
  init = list(
    log_Rec_fir = 20, logR0 = 21, logRec = rep(20, nyrs), 
    log_a50 = log(2.0), log_a95 = log(3.0), 
    logit_q = -5.0, 
    logit_h = 0.59,
    logF1 = log(0.2), logFt = rep(log(0.2), nyrs), log_sig_logF = log(0.2),log_sig_logRec = log(0.4), logtheta = log(0.05) 
  ),
  low = list(
    "log_Rec_fir" = log(10000), "logR0" = log(10000), 
    "log_a50" = log(0.00001), "log_a95" = log(1.0), 
    "logit_q" = -30.0,
    "logit_h" = -5.0,
    "logF1" = log(0.001), "log_sig_logF" = log(0.01),"log_sig_logRec" = log(0.001), "logtheta" = log(0.0000001)
  ),
  up = list(
    "log_Rec_fir" = log(1e20), "logR0" = log(1e20), 
    "log_a50" = log(10.0), "log_a95" = log(15.0), 
    "logit_q" = 30.0, 
    "logit_h" = 10.0,
    "logF1" = log(5.0), "log_sig_logF" = log(2.0),"log_sig_logRec" = log(2.0), "logtheta" = log(100.0)
  )
)

# Result checking
result_T <- run_tmb_region("Total")
result_T

fit_T <- result_T$fit
fit_T

sdrep_T <- result_T$sdrep
sdrep_T

result_T$reports

summestim_T = summary(sdrep_T)
summestim_T

TrueOrFalse_Hessian_T = sdrep_T$pdHess
TrueOrFalse_Hessian_T

max.abs.grad_T = max(abs(sdrep_T$gradient.fixed))
max.abs.grad_T
# ==============================================================================
# [Plot 3] YPR, SPR & Equilibrium Yield Curves (Reference Points)
# ==============================================================================
reports_T <- result_T$reports
years <- 2000:(2000 + nyrs - 1) 
ages <- 1:nages
YPR_vec_T <- reports_T$YPR
SPR_ratio_T <- reports_T$SPR
Yield_eq_vec_T <- reports_T$Yield_F
MSY_T <- reports_T$MSY

F_msy_idx_T <- which.max(Yield_eq_vec_T)
F_msy_T <- Fref_vec[F_msy_idx_T]

par(mfrow=c(1,2), mar=c(5,5,3,4), oma=c(0,0,0,0))


#============================================================================
# [Plot 1] Yield & CPUE Goodness of Fit (with MSY)
# ==============================================================================

par(mfrow=c(1,2), mar=c(5,5,2,3), oma=c(0,0,2,0))
# 1-1. Yield Fitting
obs_yield_T <- reports_T$yield
pred_yield_T <- reports_T$Yield_hat/1000

plot(years, obs_yield_T, ylim=c(0, max(obs_yield_T* 1.2)), 
     xlab="Year", ylab=expression(paste("MT")), 
     main="Yield - Total", 
     cex.lab=1.5, cex.axis=1.5, cex.main=2.0, pch=1)
lines(years, pred_yield_T, col="red", lwd=3)

abline(h = MSY_T, col="blue", lty=2, lwd=2)

# 1-2. CPUE Fitting
obs_cpue <- reports_T$CPUE
pred_cpue <- reports_T$predcpue

plot(years, obs_cpue, ylim=c(0, max(obs_cpue) * 1.2), 
     xlab="Year", ylab="MT/Fishing effort", 
     main="CPUE - Total", 
     cex.lab=1.5, cex.axis=1.5, cex.main=2.0, pch=1)
lines(years, pred_cpue, col="red", lwd=3)

# ==============================================================================
# [Plot 2] Age Composition Fitting
# ==============================================================================
grid_dim <- ceiling(sqrt(nyrs)) 
par(mfrow=c(grid_dim, grid_dim), oma=c(2,2,3,1), mar=c(2, 2, 1.5, 1), cex=1.0)

CAA_obs <- reports_T$CAA_obs
CAA_pred <- reports_T$CAA_pred

for(i in 1:nyrs) {
  obs_freq <- (CAA_obs[i,] / sum(CAA_obs[i,])) * SS_age_vec[i]
  pred_freq <- (CAA_pred[i,] / sum(CAA_pred[i,])) * SS_age_vec[i]
  
  y_max <- max(c(obs_freq, pred_freq), na.rm=TRUE) * 1.1
  if (is.na(y_max) | y_max == 0) y_max <- 1 
  
  plot(ages, obs_freq, type="h", ylim=c(0, y_max), main=paste(years[i]), 
       xlab="", ylab="", lwd=5, col="black") 
  lines(ages, pred_freq, col='red', lwd=3)
}
mtext("Age Composition - Total", side=3, outer=TRUE, cex=1.5, font=2)
mtext("Age", side=1, outer=TRUE, line=0.5, cex=1.2)
mtext("Frequency", side=2, outer=TRUE, line=0.5, cex=1.2)
# ---------------------------------------------------------
# 3-1. YPR (Yield Per Recruit) Curve
# ---------------------------------------------------------
# plot(Fref_vec, YPR_vec_T * 1000, type="l", col="black", lwd=3,
#      xlab="Fishing mortality", ylab="Yield Per Recruit (kg)", main="Yield Per Recruit (YPR)",
#      cex.lab=1.5, cex.axis=1.5, cex.main=1.8)
# 

# ---------------------------------------------------------
# 3-2. SPR (Spawning Potential Ratio) Curve
# ---------------------------------------------------------
par(mfrow=c(1,2))
plot(Fref_vec, SPR_ratio_T, type="l", col="black", lwd=3,
     xlab="Fishing mortality", ylab="SPR Ratio", ylim=c(0, 1.0), main="Spawning Potential Ratio (SPR)",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.8)

abline(h = c(0.3, 0.4), col=c("red", "orange"), lty=2, lwd=2)
text(x = max(Fref_vec)*0.85, y = 0.45, "40% SPR", col="orange", font=2, cex=1.2)
text(x = max(Fref_vec)*0.85, y = 0.35, "30% SPR", col="red", font=2, cex=1.2)

# ---------------------------------------------------------
# 3-3. Yield_F , F_MSY, reference points
# ---------------------------------------------------------

plot(Fref_vec, Yield_eq_vec_T / 10000, type="l", col="black", lwd=4,
     xlab="Fishing Mortality (F)", 
     ylab=expression(paste("Yield (x ", 10^4, " MT)")), 
     main="Yield_F",
     cex.lab=1.5, cex.axis=1.2, cex.main=1.8,
     ylim=c(0, max(Yield_eq_vec_T / 10000, na.rm=TRUE) * 1.2))

points(F_msy_T, MSY_T / 10000, pch=19, col="red", cex=2.5)

abline(v = F_msy_T, col="red", lty=2, lwd=2)
abline(h = MSY_T / 10000, col="blue", lty=2, lwd=2)

text(x = F_msy_T + 0.05, y = (MSY_T / 10000) * 0.4, 
     paste0("F_msy = ", round(F_msy_T, 2)), col="red", font=2, cex=1.3, pos=4)
text(x = max(Fref_vec) * 0.4, y = (MSY_T / 10000) * 1.05, 
     paste0("MSY = ", round(MSY_T / 10000, 2), " (x10^4 MT)"), col="red", font=2, cex=1.3)

#-----------------------------------------------------------------------------#
# 2. East sea
#-----------------------------------------------------------------------------#
Spawn_month_val <- 3
musig_logit_h_val <- c(0.6, 1.5)  # non-informative (0.64,0.58)
CV_yield_val <- 0.2
CV_cpue_val <- 0.1
region_settings[["East"]] <- list(
  init = list(
    log_Rec_fir = 18, logR0 = 20, logRec = rep(18, nyrs), 
    log_a50 = log(2.0), log_a95 = log(3.0), 
    logit_q = -5.0, 
    logit_h = 0.59,
    logF1 = log(0.2), logFt = rep(log(0.2), nyrs), 
    log_sig_logF = log(0.2), log_sig_logRec = log(0.4), logtheta = log(0.05) 
  ),
  low = list(
    "log_Rec_fir" = log(10000), "logR0" = log(10000), 
    "log_a50" = log(0.00001), "log_a95" = log(1.0), 
    "logit_q" = -30.0,
    "logit_h" = -5.0,
    "logF1" = log(0.001), "log_sig_logF" = log(0.01), "log_sig_logRec" = log(0.001), "logtheta" = log(0.0000001)
  ),
  up = list(
    "log_Rec_fir" = log(1e20), "logR0" = log(1e20), 
    "log_a50" = log(10.0), "log_a95" = log(15.0), 
    "logit_q" = 30.0, 
    "logit_h" = 10.0,
    "logF1" = log(5.0), "log_sig_logF" = log(2.0), "log_sig_logRec" = log(2.0), "logtheta" = log(100.0)
  )
)

# ==============================================================================
# Run model
# ==============================================================================
result_E <- run_tmb_region("East")

fit_E <- result_E$fit
fit_E

sdrep_E <- result_E$sdrep
sdrep_E

Ssummestim_E = summary(sdrep_E)
summestim_E

max.abs.grad_E = max(abs(sdrep_E$gradient.fixed))
max.abs.grad_E

# ==============================================================================
# Reference Points 
# ==============================================================================
reports_E <- result_E$reports

YPR_vec_E <- reports_E$YPR
SPR_ratio_E <- reports_E$SPR
Yield_eq_vec_E <- reports_E$Yield_F
MSY_E <- reports_E$MSY

F_msy_idx_E <- which.max(Yield_eq_vec_E)
F_msy_E <- Fref_vec[F_msy_idx_E]

# ==============================================================================
# [Plot 1] Yield & CPUE Fitting (East sea)
# ==============================================================================
par(mfrow=c(2,1), mar=c(5,5,2,3), oma=c(0,0,2,0))

# 1-1. Yield Fitting
obs_yield_E <- reports_E$yield
pred_yield_E <- reports_E$Yield_hat / 1000

plot(years, obs_yield_E, ylim=c(0, max(c(obs_yield_E, MSY_E), na.rm=TRUE) * 1.2), 
     xlab="Year", ylab=expression(paste("MT")), 
     main="Yield - East Sea", 
     cex.lab=1.5, cex.axis=1.5, cex.main=2.0, pch=1)
lines(years, pred_yield_E, col="red", lwd=3)

abline(h = MSY_E, col="blue", lty=2, lwd=2)

# 1-2. CPUE Fitting
obs_cpue_E <- reports_E$CPUE
pred_cpue_E <- reports_E$predcpue

plot(years, obs_cpue_E, ylim=c(0, max(obs_cpue_E) * 1.2), 
     xlab="Year", ylab="MT/Fishing effort", 
     main="CPUE - East Sea", 
     cex.lab=1.5, cex.axis=1.5, cex.main=2.0, pch=1)
lines(years, pred_cpue_E, col="red", lwd=3)


# ==============================================================================
# [Plot 2] Age Composition Fitting (East sea)
# ==============================================================================
grid_dim <- ceiling(sqrt(nyrs)) 
par(mfrow=c(grid_dim, grid_dim), oma=c(2,2,3,1), mar=c(2, 2, 1.5, 1), cex=1.0)

CAA_obs_E <- reports_E$CAA_obs
CAA_pred_E <- reports_E$CAA_pred

for(i in 1:nyrs) {
  obs_freq_E <- (CAA_obs_E[i,] / sum(CAA_obs_E[i,])) * SS_age_vec[i]
  pred_freq_E <- (CAA_pred_E[i,] / sum(CAA_pred_E[i,])) * SS_age_vec[i]
  
  y_max <- max(c(obs_freq_E, pred_freq_E), na.rm=TRUE) * 1.1
  if (is.na(y_max) | y_max == 0) y_max <- 1 
  
  plot(ages, obs_freq_E, type="h", ylim=c(0, y_max), main=paste(years[i]), 
       xlab="", ylab="", lwd=5, col="black") 
  lines(ages, pred_freq_E, col='red', lwd=3)
}
mtext("Age Composition - East Sea", side=3, outer=TRUE, cex=1.5, font=2)
mtext("Age", side=1, outer=TRUE, line=0.5, cex=1.2)
mtext("Frequency", side=2, outer=TRUE, line=0.5, cex=1.2)


# ==============================================================================
# [Plot 3] YPR, SPR & Equilibrium Yield Curves (East sea)
# ==============================================================================
par(mfrow=c(1,2), mar=c(5,5,3,2), oma=c(0,0,2,0))

# 3-1. YPR Curve
# plot(Fref_vec, YPR_vec_E * 1000, type="l", col="black", lwd=3,
#      xlab="Fishing mortality (F)", ylab="Yield Per Recruit", main="YPR - East Sea",
#      cex.lab=1.5, cex.axis=1.5, cex.main=1.8)
# 3-2. SPR Curve
plot(Fref_vec, SPR_ratio_E, type="l", col="black", lwd=3,
     xlab="Fishing mortality (F)", ylab="SPR Ratio", ylim=c(0, 1.0), main="SPR - East Sea",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.8)
abline(h = c(0.3, 0.4), col=c("red", "orange"), lty=2, lwd=2)
text(x = max(Fref_vec)*0.85, y = 0.45, "40% SPR", col="orange", font=2, cex=1.2)
text(x = max(Fref_vec)*0.85, y = 0.35, "30% SPR", col="red", font=2, cex=1.2)

# 3-3. Yield Curve
plot(Fref_vec, Yield_eq_vec_E / 10000, type="l", col="black", lwd=4,
     xlab="Fishing Mortality (F)", 
     ylab=expression(paste("Yield (x ", 10^4, " MT)")), 
     main="Yield Curve - East Sea",
     cex.lab=1.5, cex.axis=1.2, cex.main=1.8,
     ylim=c(0, max(Yield_eq_vec_E / 10000, na.rm=TRUE) * 1.2))

points(F_msy_E, MSY_E / 10000, pch=19, col="red", cex=2.5)
abline(v = F_msy_E, col="red", lty=2, lwd=2)
abline(h = MSY_E / 10000, col="blue", lty=2, lwd=2)

text(x = F_msy_E + 0.05, y = (MSY_E / 10000) * 0.4, 
     paste0("F_msy = ", round(F_msy_E, 2)), col="red", font=2, cex=1.3, pos=4)
text(x = max(Fref_vec) * 0.4, y = (MSY_E / 10000) * 1.05, 
     paste0("MSY = ", round(MSY_E / 10000, 2), " (x10^4 MT)"), col="red", font=2, cex=1.3)

#-----------------------------------------------------------------------------#
# 3. West sea
#-----------------------------------------------------------------------------#
CV_yield_val <- 0.1
CV_cpue_val <- 0.05
musig_logit_h_val <- c(0.64, 0.8)
Spawn_month_val <- 4
region_settings[["West"]] <- list(
  init = list(
    log_Rec_fir = 19, logR0 = 22, logRec = rep(18, nyrs), 
    log_a50 = log(2.0), log_a95 = log(3.0), 
    logit_q = -5.0, 
    logit_h = 0.59,
    logF1 = log(0.2), logFt = rep(log(0.2), nyrs), 
    log_sig_logF = log(0.2), log_sig_logRec = log(0.4), logtheta = log(0.05) 
  ),
  low = list(
    "log_Rec_fir" = log(10000), "logR0" = log(10000), 
    "log_a50" = log(0.00001), "log_a95" = log(1.0), 
    "logit_q" = -30.0,
    "logit_h" = -5.0,
    "logF1" = log(0.001), "log_sig_logF" = log(0.01), "log_sig_logRec" = log(0.001), "logtheta" = log(0.0000001)
  ),
  up = list(
    "log_Rec_fir" = log(1e20), "logR0" = log(1e20), 
    "log_a50" = log(10.0), "log_a95" = log(15.0), 
    "logit_q" = 30.0, 
    "logit_h" = 10.0,
    "logF1" = log(5.0), "log_sig_logF" = log(2.0), "log_sig_logRec" = log(2.0), "logtheta" = log(100.0)
  )
)

# ==============================================================================
# Run model
# ==============================================================================
result_W <- run_tmb_region("West")
result_W
fit_W <- result_W$fit
fit_W

sdrep_W <- result_W$sdrep
sdrep_W

summestim_W = summary(sdrep_W)
summestim_W

max.abs.grad_W = max(abs(sdrep_W$gradient.fixed))
max.abs.grad_W
# ==============================================================================
# Reference Points 
# ==============================================================================
reports_W <- result_W$reports
YPR_vec_W <- reports_W$YPR
SPR_ratio_W <- reports_W$SPR
Yield_eq_vec_W <- reports_W$Yield_F
MSY_W <- reports_W$MSY

F_msy_idx_W <- which.max(Yield_eq_vec_W)
F_msy_W <- Fref_vec[F_msy_idx_W]

# ==============================================================================
# [Plot 1] Yield & CPUE Fitting (West sea)
# ==============================================================================
par(mfrow=c(2,1), mar=c(5,5,2,3), oma=c(0,0,2,0))

# 1-1. Yield Fitting
obs_yield_W <- reports_W$yield
pred_yield_W <- reports_W$Yield_hat / 1000

plot(years, obs_yield_W, 
     xlab="Year", ylab=expression(paste("MT")), 
     main="Yield - West Sea", 
     cex.lab=1.5, cex.axis=1.5, cex.main=2.0, pch=1)
lines(years, pred_yield_W, col="red", lwd=3)

abline(h = MSY_W, col="blue", lty=2, lwd=2)

# 1-2. CPUE Fitting
obs_cpue_W <- reports_W$CPUE
pred_cpue_W <- reports_W$predcpue

plot(years, obs_cpue_W, ylim=c(0, max(pred_cpue_W) * 1.2), 
     xlab="Year", ylab="MT/Fishing effort", 
     main="CPUE - West Sea", 
     cex.lab=1.5, cex.axis=1.5, cex.main=2.0, pch=1)
lines(years, pred_cpue_W, col="red", lwd=3)


# ==============================================================================
# [Plot 2] Age Composition Fitting (West sea)
# ==============================================================================
grid_dim <- ceiling(sqrt(nyrs)) 
par(mfrow=c(grid_dim, grid_dim), oma=c(2,2,3,1), mar=c(2, 2, 1.5, 1), cex=1.0)

CAA_obs_W <- reports_W$CAA_obs
CAA_pred_W <- reports_W$CAA_pred

for(i in 1:nyrs) {
  obs_freq_W <- (CAA_obs_W[i,] / sum(CAA_obs_W[i,])) * SS_age_vec[i]
  pred_freq_W <- (CAA_pred_W[i,] / sum(CAA_pred_W[i,])) * SS_age_vec[i]
  
  y_max <- max(c(obs_freq_W, pred_freq_W), na.rm=TRUE) * 1.1
  if (is.na(y_max) | y_max == 0) y_max <- 1 
  
  plot(ages, obs_freq_W, type="h", ylim=c(0, y_max), main=paste(years[i]), 
       xlab="", ylab="", lwd=5, col="black") 
  lines(ages, pred_freq_W, col='red', lwd=3)
}
mtext("Age Composition - West Sea", side=3, outer=TRUE, cex=1.5, font=2)
mtext("Age", side=1, outer=TRUE, line=0.5, cex=1.2)
mtext("Frequency", side=2, outer=TRUE, line=0.5, cex=1.2)


# ==============================================================================
# [Plot 3] YPR, SPR & Equilibrium Yield Curves (West sea)
# ==============================================================================
par(mfrow=c(1,2), mar=c(5,5,3,2), oma=c(0,0,2,0))

# 3-1. YPR Curve
# plot(Fref_vec, YPR_vec_W * 1000, type="l", col="black", lwd=3,
#      xlab="Fishing mortality (F)", ylab="Yield Per Recruit", main="YPR - West Sea",
#      cex.lab=1.5, cex.axis=1.5, cex.main=1.8)
# 3-2. SPR Curve
plot(Fref_vec, SPR_ratio_W, type="l", col="black", lwd=3,
     xlab="Fishing mortality (F)", ylab="SPR Ratio", ylim=c(0, 1.0), main="SPR - West Sea",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.8)
abline(h = c(0.3, 0.4), col=c("red", "orange"), lty=2, lwd=2)
text(x = max(Fref_vec)*0.85, y = 0.45, "40% SPR", col="orange", font=2, cex=1.2)
text(x = max(Fref_vec)*0.85, y = 0.35, "30% SPR", col="red", font=2, cex=1.2)

# 3-3. Yield Curve
plot(Fref_vec, Yield_eq_vec_W / 10000, type="l", col="black", lwd=4,
     xlab="Fishing Mortality (F)", 
     ylab=expression(paste("Yield (x ", 10^4, " MT)")), 
     main="Yield Curve - West Sea",
     cex.lab=1.5, cex.axis=1.2, cex.main=1.8,
     ylim=c(0, max(Yield_eq_vec_W / 10000, na.rm=TRUE) * 1.2))

points(F_msy_W, MSY_W / 10000, pch=19, col="red", cex=2.5)
abline(v = F_msy_W, col="red", lty=2, lwd=2)
abline(h = MSY_W / 10000, col="blue", lty=2, lwd=2)

text(x = F_msy_W + 0.05, y = (MSY_W / 10000) * 0.4, 
     paste0("F_msy = ", round(F_msy_W, 2)), col="red", font=2, cex=1.3, pos=4)
text(x = max(Fref_vec) * 0.4, y = (MSY_W / 10000) * 1.05, 
     paste0("MSY = ", round(MSY_W / 10000, 2), " (x10^4 MT)"), col="red", font=2, cex=1.3)

#-----------------------------------------------------------------------------#
# 4. South sea
#-----------------------------------------------------------------------------#
Spawn_month_val <- 4
musig_logit_h_val <- c(0.6, 1.5)  # non-informative (0.64,0.58)
CV_yield_val <- 0.15
CV_cpue_val <- 0.15

region_settings[["South"]] <- list(
  init = list(
    log_Rec_fir = 19, logR0 = 20, logRec = rep(18, nyrs), 
    log_a50 = log(2.0), log_a95 = log(3.0), 
    logit_q = -5.0, 
    logit_h = 0.59,
    logF1 = log(0.2), logFt = rep(log(0.2), nyrs), 
    log_sig_logF = log(0.2), log_sig_logRec = log(0.4), logtheta = log(0.05) 
  ),
  low = list(
    "log_Rec_fir" = log(10000), "logR0" = log(10000), 
    "log_a50" = log(0.00001), "log_a95" = log(1.0), 
    "logit_q" = -30.0,
    "logit_h" = -5.0,
    "logF1" = log(0.001), "log_sig_logF" = log(0.01), "log_sig_logRec" = log(0.001), "logtheta" = log(0.0000001)
  ),
  up = list(
    "log_Rec_fir" = log(1e20), "logR0" = log(1e20), 
    "log_a50" = log(10.0), "log_a95" = log(15.0), 
    "logit_q" = 30.0, 
    "logit_h" = 10.0,
    "logF1" = log(8.0), "log_sig_logF" = log(2.0), "log_sig_logRec" = log(2.0), "logtheta" = log(100.0)
  )
)


# ==============================================================================
# Run model
# ==============================================================================
result_S <- run_tmb_region("South")

fit_S <- result_S$fit
fit_S

sdrep_S <- result_S$sdrep
sdrep_S

summestim_S = summary(sdrep_S)
summestim_S

max.abs.grad_S = max(abs(sdrep_S$gradient.fixed))
max.abs.grad_S
# ==============================================================================
# Reference Points 
# ==============================================================================
reports_S <- result_S$reports
YPR_vec_S <- reports_S$YPR
SPR_ratio_S <- reports_S$SPR
Yield_eq_vec_S <- reports_S$Yield_F
MSY_S <- reports_S$MSY

F_msy_idx_S <- which.max(Yield_eq_vec_S)
F_msy_S <- Fref_vec[F_msy_idx_S]

# ==============================================================================
# [Plot 1] Yield & CPUE Fitting (South sea)
# ==============================================================================
par(mfrow=c(1,2), mar=c(5,5,2,3), oma=c(0,0,2,0))

# 1-1. Yield Fitting
obs_yield_S <- reports_S$yield
pred_yield_S <- reports_S$Yield_hat / 1000

plot(years, obs_yield_S, ylim=c(0, max(c(obs_yield_S, MSY_S), na.rm=TRUE) * 1.2), 
     xlab="Year", ylab=expression(paste("MT")), 
     main="Yield - South Sea", 
     cex.lab=1.5, cex.axis=1.5, cex.main=2.0, pch=1)
lines(years, pred_yield_S, col="red", lwd=3)

abline(h = MSY_S, col="blue", lty=2, lwd=2)

# 1-2. CPUE Fitting
obs_cpue_S <- reports_S$CPUE
pred_cpue_S <- reports_S$predcpue

plot(years, obs_cpue_S, ylim=c(0, max(obs_cpue_S) * 1.2), 
     xlab="Year", ylab="MT/Fishing effort", 
     main="CPUE - South Sea", 
     cex.lab=1.5, cex.axis=1.5, cex.main=2.0, pch=1)
lines(years, pred_cpue_S, col="red", lwd=3)


# ==============================================================================
# [Plot 2] Age Composition Fitting (South sea)
# ==============================================================================
grid_dim <- ceiling(sqrt(nyrs)) 
par(mfrow=c(grid_dim, grid_dim), oma=c(2,2,3,1), mar=c(2, 2, 1.5, 1), cex=1.0)

CAA_obs_S <- reports_S$CAA_obs
CAA_pred_S <- reports_S$CAA_pred

for(i in 1:nyrs) {
  obs_freq_S <- (CAA_obs_S[i,] / sum(CAA_obs_S[i,])) * SS_age_vec[i]
  pred_freq_S <- (CAA_pred_S[i,] / sum(CAA_pred_S[i,])) * SS_age_vec[i]
  
  y_max <- max(c(obs_freq_S, pred_freq_S), na.rm=TRUE) * 1.1
  if (is.na(y_max) | y_max == 0) y_max <- 1 
  
  plot(ages, obs_freq_S, type="h", ylim=c(0, y_max), main=paste(years[i]), 
       xlab="", ylab="", lwd=5, col="black") 
  lines(ages, pred_freq_S, col='red', lwd=3)
}
mtext("Age Composition - South Sea", side=3, outer=TRUE, cex=1.5, font=2)
mtext("Age", side=1, outer=TRUE, line=0.5, cex=1.2)
mtext("Frequency", side=2, outer=TRUE, line=0.5, cex=1.2)


# ==============================================================================
# [Plot 3] YPR, SPR & Equilibrium Yield Curves (South sea)
# ==============================================================================
par(mfrow=c(1,2), mar=c(5,5,3,2), oma=c(0,0,2,0))

# 3-1. YPR Curve
plot(Fref_vec, YPR_vec_S * 1000, type="l", col="black", lwd=3,
     xlab="Fishing mortality (F)", ylab="Yield Per Recruit", main="YPR - South Sea",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.8)

# 3-2. SPR Curve
plot(Fref_vec, SPR_ratio_S, type="l", col="black", lwd=3,
     xlab="Fishing mortality (F)", ylab="SPR Ratio", ylim=c(0, 1.0), main="SPR - South Sea",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.8)
abline(h = c(0.3, 0.4), col=c("red", "orange"), lty=2, lwd=2)
text(x = max(Fref_vec)*0.85, y = 0.45, "40% SPR", col="orange", font=2, cex=1.2)
text(x = max(Fref_vec)*0.85, y = 0.35, "30% SPR", col="red", font=2, cex=1.2)

# 3-3. Yield Curve
plot(Fref_vec, Yield_eq_vec_S / 10000, type="l", col="black", lwd=4,
     xlab="Fishing Mortality (F)", 
     ylab=expression(paste("Yield (x ", 10^4, " MT)")), 
     main="Yield Curve - South Sea",
     cex.lab=1.5, cex.axis=1.2, cex.main=1.8,
     ylim=c(0, max(Yield_eq_vec_S / 10000, na.rm=TRUE) * 1.2))

points(F_msy_S, MSY_S / 10000, pch=19, col="red", cex=2.5)
abline(v = F_msy_S, col="red", lty=2, lwd=2)
abline(h = MSY_S / 10000, col="blue", lty=2, lwd=2)

text(x = F_msy_S + 0.05, y = (MSY_S / 10000) * 0.4, 
     paste0("F_msy = ", round(F_msy_S, 2)), col="red", font=2, cex=1.3, pos=4)
text(x = max(Fref_vec) * 0.4, y = (MSY_S / 10000) * 1.05, 
     paste0("MSY = ", round(MSY_S / 10000, 2), " (x10^4 MT)"), col="red", font=2, cex=1.3)


# Compare MSY
MSY_sum_regions <- MSY_E + MSY_W + MSY_S
MSY_E
MSY_W
MSY_S
MSY_sum_regions
MSY_T
# ==============================================================================
# [Plot 4] Biomass Trend Comparison (Total vs Sum of Regions)
# ==============================================================================
scale_factor <- 10000

biomass_T <- (reports_T$Biomass_temp / 1000) / scale_factor
biomass_E <- (reports_E$Biomass_temp / 1000) / scale_factor
biomass_W <- (reports_W$Biomass_temp / 1000) / scale_factor
biomass_S <- (reports_S$Biomass_temp / 1000) / scale_factor

biomass_sum_regions <- biomass_E + biomass_W + biomass_S
y_max_bio <- max(c(biomass_T, biomass_sum_regions), na.rm=TRUE) * 1.2

par(mfrow=c(1,1), mar=c(5,5,4,2)) 

plot(years, biomass_T, type="l", col="black", lwd=4,
     xlab="Year", 
     ylab=expression(paste("Biomass (x ", 10^4, " MT)")), 
     main="Biomass Trend Comparison",
     ylim=c(0, y_max_bio),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.8)

lines(years, biomass_sum_regions, col="red", lwd=4, lty=1)

legend("bottomright",   
       legend=c("Single unit", "Sum of regions (East+West+South)"), 
       col=c("black", "red"), 
       lwd=4, lty=1, cex=1.3, bty="n")
# ==============================================================================
# [Plot 5] Yield Curves Comparison (Total vs East vs West vs South)
# ==============================================================================
scale_y <- 10000

yield_T_scaled <- Yield_eq_vec_T / scale_y
yield_E_scaled <- Yield_eq_vec_E / scale_y
yield_W_scaled <- Yield_eq_vec_W / scale_y
yield_S_scaled <- Yield_eq_vec_S / scale_y

msy_T_scaled <- MSY_T / scale_y
msy_E_scaled <- MSY_E / scale_y
msy_W_scaled <- MSY_W / scale_y
msy_S_scaled <- MSY_S / scale_y

y_max_yield <- max(c(yield_T_scaled, yield_E_scaled, yield_W_scaled, yield_S_scaled), na.rm=TRUE) * 1.1

par(mfrow=c(1,1), mar=c(5,5,4,2))

plot(Fref_vec, yield_T_scaled, type="n", 
     xlab="Fishing Mortality (F)", 
     ylab=expression(paste("Equilibrium Yield (x ", 10^4, " MT)")), 
     main="Equilibrium Yield Curves by Region",
     ylim=c(0, y_max_yield),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.8)

lines(Fref_vec, yield_T_scaled, col="black", lwd=3)       
lines(Fref_vec, yield_E_scaled, col="blue", lwd=3)        
lines(Fref_vec, yield_W_scaled, col="red", lwd=3) 
lines(Fref_vec, yield_S_scaled, col="forestgreen", lwd=3)  




points(F_msy_T, msy_T_scaled, pch=19, col="black", cex=1)
points(F_msy_E, msy_E_scaled, pch=19, col="blue", cex=1)
points(F_msy_W, msy_W_scaled, pch=19, col="red", cex=1)
points(F_msy_S, msy_S_scaled, pch=19, col="forestgreen", cex=1)

abline(v = F_msy_T, col="black", lty=2, lwd=1.5)
abline(v = F_msy_E, col="blue", lty=2, lwd=1.5)
abline(v = F_msy_W, col="red", lty=2, lwd=1.5)
abline(v = F_msy_S, col="forestgreen", lty=2, lwd=1.5)

x_offset <- 0.08               
y_offset <- y_max_yield * 0.05 

text(F_msy_T + x_offset, msy_T_scaled + y_offset, sprintf("F=%.2f", F_msy_T), col="black", font=2, cex=1.2)
text(F_msy_E + x_offset, msy_E_scaled + y_offset, sprintf("F=%.2f", F_msy_E), col="blue", font=2, cex=1.2)
text(F_msy_W + x_offset, msy_W_scaled + y_offset, sprintf("F=%.2f", F_msy_W), col="red", font=2, cex=1.2)
text(F_msy_S + x_offset, msy_S_scaled + y_offset, sprintf("F=%.2f", F_msy_S), col="forestgreen", font=2, cex=1.2)
legend("topright", 
       legend=c(sprintf("Single unit"),
                sprintf("East"),
                sprintf("West"),
                sprintf("South")),
       col=c("black", "blue", "red", "forestgreen"), 
       lwd=c(4, 3, 3, 3), lty=1, cex=1.3, bty="n")

# ==============================================================================
# [Plot 5] Spawning Stock Biomass (SSB) Trend Comparison
# ==============================================================================
scale_ssb <- 10000

ssb_T <- summestim_T[rownames(summestim_T) == "SSB", "Estimate"] / scale_ssb

ssb_E <- summestim_E[rownames(summestim_E) == "SSB", "Estimate"] / scale_ssb
ssb_W <- summestim_W[rownames(summestim_W) == "SSB", "Estimate"] / scale_ssb
ssb_S <- summestim_S[rownames(summestim_S) == "SSB", "Estimate"] / scale_ssb

ssb_sum_regions <- ssb_E + ssb_W + ssb_S

y_max_ssb <- max(c(ssb_T, ssb_sum_regions), na.rm=TRUE) * 1.2

par(mfrow=c(1,1), mar=c(5,5,4,2)) 
plot(years, ssb_T, type="l", col="black", lwd=4,
     xlab="Year", 
     ylab=expression(paste("SSB (x ", 10^4, " MT)")), 
     main="Spawning Stock Biomass (SSB) Trend Comparison",
     ylim=c(0, y_max_ssb),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.8)

lines(years, ssb_sum_regions, col="red", lwd=4, lty=1)

legend("bottomright",   
       legend=c("Single unit", "Sum of Regions (East+West+South)"), 
       col=c("black", "red"), 
       lwd=4, lty=1, cex=1.3, bty="n")

# ==============================================================================
# [Plot 6] Recruitment Trend Comparison 
# ==============================================================================
scale_rec <- 1000000 

rec_T <- summestim_T[rownames(summestim_T) == "Rec", "Estimate"] / scale_rec
rec_E <- summestim_E[rownames(summestim_E) == "Rec", "Estimate"] / scale_rec
rec_W <- summestim_W[rownames(summestim_W) == "Rec", "Estimate"] / scale_rec
rec_S <- summestim_S[rownames(summestim_S) == "Rec", "Estimate"] / scale_rec

rec_sum_regions <- rec_E + rec_W + rec_S

y_max_rec <- max(c(rec_T, rec_sum_regions), na.rm=TRUE) * 1.2

par(mfrow=c(1,1), mar=c(5,5,4,2)) 
plot(years, rec_T, type="l", col="black", lwd=4,
     xlab="Year", 
     ylab=expression(paste("Recruitment (x ", 10^6, " fish)")), 
     main="Recruitment Trend Comparison",
     ylim=c(0, y_max_rec),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.8)

lines(years, rec_sum_regions, col="red", lwd=4, lty=1)

legend("bottomright",   
       legend=c("Single unit", "Sum of Regions (East+West+South)"), 
       col=c("black", "red"), 
       lwd=4, lty=1, cex=1.3, bty="n")

# ==============================================================================
# [Plot] Fishing Mortality (Ft) Trend Comparison (Total vs East vs West vs South)
# ==============================================================================
Ft_T <- exp(summestim_T[rownames(summestim_T) == "logFt", "Estimate"])
Ft_E <- exp(summestim_E[rownames(summestim_E) == "logFt", "Estimate"])
Ft_W <- exp(summestim_W[rownames(summestim_W) == "logFt", "Estimate"])
Ft_S <- exp(summestim_S[rownames(summestim_S) == "logFt", "Estimate"])

# Y축 최대값 설정 (가장 높은 F값 기준 20% 여유 공간)
y_max_F <- max(c(Ft_T, Ft_E, Ft_W, Ft_S), na.rm = TRUE) * 1.2

# 해역별 색상 세팅 (기존과 동일하게 통일)
col_single <- "black"
col_east   <- "#4575B4"
col_west   <- "#D73027"
col_south  <- "#74C476"

# 여백 설정 (좌측 Y축 이름 공간 확보)
par(mfrow = c(1, 1), mar = c(5, 6, 4, 2)) 

# -------------------------------------------------------------------
plot(years, Ft_T, type = "l", col = col_single, lwd = 4,
     xlab = "Year", 
     ylab = "Fishing Mortality (F)", 
#     main = "Fishing Mortality Trend by Region",
     ylim = c(0, y_max_F),
     las = 1, cex.lab = 1.6, cex.axis = 1.5, cex.main = 1.8)

lines(years, Ft_E, col = col_east,  lwd = 3)
lines(years, Ft_W, col = col_west,  lwd = 3)
lines(years, Ft_S, col = col_south, lwd = 3)

legend("topright", 
       legend = c("Total", "East Sea", "West Sea", "South Sea"), 
       col = c(col_single, col_east, col_west, col_south), 
       lwd = c(4, 3, 3, 3), lty = 1, cex = 1.3, bty = "n")

box(lwd = 1.5)
