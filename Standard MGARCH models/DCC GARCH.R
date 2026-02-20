# Copyright 2026 Fayçal Djebari
# Licensed under the Apache License, Version 2.0
# See the LICENSE file in the repository root for full license information.

# ====== Standard DCC-GARCH Forecasting for Multiple T0 and Epsilon ======

library(rugarch)
library(rmgarch)
library(pbapply)
library(xts)

set.seed(123)  

# ----------- 1. Load and Prepare Data -----------
setwd("/Users/faycal/Library/CloudStorage/Dropbox/R code for Network log-ARCH framework/")

log_returns <- read.csv("log_returns_data.csv")
log_returns_xts <- xts(log_returns[, -1], order.by = as.Date(log_returns[, 1]))

Y <- t(as.matrix(log_returns_xts))  

# ----------- 2. GARCH and DCC Specification -----------
spec_univ <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "std"
)

dcc_spec <- dccspec(
  uspec        = multispec(replicate(ncol(log_returns_xts), spec_univ)),
  dccOrder     = c(1, 1),
  distribution = "mvt"
)

# ----------- 3. Define Helper Function -----------
calculate_metrics <- function(fERR) {
  RMSFE <- sqrt(mean(fERR^2, na.rm = TRUE))
  MAFE  <- mean(abs(fERR), na.rm = TRUE)
  return(list(RMSFE = RMSFE, MAFE = MAFE))
}

# ----------- 4. Define Epsilon Values -----------
eps_vals <- c(
  min     = min(Y[Y != 0]^2, na.rm = TRUE),
  one_pct = quantile(Y[Y != 0]^2, 0.01, na.rm = TRUE),
  small   = 1e-6
)
cat("Epsilon values used:\n")
print(eps_vals)


time_save_dir_dcc <- "Standard MGARCH models/Computation_Metrics_dcc"
if (!dir.exists(time_save_dir_dcc)) dir.create(time_save_dir_dcc, recursive = TRUE)

# ----------- 5. Forecasting for T0 ∈ {200, 250, 300, 350, 400, 450} with epsilon = "min" -----------
T0_vec <- c(200, 250, 300, 350, 400, 450)
eps_baseline <- eps_vals["min"]

results_T0 <- lapply(T0_vec, function(T0_current) {
  T_total     <- nrow(log_returns_xts)
  n_forecasts <- T_total - T0_current
  n_assets    <- ncol(log_returns_xts)
  
  cat(sprintf("\nProcessing T0: %d (eps=min)\n", T0_current))
  
  tmp_results <- pbapply::pblapply(1:n_forecasts, function(i) {
    t_train <- T0_current + i - 1
    
    
    gc(reset = TRUE)
    step_start_time <- Sys.time()
    
    
    dcc_roll <- dccfit(dcc_spec, data = log_returns_xts[1:t_train, ])
    dcc_forecast <- dccforecast(dcc_roll, n.ahead = 1)
    
    predicted_vol <- as.numeric(sigma(dcc_forecast)[,,1])
    actual_vol    <- log(pmax(log_returns_xts[t_train + 1, ]^2, eps_baseline))
    predicted_vol <- log(pmax(predicted_vol^2, eps_baseline))
    
    forecast_error <- actual_vol - predicted_vol
    
    
    step_end_time <- Sys.time()
    mem_used <- gc()
    peak_ram_mb <- sum(mem_used[, 6])
    step_duration <- as.numeric(difftime(step_end_time, step_start_time, units = "secs"))
    
    if (!any(is.na(forecast_error)) && length(forecast_error) == n_assets) {
      return(list(error = forecast_error, time = step_duration, ram = peak_ram_mb))
    } else {
      return(list(error = rep(NA_real_, n_assets), time = step_duration, ram = peak_ram_mb))
    }
  })
  
  
  fERR_dcc <- do.call(rbind, lapply(tmp_results, `[[`, "error"))
  time_vec <- do.call(c, lapply(tmp_results, `[[`, "time"))
  ram_vec  <- do.call(c, lapply(tmp_results, `[[`, "ram"))
  
  
  write.csv(fERR_dcc,
            file = sprintf("Standard MGARCH models/fERR_dcc_T0/fERR_dcc_T0_%d_eps_min.csv", T0_current),
            row.names = FALSE)
  
  
  metrics_df <- data.frame(Forecast_Step = 1:n_forecasts, Window_Size = T0_current + (1:n_forecasts) - 1, Time_Sec = time_vec, Peak_RAM_MB = ram_vec)
  write.csv(metrics_df, file = file.path(time_save_dir_dcc, sprintf("Metrics_dcc_T0_%d_eps_min.csv", T0_current)), row.names = FALSE)
  
  metrics <- calculate_metrics(rowMeans(fERR_dcc, na.rm = TRUE))
  data.frame(T0 = T0_current, RMSFE = metrics$RMSFE, MAFE = metrics$MAFE, Avg_Time_Sec = mean(time_vec), Avg_Peak_RAM_MB = mean(ram_vec))
})

print("--- Summary for T0 varying (eps = min) ---")
print(do.call(rbind, results_T0))

# =======================================================
# SCRIPT 2: Standard DCC-GARCH Forecasts (T0=450)
# Missing Epsilons: one_pct, small
# =======================================================

# 1. SETUP: Define Data and Epsilons FIRST
# ----------------------------------------
T0_eps      <- 450
T_total     <- nrow(log_returns_xts)
n_forecasts <- T_total - T0_eps
n_assets    <- ncol(log_returns_xts)

# Calculate epsilon values from the data
vals_sq <- as.numeric(log_returns_xts^2)
vals_sq_nz <- vals_sq[vals_sq > 0]

eps_vals_2 <- list(
  "one_pct" = quantile(vals_sq_nz, 0.01),
  "small"   = 1e-6
)

eps_names_to_test <- c("one_pct", "small")

eps_file_suffix <- list(
  "one_pct" = "one_pct.1%", 
  "small"   = "small"
)


summary_eps_450 <- list()

# 2. RUN LOOP
# -----------
for (eps_name in eps_names_to_test) {
  
  eps_val <- eps_vals_2[[eps_name]]
  file_suffix <- eps_file_suffix[[eps_name]]
  
  cat(sprintf("\n--- Running DCC-GARCH for T0=%d, epsilon=%s ---\n", T0_eps, eps_name))
  
  results_list <- pbapply::pblapply(1:n_forecasts, function(i) {
    t_train <- T0_eps + i - 1
    
    # NEW: Start timer outside tryCatch to capture time even if it crashes
    gc(reset = TRUE)
    step_start_time <- Sys.time()
    
    # Wrapped your exact logic in the tryCatch
    res_error <- tryCatch({
      # 1. Fit DCC
      fit_roll <- dccfit(dcc_spec, data = log_returns_xts[1:t_train, ])
      
      # 2. Forecast
      fcast    <- dccforecast(fit_roll, n.ahead = 1)
      pred_vol <- as.numeric(sigma(fcast)[,,1])
      
      # 3. Calculate Error
      actual_vals   <- as.numeric(log_returns_xts[t_train + 1, ])
      actual_vol    <- log(pmax(actual_vals^2, eps_val))
      pred_vol_adj  <- log(pmax(pred_vol^2, eps_val))
      
      error <- actual_vol - pred_vol_adj
      
      if (!any(is.na(error)) && length(error) == n_assets) {
        error
      } else {
        rep(NA_real_, n_assets)
      }
    }, error = function(e) {
      rep(NA_real_, n_assets)
    })
    
    # NEW: Stop timer and calculate RAM
    step_end_time <- Sys.time()
    mem_used <- gc()
    peak_ram_mb <- sum(mem_used[, 6])
    step_duration <- as.numeric(difftime(step_end_time, step_start_time, units = "secs"))
    
    return(list(error = res_error, time = step_duration, ram = peak_ram_mb))
  })
  
  # NEW: Unpack lists
  fERR_dcc <- do.call(rbind, lapply(results_list, `[[`, "error"))
  time_vec <- do.call(c, lapply(results_list, `[[`, "time"))
  ram_vec  <- do.call(c, lapply(results_list, `[[`, "ram"))
  
  # 3. SAVE ERRORS
  # -------
  save_dir <- file.path("Standard MGARCH models", "fERR_dcc_T0")
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  file_name <- sprintf("fERR_dcc_T0_%d_eps_%s.csv", T0_eps, file_suffix)
  write.csv(fERR_dcc, file = file.path(save_dir, file_name), row.names = FALSE)
  cat(sprintf("Saved Errors: %s\n", file_name))
  
  # 4. SAVE METRICS
  # -------
  metrics_df <- data.frame(Forecast_Step = 1:n_forecasts, Window_Size = T0_eps + (1:n_forecasts) - 1, Time_Sec = time_vec, Peak_RAM_MB = ram_vec)
  metrics_file_name <- sprintf("Metrics_dcc_T0_%d_eps_%s.csv", T0_eps, file_suffix)
  write.csv(metrics_df, file = file.path(time_save_dir_dcc, metrics_file_name), row.names = FALSE)
  cat(sprintf("Saved Metrics: %s\n", metrics_file_name))
  
  # Calculate summary for final print
  metrics <- calculate_metrics(rowMeans(fERR_dcc, na.rm = TRUE))
  summary_eps_450[[eps_name]] <- data.frame(
    Epsilon = eps_name, 
    RMSFE = metrics$RMSFE, 
    MAFE = metrics$MAFE, 
    Avg_Time_Sec = mean(time_vec, na.rm=TRUE), 
    Avg_Peak_RAM_MB = mean(ram_vec, na.rm=TRUE)
  )
}

print("--- Summary for T0 = 450 (missing epsilons) ---")
print(do.call(rbind, summary_eps_450))




