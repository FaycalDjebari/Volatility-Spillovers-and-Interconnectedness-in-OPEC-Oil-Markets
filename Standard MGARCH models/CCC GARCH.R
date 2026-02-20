# Copyright 2026 Fayçal Djebari
# Licensed under the Apache License, Version 2.0
# See the LICENSE file in the repository root for full license information.

# ====== Standard CCC-GARCH Forecasting for Multiple T0 and Epsilon ======


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
# ----------- 2. GARCH and CCC Specification -----------
spec_univ <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model     = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "std"
)

ccc_spec <- dccspec(
  uspec        = multispec(replicate(ncol(log_returns_xts), spec_univ)),
  dccOrder     = c(0, 0),  # CCC-GARCH specification
  distribution = "mvt"
)

# ----------- 3. Helper Function -----------
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

# ----------- Directories for saving Computational Metrics -----------
time_save_dir <- "Standard MGARCH models/Computation_Metrics_ccc"
if (!dir.exists(time_save_dir)) dir.create(time_save_dir, recursive = TRUE)

# ----------- 5. Forecasting for T0 ∈ {200, 250, 300, 350, 400, 450} with epsilon = min -----------
T0_vec <- c(200, 250, 300, 350, 400, 450)
eps_baseline <- eps_vals["min"]

results_T0 <- lapply(T0_vec, function(T0_current) {
  T_total     <- nrow(log_returns_xts)
  n_forecasts <- T_total - T0_current
  n_assets    <- ncol(log_returns_xts)
  
  cat(sprintf("\nProcessing T0: %d\n", T0_current))
  
  tmp_results <- pbapply::pblapply(1:n_forecasts, function(i) {
    t_train <- T0_current + i - 1
    
    gc(reset = TRUE) # Reset GC to measure peak memory for this specific step
    step_start_time <- Sys.time()
    
    # Core CCC logic
    ccc_roll <- dccfit(ccc_spec, data = log_returns_xts[1:t_train, ])
    covmat   <- rcov(ccc_roll)[,,1]
    predicted_vol <- sqrt(diag(covmat))
    actual_vol    <- log(pmax(log_returns_xts[t_train + 1, ]^2, eps_baseline))
    predicted_vol <- log(pmax(predicted_vol^2, eps_baseline))
    forecast_error <- actual_vol - predicted_vol
    
    # Stop timer and memory tracker
    step_end_time <- Sys.time()
    mem_used <- gc()
    peak_ram_mb <- sum(mem_used[, 6]) # Max memory used in MB
    step_duration <- as.numeric(difftime(step_end_time, step_start_time, units = "secs"))
    
    if (!any(is.na(forecast_error)) && length(forecast_error) == n_assets) {
      return(list(error = forecast_error, time = step_duration, ram = peak_ram_mb))
    } else {
      return(list(error = rep(NA_real_, n_assets), time = step_duration, ram = peak_ram_mb))
    }
  })
  
  # Unpack results
  fERR_ccc <- do.call(rbind, lapply(tmp_results, `[[`, "error"))
  time_vec <- do.call(c, lapply(tmp_results, `[[`, "time"))
  ram_vec  <- do.call(c, lapply(tmp_results, `[[`, "ram"))
  
  # Save forecast error matrix
  write.csv(fERR_ccc,
            file = sprintf("Standard MGARCH models/fERR_ccc_T0/fERR_ccc_T0_%d_eps_min.csv", T0_current),
            row.names = FALSE)
  
  # Save Time and RAM matrix
  metrics_df <- data.frame(Forecast_Step = 1:n_forecasts, Window_Size = T0_current + (1:n_forecasts) - 1, Time_Sec = time_vec, Peak_RAM_MB = ram_vec)
  write.csv(metrics_df, file = file.path(time_save_dir, sprintf("Metrics_ccc_T0_%d_eps_min.csv", T0_current)), row.names = FALSE)
  
  metrics <- calculate_metrics(rowMeans(fERR_ccc, na.rm = TRUE))
  data.frame(T0 = T0_current, RMSFE = metrics$RMSFE, MAFE = metrics$MAFE, Avg_Time_Sec = mean(time_vec), Avg_Peak_RAM_MB = mean(ram_vec))
})

print("Summary for T0 testing (eps = min):")
print(do.call(rbind, results_T0))



# ----------- 7. Setup Parameters for T0 = 450 (excluding "min") -----------
T0_eps_450  <- 450
T_total     <- nrow(log_returns_xts)
n_forecasts <- T_total - T0_eps_450
n_assets    <- ncol(log_returns_xts)

# Ensure the metrics directory variable exists in this session
time_save_dir <- "Standard MGARCH models/Computation_Metrics_ccc"
if (!dir.exists(time_save_dir)) dir.create(time_save_dir, recursive = TRUE)

eps_names_to_test_450 <- c("one_pct.1%", "small")

results_eps_450 <- lapply(eps_names_to_test_450, function(eps_name) {
  eps <- eps_vals[[eps_name]]
  cat(sprintf("\nProcessing epsilon: %s (Value: %e) for T0 = 450\n", eps_name, eps))
  
  tmp_results <- pbapply::pblapply(1:n_forecasts, function(i) {
    t_train <- T0_eps_450 + i - 1
    
    gc(reset = TRUE)
    step_start_time <- Sys.time()
    
    ccc_roll <- dccfit(ccc_spec, data = log_returns_xts[1:t_train, ])
    covmat   <- rcov(ccc_roll)[,,1]
    predicted_vol <- sqrt(diag(covmat))
    actual_vals   <- as.numeric(log_returns_xts[t_train + 1, ])
    actual_vol    <- log(pmax(actual_vals^2, eps))
    predicted_vol <- log(pmax(predicted_vol^2, eps))
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
  
  fERR_ccc <- do.call(rbind, lapply(tmp_results, `[[`, "error"))
  time_vec <- do.call(c, lapply(tmp_results, `[[`, "time"))
  ram_vec  <- do.call(c, lapply(tmp_results, `[[`, "ram"))
  
  save_dir <- "Standard MGARCH models/fERR_ccc_T0"
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  # Updated the map to expect "one_pct.1%" and save the file cleanly as "one_pct"
  suffix_map <- list("one_pct.1%" = "one_pct", "small" = "small")
  file_suffix <- suffix_map[[eps_name]]
  
  full_path <- file.path(save_dir, sprintf("fERR_ccc_T0_%d_eps_%s.csv", T0_eps_450, file_suffix))
  write.csv(fERR_ccc, file = full_path, row.names = FALSE)
  
  metrics_df <- data.frame(Forecast_Step = 1:n_forecasts, Window_Size = T0_eps_450 + (1:n_forecasts) - 1, Time_Sec = time_vec, Peak_RAM_MB = ram_vec)
  write.csv(metrics_df, file = file.path(time_save_dir, sprintf("Metrics_ccc_T0_450_eps_%s.csv", file_suffix)), row.names = FALSE)
  
  metrics <- calculate_metrics(rowMeans(fERR_ccc, na.rm = TRUE))
  data.frame(Epsilon = eps_name, RMSFE = metrics$RMSFE, MAFE = metrics$MAFE, Avg_Time_Sec = mean(time_vec), Avg_Peak_RAM_MB = mean(ram_vec))
})

print(do.call(rbind, results_eps_450))

