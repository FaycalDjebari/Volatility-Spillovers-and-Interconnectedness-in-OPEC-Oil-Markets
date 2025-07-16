# ====== Standard GO-GARCH Forecasting for Multiple T0 and Epsilon ======

library(rmgarch)
library(rugarch)
library(pbapply)
library(xts)

set.seed(123)  # Reproducibility

# ----------- 1. Load and Prepare Data -----------
setwd("/Users/faycal/Library/CloudStorage/Dropbox/R code for Network log-ARCH framework/")

log_returns <- read.csv("log_returns_data.csv")
Y <- t(as.matrix(log_returns[, 2:7]))  # (n_assets x T)

start_date <- as.Date("1983-02-01")
date_index <- seq(from = start_date, by = "months", length.out = ncol(Y))
Y_xts <- xts(t(Y), order.by = date_index)

# ----------- 2. Specification -----------
spec_univ <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0)),
  distribution.model = "std"
)

spec_go <- gogarchspec(
  multispec(replicate(nrow(Y), spec_univ)),
  distribution.model = "mvnowm"
)

# ----------- 3. Helper Function -----------
calculate_metrics <- function(fERR) {
  RMSFE <- sqrt(mean(fERR^2, na.rm = TRUE))
  MAFE  <- mean(abs(fERR), na.rm = TRUE)
  return(list(RMSFE = RMSFE, MAFE = MAFE))
}

# ----------- 4. Define Epsilon Values and T0s -----------
eps_vals <- c(
  min     = min(Y[Y != 0]^2, na.rm = TRUE),
  one_pct = quantile(Y[Y != 0]^2, 0.01, na.rm = TRUE),
  small   = 1e-6
)

cat("Epsilon values used:\n")
print(eps_vals)

T0_vec <- c(200, 250, 300, 350)

# ----------- 5. Forecasting for T0_vec with epsilon = min -----------
all_results_T0 <- list()
eps_baseline <- eps_vals["min"]

for (T0_current in T0_vec) {
  cat("\n=== GO-GARCH: T0 =", T0_current, "| epsilon = min ===\n")
  T_total     <- ncol(Y)
  n_forecasts <- T_total - T0_current
  n_assets    <- nrow(Y)
  
  fERR_go <- matrix(NA_real_, nrow = n_forecasts, ncol = n_assets)
  
  tmp_results <- pbapply::pblapply(1:n_forecasts, function(i) {
    t_train <- T0_current + i - 1
    go_roll <- gogarchfit(spec_go, data = t(Y[, 1:t_train]))
    go_forecast <- gogarchforecast(go_roll, n.ahead = 1)
    
    predicted_vol <- as.numeric(sigma(go_forecast)[,,1])
    actual_vol    <- log(pmax(Y[, t_train + 1]^2, eps_baseline))
    predicted_vol <- log(pmax(predicted_vol^2, eps_baseline))
    
    forecast_error <- actual_vol - predicted_vol
    
    if (!any(is.na(forecast_error)) && length(forecast_error) == n_assets) {
      return(forecast_error)
    } else {
      return(rep(NA_real_, n_assets))
    }
  })
  
  fERR_go <- do.call(rbind, tmp_results)
  
  # Save forecast error matrix
  write.csv(fERR_go,
            file = sprintf("Standard MGARCH models/fERR_go_T0/fERR_go_T0_%d_eps_min.csv", T0_current),
            row.names = FALSE)
  
  metrics <- calculate_metrics(rowMeans(fERR_go, na.rm = TRUE))
  all_results_T0[[as.character(T0_current)]] <- data.frame(T0 = T0_current, RMSFE = metrics$RMSFE, MAFE = metrics$MAFE)
}

# ----------- 6. Sensitivity Analysis to epsilon at T0 = 300 (excluding "min") -----------
T0_eps      <- 300
T_total     <- ncol(Y)
n_forecasts <- T_total - T0_eps
n_assets    <- nrow(Y)

# Exclude "min" to avoid duplication
eps_names_to_test <- setdiff(names(eps_vals), "min")

results_eps <- lapply(eps_names_to_test, function(eps_name) {
  eps <- eps_vals[eps_name]
  
  fERR_go <- matrix(NA_real_, nrow = n_forecasts, ncol = n_assets)
  
  tmp_results <- pbapply::pblapply(1:n_forecasts, function(i) {
    t_train <- T0_eps + i - 1
    go_roll <- gogarchfit(spec_go, data = t(Y[, 1:t_train]))
    go_forecast <- gogarchforecast(go_roll, n.ahead = 1)
    
    predicted_vol <- as.numeric(sigma(go_forecast)[,,1])
    actual_vol    <- log(pmax(Y[, t_train + 1]^2, eps))
    predicted_vol <- log(pmax(predicted_vol^2, eps))
    
    forecast_error <- actual_vol - predicted_vol
    
    if (!any(is.na(forecast_error)) && length(forecast_error) == n_assets) {
      return(forecast_error)
    } else {
      return(rep(NA_real_, n_assets))
    }
  })
  
  fERR_go <- do.call(rbind, tmp_results)
  
  # Save forecast error matrix for each epsilon ≠ "min"
  write.csv(fERR_go,
            file = sprintf("Standard MGARCH models/fERR_go_T0/fERR_go_T0_300_eps_%s.csv", eps_name),
            row.names = FALSE)
  
  metrics <- calculate_metrics(rowMeans(fERR_go, na.rm = TRUE))
  data.frame(Epsilon = eps_name, RMSFE = metrics$RMSFE, MAFE = metrics$MAFE)
})





