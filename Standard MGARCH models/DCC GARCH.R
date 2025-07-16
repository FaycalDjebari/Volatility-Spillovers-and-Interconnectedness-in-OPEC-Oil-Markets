# ====== Standard DCC-GARCH Forecasting for Multiple T0 and Epsilon ======


library(rugarch)
library(rmgarch)
library(pbapply)
library(xts)

set.seed(123)  # Ensure reproducibility

# ----------- 1. Load and Prepare Data -----------
setwd("/Users/faycal/Library/CloudStorage/Dropbox/R code for Network log-ARCH framework/")
log_returns <- read.csv("log_returns_data.csv")
log_returns_xts <- xts(log_returns[, -1], order.by = as.Date(log_returns[, 1]))

Y <- t(as.matrix(log_returns_xts))  # For compatibility with GMM code

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

# ----------- 5. Forecasting for T0 ∈ {200, 250, 300, 350} with epsilon = "min" -----------
T0_vec <- c(200, 250, 300, 350)
eps_baseline <- eps_vals["min"]

results_T0 <- lapply(T0_vec, function(T0_current) {
  T_total     <- nrow(log_returns_xts)
  n_forecasts <- T_total - T0_current
  n_assets    <- ncol(log_returns_xts)
  
  fERR_dcc <- matrix(NA_real_, nrow = n_forecasts, ncol = n_assets)
  
  tmp_results <- pbapply::pblapply(1:n_forecasts, function(i) {
    t_train <- T0_current + i - 1
    dcc_roll <- dccfit(dcc_spec, data = log_returns_xts[1:t_train, ])
    dcc_forecast <- dccforecast(dcc_roll, n.ahead = 1)
    
    predicted_vol <- as.numeric(sigma(dcc_forecast)[,,1])
    actual_vol    <- log(pmax(log_returns_xts[t_train + 1, ]^2, eps_baseline))
    predicted_vol <- log(pmax(predicted_vol^2, eps_baseline))
    
    forecast_error <- actual_vol - predicted_vol
    
    if (!any(is.na(forecast_error)) && length(forecast_error) == n_assets) {
      return(forecast_error)
    } else {
      return(rep(NA_real_, n_assets))
    }
  })
  
  fERR_dcc <- do.call(rbind, tmp_results)
  
  # Save forecast error matrix
  write.csv(fERR_dcc,
            file = sprintf("Standard MGARCH models/fERR_dcc_T0/fERR_dcc_T0_%d_eps_min.csv", T0_current),
            row.names = FALSE
  )
  
  metrics <- calculate_metrics(rowMeans(fERR_dcc, na.rm = TRUE))
  data.frame(T0 = T0_current, RMSFE = metrics$RMSFE, MAFE = metrics$MAFE)
})



# ----------- 6. Sensitivity to Epsilon at T0 = 300 -----------
T0_eps      <- 300
T_total     <- nrow(log_returns_xts)
n_forecasts <- T_total - T0_eps
n_assets    <- ncol(log_returns_xts)

results_eps <- lapply(names(eps_vals), function(eps_name) {
  eps <- eps_vals[eps_name]
  fERR_dcc <- matrix(NA_real_, nrow = n_forecasts, ncol = n_assets)
  
  tmp_results <- pbapply::pblapply(1:n_forecasts, function(i) {
    t_train <- T0_eps + i - 1
    dcc_roll <- dccfit(dcc_spec, data = log_returns_xts[1:t_train, ])
    dcc_forecast <- dccforecast(dcc_roll, n.ahead = 1)
    
    predicted_vol <- as.numeric(sigma(dcc_forecast)[,,1])
    actual_vol    <- log(pmax(log_returns_xts[t_train + 1, ]^2, eps))
    predicted_vol <- log(pmax(predicted_vol^2, eps))
    
    forecast_error <- actual_vol - predicted_vol
    
    if (!any(is.na(forecast_error)) && length(forecast_error) == n_assets) {
      return(forecast_error)
    } else {
      return(rep(NA_real_, n_assets))
    }
  })
  
  fERR_dcc <- do.call(rbind, tmp_results)
  
  # Save forecast error matrix
  write.csv(fERR_dcc,
            file = sprintf("Standard MGARCH models/fERR_dcc_T0/fERR_dcc_T0_300_eps_%s.csv", eps_name),
            row.names = FALSE
  )
  
  metrics <- calculate_metrics(rowMeans(fERR_dcc, na.rm = TRUE))
  data.frame(Epsilon = eps_name, RMSFE = metrics$RMSFE, MAFE = metrics$MAFE)
})






