# Copyright 2026 Fayçal Djebari
# Licensed under the Apache License, Version 2.0
# See the LICENSE file in the repository root for full license information.

# ==============================================================================
# Forecasting Oil Volatility through Network Models with GARCH-Informed
# Correlation Weights
# Djebari, Mehidi, Mazouz and Otto
#
# Standard multivariate GARCH benchmarks: CCC-, DCC- and GO-GARCH, re-estimated
# at every forecast origin on the pre-whitened returns. Run this script before
# Network_based_log_ARCH_Framework.R, which reads its output.
#
# For each model, training window T0 and zero adjustment eps it writes, to
# "Standard MGARCH models (whitened)/":
#   fERR_<m>_T0/           errors of the forecast ln(sigma^2)
#   fERR_<m>_T0_centred/   errors of the centred forecast ln(sigma^2) + c, where
#                          c is the training-window mean of
#                          ln max(e^2, eps) - ln(sigma^2), re-estimated at each
#                          origin (Section 6.3.1)
#   Computation_Metrics_<m>/  computing time per origin (Table D.1)
#
#  * whitening: one AR(1) + ZCA whitening per training window, as in the
#    network script
#  * eps: from the whitened panel of T0 = 450, as in the network script
#  * origins: t_train = T0 + i - 1, i = 1, ..., 502 - T0
#  * error: ln max(y^2, eps) - ln max(sigma_hat^2, eps)
#  * each origin is fitted once and serves every eps value
# ==============================================================================

library(rugarch)
library(rmgarch)
library(pbapply)
library(xts)

set.seed(123)
# setwd("path/to/the/replication/folder")

log_returns <- read.csv("log_returns_data.csv")
Y_raw <- t(as.matrix(log_returns[, 2:7]))


# ==============================================================================
# 1. Pre-whitening (as in the network script)
# ==============================================================================

# Pre-whitening, paper Eqs. (35)-(36):
#   1. AR(1) mean model per country, estimated on the training window
#   2. residuals u_t = y_t - mu - phi y_{t-1}
#   3. covariance of the training-window residuals
#   4. symmetric (ZCA) whitening matrix Sigma^{-1/2} = V D^{-1/2} V'
#   5. e_t = Sigma^{-1/2} u_t for all months, with the training-window matrix
fit_whitening <- function(Y_in, Y, mean_model = c("AR1", "constant"), eigen_tol = 1e-8) {

  mean_model <- match.arg(mean_model)
  Y_in <- as.matrix(Y_in)
  Y    <- as.matrix(Y)
  N    <- nrow(Y_in)
  T_in <- ncol(Y_in)
  T    <- ncol(Y)
  stopifnot(nrow(Y) == N)

  mu_hat  <- numeric(N)
  phi_hat <- numeric(N)
  if (mean_model == "AR1") {
    for (i in seq_len(N)) {
      y <- Y_in[i, 2:T_in]
      x <- Y_in[i, 1:(T_in - 1)]
      fit_i <- lm(y ~ x)
      mu_hat[i]  <- coef(fit_i)[1]
      phi_hat[i] <- coef(fit_i)[2]
    }
  } else {
    mu_hat  <- rowMeans(Y_in)
    phi_hat <- rep(0, N)
  }

  calc_mean_residuals <- function(Ydata) {
    TT <- ncol(Ydata)
    U <- matrix(NA_real_, nrow = N, ncol = TT)
    for (tt in 2:TT) U[, tt] <- Ydata[, tt] - mu_hat - phi_hat * Ydata[, tt - 1]
    U
  }
  U_in <- calc_mean_residuals(Y_in)
  U    <- calc_mean_residuals(Y)

  Sigma_hat <- cov(t(U_in[, 2:T_in, drop = FALSE]))
  ee <- eigen(Sigma_hat, symmetric = TRUE)
  d  <- ee$values
  V  <- ee$vectors
  d_safe <- pmax(d, eigen_tol * max(d))
  Sigma_inv_half <- V %*% diag(1 / sqrt(d_safe), N) %*% t(V)

  E <- matrix(NA_real_, nrow = N, ncol = T)
  E[, 2:T] <- Sigma_inv_half %*% U[, 2:T, drop = FALSE]
  E_in <- matrix(NA_real_, nrow = N, ncol = T_in)
  E_in[, 2:T_in] <- Sigma_inv_half %*% U_in[, 2:T_in, drop = FALSE]

  list(coefficients = list(mu = mu_hat, phi = phi_hat),
       Sigma = Sigma_hat,
       whitening_matrix = Sigma_inv_half,
       mean_residuals = U,
       residuals = E,
       residuals_in = E_in,
       correlation_before = cor(t(U_in[, 2:T_in, drop = FALSE])),
       correlation_after  = cor(t(E_in[, 2:T_in, drop = FALSE])),
       eigenvalues = d)
}

# Whitened panel for training window T0. T0 counts whitened months: whitened
# column k is log-return column k + 1, so the whitening for window T0 uses
# log-return columns 1:(T0 + 1), the data available at the end of the window.
.whiten_cache <- new.env()
Y_w <- function(T0, Y_logret = Y_raw, mean_model = "AR1") {
  key <- paste(T0, mean_model)
  if (!exists(key, envir = .whiten_cache, inherits = FALSE)) {
    fit <- fit_whitening(Y_in = Y_logret[, 1:(T0 + 1), drop = FALSE],
                         Y = Y_logret, mean_model = mean_model)
    E <- fit$residuals[, -1, drop = FALSE]
    dimnames(E) <- list(rownames(Y_logret), colnames(Y_logret)[-1])
    attr(E, "whitening_fit") <- fit
    assign(key, E, envir = .whiten_cache)
  }
  get(key, envir = .whiten_cache)
}

# First month of the whitened panel (the first log-return is February 1983)
DATA_START <- "1983-03-01"


# ==============================================================================
# 2. Rolling-window engine
# ==============================================================================

out_root   <- "Standard MGARCH models (whitened)"
eps_suffix <- c(min = "min", one_pct = "one_pct.1%", small = "small")

E450 <- Y_w(450)
eps_vals <- c(
  min     = min(E450[E450 != 0]^2),
  one_pct = as.numeric(quantile(E450[E450 != 0]^2, 0.01)),
  small   = 1e-6
)
print(eps_vals)

to_xts <- function(E) {
  xts(t(E), order.by = seq(as.Date(DATA_START), by = "months", length.out = ncol(E)))
}

# fit_fc(X) returns list(sigma_fc = one-step forecast of the conditional
# standard deviations, sigma_in = in-sample conditional standard deviations,
# resid = in-sample residuals, ram = size of the fitted object in MB)
run_rolling <- function(model, fit_fc, T0_vec = c(200, 250, 300, 350, 400, 450), all_eps_T0 = 450) {

  dir_raw <- file.path(out_root, sprintf("fERR_%s_T0", model))
  dir_cen <- file.path(out_root, sprintf("fERR_%s_T0_centred", model))
  dir_met <- file.path(out_root, sprintf("Computation_Metrics_%s", model))
  for (d in c(dir_raw, dir_cen, dir_met)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

  summary_rows <- list()
  for (T0 in T0_vec) {
    E <- Y_w(T0)
    attr(E, "whitening_fit") <- NULL
    X <- to_xts(E)
    N <- nrow(E); n_fc <- ncol(E) - T0
    eps_set <- if (T0 == all_eps_T0) names(eps_vals) else "min"
    cat(sprintf("\n=== %s-GARCH: T0 = %d, %d origins ===\n", toupper(model), T0, n_fc))

    steps <- pbapply::pblapply(seq_len(n_fc), function(i) {
      t_train <- T0 + i - 1
      t_start <- Sys.time()
      out <- tryCatch(fit_fc(X[1:t_train, ]), error = function(e) NULL)
      dur <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
      actual <- as.numeric(E[, t_train + 1])
      errs <- lapply(eps_set, function(k) {
        ep <- eps_vals[[k]]
        if (is.null(out) || any(!is.finite(out$sigma_fc)) || length(out$sigma_fc) != N)
          return(list(raw = rep(NA_real_, N), cen = rep(NA_real_, N)))
        raw <- log(pmax(actual^2, ep)) - log(pmax(out$sigma_fc^2, ep))
        c_k <- colMeans(log(pmax(out$resid^2, ep)) - log(out$sigma_in^2), na.rm = TRUE)
        list(raw = raw, cen = raw - c_k)
      })
      names(errs) <- eps_set
      list(errs = errs, time = dur, ram = if (is.null(out)) NA_real_ else out$ram)
    })

    time_vec <- vapply(steps, `[[`, numeric(1), "time")
    ram_vec  <- vapply(steps, `[[`, numeric(1), "ram")
    for (k in eps_set) {
      f_raw <- do.call(rbind, lapply(steps, function(s) s$errs[[k]]$raw))
      f_cen <- do.call(rbind, lapply(steps, function(s) s$errs[[k]]$cen))
      colnames(f_raw) <- colnames(f_cen) <- rownames(E)
      fname <- sprintf("fERR_%s_T0_%d_eps_%s.csv", model, T0, eps_suffix[[k]])
      write.csv(f_raw, file.path(dir_raw, fname), row.names = FALSE)
      write.csv(f_cen, file.path(dir_cen, fname), row.names = FALSE)
      write.csv(data.frame(Forecast_Step = seq_len(n_fc), Window_Size = T0 + seq_len(n_fc) - 1,
                           Time_Sec = time_vec, Peak_RAM_MB = ram_vec),
                file.path(dir_met, sprintf("Metrics_%s_T0_%d_eps_%s.csv", model, T0, eps_suffix[[k]])),
                row.names = FALSE)
      summary_rows[[paste(T0, k)]] <- data.frame(
        T0 = T0, eps = k, n_failed = sum(!complete.cases(f_raw)),
        RMSFE = sqrt(mean(f_raw^2, na.rm = TRUE)), RMSFE_centred = sqrt(mean(f_cen^2, na.rm = TRUE)),
        Bias = mean(f_raw, na.rm = TRUE), Bias_centred = mean(f_cen, na.rm = TRUE),
        Avg_Time_Sec = mean(time_vec, na.rm = TRUE))
    }
  }
  res <- do.call(rbind, summary_rows); rownames(res) <- NULL
  print(res, digits = 4)
  invisible(res)
}


# ==============================================================================
# 3. Specifications and one-step forecasts
# ==============================================================================

spec_univ <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "std"
)
N_assets <- nrow(Y_raw)
ccc_spec <- dccspec(uspec = multispec(replicate(N_assets, spec_univ)), dccOrder = c(0, 0), distribution = "mvt")
dcc_spec <- dccspec(uspec = multispec(replicate(N_assets, spec_univ)), dccOrder = c(1, 1), distribution = "mvt")
spec_go  <- gogarchspec(multispec(replicate(N_assets, spec_univ)), distribution.model = "mvnorm")

# CCC: dccforecast() does not accept dccOrder = c(0, 0). With a constant
# correlation matrix and a marginal target, the forecast is the univariate
# recursion h_{i,T+1} = omega_i + alpha_i e_{i,T}^2 + beta_i h_{i,T}.
ccc_forecast_sigma <- function(fit, n_assets, asset_names) {
  cf  <- coef(fit)
  sig <- sigma(fit)
  res <- residuals(fit)
  Tt  <- nrow(sig)
  vapply(seq_len(n_assets), function(k) {
    nm  <- asset_names[k]
    sqrt(unname(cf[paste0("[", nm, "].omega")] +
                cf[paste0("[", nm, "].alpha1")] * as.numeric(res[Tt, k])^2 +
                cf[paste0("[", nm, "].beta1")]  * as.numeric(sig[Tt, k])^2))
  }, numeric(1))
}

fit_fc_ccc <- function(X) {
  fit <- dccfit(ccc_spec, data = X)
  list(sigma_fc = ccc_forecast_sigma(fit, ncol(X), colnames(X)),
       sigma_in = as.matrix(sigma(fit)),
       resid    = as.matrix(residuals(fit)),
       ram      = as.numeric(object.size(fit)) / 1024^2)
}

fit_fc_dcc <- function(X) {
  fit <- dccfit(dcc_spec, data = X)
  fc  <- dccforecast(fit, n.ahead = 1)
  list(sigma_fc = as.numeric(sigma(fc)[, , 1]),
       sigma_in = as.matrix(sigma(fit)),
       resid    = as.matrix(residuals(fit)),
       ram      = as.numeric(object.size(fit)) / 1024^2)
}

# GO: the return covariance comes from rcov(); sigma() of a GO-GARCH object
# returns the volatilities of the latent factors.
fit_fc_go <- function(X) {
  fit <- gogarchfit(spec_go, data = X)
  fc  <- gogarchforecast(fit, n.ahead = 1)
  cv  <- rcov(fc); if (is.list(cv)) cv <- cv[[1]]
  ci  <- rcov(fit)
  res <- tryCatch(as.matrix(residuals(fit)), error = function(e) as.matrix(X))
  list(sigma_fc = sqrt(diag(as.matrix(cv[, , 1]))),
       sigma_in = sqrt(t(apply(ci, 3, diag))),
       resid    = res,
       ram      = as.numeric(object.size(fit)) / 1024^2)
}


# ==============================================================================
# 4. Run
# ==============================================================================

res_ccc <- run_rolling("ccc", fit_fc_ccc)
res_dcc <- run_rolling("dcc", fit_fc_dcc)
res_go  <- run_rolling("go",  fit_fc_go)
