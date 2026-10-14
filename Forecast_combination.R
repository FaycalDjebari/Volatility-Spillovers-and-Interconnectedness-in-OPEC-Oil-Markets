# Copyright 2026 Fayçal Djebari
# Licensed under the Apache License, Version 2.0
# ==============================================================================
# Forecast_combination.R: combined forecasts on the log-squared scale (Table 7)
#
# Runs in a FRESH R session. It needs two sets of files, both written by the
# existing scripts:
#   Results/network_fERR_T0_450.rds   saved by Network_based_log_ARCH_Framework.R
#   Standard MGARCH models (whitened)/fERR_<m>_T0_centred/...  written by
#                                              MGARCH_benchmarks.R
# If it is still open, the session of the main script can be used as well.
#
# Combination schemes, as in Mattera and Otto (2024, Section 4.3), applied
# country by country:
#   Average : equal weights
#   MinVar  : weights minimising the variance of the combined forecast error,
#             w = argmin w' Cov(e) w, w >= 0, sum(w) = 1
#   COLS    : constrained OLS, actual = a + F w + error, w >= 0, sum(w) = 1
# NO LOOK-AHEAD: at each origin h the weights use only the forecast errors of
# origins 1, ..., h-1 (expanding window); the first BURN origins use equal
# weights. Results are reported on all 52 months and on the months after the
# burn-in, where the schemes differ.
#
# Model sets:
#   Network      : the seven network log-ARCH specifications
#   Net + MGARCH : the seven network specifications, centred DCC and GO-GARCH
#   Net.GO + DCC : Net. GO-GARCH and the centred DCC-GARCH only
# CCC is not added because its forecasts coincide with DCC.
# ==============================================================================

# setwd("path/to/the/replication/folder")
T0_c   <- 450
BURN   <- 12

# ---- inputs -----------------------------------------------------------------
f_net <- file.path("Results", "network_fERR_T0_450.rds")
if (!file.exists(f_net)) stop("Run the main network script once: it now saves ", f_net)
net <- readRDS(f_net)
Y <- net$Y; eps_c <- net$eps_vals[["one_pct"]]

read_centred <- function(model, T0, eps_tag = "one_pct.1%") {
  f <- file.path("Standard MGARCH models (whitened)", sprintf("fERR_%s_T0_centred", model),
                 sprintf("fERR_%s_T0_%d_eps_%s.csv", model, T0, eps_tag))
  if (!file.exists(f)) stop("Missing centred MGARCH file: ", f)
  as.matrix(sapply(read.csv(f), as.numeric))
}
reconstruct_actual_log <- function(Y, T0, eps) {          # as in the main script
  m <- matrix(NA, ncol(Y) - T0, nrow(Y))
  for (i in seq_len(nrow(m))) m[i, ] <- log(pmax(Y[, T0 + i]^2, eps))
  colnames(m) <- rownames(Y); m
}
panel_robust_dm_test <- function(e1, e2, lag_trunc = NULL) {   # as in the main script
  d <- rowMeans(e1^2 - e2^2, na.rm = TRUE); d <- d[!is.na(d)]
  Te <- length(d)
  if (Te < 5) return(list(p.value = NA, mean_diff = NA, T_eff = Te))
  if (is.null(lag_trunc)) lag_trunc <- floor(1.2 * Te^(1/3))
  db <- mean(d); c0 <- d - db; hv <- mean(c0^2)
  if (lag_trunc > 0) for (k in 1:lag_trunc)
    hv <- hv + 2 * (1 - k/(lag_trunc+1)) * mean(c0[(k+1):Te] * c0[1:(Te-k)])
  if (!is.finite(hv) || hv <= 0) return(list(p.value = NA, mean_diff = db, T_eff = Te))
  s <- db / sqrt(hv / Te)
  list(p.value = 2*pnorm(-abs(s)), mean_diff = db, T_eff = Te)
}

net_names <- c("Euclidean", "Correlation", "Piccolo", "Net_CCC", "Net_DCC", "Net_GO", "Net_DY")
E_all <- net$fERR[net_names]
E_all$DCC_centred <- read_centred("dcc", T0_c)
E_all$GO_centred  <- read_centred("go",  T0_c)

actual <- reconstruct_actual_log(Y, T0_c, eps_c)              # H x N
H <- nrow(actual); N <- ncol(actual)
stopifnot(all(sapply(E_all, function(e) all(dim(e) == c(H, N)))))
F_all <- lapply(E_all, function(e) actual - e)                 # forecasts

# ---- least squares on the simplex: min ||A w - b||^2, w >= 0, sum(w) = 1 ----
proj_simplex <- function(v) {                                  # Duchi et al. (2008)
  u <- sort(v, decreasing = TRUE); cs <- cumsum(u)
  k <- max(which(u - (cs - 1) / seq_along(u) > 0))
  pmax(v - (cs[k] - 1) / k, 0)
}
simplex_ls <- function(A, b, iter = 5000, tol = 1e-12) {
  K <- ncol(A); w <- rep(1 / K, K)
  AtA <- crossprod(A); Atb <- crossprod(A, b)
  L <- max(eigen(AtA, symmetric = TRUE, only.values = TRUE)$values)
  if (!is.finite(L) || L <= 0) return(w)
  z <- w; tk <- 1
  for (it in 1:iter) {                                         # accelerated projected gradient
    w_new <- proj_simplex(z - (AtA %*% z - Atb) / L)
    tk1 <- (1 + sqrt(1 + 4 * tk^2)) / 2
    z <- w_new + ((tk - 1) / tk1) * (w_new - w)
    if (sum(abs(w_new - w)) < tol) { w <- w_new; break }
    w <- w_new; tk <- tk1
  }
  as.numeric(w)
}

combine <- function(set, scheme) {
  Fm <- F_all[set]; Em <- E_all[set]; K <- length(set)
  fc <- matrix(NA, H, N); wsum <- matrix(0, N, K, dimnames = list(colnames(actual), set))
  for (i in 1:N) for (h in 1:H) {
    f_now <- sapply(Fm, function(f) f[h, i])
    if (scheme == "Average" || h <= BURN) {
      w <- rep(1 / K, K); a <- 0
    } else {
      past <- 1:(h - 1)
      Ep <- sapply(Em, function(e) e[past, i])                 # (h-1) x K errors
      Fp <- sapply(Fm, function(f) f[past, i]); yp <- actual[past, i]
      if (scheme == "MinVar") {
        w <- simplex_ls(scale(Ep, scale = FALSE), rep(0, length(past))); a <- 0
      } else {                                                 # COLS with intercept
        Fc <- scale(Fp, scale = FALSE); yc <- yp - mean(yp)
        w <- simplex_ls(Fc, yc); a <- mean(yp) - sum(colMeans(Fp) * w)
      }
    }
    fc[h, i] <- a + sum(w * f_now)
    if (h == H) wsum[i, ] <- w
  }
  list(err = actual - fc, w_last = wsum)
}

sets <- list(Network = net_names,
             `Net + MGARCH` = c(net_names, "DCC_centred", "GO_centred"),
             `Net.GO + DCC` = c("Net_GO", "DCC_centred"))
schemes <- c("Average", "MinVar", "COLS")
comb <- list(); wlast <- list()
for (s in names(sets)) for (sc in schemes) {
  r <- combine(sets[[s]], sc)
  comb[[paste(s, sc, sep = " | ")]] <- r$err; wlast[[paste(s, sc, sep = " | ")]] <- r$w_last
}

# ---- evaluation --------------------------------------------------------------
loss_row <- function(e, rows = 1:H) {
  e <- e[rows, , drop = FALSE]
  c(RMSFE = sqrt(mean(e^2, na.rm = TRUE)), MAFE = mean(abs(e), na.rm = TRUE),
    Bias = mean(e, na.rm = TRUE))
}
all_models <- c(E_all, comb)
eval_tab <- function(rows) {
  t(sapply(names(all_models), function(nm) {
    e <- all_models[[nm]][rows, , drop = FALSE]
    p_go  <- if (nm == "Net_GO") NA else
      panel_robust_dm_test(e, E_all$Net_GO[rows, , drop = FALSE])$p.value
    p_dcc <- if (nm == "DCC_centred") NA else
      panel_robust_dm_test(e, E_all$DCC_centred[rows, , drop = FALSE])$p.value
    c(loss_row(all_models[[nm]], rows), DM_vs_NetGO = p_go, DM_vs_DCCc = p_dcc)
  }))
}
tab_all  <- eval_tab(1:H)
tab_post <- eval_tab((BURN + 1):H)
cat("\n=== TABLE C1: combined forecasts, log scale, T0 = 450, eps = 1%, all", H, "months ===\n")
cat("DM: panel-robust, squared errors; p-values against Net_GO and against the centred DCC.\n")
print(round(tab_all, 4))
cat("\n=== TABLE C2: same, months", BURN + 1, "to", H, "(after the equal-weight burn-in) ===\n")
print(round(tab_post, 4))

cat("\n--- Combination weights at the last origin (country by country) ---\n")
for (nm in names(wlast)) if (!grepl("Average", nm)) { cat("\n", nm, "\n"); print(round(wlast[[nm]], 3)) }

# ---- MCS over individual models and combinations (squared error, all months) ----
suppressMessages(library(MCS))
loss_mat <- sapply(all_models, function(e) rowMeans(e^2, na.rm = TRUE))
colnames(loss_mat) <- make.names(colnames(loss_mat))
set.seed(123)
cat("\n--- MCS (Tmax, B = 5000), squared error, eps = 1%, all models and combinations ---\n")
mcs_comb <- MCS::MCSprocedure(loss_mat, alpha = 0.05, B = 5000, statistic = "Tmax")

out_dir_c <- "Results"
if (!dir.exists(out_dir_c)) dir.create(out_dir_c, recursive = TRUE)
write.csv(tab_all,  file.path(out_dir_c, "Table7_combination_all_months.csv"))
write.csv(tab_post, file.path(out_dir_c, "Table7_combination_after_burnin.csv"))
cat("\nWritten to", normalizePath(out_dir_c), "\n")
