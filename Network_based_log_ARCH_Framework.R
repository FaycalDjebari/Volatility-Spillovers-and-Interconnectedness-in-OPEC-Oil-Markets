# Copyright 2026 Fayçal Djebari
# Licensed under the Apache License, Version 2.0
# See the LICENSE file in the repository root for full license information.

# ==============================================================================
# Forecasting Oil Volatility through Network Models with GARCH-Informed
# Correlation Weights
# Djebari, Mehidi, Mazouz and Otto
#
# Replication script for the network log-ARCH models.
#
#  1. Packages and input files
#  2. GMM estimator of the network log-ARCH model, forecast intercept (Eq. 17)
#  3. Distance measures and pre-whitening (Section 5.3)
#  4. Data and descriptive statistics (Tables 1-2)
#  5. Full-sample GARCH models (Tables 3, 4, C.1) and weight matrices
#  6. Weight matrices for each training window
#  7. Network figures (Figures 2, 3, F.1, F.2)
#  8. GMM estimates at T0 = 450 (Table 5)
#  9. Rolling one-step forecasts (Tables 6 and E.1)
# 10. Forecast evaluation: bootstrap intervals, DM tests, MCS (Tables 6, E.1, E.2)
# 11. Weight stability, subperiods and structural breaks (Table 8)
# 12. Economic value (Table 9)
# 13. Simulation: estimation cost (Tables D.1 and D.2)
# 14. Simulation: forecast accuracy (Table D.3)
# 15. Extended panel with Iraq and Kuwait (Section 6.3.7, Appendix G)
#
# Inputs
#  * log_returns_data.csv    monthly log-returns of the six producers
#  * log_returns_data_8.csv  the same with Iraq and Kuwait, from 1998
#  * "Standard MGARCH models (whitened)/": forecast errors and timings of the
#    CCC-, DCC- and GO-GARCH benchmarks, written by MGARCH_benchmarks.R
#  * W_go_cached_logret_v2*.rds: GO-GARCH fits, saved on the first run and
#    reused afterwards. The ICA step of gogarchfit() is not fully controlled
#    by set.seed(), so a new run can give slightly different GO-GARCH weights.
#
# Conventions
#  * The volatility proxy is ln(max(y^2, eps)). The same eps is used for
#    estimation and evaluation. Baseline: eps at the 1% quantile of the
#    non-zero squared whitened values; robustness: the smallest non-zero
#    value and 1e-6 (Table E.1).
#  * Weight matrices are built from the log-returns. All volatility models are
#    estimated and evaluated on the pre-whitened returns.
#  * CCC- and DCC-GARCH give identical marginal variance forecasts and are
#    reported as one benchmark.
# ==============================================================================


# ==============================================================================
# 1. Packages and input files
# ==============================================================================

library(MASS)
library(rugarch)
library(rmgarch)
library(xts)
library(vars)
library(MCS)
library(strucchange)
library(igraph)
library(wesanderson)
library(classInt)
library(R.utils)

# setwd("path/to/the/replication/folder")


log_returns <- read.csv("log_returns_data.csv")

out_dir <- "Results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# ==============================================================================
# 2. GMM estimator of the network log-ARCH model
# ==============================================================================

GMM_SDPD_2SLS_ARCH_ind_timelags <- function(Y, X, W, info){
  
  ksy = info$ksy  # spatial expansion order of y
  ksx = info$ksx  # spatial expansion order of x
  
  if(is.null(Y)){
    stop("Y is missing")
  }
  if(is.null(W)){
    stop("W is missing")
  }
  if(is.null(info)){
    stop("info is missing")
  }
  
  dimW <- dim(W)
  n    <- dimW[1]
  if(length(dimW) == 2){
    p <- 1
    new.W <- array(, dim = c(n,n,p))
    new.W[,,1] <- W
    W <- new.W
  } else {
    p <- dimW[3]
  }
  
  if(length(dimW) > 3 |  dimW[1] !=  dimW[2] | dimW[2] !=  n){
    stop("Check W matrix (W must be of dimension n x n x p)")
  }
  
  if(dim(Y)[1] != n | length(dim(Y)) != 2){
    stop("Y should be a matrix of dimension n x T")
  }
  
  t <- dim(Y)[2] - 1
  
  
  
  s  <- t - 1
  nt <- n * t
  ns <- n * s
  
  yt   <- as.vector(Y[, 2:(t+1)])
  ytlv <- as.vector(Y[, 1:(t)]) # ytl vector
  ytl  <- array(0, dim = c(nt,n)) # ytl matrix to get individual temporal coefficients
  ysl  <- array(0, dim = c(nt,p))
  ystl <- array(0, dim = c(nt,p))
  
  for(i in 1:n){
    ytl[seq(1, nt, by = n) + i - 1, i] <- Y[i, 1:(t)]
  }
  
  for(i in 1:t){
    for(j in 1:p){
      ysl[(1+(i-1)*n):(i*n), j]  <- W[,,j] %*% yt[(1+(i-1)*n):(i*n)];
      ystl[(1+(i-1)*n):(i*n), j] <- W[,,j] %*% ytlv[(1+(i-1)*n):(i*n)];
    }
  }
  
  if(info$stl + info$tl == 2){
    xw <- cbind(ytl, ystl);
  } else if (info$stl + info$tl == 1){
    if(info$stl == 1){
      xw <- ystl
    } else {
      xw <- ytl
    }
  } else if (info$stl + info$tl == 0){
    stop("No spatial and no temporal lag given")
  } else {
    stop("Double-Check stl & tl # in Info structure")
  }
  
  X_stacked <- NULL
  W1hx <- NULL
  MHX <- NULL
  if(!is.null(X)){
    for(i in 1:t){
      X_stacked <- rbind(X_stacked, X[,i+1,])
    }
    
    if(dim(X)[3] == 1){
      X_stacked <- array(X_stacked, dim = c(length(X_stacked), 1))
    }
  }
  
  xs  <- X_stacked;
  xt  <- cbind(xw, xs);
  zt  <- cbind(ysl, xt);
  
  kz  <- dim(zt)[2]
  kx  <- dim(xt)[2]
  kxs <- dim(xs)[2]
  kxw <- dim(xw)[2]
  
  c <- sqrt((t-(1:s))/(t-(1:s)+1)); 
  
  F <- diag(t)
  F <- F[, 1:(t-1)]
  for(i in 1:(t-1)){
    F[(i+1):t, i] <- -1/(t-i);
    F[, i]        <- c[i] * F[,i];
  }
  
  hyt  <- array(yt, dim = c(n,t));
  hyt  <- hyt %*% F;
  hyt  <- as.vector(hyt);
  hytlv <- array(ytlv, dim = c(n,t));
  hytlv <- hytlv %*% F;
  hytlv <- as.vector(hytlv);
  
  hytl  <- array(0, dim = c(ns,n)) # hytl matrix to get individual temporal coefficients
  for(i in 1:n){
    hytl[seq(1, ns, by = n) + i - 1, i] <- array(hytlv, dim = c(n,s))[i, ]
  }
  
  hysl <- array(ysl, dim = c(n,t,p));
  hysltemp <- array(0, dim = c(n,t-1,p));
  for(i in 1:p){
    hysltemp[,,i] <- hysl[,,i] %*% F;
  }
  hysl <- array(hysltemp, dim = c(ns,p));
  
  hystl <- array(ystl, dim = c(n,t,p));
  hystltemp <- array(0, dim = c(n,t-1,p));
  for(i in 1:p){
    hystltemp[,,i] <- hystl[,,i] %*% F;
  }
  hystl <- array(hystltemp, dim = c(ns,p));
  
  if(!is.null(X_stacked)){
    kx <- dim(X_stacked)[2]
    hx <- array(X_stacked, dim = c(n,t,kx));
    hxtemp <- array(0, dim = c(n,t-1,kx));
    for(i in 1:kx){
      hxtemp[,,i] <- hx[,,i] %*% F;
    }
    hx <- array(hxtemp, dim = c(ns,kx));
  } else {
    hx <- NULL
  }
  
  if(info$stl + info$tl == 2){
    hxw <- cbind(hytl, hystl);
  } else if(info$stl + info$tl == 1){
    if(info$stl == 1){
      hxw <- hystl
    } else {
      hxw <- hytl
    } 
  } else if(info$stl + info$tl == 0){
    stop("no spatial and temporal lag given, check suitability")
  } else {
    stop("Doube-check info$stl and info$tl")
  }
  
  pyid     <- array(0, dim = c(ksy,2));
  pa       <- 1
  pb       <- p
  pyid[1,] <- c(pa, pb);
  for(k in 2:ksy){
    pa <- pa + p^(k-1)
    pb <- pb + p^k
    pyid[k,] <- c(pa, pb);
  }
  
  WY <- array(0, dim = c(nt, pyid[ksy,2]));
  WY[, 1:p] = ystl
  for(i in 1:t){
    for(k in 1:(ksy-1)){
      for(j in 1:p){
        WY[(1+(i-1)*n):(i*n), (pyid[k,2] + 1 + (j-1)*p^k):(pyid[k,2]+j*p^k)] <- W[,,j] %*% WY[(1+(i-1)*n):(i*n), pyid[k,1]:pyid[k,2]];
      }
    }
  }
  
  if(!is.null(X_stacked)){
    kx <- dim(X_stacked)[2]
    W1hx <- array(0, dim = c(ns,kx*p))
    
    for(i in 1:s){
      for(j in 1:p){
        W1hx[(1+(i-1)*n):(i*n), (1+(j-1)*kx):(j*kx)] <- W[,,j] %*% hx[(1+(i-1)*n):(i*n),];
      }
    }
    
    pxid <- array(0, dim = c(ksx,2))
    pa       <- 1
    pb       <- p*kx
    pxid[1,] <- c(pa, pb);
    for(k in 2:ksx){
      pa <- pa + p^(k-1)*kx
      pb <- pb + p^k*kx
      pxid[k,] <- c(pa, pb);
    }
    
    WHX <- array(0, dim = c(ns, pxid[ksx,2]*kx))
    WHX[,1:(p*kx)] <- W1hx;
    for(i in 1:s){
      for(k in 1:(ksx-1)){
        for(j in 1:p){
          WHX[(1+(i-1)*n):(i*n), (pxid[k,2]+1+(j-1)*p^k*kx):(pxid[k,2]+j*p^k*kx)] <- W[,,j] %*% WHX[(1+(i-1)*n):(i*n), pxid[k,1]:pxid[k,2]];
        }
      }
    }
  }
  
  
  ## Following is the IV without interaction term
  
  pyid <- array(0, dim = c(ksy,2));
  pa   <- 1;
  pb   <- p;
  pyid[1,] <- c(pa, pb);
  for(k in 2:ksy){
    pa <- pa + p
    pb <- pb + p
    pyid[k,] <- c(pa, pb)
  }
  
  MY <- array(0, dim = c(nt, pyid[ksy,2]))
  MY[,1:p] <- ystl;
  
  for(i in 1:t){
    for(k in 1:(ksy-1)){
      for(j in 1:p){
        MY[(1+(i-1)*n):(i*n), (pyid[k,2]+j):(pyid[k,2]+j)] <- W[,,j] %*% MY[(1+(i-1)*n):(i*n), (pyid[k,1]-1+j):(pyid[k,1]-1+j)];
      }
    }
  }
  
  if(!is.null(X_stacked)){
    
    kx <- dim(X_stacked)[2]
    pxid <- array(0, dim = c(ksx,2))
    pa <- 1
    pb <- p*kx
    pxid[1,] <- c(pa, pb)
    for(k in 2:ksx){
      pa <- pa + p*kx
      pb <- pb + p*kx
      pxid[k, ] <- c(pa, pb)
    }
    
    MHX <- array(0, dim = c(ns,pyid[ksx,2]*kx));
    MHX[,1:(p*kx)] <- W1hx;
    for(i in 1:s){
      for(k in 1:(ksx-1)){
        for(j in 1:p){
          MHX[(1+(i-1)*n):(i*n),(pxid[k,2]+1+(j-1)*kx):(pxid[k,2]+j*kx)] <- W[,,j] %*% MHX[(1+(i-1)*n):(i*n),(pxid[k,1]+(j-1)*kx):(pxid[k,1]-1+j*kx)];
        }
      }
    }
    
  }
  
  Qw     <- cbind(ytl, WY);
  Qw_alt <- cbind(ytl, MY);
  
  if(ksx == 0){
    Qs     <- NULL
    Qs_alt <- NULL
    qs     <- 0
    qs_alt <- 0
  } else {
    Qs     <- cbind(hx, W1hx)
    Qs_alt <- cbind(hx, MHX)
    qs     <- dim(Qs)[2]
    qs_slt <- dim(Qs_alt)[2]
  }
  
  Qw     <- Qw[1:ns,]
  Qw_alt <- Qw_alt[1:ns,]
  qw     <- dim(Qw)[2]
  qw_alt <- dim(Qw_alt)[2]
  
  Q      <- cbind(Qw, Qs)
  Q_alt  <- cbind(Qw_alt, Qs_alt)
  kq     <- dim(Q)[2]
  kq_alt <- dim(Q_alt)[2]
  
  hxs <- hx
  hxt <- cbind(hxw, hxs)
  hzt <- cbind(hysl, hxt)
  
  Qhz <- array(0, dim = c(kq, kz))
  QQ  <- array(0, dim = c(kq, kq))
  Qhy <- array(0, dim = c(kq, 1))
  
  Qhz_alt <- array(0, dim = c(kq_alt, kz))
  QQ_alt  <- array(0, dim = c(kq_alt, kq_alt))
  Qhy_alt <- array(0, dim = c(kq_alt, 1))
  
  Jn <- diag(n) - array(1/n, dim = c(n,n));
  
  for(i in 1:s){
    if(info$ted == 1){
      
      Qhz <- Qhz + t(Q[(1+(i-1)*n):(i*n),]) %*% Jn %*% hzt[(1+(i-1)*n):(i*n),];
      QQ  <- QQ  + t(Q[(1+(i-1)*n):(i*n),]) %*% Jn %*% Q[(1+(i-1)*n):(i*n),];
      Qhy <- Qhy + t(Q[(1+(i-1)*n):(i*n),]) %*% Jn %*% hyt[(1+(i-1)*n):(i*n)];
      
      Qhz_alt <- Qhz_alt + t(Q_alt[(1+(i-1)*n):(i*n),]) %*% Jn %*% hzt[(1+(i-1)*n):(i*n),];
      QQ_alt  <- QQ_alt  + t(Q_alt[(1+(i-1)*n):(i*n),]) %*% Jn %*% Q_alt[(1+(i-1)*n):(i*n),];
      Qhy_alt <- Qhy_alt + t(Q_alt[(1+(i-1)*n):(i*n),]) %*% Jn %*% hyt[(1+(i-1)*n):(i*n)];
      
    } else {
      
      Qhz <- Qhz + t(Q[(1+(i-1)*n):(i*n),]) %*% hzt[(1+(i-1)*n):(i*n),];
      QQ  <- QQ  + t(Q[(1+(i-1)*n):(i*n),]) %*% Q[(1+(i-1)*n):(i*n),];
      Qhy <- Qhy + t(Q[(1+(i-1)*n):(i*n),]) %*% hyt[(1+(i-1)*n):(i*n)];
      
      Qhz_alt <- Qhz_alt + t(Q_alt[(1+(i-1)*n):(i*n),]) %*% hzt[(1+(i-1)*n):(i*n),];
      QQ_alt  <- QQ_alt  + t(Q_alt[(1+(i-1)*n):(i*n),]) %*% Q_alt[(1+(i-1)*n):(i*n),];
      Qhy_alt <- Qhy_alt + t(Q_alt[(1+(i-1)*n):(i*n),]) %*% hyt[(1+(i-1)*n):(i*n)];
      
    }
  }
  
  theta <- mldivide(t(Qhz) %*% ginv(QQ) %*% Qhz, t(Qhz) %*% ginv(QQ) %*% Qhy);
  theta_alt <- mldivide(t(Qhz_alt) %*% ginv(QQ_alt) %*% Qhz_alt, t(Qhz_alt) %*% ginv(QQ_alt) %*% Qhy_alt);
  
  e <- hyt - hzt %*% theta;
  e_alt <- hyt - hzt %*% theta_alt;
  if(info$ted == 1){
    for(i in 1:s){
      e[(1+(i-1)*n):(i*n)] <- Jn %*% e[(1+(i-1)*n):(i*n)];
      e_alt[(1+(i-1)*n):(i*n)] <- Jn %*% e_alt[(1+(i-1)*n):(i*n)];
    }
  }
  
  sigma2 <- mean((e-mean(e))^2);
  sigma4 <- mean((e-mean(e))^4);
  
  sigma2_alt <- mean((e_alt-mean(e_alt))^2);
  sigma4_alt <- mean((e_alt-mean(e_alt))^4);
  
  lambda <- theta[1:p];
  delta  <- theta[(p+1):kz];
  
  lambda_alt <- theta_alt[1:p];
  delta_alt  <- theta_alt[p+1:kz];
  
  lambdaW     <- array(0, dim = c(n,n));
  lambdaW_alt <- array(0, dim = c(n,n));
  for(j in 1:p){
    lambdaW     <- lambdaW     + lambda[j] * W[,,j];
    lambdaW_alt <- lambdaW_alt + lambda_alt[j] * W[,,j];
  }
  Sn <- diag(n) - lambdaW;
  Sn_alt <- diag(n) - lambdaW_alt;
  
  DSiD <- 1/(sigma2*ns) * t(Qhz) %*% ginv(QQ) %*% Qhz;
  DSiD_alt <- 1/(sigma2_alt*ns) * t(Qhz_alt) %*% ginv(QQ_alt) %*% Qhz_alt;
  
  SIG <- tryCatch(1/ns * solve(DSiD), error = function(e){cat("DSiD not invertible \n"); return(array(NA, dim = dim(DSiD)))})
  
  std <- sqrt(abs(diag(SIG)));
  tstat <- theta/std;
  
  SIG_alt <- tryCatch(1/ns * solve(DSiD_alt), error = function(e){cat("DSiD_alt not invertible \n"); return(array(NA, dim = dim(DSiD_alt)))})
  std_alt <- sqrt(abs(diag(SIG_alt)));
  tstat_alt <- theta_alt/std_alt;
  
  # Forecast intercept, Mattera and Otto (2024, Eq. 9), paper Eq. (17): the
  # orthonormal transformation removes phi0, which is recovered for each unit as
  # the time average of (I - rho W) y*_t - Gamma y*_{t-1} over the training window.
  phi0 <- NULL; smear <- NULL
  if (info$tl == 1 && info$stl == 0 && is.null(X)) {
    gam  <- as.numeric(theta[(p + 1):(p + n)])
    phi0 <- rowMeans(Sn %*% Y[, 2:(t + 1), drop = FALSE] - gam * Y[, 1:t, drop = FALSE])
    # In-sample Duan (1983) smearing factor of each unit, paper Eq. (37).
    yhat  <- solve(Sn, gam * Y[, 1:t, drop = FALSE] + phi0)
    smear <- rowMeans(exp(Y[, 2:(t + 1), drop = FALSE] - yhat))
  }

  results <- list(theta = theta, std = std, SIG = SIG, tstat = tstat, sigma2 = sigma2,
                  theta_alt = theta_alt, std_alt = std_alt, SIG_alt = SIG_alt, tstat_alt = tstat_alt, sigma2_alt = sigma2_alt,
                  e = array(e, dim = c(n, t)), e_alt = array(e_alt, dim = c(n, t)), hyt = array(hyt, dim = c(n, t)),
                  phi0 = phi0, smear = smear)
  
  return(results)
  
}


mldivide <- function(A, b){
  return(ginv(t(A) %*% A) %*% t(A) %*% b)
}


mldivide <- function(A, b){
  return(ginv(t(A) %*% A) %*% t(A) %*% b)
}

# One-step forecast, paper Eq. (18): (I - rho W)^{-1} (Gamma y*_T + phi0)
forecast_gmm <- function(est_obj, W, prev_vol) {
  rho   <- est_obj$theta[1]
  gamma <- est_obj$theta[-1]
  n     <- nrow(W)
  if (is.null(est_obj$phi0)) stop("est_obj$phi0 missing")
  solve(diag(n) - rho * W) %*% (gamma * prev_vol + est_obj$phi0)
}


# ==============================================================================
# 3. Distance measures and pre-whitening
# ==============================================================================

# Euclidean distance between the rows of Y (as TSclust::diss.EUCL)
diss_euclidean <- function(Y) {
  n <- nrow(Y)
  D <- matrix(0, n, n)
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    d <- sqrt(sum((Y[i, ] - Y[j, ])^2))
    D[i, j] <- d; D[j, i] <- d
  }
  D
}

# Correlation distance sqrt(2 (1 - r)) (as TSclust::diss.COR)
diss_correlation <- function(Y) {
  R <- cor(t(Y))
  sqrt(2 * (1 - R))
}

# Piccolo (1990) distance between AR coefficient vectors, Yule-Walker with
# AIC-selected order, zero-padded to a common length
diss_ar_pic_pair <- function(x, y, permissive = TRUE) {
  ar_coefs <- function(v) {
    fit <- ar(v)
    if (fit$order < 1 && permissive) fit <- ar(v, aic = FALSE, order.max = 1)
    if (fit$order < 1) stop("Could not find a valid AR order for the series")
    fit$ar
  }
  PIx <- ar_coefs(x)
  PIy <- ar_coefs(y)
  k <- max(length(PIx), length(PIy))
  PRIMAx <- c(PIx, rep(0, k - length(PIx)))
  PRIMAy <- c(PIy, rep(0, k - length(PIy)))
  as.numeric(dist(rbind(PRIMAx, PRIMAy)))
}

# Piccolo distance on the log-squared returns (log-ARCH coefficients), with
# zero returns floored at the smallest non-zero squared return
diss_piccolo <- function(Y) {
  n <- nrow(Y)
  D <- matrix(0, n, n)
  eps_pic <- min(Y[Y != 0]^2)
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    d <- diss_ar_pic_pair(log(pmax(Y[i, ]^2, eps_pic)), log(pmax(Y[j, ]^2, eps_pic)))
    D[i, j] <- d; D[j, i] <- d
  }
  D
}

# Generalised FEVD of Pesaran and Shin (1998), rows normalised to one (Eq. 33).
# Rows are receivers, columns transmitters.
compute_gfevd <- function(var_fit, H = 10) {
  N <- var_fit$K
  Sigma <- summary(var_fit)$covres
  Phi <- vars::Phi(var_fit, nstep = H - 1)
  sigma_diag <- diag(Sigma)
  theta <- matrix(0, N, N)
  for (i in 1:N) {
    denom <- 0
    for (h in 1:H) denom <- denom + (Phi[i, , h] %*% Sigma %*% Phi[i, , h])
    for (j in 1:N) {
      num <- 0
      for (h in 1:H) num <- num + (Phi[i, , h] %*% Sigma[, j])^2
      theta[i, j] <- (1 / sigma_diag[j]) * num / denom
    }
  }
  theta / rowSums(theta)
}

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
# 4. Data and descriptive statistics
# ==============================================================================

plot(as.Date(log_returns$Date), log_returns$Algeria, type = "l")
lines(as.Date(log_returns$Date), log_returns$Iran, col = "red")
lines(as.Date(log_returns$Date), log_returns$Libya, col = "blue")
lines(as.Date(log_returns$Date), log_returns$Nigeria, col = "darkgreen")
lines(as.Date(log_returns$Date), log_returns$Saudi_arabia, col = "orange")
lines(as.Date(log_returns$Date), log_returns$UAE, col = "purple")

cor(log_returns[, 2:7])                       # Table 2

Y <- t(as.matrix(log_returns[, 2:7]))

#   Y_raw   : log-returns, 503 months (February 1983 to December 2024)
#   Y_W     : log-returns aligned with the whitened panel, used to build W
#   Y_w(T0) : whitened panel, 502 months, whitening estimated on window T0
# Y stays equal to the log-returns until Section 8, where it is set to the
# whitened panel.
Y_raw <- Y
Y_W   <- Y_raw[, -1, drop = FALSE]
T_in  <- 450
Y_wh  <- Y_w(T_in)
wh_fit <- attr(Y_wh, "whitening_fit")
cat("Correlation of the AR(1) residuals, training window, before whitening:\n")
print(round(wh_fit$correlation_before, 3))
cat("After whitening:\n")
print(round(wh_fit$correlation_after, 3))
attr(Y_wh, "whitening_fit") <- NULL

n <- dim(Y)[1]
t <- dim(Y)[2]


# ==============================================================================
# 5. Full-sample GARCH models and weight matrices
# ==============================================================================

start_date <- as.Date("1983-02-01")
date_index <- seq(from = start_date, by = "months", length.out = ncol(Y))
Y_xts <- xts(t(Y), order.by = date_index)

spec_univ <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0)),
  distribution.model = "std"
)

# ---- Univariate GARCH(1,1), Student-t (Table C.1) ----------------------------
fit_list <- list()
results_df <- data.frame(Variable = character(), Mu = numeric(), Omega = numeric(),
                         Alpha1 = numeric(), Beta1 = numeric(), Shape = numeric(),
                         LogLik = numeric(), stringsAsFactors = FALSE)
pval_df <- data.frame(Variable = character(), Mu = numeric(), Omega = numeric(),
                      Alpha1 = numeric(), Beta1 = numeric(), Shape = numeric(),
                      stringsAsFactors = FALSE)

for (i in 1:ncol(Y_xts)) {
  series_name <- colnames(Y_xts)[i]
  cat("Fitting series:", series_name, "\n")
  fit <- ugarchfit(spec = spec_univ, data = Y_xts[, i])
  fit_list[[series_name]] <- fit
  coef_vals  <- coef(fit)
  coef_table <- as.data.frame(fit@fit$matcoef)
  results_df <- rbind(results_df, data.frame(
    Variable = series_name, Mu = coef_vals["mu"], Omega = coef_vals["omega"],
    Alpha1 = coef_vals["alpha1"], Beta1 = coef_vals["beta1"], Shape = coef_vals["shape"],
    LogLik = likelihood(fit)))
  pval_df <- rbind(pval_df, data.frame(
    Variable = series_name, Mu = coef_table["mu", "Pr(>|t|)"],
    Omega = coef_table["omega", "Pr(>|t|)"], Alpha1 = coef_table["alpha1", "Pr(>|t|)"],
    Beta1 = coef_table["beta1", "Pr(>|t|)"], Shape = coef_table["shape", "Pr(>|t|)"]))
}

combined_df <- data.frame(Variable = results_df$Variable)
for (param in c("Mu", "Omega", "Alpha1", "Beta1", "Shape")) {
  combined_df[[param]] <- sprintf("%.4f (%.4f)", results_df[[param]], pval_df[[param]])
}
combined_df$LogLik <- round(results_df$LogLik, 4)
print(combined_df)

# ---- CCC-GARCH (Table 3) and CCC-based weights, Eqs. (27)-(28) ---------------
spec_ccc <- dccspec(uspec = multispec(replicate(ncol(Y_xts), spec_univ)),
                    dccOrder = c(0, 0), distribution = "mvt")
ccc_fit <- dccfit(spec_ccc, data = Y_xts)
print(ccc_fit)

ccc_corr <- rcor(ccc_fit)[, , 1]
distance_matrix_ccc <- sqrt(2 * (1 - ccc_corr))
diag(distance_matrix_ccc) <- 0
Wmat_ccc <- 1 / (distance_matrix_ccc^2 + 1)
diag(Wmat_ccc) <- 0
Wmat_ccc <- Wmat_ccc / max(eigen(Wmat_ccc, symmetric = TRUE, only.values = TRUE)$values)
print(Wmat_ccc)

# ---- DCC-GARCH (Table 3) and DCC-based weights, Eqs. (29)-(30) ---------------
spec_dcc <- dccspec(uspec = multispec(replicate(ncol(Y_xts), spec_univ)),
                    dccOrder = c(1, 1), distribution = "mvt")
dcc_fit <- dccfit(spec_dcc, data = Y_xts)
print(dcc_fit)

corr_array_dcc <- aperm(rcor(dcc_fit), perm = c(3, 1, 2))
avg_corr_dcc <- apply(corr_array_dcc, c(2, 3), mean, na.rm = TRUE)
distance_matrix_dcc <- sqrt(2 * (1 - avg_corr_dcc))
Wmat_dcc <- 1 / (distance_matrix_dcc^2 + 1)
diag(Wmat_dcc) <- 0
Wmat_dcc <- (Wmat_dcc + t(Wmat_dcc)) / 2
Wmat_dcc <- Wmat_dcc / max(eigen(Wmat_dcc, only.values = TRUE)$values)
print(Wmat_dcc)

# ---- GO-GARCH (Table 4) and GO-based weights, Eqs. (31)-(32) -----------------
# The fit is stored on first estimation and reloaded afterwards (see header).
spec_go <- gogarchspec(multispec(replicate(nrow(Y), spec_univ)), distribution.model = "mvnorm")
go_cache <- "W_go_cached_logret_v2.rds"
if (file.exists(go_cache)) {
  go_fit <- readRDS(go_cache)$go_fit
  cat("GO-GARCH: stored fit loaded (logLik =", round(likelihood(go_fit), 4), ")\n")
} else {
  set.seed(123)
  go_fit <- gogarchfit(spec_go, data = t(Y))
  cat("GO-GARCH: estimated (logLik =", round(likelihood(go_fit), 4), ")\n")
}
print(go_fit)

loglik <- likelihood(go_fit)
n_obs  <- ncol(Y)
k_par  <- 3 * nrow(Y)
cat("AIC:", (-2 * loglik + 2 * k_par) / n_obs, "\n")
cat("BIC:", (-2 * loglik + k_par * log(n_obs)) / n_obs, "\n")

T_go <- dim(rcov(go_fit))[3]
R_go_list <- lapply(1:T_go, function(tt) {
  V_t <- rcov(go_fit)[, , tt]
  D_inv <- diag(1 / sqrt(diag(V_t)))
  D_inv %*% V_t %*% D_inv
})
R_go_avg <- Reduce("+", R_go_list) / T_go
distance_go <- sqrt(2 * (1 - R_go_avg))
W_go <- 1 / (distance_go^2 + 1)
diag(W_go) <- 0
W_go <- (W_go + t(W_go)) / 2
W_go <- W_go / max(eigen(W_go, only.values = TRUE)$values)
if (!file.exists(go_cache)) saveRDS(list(go_fit = go_fit, W_go = W_go), go_cache)
print(W_go)

# ---- Connectedness-based weights, Eqs. (33)-(34) -----------------------------
# Directional: rows are receivers, columns transmitters, no symmetrisation.
# vars:: is written out because other packages also export VAR().
eps_dy_full   <- min(Y[Y != 0]^2)
Ystar_dy_full <- t(log(Y^2 + eps_dy_full))
lag_sel_dy    <- vars::VARselect(Ystar_dy_full, lag.max = 4, type = "const")$selection["AIC(n)"]
p_dy          <- max(as.integer(lag_sel_dy), 1)
var_fit_dy    <- vars::VAR(Ystar_dy_full, p = p_dy, type = "const")
cat("Connectedness VAR lag order (AIC):", p_dy, "\n")

connectedness_dy <- compute_gfevd(var_fit_dy, H = 10)
dimnames(connectedness_dy) <- list(rownames(Y), rownames(Y))
Wmat_dy <- connectedness_dy
diag(Wmat_dy) <- 0
Wmat_dy <- Wmat_dy / max(Mod(eigen(Wmat_dy, only.values = TRUE)$values))
print(round(Wmat_dy, 4))

from_others_dy <- rowSums(connectedness_dy) - diag(connectedness_dy)
to_others_dy   <- colSums(connectedness_dy) - diag(connectedness_dy)
print(round(data.frame(FROM = from_others_dy, TO = to_others_dy,
                       NET = to_others_dy - from_others_dy), 4))

# ---- Euclidean, correlation and Piccolo weights, Eqs. (19)-(26) --------------
Wmat1 <- diss_euclidean(Y)
Wmat1 <- 1 / Wmat1
diag(Wmat1) <- 0
Wmat1 <- Wmat1 / max(eigen(Wmat1, only.values = TRUE)$values)

Wmat2 <- diss_correlation(Y)
Wmat2 <- 1 / Wmat2
diag(Wmat2) <- 0
Wmat2 <- Wmat2 / max(eigen(Wmat2, only.values = TRUE)$values)

Wmat3 <- diss_piccolo(Y)
Wmat3 <- 1 / Wmat3
diag(Wmat3) <- 0
Wmat3 <- Wmat3 / max(eigen(Wmat3, only.values = TRUE)$values)


# ==============================================================================
# 6. Weight matrices for each training window
# ==============================================================================

generate_weights_for_T0 <- function(Y, T0, start_date = DATA_START,
                                    go_cache_file = sprintf("W_go_cached_logret_v2_T0_%d.rds", T0)) {

  cat(paste0("\n>>> Weight matrices for training window T0 = ", T0, "\n"))
  Y_train <- Y[, 1:T0]
  date_index_sub <- seq(from = as.Date(start_date), by = "months", length.out = T0)
  Y_xts_sub <- xts(t(Y_train), order.by = date_index_sub)

  spec_univ <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0)),
    distribution.model = "std"
  )
  mspec <- multispec(replicate(ncol(Y_xts_sub), spec_univ))

  # CCC
  spec_ccc <- dccspec(uspec = mspec, dccOrder = c(0, 0), distribution = "mvt")
  ccc_fit <- tryCatch(dccfit(spec_ccc, data = Y_xts_sub), error = function(e) NULL)
  if (!is.null(ccc_fit)) {
    ccc_corr <- rcor(ccc_fit)[, , 1]
    dist_ccc <- sqrt(2 * (1 - ccc_corr))
    W_ccc <- 1 / (dist_ccc^2 + 1)
    diag(W_ccc) <- 0
    W_ccc <- W_ccc / max(eigen(W_ccc, symmetric = TRUE, only.values = TRUE)$values)
  } else { W_ccc <- NULL }

  # DCC, time-averaged correlations
  spec_dcc <- dccspec(uspec = mspec, dccOrder = c(1, 1), distribution = "mvt")
  dcc_fit <- tryCatch(dccfit(spec_dcc, data = Y_xts_sub), error = function(e) NULL)
  if (!is.null(dcc_fit)) {
    avg_corr_dcc <- apply(rcor(dcc_fit), c(1, 2), mean)
    dist_dcc <- sqrt(2 * (1 - avg_corr_dcc))
    W_dcc <- 1 / (dist_dcc^2 + 1)
    diag(W_dcc) <- 0
    W_dcc <- W_dcc / max(eigen(W_dcc, symmetric = TRUE, only.values = TRUE)$values)
  } else { W_dcc <- NULL }

  # GO, time-averaged correlations; stored per window (see header)
  spec_go <- gogarchspec(mspec, distribution.model = "mvnorm")
  if (file.exists(go_cache_file)) {
    W_go <- readRDS(go_cache_file)
  } else {
    set.seed(123)
    go_fit <- tryCatch(gogarchfit(spec_go, data = Y_xts_sub), error = function(e) NULL)
    if (!is.null(go_fit)) {
      rcov_go <- rcov(go_fit)
      T_sub <- dim(rcov_go)[3]
      R_go_avg <- matrix(0, nrow = nrow(Y_train), ncol = nrow(Y_train))
      for (tt in 1:T_sub) {
        V_t <- rcov_go[, , tt]
        D_inv <- diag(1 / sqrt(diag(V_t)))
        R_go_avg <- R_go_avg + (D_inv %*% V_t %*% D_inv)
      }
      R_go_avg <- R_go_avg / T_sub
      dist_go <- sqrt(2 * (1 - R_go_avg))
      W_go <- 1 / (dist_go^2 + 1)
      diag(W_go) <- 0
      W_go <- W_go / max(eigen(W_go, symmetric = TRUE, only.values = TRUE)$values)
      saveRDS(W_go, go_cache_file)
    } else { W_go <- NULL }
  }

  # Connectedness, VAR on log-squared returns with AIC lag order
  dy_fit <- tryCatch({
    eps_local_dy <- min(Y_train[Y_train != 0]^2)
    Ystar_dy <- t(log(Y_train^2 + eps_local_dy))
    n_assets_dy <- ncol(Ystar_dy)
    lag_sel <- vars::VARselect(Ystar_dy,
                               lag.max = min(4, floor(nrow(Ystar_dy) / (n_assets_dy + 1)) - 1),
                               type = "const")$selection["AIC(n)"]
    p_local <- max(as.integer(lag_sel), 1)
    vars::VAR(Ystar_dy, p = p_local, type = "const")
  }, error = function(e) {
    message("   Connectedness VAR failed at T0 = ", T0, ": ", conditionMessage(e))
    NULL
  })
  if (!is.null(dy_fit)) {
    W_dy <- compute_gfevd(dy_fit, H = 10)
    diag(W_dy) <- 0
    W_dy <- W_dy / max(Mod(eigen(W_dy, only.values = TRUE)$values))
  } else { W_dy <- NULL }

  # Euclidean, correlation and Piccolo
  W_euc <- 1 / (diss_euclidean(Y_train) + 1e-10)
  diag(W_euc) <- 0
  W_euc <- W_euc / max(eigen(W_euc, symmetric = TRUE, only.values = TRUE)$values)

  W_cor <- 1 / (diss_correlation(Y_train) + 1e-10)
  diag(W_cor) <- 0
  W_cor <- W_cor / max(eigen(W_cor, symmetric = TRUE, only.values = TRUE)$values)

  W_pic <- 1 / (diss_piccolo(Y_train) + 1e-10)
  diag(W_pic) <- 0
  W_pic <- W_pic / max(eigen(W_pic, symmetric = TRUE, only.values = TRUE)$values)

  list(CCC = W_ccc, DCC = W_dcc, GO = W_go, DY = W_dy,
       Euclidean = W_euc, Correlation = W_cor, Piccolo = W_pic)
}

T0_values <- c(200, 250, 300, 350, 400, 450)
all_weights_storage <- list()

for (T0 in T0_values) {
  weights <- generate_weights_for_T0(Y_W, T0)   # W from the log-returns of the window
  all_weights_storage[[as.character(T0)]] <- weights
  assign(paste0("Wmat_CCC_", T0), weights$CCC)
  assign(paste0("Wmat_DCC_", T0), weights$DCC)
  assign(paste0("Wmat_GO_", T0),  weights$GO)
  assign(paste0("Wmat_DY_", T0),  weights$DY)
  assign(paste0("Wmat1_", T0),    weights$Euclidean)
  assign(paste0("Wmat2_", T0),    weights$Correlation)
  assign(paste0("Wmat3_", T0),    weights$Piccolo)
}


# ==============================================================================
# 7. Network figures
# ==============================================================================

country_names  <- c("Algeria", "Iran", "Libya", "Nigeria", "Saudi Arabia", "UAE")
country_colors <- wes_palette("Zissou1", length(country_names), type = "continuous")

get_network_layout <- function(Wmat) {
  Wmat <- as.matrix(Wmat); colnames(Wmat) <- country_names
  Wsym <- (Wmat + t(Wmat)) / 2
  dmat <- 1 / (Wsym + 1e-6); diag(dmat) <- 0
  layout_with_mds(graph_from_adjacency_matrix(Wsym, mode = "undirected", weighted = TRUE),
                  dist = dmat)
}

# Edge shading: quantile classes of the weight (light grey = weakest, dark grey
# = strongest). Width is proportional to the weight.
plot_network <- function(Wmat) {
  Wmat <- as.matrix(Wmat); dimnames(Wmat) <- list(country_names, country_names)
  g <- graph_from_adjacency_matrix((Wmat + t(Wmat)) / 2, mode = "undirected", weighted = TRUE)
  V(g)$color <- country_colors; V(g)$size <- 35
  V(g)$label <- NA; V(g)$frame.color <- "white"
  cl <- classIntervals(E(g)$weight, n = 8, style = "quantile")
  E(g)$color <- findColours(cl, pal = gray(seq(0.9, 0.1, length = 9)))
  E(g)$width <- E(g)$weight * 6
  plot(g, layout = get_network_layout(Wmat), main = NULL, vertex.label = NA)
}

plot_all_networks <- function(Wlist) {
  n <- length(Wlist); ncol_grid <- 3
  nrow_plots <- ceiling(n / ncol_grid)
  slots <- c(seq_len(n), rep(0, nrow_plots * ncol_grid - n))
  if (n %% ncol_grid == 1) {
    s <- (nrow_plots - 1) * ncol_grid + 1; slots[s:(s + 2)] <- c(0, n, 0)
  }
  layout(mat = rbind(matrix(slots, nrow = nrow_plots, byrow = TRUE), rep(n + 1, ncol_grid)),
         heights = c(rep(4, nrow_plots), 1))
  par(mar = c(1, 1, 1, 1))
  for (i in seq_along(Wlist)) {
    plot_network(Wlist[[i]])
    mtext(paste0("(", letters[i], ")"), side = 1, line = -1.5, adj = 0.5, cex = 1, font = 2)
  }
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", legend = country_names, ncol = 3, pch = 21, pt.bg = country_colors,
         pt.cex = 1.8, cex = 0.95, x.intersp = 0.8, inset = c(0, 0.05), bty = "n")
}

# Figure 2: the six symmetric weight matrices, full sample
plot_all_networks(list(Wmat1, Wmat2, Wmat3, Wmat_ccc, Wmat_dcc, W_go))

# Figure F.1: the same for each training window
for (T0 in T0_values) {
  plot_all_networks(list(get(paste0("Wmat1_", T0)), get(paste0("Wmat2_", T0)),
                         get(paste0("Wmat3_", T0)), get(paste0("Wmat_CCC_", T0)),
                         get(paste0("Wmat_DCC_", T0)), get(paste0("Wmat_GO_", T0))))
  mtext(sprintf("T0 = %d", T0), side = 3, line = -1.5, outer = TRUE, cex = 1.1, font = 2)
}

# Connectedness network (Figure 3): MDS layout with a small repulsion step so
# that nearly identical dissimilarities do not overlap nodes. Only the dominant
# direction of each pair is drawn, from W' - W.
spread_layout <- function(lay, min_dist = 0.28, iter = 200) {
  lay <- scale(lay, center = TRUE, scale = FALSE)
  lay <- lay / max(abs(lay))
  n <- nrow(lay)
  for (it in seq_len(iter)) {
    moved <- FALSE
    for (i in 1:(n - 1)) for (j in (i + 1):n) {
      d    <- lay[j, ] - lay[i, ]
      dist <- sqrt(sum(d^2))
      if (dist < min_dist) {
        if (dist < 1e-8) {
          d    <- c(runif(1, -1, 1), runif(1, -1, 1))
          dist <- sqrt(sum(d^2))
        }
        push     <- (min_dist - dist) / 2 * (d / dist)
        lay[i, ] <- lay[i, ] - push
        lay[j, ] <- lay[j, ] + push
        moved    <- TRUE
      }
    }
    if (!moved) break
  }
  lay
}

plot_dy_netflow <- function(Wdy, max_ref = NULL, max_width = 3.2, vsize = 34, lay = NULL) {
  Wdy <- as.matrix(Wdy); dimnames(Wdy) <- list(country_names, country_names)
  net <- colSums(Wdy) - diag(Wdy) - (rowSums(Wdy) - diag(Wdy))   # TO - FROM
  frame_col <- ifelse(net > 0, "black", "grey80")                 # net transmitters dark
  # rows of Wdy are receivers: (W' - W)[i, j] > 0 means i transmits more to j,
  # and igraph draws an edge i -> j for a positive entry [i, j]
  W_net <- t(Wdy) - Wdy; W_net[W_net < 0] <- 0
  g <- graph_from_adjacency_matrix(W_net, mode = "directed", weighted = TRUE)
  V(g)$color <- country_colors
  V(g)$size  <- vsize
  V(g)$label <- NA
  V(g)$frame.color <- frame_col
  if (is.null(max_ref)) max_ref <- max(E(g)$weight)
  E(g)$width <- 0.6 + (E(g)$weight / max_ref) * (max_width - 0.6)
  E(g)$color <- gray(0.75 - 0.40 * (E(g)$weight / max_ref))
  if (is.null(lay)) lay <- spread_layout(get_network_layout(Wdy))
  plot(g, layout = lay, vertex.label = NA, main = NULL, vertex.frame.width = 2,
       edge.arrow.size = 0.35, edge.arrow.width = 0.9, edge.curved = 0.12)
}

# One full-sample layout for all panels, so differences reflect the weights only
lay_fixed <- spread_layout(get_network_layout(Wmat_dy), min_dist = 0.28)
W_dy_list <- lapply(T0_values, function(t0) get(paste0("Wmat_DY_", t0)))
names(W_dy_list) <- paste0("T0 = ", T0_values)
global_max <- max(sapply(W_dy_list, function(W) { M <- t(W) - W; max(M[M > 0]) }))

# Figure 3
pdf("Figure_5.pdf", width = 6, height = 6)
par(mfrow = c(1, 1), mar = c(5, 2, 2, 2))
plot_dy_netflow(Wmat_dy, lay = lay_fixed)
legend("bottom", legend = country_names, ncol = 3, pch = 21, pt.bg = country_colors,
       pt.cex = 1.6, cex = 0.9, bty = "n", inset = c(0, -0.18), xpd = TRUE)
dev.off()

# Figure F.2
pdf("Figure_6.pdf", width = 9, height = 7)
layout(rbind(matrix(1:6, nrow = 2, byrow = TRUE), rep(7, 3)), heights = c(4, 4, 1))
par(mar = c(1.5, 1, 1, 1))
for (i in seq_along(W_dy_list)) {
  plot_dy_netflow(W_dy_list[[i]], max_ref = global_max, vsize = 30, lay = lay_fixed)
  mtext(sprintf("(%s) %s", letters[i], names(W_dy_list)[i]),
        side = 1, line = -0.5, adj = 0.5, cex = 0.9, font = 2)
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", legend = country_names, ncol = 3, pch = 21, pt.bg = country_colors,
       pt.cex = 1.8, cex = 0.95, x.intersp = 0.8, bty = "n")
dev.off()


# ==============================================================================
# 8. GMM estimates at T0 = 450 (Table 5)
# ==============================================================================

# From here on, the network models and the MGARCH benchmarks use the whitened
# panel. The weight matrices above were built from the log-returns.
Y <- Y_wh
cat("Whitened panel:", nrow(Y), "series x", ncol(Y), "months\n")

T_in <- 450
Y_in <- Y[, 1:T_in]

weight_matrices <- list(
  Net_CCC     = Wmat_CCC_450,
  Net_DCC     = Wmat_DCC_450,
  Net_GO      = Wmat_GO_450,
  Net_DY      = Wmat_DY_450,
  Euclidean   = Wmat1_450,
  Correlation = Wmat2_450,
  Piccolo     = Wmat3_450
)

for (name in names(weight_matrices)) {
  min_val <- min(Y_in[Y_in != 0]^2)                 # training window only
  model_output <- GMM_SDPD_2SLS_ARCH_ind_timelags(
    Y = log(pmax(Y_in^2, min_val)), X = NULL, W = weight_matrices[[name]],
    info = list(ksy = 10, ksx = 10, stl = 0, tl = 1, ted = 0))
  assign(name, model_output)
}

for (name in names(weight_matrices)) {
  res <- get(name)
  cat("\n=========== Results for", name, "===========\n")
  cat("Theta:\n");           print(res$theta)
  cat("Standard errors:\n"); print(res$std)
  cat("t-statistics:\n");    print(round(res$tstat, 3))
  cat("Sigma2:\n");          print(res$sigma2)
  cat("Largest eigenvalue modulus of rho W:",
      round(max(Mod(eigen(res$theta[1] * weight_matrices[[name]], only.values = TRUE)$values)), 4), "\n")
}


# ==============================================================================
# 9. Rolling one-step forecasts
# ==============================================================================

# Forecast errors of the MGARCH benchmarks, written by MGARCH_benchmarks.R
forecast_subfolder <- "Standard MGARCH models (whitened)"
models <- c("ccc", "dcc", "go")
eps_filenames <- list(min = "min", one_pct = "one_pct.1%", small = "small")
T0_all_eps  <- 450
T0_min_only <- c(200, 250, 300, 350, 400)

read_ferr_matrix <- function(model, T0, eps) {
  path <- file.path(forecast_subfolder, paste0("fERR_", model, "_T0"),
                    sprintf("fERR_%s_T0_%d_eps_%s.csv", model, T0, eps_filenames[[eps]]))
  if (file.exists(path)) {
    return(as.matrix(sapply(read.csv(path), as.numeric)))
  } else {
    message("File not found: ", path)
    return(NULL)
  }
}

standard_fERR_list <- list()
for (model in models) {
  model_upper <- toupper(model)
  standard_fERR_list[[model_upper]] <- list()
  T0_chr <- as.character(T0_all_eps)
  standard_fERR_list[[model_upper]][[T0_chr]] <- list()
  for (eps in names(eps_filenames)) {
    mat <- read_ferr_matrix(model, T0_all_eps, eps)
    if (!is.null(mat)) standard_fERR_list[[model_upper]][[T0_chr]][[eps]] <- mat
  }
  for (T0 in T0_min_only) {
    T0_chr <- as.character(T0)
    standard_fERR_list[[model_upper]][[T0_chr]] <- list()
    mat <- read_ferr_matrix(model, T0, "min")
    if (!is.null(mat)) standard_fERR_list[[model_upper]][[T0_chr]][["min"]] <- mat
  }
}

set.seed(123)

calculate_metrics <- function(err_mat) {
  RMSFE <- apply(err_mat, 2, function(x) sqrt(mean(x^2, na.rm = TRUE)))
  MAFE  <- apply(err_mat, 2, function(x) mean(abs(x), na.rm = TRUE))
  list(RMSFE = RMSFE, MAFE = MAFE)
}

# Bootstrap intervals for the pooled RMSFE and MAFE, resampling forecast dates
# and keeping the countries together
bootstrap_CI <- function(err_mat, R = 2000, seed = 123) {
  set.seed(seed)
  n_fore <- nrow(err_mat)
  rms_b <- numeric(R)
  maf_b <- numeric(R)
  for (b in seq_len(R)) {
    idx  <- sample(n_fore, n_fore, replace = TRUE)
    errs <- as.vector(err_mat[idx, , drop = FALSE])
    rms_b[b] <- sqrt(mean(errs^2, na.rm = TRUE))
    maf_b[b] <- mean(abs(errs), na.rm = TRUE)
  }
  list(RMSFE_CI = quantile(rms_b, c(0.025, 0.975), na.rm = TRUE),
       MAFE_CI  = quantile(maf_b, c(0.025, 0.975), na.rm = TRUE))
}

# GMM fits of the seven specifications on the training window T0
get_models_for_T0 <- function(T0, Y, eps_baseline) {
  w_euc <- get(paste0("Wmat1_", T0))
  w_cor <- get(paste0("Wmat2_", T0))
  w_pic <- get(paste0("Wmat3_", T0))
  w_ccc <- get(paste0("Wmat_CCC_", T0))
  w_dcc <- get(paste0("Wmat_DCC_", T0))
  w_go  <- get(paste0("Wmat_GO_", T0))
  w_dy  <- get(paste0("Wmat_DY_", T0))
  Y_star <- log(pmax(Y[, 1:T0]^2, eps_baseline))

  fit_gmm <- function(W_mat, name) {
    start_time <- Sys.time()
    obj <- tryCatch(
      GMM_SDPD_2SLS_ARCH_ind_timelags(Y = Y_star, X = NULL, W = W_mat,
                                      info = list(ksy = 10, ksx = 10, stl = 0, tl = 1, ted = 0)),
      error = function(e) { message("Estimation failed for T0 = ", T0, ", ", name); NULL })
    cat(sprintf("   -> %s estimation: %.2f sec\n", name,
                as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
    obj
  }

  list(
    Euclidean   = list(obj = fit_gmm(w_euc, "Euclidean"),   W = w_euc),
    Correlation = list(obj = fit_gmm(w_cor, "Correlation"), W = w_cor),
    Piccolo     = list(obj = fit_gmm(w_pic, "Piccolo"),     W = w_pic),
    Net_CCC     = list(obj = fit_gmm(w_ccc, "Net_CCC"),     W = w_ccc),
    Net_DCC     = list(obj = fit_gmm(w_dcc, "Net_DCC"),     W = w_dcc),
    Net_GO      = list(obj = fit_gmm(w_go,  "Net_GO"),      W = w_go),
    Net_DY      = list(obj = fit_gmm(w_dy,  "Net_DY"),      W = w_dy)
  )
}

# Rolling one-step forecasts with W and parameters fixed at T0. The MGARCH
# errors (re-estimated at every origin) are appended through standard_errors.
run_forecasts <- function(Y, models, T0, eps, standard_errors = NULL) {
  T_total     <- ncol(Y)
  n_forecasts <- T_total - T0
  n_series    <- nrow(Y)
  fERR <- lapply(models, function(x) matrix(NA, nrow = n_forecasts, ncol = n_series))
  time_mat <- matrix(0, nrow = n_forecasts, ncol = length(models))
  colnames(time_mat) <- names(models)

  pb <- txtProgressBar(min = 1, max = n_forecasts, style = 3)
  for (i in seq_len(n_forecasts)) {
    setTxtProgressBar(pb, i)
    t_train  <- T0 + i - 1
    prev_vol <- log(pmax(Y[, t_train]^2, eps))
    actual   <- log(pmax(Y[, t_train + 1]^2, eps))
    for (j in seq_along(names(models))) {
      name <- names(models)[j]
      est  <- models[[name]]$obj
      Wmat <- models[[name]]$W
      if (!is.null(est)) {
        t_start <- Sys.time()
        fcast <- forecast_gmm(est, Wmat, prev_vol)
        fERR[[name]][i, ] <- actual - fcast
        time_mat[i, j] <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
      }
    }
  }
  close(pb)
  # in-sample smearing factors travel with the error matrices (used by predicted_var)
  for (name in names(models)) if (!is.null(models[[name]]$obj$smear))
    attr(fERR[[name]], "smear") <- models[[name]]$obj$smear

  if (!is.null(standard_errors)) {
    for (model in names(standard_errors)) fERR[[model]] <- standard_errors[[model]]
  }
  metrics <- lapply(fERR, calculate_metrics)
  list(fERR = fERR,
       RMSFE = sapply(metrics, function(m) mean(m$RMSFE, na.rm = TRUE)),
       MAFE  = sapply(metrics, function(m) mean(m$MAFE,  na.rm = TRUE)),
       avg_time = colMeans(time_mat, na.rm = TRUE))
}

eps_vals <- c(
  min     = min(Y[Y != 0]^2),
  one_pct = as.numeric(quantile(Y[Y != 0]^2, 0.01)),
  small   = 1e-6
)

results <- list()

# ---- T0 = 450, all three zero adjustments ------------------------------------
T0 <- 450
results[[as.character(T0)]] <- list()
for (eps_name in names(eps_vals)) {
  cat("T0 =", T0, ", eps =", eps_name, "\n")
  current_models <- get_models_for_T0(T0, Y, eps_vals[[eps_name]])
  garch_models <- list()
  for (model in c("CCC", "DCC", "GO")) {
    ferr_mat <- standard_fERR_list[[model]][[as.character(T0)]][[eps_name]]
    if (!is.null(ferr_mat)) garch_models[[model]] <- ferr_mat
  }
  res <- run_forecasts(Y, current_models, T0, eps_vals[[eps_name]], standard_errors = garch_models)
  results[[as.character(T0)]][[eps_name]] <- res
}

# ---- Shorter training windows, eps = smallest non-zero value -----------------
# Each window has its own whitening and weight matrices.
eps_name <- "min"
for (T0 in c(200, 250, 300, 350, 400, 450)) {
  cat("\nT0 =", T0, ", eps = min\n")
  current_models <- get_models_for_T0(T0, Y_w(T0), eps_vals[[eps_name]])
  garch_models <- list()
  for (model in c("CCC", "DCC", "GO")) {
    ferr_mat <- standard_fERR_list[[model]][[as.character(T0)]][[eps_name]]
    if (!is.null(ferr_mat)) garch_models[[model]] <- ferr_mat
  }
  res <- run_forecasts(Y_w(T0), current_models, T0, eps_vals[[eps_name]], standard_errors = garch_models)
  if (is.null(results[[as.character(T0)]])) results[[as.character(T0)]] <- list()
  results[[as.character(T0)]][[eps_name]] <- res
}

# ---- Parameters re-estimated at every origin, W fixed (Section 6.3.2) --------
run_forecasts_recursive <- function(Y, weight_list, T0, eps) {
  n_forecasts <- ncol(Y) - T0
  fERR <- lapply(weight_list, function(x) matrix(NA, nrow = n_forecasts, ncol = nrow(Y)))
  pb <- txtProgressBar(min = 1, max = n_forecasts, style = 3)
  for (i in seq_len(n_forecasts)) {
    setTxtProgressBar(pb, i)
    t_train  <- T0 + i - 1
    Y_star   <- log(pmax(Y[, 1:t_train]^2, eps))
    prev_vol <- log(pmax(Y[, t_train]^2, eps))
    actual   <- log(pmax(Y[, t_train + 1]^2, eps))
    for (name in names(weight_list)) {
      Wmat <- weight_list[[name]]
      est <- tryCatch(GMM_SDPD_2SLS_ARCH_ind_timelags(
        Y = Y_star, X = NULL, W = Wmat,
        info = list(ksy = 10, ksx = 10, stl = 0, tl = 1, ted = 0)), error = function(e) NULL)
      if (!is.null(est)) fERR[[name]][i, ] <- actual - forecast_gmm(est, Wmat, prev_vol)
    }
  }
  close(pb)
  metrics <- lapply(fERR, calculate_metrics)
  list(fERR = fERR, RMSFE = sapply(metrics, function(m) mean(m$RMSFE, na.rm = TRUE)))
}

weight_list_450 <- list(Euclidean = Wmat1_450, Correlation = Wmat2_450, Piccolo = Wmat3_450,
                        Net_CCC = Wmat_CCC_450, Net_DCC = Wmat_DCC_450,
                        Net_GO = Wmat_GO_450, Net_DY = Wmat_DY_450)
res_recursive_450 <- run_forecasts_recursive(Y, weight_list_450, T0 = 450, eps = eps_vals[["min"]])
cmp_recursive <- rbind(
  fixed     = results[["450"]][["min"]]$RMSFE[names(weight_list_450)],
  recursive = res_recursive_450$RMSFE)
cat("\nAverage country RMSFE, parameters fixed vs re-estimated (T0 = 450, eps = min):\n")
print(round(rbind(cmp_recursive, change_pct = 100 * (cmp_recursive[2, ] / cmp_recursive[1, ] - 1)), 4))

# ---- Table E.1 ---------------------------------------------------------------
pooled_rmsfe <- function(e) sqrt(mean(e^2, na.rm = TRUE))
pooled_mafe  <- function(e) mean(abs(e), na.rm = TRUE)

tabE1_rows <- c(Euclidean = "Euclidean-Diss", Correlation = "Correlation-Diss",
                Piccolo = "Piccolo-Diss", Net_CCC = "Net CCC-GARCH",
                Net_DCC = "Net DCC-GARCH", Net_GO = "Net GO-GARCH",
                Net_DY = "Net Connectedness")
tabE1 <- t(sapply(names(tabE1_rows), function(m) {
  base <- results[["450"]][["one_pct"]]$fERR[[m]]
  win  <- sapply(c(200, 250, 300, 350, 400), function(t0)
    pooled_rmsfe(results[[as.character(t0)]][["min"]]$fERR[[m]]))
  epsc <- sapply(c("min", "one_pct", "small"), function(k)
    pooled_rmsfe(results[["450"]][[k]]$fERR[[m]]))
  c(RMSFE = pooled_rmsfe(base), MAFE = pooled_mafe(base),
    setNames(win, paste0("T0_", c(200, 250, 300, 350, 400))),
    setNames(epsc, paste0("eps_", c("min", "one_pct", "small"))))
}))
rownames(tabE1) <- tabE1_rows
cat("\n=== Table E.1 ===\n")
print(round(tabE1, 4))
write.csv(tabE1, file.path(out_dir, "TableE1_sensitivity.csv"))

# Network forecast errors at T0 = 450, used by the forecast-combination script
saveRDS(list(fERR = results[["450"]][["one_pct"]]$fERR, fERR_min = results[["450"]][["min"]]$fERR,
             Y = Y, eps_vals = eps_vals),
        file.path(out_dir, "network_fERR_T0_450.rds"))


# ==============================================================================
# 10. Forecast evaluation (Tables 6, E.1 and E.2)
# ==============================================================================

eval_T0  <- 450
eval_eps <- "one_pct"
all_ferr    <- results[[as.character(eval_T0)]][[eval_eps]]$fERR
model_names <- names(all_ferr)
eps_used    <- eps_vals[[eval_eps]]
country_names <- gsub(".", " ", rownames(Y), fixed = TRUE)
net_names <- c("Euclidean", "Correlation", "Piccolo", "Net_CCC", "Net_DCC", "Net_GO", "Net_DY")
benchmark_name <- "Net_GO"

cat("CCC and DCC forecast errors identical:",
    isTRUE(all.equal(all_ferr$CCC, all_ferr$DCC, check.attributes = FALSE)), "\n")

# ---- Bootstrap intervals (Table 6) -------------------------------------------
ci_list <- lapply(all_ferr, function(mat) bootstrap_CI(mat, R = 2000, seed = 123))
print(ci_list)

# ---- Heteroskedasticity-robust DM test of Harvey et al. (2024), Eqs. (44)-(47)
lrvarnw_custom <- function(data, lags = NULL) {
  data <- as.matrix(data)
  T_len <- nrow(data)
  if (is.null(lags)) lags <- min(floor(1.2 * T_len^(1/3)), T_len)
  w <- (lags + 1 - (0:lags)) / (lags + 1)
  lrv <- (t(data) %*% data) / T_len
  if (lags > 0) {
    for (i in 1:lags) {
      Gammai <- (t(data[(i + 1):T_len, , drop = FALSE]) %*% data[1:(T_len - i), , drop = FALSE]) / T_len
      lrv <- lrv + w[i + 1] * (Gammai + t(Gammai))
    }
  }
  lrv
}

# leave-k-out kernel estimate of the spot variance
spv_lko_custom <- function(tgrid, x, t, h, k) {
  m <- length(t); n <- length(tgrid)
  kh <- function(val) dnorm(val / h) / h
  ktth <- kh(outer(t, tgrid, "-"))
  ktth[abs(row(ktth) - col(ktth)) < k] <- 0
  col_sums <- colSums(ktth)
  col_sums[col_sums == 0] <- 1e-10
  w <- ktth / matrix(col_sums, nrow = m, ncol = n, byrow = TRUE)
  colSums(w * matrix(x, nrow = m, ncol = n))
}

# bandwidth by cross-validation
spv_cvbw_lko_custom <- function(x, t, k, h_range) {
  cv <- numeric(length(h_range))
  for (i in seq_along(h_range)) {
    xhat <- spv_lko_custom(t, x, t, h_range[i], k)
    cv[i] <- mean((x - xhat)^2, na.rm = TRUE)
  }
  h_range[which.min(cv)]
}

# Returns the standard DM p-value and the two Harvey et al. statistics
# (DM*, used in the paper, and DM**)
harvey_dm_test <- function(d, h_range = NULL, k = 20) {
  d <- as.numeric(as.vector(d))
  n <- length(d)
  t <- (1:n) / n
  if (is.null(h_range)) h_range <- seq(0.05, 0.3, length.out = 50)
  h_opt <- tryCatch(spv_cvbw_lko_custom(x = d^2, t = t, k = k, h_range = h_range),
                    error = function(e) 0.15)
  kh <- function(val) dnorm(val / h_opt) / h_opt
  ktth <- kh(outer(t, t, "-"))
  w  <- ktth / matrix(colSums(ktth), nrow = n, ncol = n, byrow = TRUE)
  sv <- colSums(w * matrix(d^2, nrow = n, ncol = n))
  sigma <- sqrt(sv)
  sigma[sigma < 1e-8] <- 1e-8
  lags <- min(floor(1.2 * n^(1/3)), n)
  dv <- d / sigma
  x_term <- 1 / sigma
  stat1 <- sqrt(n) * mean(dv) / sqrt(lrvarnw_custom(as.matrix(dv), lags))
  stat2 <- sqrt(n) * mean(dv * x_term) / sqrt(lrvarnw_custom(as.matrix(dv * x_term), lags))
  stat_std <- sqrt(n) * mean(d) / sqrt(lrvarnw_custom(as.matrix(d), lags))
  list(Pval_Standard = 2 * (1 - pnorm(abs(stat_std))),
       Pval_Harvey_1 = 2 * (1 - pnorm(abs(stat1))),
       Pval_Harvey_2 = 2 * (1 - pnorm(abs(stat2))))
}

# Panel-robust DM test, Eq. (43): losses averaged across countries at each
# date, Bartlett HAC variance with K = floor(1.2 T^(1/3)) lags.
# panel_robust_dm_test() works on forecast errors (squared-error loss),
# panel_robust_dm_loss() on any pair of loss matrices (e.g. QLIKE).
panel_robust_dm_loss <- function(L1, L2, lag_trunc = NULL) {
  d <- rowMeans(L1 - L2, na.rm = TRUE); d <- d[is.finite(d)]
  Te <- length(d)
  if (Te < 5) return(list(p.value = NA, mean_diff = NA, T_eff = Te))
  if (is.null(lag_trunc)) lag_trunc <- floor(1.2 * Te^(1/3))
  db <- mean(d); c0 <- d - db; hv <- mean(c0^2)
  if (lag_trunc > 0) for (k in 1:lag_trunc)
    hv <- hv + 2 * (1 - k / (lag_trunc + 1)) * mean(c0[(k + 1):Te] * c0[1:(Te - k)])
  if (!is.finite(hv) || hv <= 0) return(list(p.value = NA, mean_diff = db, T_eff = Te))
  list(p.value = 2 * pnorm(-abs(db / sqrt(hv / Te))), mean_diff = db, T_eff = Te)
}
panel_robust_dm_test <- function(e1, e2, lag_trunc = NULL) panel_robust_dm_loss(e1^2, e2^2, lag_trunc)

# Standard and Harvey DM p-values on pooled squared-error differentials
dm_pair <- function(e1, e2) {
  d <- na.omit(as.vector(e1)^2 - as.vector(e2)^2)
  if (sd(d) < 1e-10) return(c(std = 1, harvey = 1))
  r <- tryCatch(harvey_dm_test(d = d, h_range = seq(0.04, 0.3, length.out = 50), k = 20),
                error = function(e) NULL)
  if (is.null(r)) c(std = NA, harvey = NA) else c(std = r$Pval_Standard, harvey = r$Pval_Harvey_1)
}

# ---- Table 6: losses and DM tests against Net. GO-GARCH ----------------------
loss_row <- function(e) c(RMSFE = sqrt(mean(e^2, na.rm = TRUE)), MAFE = mean(abs(e), na.rm = TRUE),
                          MedAE = median(abs(e), na.rm = TRUE), Bias = mean(e, na.rm = TRUE))

tab6 <- do.call(rbind, lapply(net_names, function(nm) {
  e  <- all_ferr[[nm]]
  ci <- ci_list[[nm]]$RMSFE_CI
  if (nm == benchmark_name) {
    p <- c(std = NA, harvey = NA); pr <- NA
  } else {
    p  <- dm_pair(e, all_ferr[[benchmark_name]])
    pr <- panel_robust_dm_test(e, all_ferr[[benchmark_name]])$p.value
  }
  data.frame(Model = nm, t(loss_row(e)), CI_low = ci[1], CI_high = ci[2],
             DM_Standard = p["std"], DM_Harvey = p["harvey"], DM_Panel = pr, row.names = NULL)
}))

# MCS, Tmax with B = 5000, eps at the smallest non-zero value, on all
# forecasts of this run
loss_sq_all_models <- sapply(results[["450"]][["min"]]$fERR,
                             function(err_mat) rowMeans(err_mat^2, na.rm = TRUE))
set.seed(123)
mcs_result_all <- MCS::MCSprocedure(loss_sq_all_models, alpha = 0.05, B = 5000, statistic = "Tmax")
print(mcs_result_all)

cat("\n=== Table 6 (T0 = 450, eps = one_pct) ===\n")
print(tab6, digits = 4, row.names = FALSE)
write.csv(tab6, file.path(out_dir, "Table6_network_specifications.csv"), row.names = FALSE)

# Directional accuracy with the Pesaran and Timmermann (1992) test
reconstruct_actual_log <- function(Y, T0, eps) {
  m <- matrix(NA, ncol(Y) - T0, nrow(Y))
  for (i in seq_len(nrow(m))) m[i, ] <- log(pmax(Y[, T0 + i]^2, eps))
  colnames(m) <- rownames(Y); m
}
reconstruct_prev_log <- function(Y, T0, eps) {
  m <- matrix(NA, ncol(Y) - T0, nrow(Y))
  for (i in seq_len(nrow(m))) m[i, ] <- log(pmax(Y[, T0 + i - 1]^2, eps))
  colnames(m) <- rownames(Y); m
}

pt_test <- function(ad, pd) {
  ok <- !is.na(ad) & !is.na(pd); ad <- ad[ok]; pd <- pd[ok]; n <- length(ad)
  if (n < 5) return(list(hit = NA, exp = NA, stat = NA, p = NA))
  ph <- mean(ad == pd); py <- mean(ad == 1); px <- mean(pd == 1)
  ps <- py * px + (1 - py) * (1 - px)
  vh <- ps * (1 - ps) / n
  vs <- ((2 * py - 1)^2 * px * (1 - px) + (2 * px - 1)^2 * py * (1 - py) +
           4 * px * py * (1 - px) * (1 - py) / n) / n
  den <- vh - vs
  st  <- if (is.finite(den) && den > 0) (ph - ps) / sqrt(den) else NA
  list(hit = ph, exp = ps, stat = st, p = if (is.na(st)) NA else 2 * pnorm(-abs(st)))
}

directional_accuracy <- function(fl, Y, T0, eps) {
  al <- reconstruct_actual_log(Y, T0, eps); pl <- reconstruct_prev_log(Y, T0, eps)
  ad <- sign(al - pl)
  do.call(rbind, lapply(names(fl), function(nm) {
    r <- pt_test(as.vector(ad), as.vector(sign((al - fl[[nm]]) - pl)))
    data.frame(Model = nm, HitRate = r$hit, Expected = r$exp, Stat = r$stat, PValue = r$p)
  }))
}
tab_dir <- directional_accuracy(all_ferr[net_names], Y, eval_T0, eps_used)
cat("\n=== Directional accuracy ===\n")
print(tab_dir, digits = 4, row.names = FALSE)

# ---- Table 6 (MGARCH part): network specifications and MGARCH benchmarks on both scales ---
# Log-squared scale: centred benchmark forecasts ln(sigma^2) + training-window
# mean of ln max(e^2, eps) - ln(sigma^2), written by MGARCH_benchmarks.R.
read_centred <- function(model, T0, eps) {
  f <- file.path(forecast_subfolder, sprintf("fERR_%s_T0_centred", model),
                 sprintf("fERR_%s_T0_%d_eps_%s.csv", model, T0, eps_filenames[[eps]]))
  if (file.exists(f)) as.matrix(sapply(read.csv(f), as.numeric)) else NULL
}
dcc_c <- read_centred("dcc", eval_T0, "one_pct")
go_c  <- read_centred("go",  eval_T0, "one_pct")
if (is.null(dcc_c) || is.null(go_c)) stop("Centred MGARCH forecast errors not found: run MGARCH_benchmarks.R first.")

# Variance scale: QLIKE and MSE, Eq. (40). The network forecasts are mapped by
# the in-sample smearing factor, Eq. (37). The MGARCH models forecast the
# variance directly.
LOGSPACE_MODELS <- net_names

# Smearing factor estimated recursively from past forecast errors (robustness).
# For the first min_obs origins it is the Gaussian value exp(-E[ln chi2_1]).
recursive_smearing <- function(fERR, min_obs = 10, init = exp(-(digamma(0.5) + log(2)))) {
  vapply(seq_len(nrow(fERR)), function(i)
    if (i <= min_obs) init else mean(exp(as.vector(fERR[1:(i - 1), , drop = FALSE])), na.rm = TRUE),
    numeric(1))
}
recursive_error_var <- function(fERR, min_obs = 10) {
  vapply(seq_len(nrow(fERR)), function(i)
    if (i <= min_obs) 0 else var(as.vector(fERR[1:(i - 1), , drop = FALSE]), na.rm = TRUE),
    numeric(1))
}

# correction = "duan": in-sample smearing attached by run_forecasts(), or the
# recursive factor when none is attached; "lognormal": exp(s2 / 2)
predicted_var <- function(fERR, actual_log, model_name, bias_correct = TRUE, correction = "duan") {
  pl <- actual_log - fERR
  if (bias_correct && model_name %in% LOGSPACE_MODELS) {
    sm <- attr(fERR, "smear")
    if (correction == "duan" && !is.null(sm)) {
      pl <- pl + matrix(log(sm), nrow(pl), ncol(pl), byrow = TRUE)
    } else {
      adj <- if (correction == "duan") log(recursive_smearing(fERR)) else recursive_error_var(fERR) / 2
      pl <- pl + matrix(adj, nrow(pl), ncol(pl))
    }
  }
  pv <- exp(pl); pv[pv <= 0] <- NA; pv
}
compute_qlike <- function(av, pv) { r <- av / pv; r[r <= 0] <- NA; r - log(r) - 1 }

al <- reconstruct_actual_log(Y, eval_T0, eps_used); av <- exp(al)
qlike_of <- function(e, nm) compute_qlike(av, predicted_var(e, al, nm))
mse_of   <- function(e, nm) (av - predicted_var(e, al, nm))^2

tab7 <- do.call(rbind, lapply(c(net_names, "DCC", "GO"), function(nm) {
  e_log <- if (nm == "DCC") dcc_c else if (nm == "GO") go_c else all_ferr[[nm]]
  dm_log <- if (nm == "DCC") NA else panel_robust_dm_test(e_log, dcc_c)$p.value
  dm_var <- if (nm %in% c("DCC", "GO")) NA else
    panel_robust_dm_loss(qlike_of(all_ferr[[nm]], nm), qlike_of(all_ferr$DCC, "DCC"))$p.value
  data.frame(Model = nm, t(loss_row(e_log)), DM_log_vs_DCC = dm_log,
             QLIKE = mean(qlike_of(all_ferr[[nm]], nm), na.rm = TRUE),
             MSE_var = mean(mse_of(all_ferr[[nm]], nm), na.rm = TRUE),
             DM_QLIKE_vs_DCC = dm_var, row.names = NULL)
}))
tab7 <- rbind(tab7, data.frame(Model = "DCC_uncentred", t(loss_row(all_ferr$DCC)),
                               DM_log_vs_DCC = NA, QLIKE = NA, MSE_var = NA, DM_QLIKE_vs_DCC = NA))
cat("\n=== Table 6, MGARCH benchmarks (T0 = 450, eps = one_pct) ===\n")
print(tab7, digits = 4, row.names = FALSE)
write.csv(tab7, file.path(out_dir, "Table7_benchmarks_both_scales.csv"), row.names = FALSE)

# Benchmarks tested against Net. GO-GARCH (text of Section 6.3.3)
for (nm in c("DCC", "GO")) {
  e  <- if (nm == "DCC") dcc_c else go_c
  p  <- dm_pair(e, all_ferr[[benchmark_name]])
  pr <- panel_robust_dm_test(e, all_ferr[[benchmark_name]])
  ci <- bootstrap_CI(e, R = 2000, seed = 123)
  cat(sprintf("Centred %s vs Net_GO: standard p = %.3f, Harvey p = %.3f, panel-robust p = %.3f, RMSFE CI [%.3f, %.3f]\n",
              nm, p["std"], p["harvey"], pr$p.value, ci$RMSFE_CI[1], ci$RMSFE_CI[2]))
}

# MCS on the log scale with the centred benchmarks, and under QLIKE and MSE
cen_min <- list(DCC_centred = read_centred("dcc", eval_T0, "min"),
                GO_centred  = read_centred("go",  eval_T0, "min"))
fm <- results[[as.character(eval_T0)]][["min"]]$fERR
set.seed(123)
mcs_set <- c("Euclidean", "Correlation", "Net_CCC", "Net_DCC", "Net_GO", "Net_DY")
print(MCS::MCSprocedure(sapply(c(fm[mcs_set], cen_min), function(e) rowMeans(e^2, na.rm = TRUE)),
                        alpha = 0.05, B = 5000, statistic = "Tmax"))

al_min <- reconstruct_actual_log(Y, eval_T0, eps_vals[["min"]]); av_min <- exp(al_min)
keep <- c("Net_CCC", "Net_DCC", "Net_GO", "Net_DY", "CCC", "DCC")   # heuristic weights give unbounded QLIKE
loss_q  <- sapply(keep, function(nm)
  rowMeans(compute_qlike(av_min, predicted_var(fm[[nm]], al_min, nm)), na.rm = TRUE))
loss_mv <- sapply(keep, function(nm)
  rowMeans((av_min - predicted_var(fm[[nm]], al_min, nm))^2, na.rm = TRUE))
set.seed(123); print(MCS::MCSprocedure(loss_q,  alpha = 0.05, B = 5000, statistic = "Tmax"))
set.seed(123); print(MCS::MCSprocedure(loss_mv, alpha = 0.05, B = 5000, statistic = "Tmax"))

# ---- Merged Table 6 (landscape): both scales, benchmark = centred DCC-GARCH --
# Log scale: squared-error loss on the log-squared returns, with the MGARCH
# forecasts centred. Variance scale: QLIKE with the in-sample smearing factor.
# Standard, Harvey and panel-robust DM tests against the (centred) DCC-GARCH.
dm_pair_loss <- function(L1, L2) {
  d <- as.vector(L1 - L2); d <- d[is.finite(d)]
  if (sd(d) < 1e-10) return(c(std = 1, harvey = 1))
  r <- tryCatch(harvey_dm_test(d = d, h_range = seq(0.04, 0.3, length.out = 50), k = 20),
                error = function(e) NULL)
  if (is.null(r)) c(std = NA, harvey = NA) else c(std = r$Pval_Standard, harvey = r$Pval_Harvey_1)
}
mcs_pvals <- function(m) {
  sh <- m@show
  setNames(as.numeric(sh[, grep("^MCS", colnames(sh))[1]]), rownames(sh))
}

merged_models <- c(net_names, "DCC", "GO")
e_log_of <- function(nm) if (nm == "DCC") dcc_c else if (nm == "GO") go_c else all_ferr[[nm]]
q_of     <- function(nm) qlike_of(all_ferr[[nm]], nm)
q_dcc    <- q_of("DCC")

# MCS on the log scale (eps at the smallest non-zero value, as in Table 6),
# all seven network specifications and the two centred benchmarks
mcs_in <- sapply(c(fm[net_names], cen_min), function(e) rowMeans(e^2, na.rm = TRUE))
colnames(mcs_in) <- c(net_names, "DCC", "GO")
set.seed(123)
mcs_merged <- MCS::MCSprocedure(mcs_in, alpha = 0.05, B = 5000, statistic = "Tmax", verbose = FALSE)
p_mcs <- mcs_pvals(mcs_merged)

tab6m <- do.call(rbind, lapply(merged_models, function(nm) {
  e  <- e_log_of(nm)
  ci <- bootstrap_CI(e, R = 2000, seed = 123)$RMSFE_CI
  if (nm == "DCC") {
    pl <- c(std = NA, harvey = NA); prl <- NA; pq <- c(std = NA, harvey = NA); prq <- NA
  } else {
    pl  <- dm_pair(e, dcc_c)
    prl <- panel_robust_dm_test(e, dcc_c)$p.value
    pq  <- dm_pair_loss(q_of(nm), q_dcc)
    prq <- panel_robust_dm_loss(q_of(nm), q_dcc)$p.value
  }
  data.frame(Model = nm, t(loss_row(e)), CI_low = ci[1], CI_high = ci[2],
             DMlog_Standard = pl["std"], DMlog_Harvey = pl["harvey"], DMlog_Panel = prl,
             MCS_p = if (nm %in% names(p_mcs)) p_mcs[[nm]] else 0,
             QLIKE = mean(q_of(nm), na.rm = TRUE), MSE_var = mean(mse_of(all_ferr[[nm]], nm), na.rm = TRUE),
             DMq_Standard = pq["std"], DMq_Harvey = pq["harvey"], DMq_Panel = prq, row.names = NULL)
}))
cat("\n=== Merged Table 6 (landscape), benchmark: centred DCC-GARCH ===\n")
print(tab6m, digits = 4, row.names = FALSE)
cat("MCS_p = 0 marks a model eliminated from the MCS.\n")
write.csv(tab6m, file.path(out_dir, "Table6_merged_both_scales.csv"), row.names = FALSE)

# Within-network comparison against Net. GO-GARCH (former Table 6, now in the
# appendix): tab6 above, written to Table6_network_specifications.csv

# ---- Table E.1: centred MGARCH rows -------------------------------------------
tabE1_mg <- t(sapply(c("dcc", "go"), function(m) {
  base <- read_centred(m, 450, "one_pct")
  win  <- sapply(c(200, 250, 300, 350, 400), function(t0) {
    e <- read_centred(m, t0, "min"); if (is.null(e)) NA else pooled_rmsfe(e) })
  epsc <- sapply(c("min", "one_pct", "small"), function(k) {
    e <- read_centred(m, 450, k); if (is.null(e)) NA else pooled_rmsfe(e) })
  c(RMSFE = pooled_rmsfe(base), MAFE = pooled_mafe(base),
    setNames(win, paste0("T0_", c(200, 250, 300, 350, 400))),
    setNames(epsc, paste0("eps_", c("min", "one_pct", "small"))))
}))
rownames(tabE1_mg) <- c("Std. CCC / DCC-GARCH, centred", "Std. GO-GARCH, centred")
tabE1_full <- rbind(tabE1, tabE1_mg)
cat("\n=== Table E.1 with the centred MGARCH benchmarks ===\n")
print(round(tabE1_full, 4))
write.csv(tabE1_full, file.path(out_dir, "TableE1_sensitivity_with_MGARCH.csv"))

# Retransformation: in-sample smearing, recursive smearing and lognormal correction
print(round(t(sapply(net_names, function(nm) {
  e  <- all_ferr[[nm]]; e_rec <- e; attr(e_rec, "smear") <- NULL
  c(QLIKE_insample  = mean(compute_qlike(av, predicted_var(e, al, nm)), na.rm = TRUE),
    QLIKE_recursive = mean(compute_qlike(av, predicted_var(e_rec, al, nm)), na.rm = TRUE),
    QLIKE_lognormal = mean(compute_qlike(av, predicted_var(e_rec, al, nm, correction = "lognormal")), na.rm = TRUE))
})), 4))

# Calibration (mean ratio of realized to forecast variance) and QLIKE by
# quartile of the cross-sectional mean of the squared whitened returns
print(round(t(sapply(c(net_names, "DCC", "GO"), function(nm) {
  pv <- predicted_var(all_ferr[[nm]], al, nm); k <- mean(av / pv, na.rm = TRUE)
  c(mean_ratio = k, QLIKE = mean(compute_qlike(av, pv), na.rm = TRUE),
    QLIKE_rescaled = mean(compute_qlike(av, pv * k), na.rm = TRUE))
})), 4))
q_vol <- cut(rowMeans(av), quantile(rowMeans(av), 0:4 / 4), include.lowest = TRUE,
             labels = c("Q1 calm", "Q2", "Q3", "Q4 turbulent"))
print(round(sapply(c(net_names, "DCC", "GO"), function(nm)
  tapply(rowMeans(compute_qlike(av, predicted_var(all_ferr[[nm]], al, nm)), na.rm = TRUE),
         q_vol, mean, na.rm = TRUE)), 3))


# ==============================================================================
# 11. Weight stability, subperiods and structural breaks (Section 6.3.5)
# ==============================================================================

weight_types <- c("CCC", "DCC", "GO", "DY", "Euclidean", "Correlation", "Piccolo")

# Table 8: variation of the 30 off-diagonal weights across the six windows
weight_variation <- function(storage, wtype, T0s) {
  Ws <- lapply(as.character(T0s), function(k) storage[[k]][[wtype]])
  Ws <- Ws[!sapply(Ws, is.null)]
  if (length(Ws) < 2) return(NULL)
  n   <- nrow(Ws[[1]])
  off <- which(row(diag(n)) != col(diag(n)))
  M   <- sapply(Ws, function(W) W[off])
  cell_sd   <- apply(M, 1, sd)
  cell_mean <- apply(M, 1, mean)
  cv        <- cell_sd / cell_mean
  rel_change <- sapply(2:ncol(M), function(j)
    norm(as.matrix(M[, j] - M[, j - 1]), "F") / norm(as.matrix(M[, j - 1]), "F"))
  data.frame(Weights = wtype, Mean_weight = mean(cell_mean), Mean_SD = mean(cell_sd),
             Mean_CV = mean(cv), Max_CV = max(cv), Max_rel_change = max(rel_change),
             Corr_first_last = cor(M[, 1], M[, ncol(M)]))
}

tab9 <- do.call(rbind, lapply(weight_types, function(w)
  weight_variation(all_weights_storage, w, T0_values)))
cat("\n=== Table 8 ===\n")
print(tab9, digits = 4, row.names = FALSE)
write.csv(tab9, file.path(out_dir, "Table9_weight_stability.csv"), row.names = FALSE)

# Subperiods: evaluation period split at February 2022
dates_oos  <- seq(as.Date(DATA_START), by = "months", length.out = ncol(Y))[(eval_T0 + 1):ncol(Y)]
idx_early  <- which(dates_oos <  as.Date("2022-02-01"))
idx_late   <- which(dates_oos >= as.Date("2022-02-01"))
tab_sub <- do.call(rbind, lapply(model_names, function(nm) {
  e <- all_ferr[[nm]]
  f <- function(ix) c(RMSFE = sqrt(mean(e[ix, ]^2, na.rm = TRUE)),
                      MAFE  = mean(abs(e[ix, ]), na.rm = TRUE),
                      Bias  = mean(e[ix, ], na.rm = TRUE))
  a <- f(idx_early); b <- f(idx_late)
  data.frame(Model = nm, RMSFE_pre = a["RMSFE"], RMSFE_post = b["RMSFE"],
             RMSFE_ratio = b["RMSFE"] / a["RMSFE"], MAFE_pre = a["MAFE"], MAFE_post = b["MAFE"],
             Bias_pre = a["Bias"], Bias_post = b["Bias"], row.names = NULL)
}))
cat("\n=== Subperiods (split at February 2022) ===\n")
print(tab_sub, digits = 4, row.names = FALSE)
write.csv(tab_sub, file.path(out_dir, "Subperiod_split.csv"), row.names = FALSE)

# Bai and Perron (2003) breaks and sup-F test on the whitened log-squared
# series, minimum segment of 15% (75 months) and 8% (40 months)
Ystar_full <- t(log(pmax(Y^2, eps_vals[["min"]])))
dates_full <- seq(as.Date(DATA_START), by = "months", length.out = nrow(Ystar_full))

tab_breaks <- do.call(rbind, lapply(seq_along(country_names), function(i) {
  y  <- Ystar_full[, i]
  bp <- tryCatch(breakpoints(y ~ 1, h = 0.15), error = function(e) NULL)
  sf <- tryCatch(sctest(y ~ 1, type = "supF"),  error = function(e) NULL)
  if (is.null(bp) || all(is.na(bp$breakpoints))) {
    brk_dates <- "none"; n_brk <- 0
  } else {
    n_brk <- length(bp$breakpoints)
    brk_dates <- paste(format(dates_full[bp$breakpoints], "%Y-%m"), collapse = "; ")
  }
  data.frame(Country = country_names[i], N_breaks = n_brk, Break_dates = brk_dates,
             supF_stat = if (is.null(sf)) NA else unname(sf$statistic),
             supF_pval = if (is.null(sf)) NA else sf$p.value, row.names = NULL)
}))
cat("\n=== Structural breaks, minimum segment 15% ===\n")
print(tab_breaks, digits = 4, row.names = FALSE)
write.csv(tab_breaks, file.path(out_dir, "Structural_breaks.csv"), row.names = FALSE)

for (i in seq_along(country_names)) {
  bp <- tryCatch(breakpoints(Ystar_full[, i] ~ 1, h = 0.08), error = function(e) NULL)
  if (!is.null(bp) && !all(is.na(bp$breakpoints)))
    cat(sprintf("%-13s h = 0.08 breaks: %s\n", country_names[i],
                paste(format(dates_full[bp$breakpoints], "%Y-%m"), collapse = "; ")))
  else cat(sprintf("%-13s h = 0.08 breaks: none\n", country_names[i]))
}


# ==============================================================================
# 12. Economic value (Table 9)
# ==============================================================================

# Both exercises use the log-returns. The variance forecast h_{t,j} of each
# whitened series is mapped to a return variance with the training-window
# whitening matrix, v_{t,i} = sum_j (Sigma^{1/2})_{ij}^2 h_{t,j}, and the
# conditional mean is the AR(1) forecast of Eq. (35).
pred_var_list <- lapply(model_names, function(nm)
  predicted_var(all_ferr[[nm]], al, nm, bias_correct = TRUE))
names(pred_var_list) <- model_names

wh_450 <- fit_whitening(Y_in = Y_raw[, 1:(eval_T0 + 1)], Y = Y_raw, mean_model = "AR1")
ee_450 <- eigen(wh_450$Sigma, symmetric = TRUE)
S_half <- ee_450$vectors %*% diag(sqrt(pmax(ee_450$values, 0)), nrow(Y_raw)) %*% t(ee_450$vectors)
stopifnot(max(abs(wh_450$residuals[, -1] - Y), na.rm = TRUE) < 1e-8)

oos_cols <- (eval_T0 + 2):ncol(Y_raw)
R_oos <- t(Y_raw[, oos_cols, drop = FALSE]); colnames(R_oos) <- country_names
M_oos <- t(wh_450$coefficients$mu + wh_450$coefficients$phi * Y_raw[, oos_cols - 1, drop = FALSE])
pred_var_list <- lapply(pred_var_list, function(pv) pv %*% t(S_half^2))
stopifnot(nrow(R_oos) == nrow(pred_var_list[[1]]))

# ---- Volatility timing, Fleming, Kirby and Ostdiek (2001) --------------------
inv_var_weights <- function(pv) {
  w <- 1 / pv
  w[!is.finite(w)] <- NA
  w / rowSums(w, na.rm = TRUE)
}

# Annualised fee (basis points) that an investor with quadratic utility and
# risk aversion gamma pays to hold portfolio r1 instead of r0
performance_fee <- function(r1, r0, gamma = 5, periods = 12) {
  u <- function(r, d) mean((1 + r - d) - (gamma / (2 * (1 + gamma))) * (1 + r - d)^2)
  f <- tryCatch(uniroot(function(d) u(r1, d) - u(r0, 0), interval = c(-0.05, 0.05),
                        extendInt = "yes", tol = 1e-14)$root, error = function(e) NA_real_)
  f * periods * 10000
}

sharpe_boot_ci <- function(r, R = 5000, seed = 123) {
  set.seed(seed)
  r <- r[!is.na(r)]; n <- length(r)
  s <- replicate(R, { i <- sample(n, n, replace = TRUE); mean(r[i]) / sd(r[i]) * sqrt(12) })
  quantile(s, c(0.025, 0.975), na.rm = TRUE)
}

port_r <- lapply(pred_var_list, function(pv) rowSums(inv_var_weights(pv) * R_oos))
names(port_r) <- model_names
ew_r <- rowMeans(R_oos)

bench_vt <- "DCC"
tab_vt <- do.call(rbind, lapply(model_names, function(nm) {
  r  <- port_r[[nm]]
  ci <- sharpe_boot_ci(r)
  data.frame(Model = nm, Mean_ann = mean(r, na.rm = TRUE) * 12,
             SD_ann = sd(r, na.rm = TRUE) * sqrt(12),
             Sharpe_ann = mean(r, na.rm = TRUE) / sd(r, na.rm = TRUE) * sqrt(12),
             Sharpe_lo = ci[1], Sharpe_hi = ci[2],
             Fee_vs_EW_bp = performance_fee(r, ew_r),
             Fee_vs_DCC_bp = performance_fee(r, port_r[[bench_vt]]), row.names = NULL)
}))
ci_ew <- sharpe_boot_ci(ew_r)
tab_vt <- rbind(tab_vt, data.frame(Model = "Equally weighted", Mean_ann = mean(ew_r) * 12,
                                   SD_ann = sd(ew_r) * sqrt(12),
                                   Sharpe_ann = mean(ew_r) / sd(ew_r) * sqrt(12),
                                   Sharpe_lo = ci_ew[1], Sharpe_hi = ci_ew[2], Fee_vs_EW_bp = 0,
                                   Fee_vs_DCC_bp = performance_fee(ew_r, port_r[[bench_vt]])))
cat("\n=== Table 9: volatility timing ===\n")
print(tab_vt, digits = 4, row.names = FALSE)

# ---- Value at Risk, Kupiec (1995) and Christoffersen (1998) ------------------
kupiec_test <- function(hits, alpha) {
  n <- length(hits); x <- sum(hits); pi_hat <- x / n
  xl <- function(k, q) if (k == 0) 0 else k * log(q)
  LR <- -2 * (xl(n - x, 1 - alpha) + xl(x, alpha) - xl(n - x, 1 - pi_hat) - xl(x, pi_hat))
  list(hits = x, rate = pi_hat, LR = LR, p = 1 - pchisq(LR, 1))
}

# hits is an (origins x countries) matrix: transitions are counted within each
# country and pooled
christoffersen_test <- function(hits, alpha) {
  hits <- as.matrix(hits)
  tr <- function(a, b) sum(sapply(seq_len(ncol(hits)), function(j) {
    h <- hits[, j]; h <- h[!is.na(h)]; sum(h[-length(h)] == a & h[-1] == b) }))
  n00 <- tr(0, 0); n01 <- tr(0, 1); n10 <- tr(1, 0); n11 <- tr(1, 1)
  uc <- kupiec_test(as.vector(hits[!is.na(hits)]), alpha)
  if (is.na(uc$LR)) return(list(LR_cc = NA, p_cc = NA))
  xlogy <- function(x, y) ifelse(x == 0, 0, x * log(y))
  pi01 <- n01 / (n00 + n01); pi11 <- if ((n10 + n11) > 0) n11 / (n10 + n11) else 0
  pi_h <- (n01 + n11) / (n00 + n01 + n10 + n11)
  LR_ind <- -2 * (xlogy(n00 + n10, 1 - pi_h) + xlogy(n01 + n11, pi_h) -
                    xlogy(n00, 1 - pi01) - xlogy(n01, pi01) -
                    xlogy(n10, 1 - pi11) - xlogy(n11, pi11))
  list(LR_cc = uc$LR + LR_ind, p_cc = 1 - pchisq(uc$LR + LR_ind, 2))
}

nu  <- 6
q_t <- function(a, nu) qt(a, df = nu) * sqrt((nu - 2) / nu)   # standardised Student-t

tab_var <- do.call(rbind, lapply(model_names, function(nm) {
  pv <- pred_var_list[[nm]]
  do.call(rbind, lapply(c(0.05, 0.01), function(a) {
    VaR   <- M_oos + q_t(a, nu) * sqrt(pv)
    hit_m <- (R_oos < VaR) * 1
    hits  <- as.vector(hit_m); hits <- hits[!is.na(hits)]
    uc <- kupiec_test(hits, a); cc <- christoffersen_test(hit_m, a)
    data.frame(Model = nm, Level = paste0(a * 100, "%"), Mean_VaR = mean(VaR, na.rm = TRUE),
               Exceed_rate = uc$rate, Expected = a, N_hits = uc$hits,
               Kupiec_p = uc$p, Christoffersen_p = cc$p_cc, row.names = NULL)
  }))
}))
cat("\n=== Table 9: Value-at-Risk backtests ===\n")
print(tab_var, digits = 4, row.names = FALSE)

write.csv(tab_vt,  file.path(out_dir, "Table10_volatility_timing.csv"), row.names = FALSE)
write.csv(tab_var, file.path(out_dir, "Table10_VaR_backtests.csv"),     row.names = FALSE)


# ==============================================================================
# 13. Simulation: estimation cost (Tables D.1 and D.2)
# ==============================================================================

# DGP parameters seeded by N, so a given N always gives the same true network
make_dgp_params <- function(N, seed_base = 7000) {
  set.seed(seed_base + N)
  A <- matrix(runif(N * N, 0, 1), N, N)
  A <- (A + t(A)) / 2
  diag(A) <- 0
  W_true <- A / max(eigen(A, only.values = TRUE, symmetric = TRUE)$values)
  list(W = W_true, gamma = runif(N, -0.1, 0.3))
}

simulate_network_logarch <- function(N, Tt, rho, gamma_vec, W, phi0 = rep(0, N),
                                     burn_in = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Tt_full <- Tt + burn_in
  IminrhoW_inv <- solve(diag(N) - rho * W)
  Ystar <- matrix(0, N, Tt_full)
  Y     <- matrix(0, N, Tt_full)
  for (t in 1:Tt_full) {
    ylag <- if (t == 1) rep(0, N) else Ystar[, t - 1]
    u_t  <- log(rnorm(N)^2) - (-1.2704)
    Ystar[, t] <- IminrhoW_inv %*% (phi0 + gamma_vec * ylag + u_t)
    sign_t     <- sample(c(-1, 1), N, replace = TRUE)
    Y[, t]     <- sign_t * exp(Ystar[, t] / 2)
  }
  keep <- (burn_in + 1):Tt_full
  list(Y = Y[, keep, drop = FALSE], Ystar = Ystar[, keep, drop = FALSE])
}

make_xts_matrix <- function(Y) {
  mat <- t(Y)
  colnames(mat) <- paste0("A", seq_len(nrow(Y)))
  xts(mat, order.by = seq(as.Date("2000-01-01"), by = "day", length.out = ncol(Y)))
}

fit_network_gmm_safe <- function(Y, W) {
  eps <- min(Y[Y != 0]^2)
  t0 <- Sys.time()
  out <- tryCatch(
    GMM_SDPD_2SLS_ARCH_ind_timelags(Y = log(pmax(Y^2, eps)), X = NULL, W = W,
                                    info = list(ksy = 10, ksx = 10, stl = 0, tl = 1, ted = 0)),
    error = function(e) NULL)
  list(fit = out, elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
       timed_out = is.null(out))
}

# All MGARCH variants use sGARCH(1,1) marginals with zero mean
fit_mgarch_safe <- function(Y, type = c("CCC", "DCC", "GO"),
                            timeout_ccc = 300, timeout_dcc = 900, timeout_go = 300) {
  type <- match.arg(type)
  N <- nrow(Y)
  Y_xts <- make_xts_matrix(Y)
  spec_u0 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                        mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                        distribution.model = "std")
  this_timeout <- switch(type, CCC = timeout_ccc, DCC = timeout_dcc, GO = timeout_go)
  t0 <- Sys.time()
  err_msg <- NULL
  out <- tryCatch({
    withTimeout({
      mspec <- multispec(replicate(N, spec_u0))
      if (type == "CCC") {
        dccfit(dccspec(uspec = mspec, dccOrder = c(0, 0), distribution = "mvt"), data = Y_xts)
      } else if (type == "DCC") {
        dccfit(dccspec(uspec = mspec, dccOrder = c(1, 1), distribution = "mvt"), data = Y_xts)
      } else {
        gogarchfit(gogarchspec(mspec, distribution.model = "mvnorm"), data = Y_xts)
      }
    }, timeout = this_timeout, cpu = this_timeout, elapsed = this_timeout, onTimeout = "silent")
  }, error = function(e) { err_msg <<- conditionMessage(e); NULL })
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(out) && !is.null(err_msg)) cat(sprintf("    [%s] error: %s\n", type, err_msg))
  list(fit = out, elapsed_sec = elapsed, timed_out = is.null(out),
       is_genuine_timeout = is.null(out) && is.null(err_msg) && elapsed >= 0.95 * this_timeout)
}

make_progress_tracker <- function(total_steps) {
  pb <- txtProgressBar(min = 0, max = max(total_steps, 1), style = 3)
  count <- 0
  list(tick = function() { count <<- count + 1; setTxtProgressBar(pb, count) },
       close = function() close(pb))
}

run_full_simulation_study <- function(N_vals = c(20, 50, 100), Tt = 500, R_vals = c(15, 10, 5),
                                      rho_true = 0.5,
                                      checkpoint_file = "simulation_study_full_v8.rds") {
  stopifnot(length(N_vals) == length(R_vals))
  results <- if (file.exists(checkpoint_file)) readRDS(checkpoint_file) else list()
  pb <- make_progress_tracker(sum(R_vals))
  already_done <- sum(sapply(N_vals, function(N) sum(grepl(paste0("^N", N, "_r"), names(results)))))
  for (i in seq_len(already_done)) pb$tick()
  for (k in seq_along(N_vals)) {
    N <- N_vals[k]; R <- R_vals[k]
    dgp <- make_dgp_params(N)
    for (r in 1:R) {
      key <- paste0("N", N, "_r", r)
      if (!is.null(results[[key]])) next
      sim <- simulate_network_logarch(N, Tt, rho_true, dgp$gamma, dgp$W, seed = 10000 * N + r)
      res_net <- fit_network_gmm_safe(sim$Y, dgp$W)
      res_ccc <- fit_mgarch_safe(sim$Y, "CCC")
      res_dcc <- fit_mgarch_safe(sim$Y, "DCC")
      res_go  <- fit_mgarch_safe(sim$Y, "GO")
      cat(sprintf("  [N = %d, rep %d/%d] Net %.3fs | CCC %.2fs | DCC %.2fs | GO %.2fs\n",
                  N, r, R, res_net$elapsed_sec, res_ccc$elapsed_sec,
                  res_dcc$elapsed_sec, res_go$elapsed_sec))
      results[[key]] <- list(
        N = N, rep = r, rho_true = rho_true,
        rho_hat = if (!is.null(res_net$fit)) res_net$fit$theta[1] else NA,
        time_net = res_net$elapsed_sec,
        time_ccc = res_ccc$elapsed_sec, ccc_failed = res_ccc$timed_out, ccc_timeout = res_ccc$is_genuine_timeout,
        time_dcc = res_dcc$elapsed_sec, dcc_failed = res_dcc$timed_out, dcc_timeout = res_dcc$is_genuine_timeout,
        time_go  = res_go$elapsed_sec,  go_failed  = res_go$timed_out,  go_timeout  = res_go$is_genuine_timeout)
      saveRDS(results, checkpoint_file)
      pb$tick()
    }
  }
  pb$close(); results
}

summarize_scaling_study <- function(results) {
  N_vals <- sort(unique(sapply(results, function(x) x$N)))
  do.call(rbind, lapply(N_vals, function(N) {
    sub <- results[sapply(results, function(x) x$N == N)]
    rho_hat <- sapply(sub, function(x) x$rho_hat); rho_true <- sub[[1]]$rho_true
    tm <- function(f) mean(sapply(sub, function(x) x[[f]]), na.rm = TRUE)
    data.frame(N = N, R = length(sub),
               Net_mean = tm("time_net"), CCC_mean = tm("time_ccc"), DCC_mean = tm("time_dcc"),
               GO_mean = tm("time_go"),
               CCC_fail = mean(sapply(sub, function(x) x$ccc_failed)),
               DCC_fail = mean(sapply(sub, function(x) x$dcc_failed)),
               GO_fail  = mean(sapply(sub, function(x) x$go_failed)),
               Speedup_vs_DCC = tm("time_dcc") / tm("time_net"),
               Speedup_vs_GO  = tm("time_go")  / tm("time_net"),
               rho_mean = mean(rho_hat, na.rm = TRUE),
               rho_bias = mean(rho_hat, na.rm = TRUE) - rho_true,
               rho_sd   = sd(rho_hat, na.rm = TRUE),
               rho_rmse = sqrt(mean((rho_hat - rho_true)^2, na.rm = TRUE)))
  }))
}

scaling_results <- run_full_simulation_study(N_vals = c(20, 50, 100), Tt = 500, R_vals = c(15, 10, 5))
df_scaling <- summarize_scaling_study(scaling_results)
cat("\n=== Table D.2 ===\n")
print(df_scaling, digits = 4, row.names = FALSE)
write.csv(df_scaling, file.path(out_dir, "TableD2_estimation_cost.csv"), row.names = FALSE)

# Table D.1: per-origin computing time in the application. The network
# forecast uses parameters fixed at T0; the benchmarks are re-estimated at each
# origin (timings written by MGARCH_benchmarks.R).
net_sec <- mean(results[["450"]][["min"]]$avg_time, na.rm = TRUE)
tabD1 <- do.call(rbind, lapply(c("ccc", "dcc", "go"), function(m) {
  f <- file.path(forecast_subfolder, sprintf("Computation_Metrics_%s", m),
                 sprintf("Metrics_%s_T0_450_eps_min.csv", m))
  if (!file.exists(f)) return(NULL)
  data.frame(Model = toupper(m), Seconds_per_origin = mean(read.csv(f)$Time_Sec, na.rm = TRUE))
}))
tabD1 <- rbind(tabD1, data.frame(Model = "Network log-ARCH", Seconds_per_origin = net_sec))
cat("\n=== Table D.1 ===\n")
print(tabD1, digits = 4, row.names = FALSE)
write.csv(tabD1, file.path(out_dir, "TableD1_per_origin_time.csv"), row.names = FALSE)


# ==============================================================================
# 14. Simulation: forecast accuracy (Table D.3)
# ==============================================================================

# Two designs, each run through the pipeline of the application:
#   (a) DCC-GARCH DGP, the benchmark is correctly specified (whitened data)
#   (b) network log-ARCH DGP, the network model is correctly specified
# In each replication: simulate 450 + 52 + 1 months, whiten on the training
# window, build W from the raw returns of the training window (bounded Pearson
# weights for Net_Corr, CCC-based weights for Net_CCC), estimate by GMM and
# forecast with the Eq. (17) intercept, fit DCC- and GO-GARCH, and produce 52
# rolling one-step forecasts with parameters fixed. QLIKE and MSE are computed
# against the realised squared values; the oracle row uses the true conditional
# variance. Runtime is several hours, checkpointed after each replication.

sim_dir <- "Simulation accuracy"
if (!dir.exists(sim_dir)) dir.create(sim_dir, recursive = TRUE)
LOGSPACE_MODELS <- union(LOGSPACE_MODELS, c("Net_Corr", "Net_CCC"))

# DCC-GARCH(1,1) with a one-factor unconditional correlation (loadings 0.75-0.95,
# i.e. correlations of about 0.55-0.90, comparable to the OPEC panel) and
# standardised multivariate t(8) innovations. Parameters are seeded by N only.
make_dcc_params <- function(N, seed_base = 8100) {
  set.seed(seed_base + N)
  lam   <- runif(N, 0.75, 0.95)
  Rbar  <- tcrossprod(lam); diag(Rbar) <- 1
  alpha <- runif(N, 0.05, 0.12); beta <- runif(N, 0.80, 0.88)
  s2    <- runif(N, 0.08, 0.12)^2                         # monthly oil-return variance
  list(Rbar = Rbar, alpha = alpha, beta = beta, omega = s2 * (1 - alpha - beta),
       a = 0.04, b = 0.93, nu = 8)
}

simulate_dcc <- function(N, Tt, p, burn_in = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  TT <- Tt + burn_in
  h  <- p$omega / (1 - p$alpha - p$beta); Q <- p$Rbar; z_prev <- rep(0, N)
  Y <- H <- matrix(0, N, TT); Sig <- vector("list", TT)
  for (t in seq_len(TT)) {
    if (t > 1) {
      h <- p$omega + p$alpha * Y[, t - 1]^2 + p$beta * h
      Q <- (1 - p$a - p$b) * p$Rbar + p$a * tcrossprod(z_prev) + p$b * Q
    }
    R <- Q / sqrt(tcrossprod(diag(Q)))
    g <- rchisq(1, p$nu) / (p$nu - 2)                      # unit-variance multivariate t
    z <- as.numeric(t(chol(R)) %*% rnorm(N)) / sqrt(g)
    Y[, t] <- sqrt(h) * z; H[, t] <- h; z_prev <- z
    Sig[[t]] <- R * sqrt(tcrossprod(h))                     # full conditional covariance
  }
  keep <- (burn_in + 1):TT
  list(Y = Y[, keep, drop = FALSE], h = H[, keep, drop = FALSE], Sig = Sig[keep])
}

# Network log-ARCH DGP with a common factor in the Gaussian innovations, so that
# the returns are correlated and W can be estimated from them. The true conditional
# variance is E[y_t^2 | past] = exp(m_t) * kappa, with m_t = S^{-1}(phi0 + Gamma y*_{t-1})
# and kappa_i = E exp([S^{-1} u_t]_i), computed once by Monte Carlo.
make_net_params <- function(N, seed_base = 7000) {
  set.seed(seed_base + N)
  gamma <- runif(N, -0.1, 0.3); rho <- 0.5; lam <- runif(N, 0.75, 0.95)
  # True W: the bounded transform of the return correlations. With
  # y = sign(z) exp(y*/2) and z Gaussian with correlation Rz,
  # corr(y) = (2/pi) asin(Rz) up to the volatility scaling.
  Rz <- tcrossprod(lam); diag(Rz) <- 1
  Ry <- (2 / pi) * asin(Rz)
  W  <- spec_norm(1 / (2 * (1 - Ry) + 1))
  phi0  <- as.numeric((diag(N) - rho * W - diag(gamma)) %*% rep(log(0.1^2), N))  # E[y*] = ln(0.01)
  Sinv  <- solve(diag(N) - rho * W)
  draw_u <- function(n) { f <- rnorm(n); z <- outer(f, lam) + matrix(rnorm(n * N), n, N) %*% diag(sqrt(1 - lam^2), N)
                          log(z^2) + 1.2704 }
  kappa <- colMeans(exp(draw_u(2e5) %*% t(Sinv)))
  list(W = W, gamma = gamma, rho = rho, lam = lam, phi0 = phi0, Sinv = Sinv, kappa = kappa)
}

simulate_network <- function(N, Tt, p, burn_in = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  TT <- Tt + burn_in
  Ys <- Y <- H <- matrix(0, N, TT)
  for (t in seq_len(TT)) {
    ylag <- if (t == 1) as.numeric(rep(log(0.01), N)) else Ys[, t - 1]
    f <- rnorm(1); z <- p$lam * f + sqrt(1 - p$lam^2) * rnorm(N)
    m <- as.numeric(p$Sinv %*% (p$phi0 + p$gamma * ylag))
    Ys[, t] <- m + as.numeric(p$Sinv %*% (log(z^2) + 1.2704))
    Y[, t]  <- sign(z) * exp(Ys[, t] / 2)
    H[, t]  <- exp(m) * p$kappa
  }
  keep <- (burn_in + 1):TT
  list(Y = Y[, keep, drop = FALSE], h = H[, keep, drop = FALSE])
}

# Benchmarks: rolling one-step forecasts with fixed parameters
spec_u <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                     mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                     distribution.model = "std")
spec_u_mean <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                          mean.model = list(armaOrder = c(0, 0)), distribution.model = "std")  # mean model as in the W construction

to_xts <- function(X) {
  m <- t(X); colnames(m) <- paste0("A", seq_len(nrow(X)))
  xts(m, order.by = seq(as.Date("1980-01-01"), by = "months", length.out = ncol(X)))
}

timed <- function(expr, timeout) {
  t0 <- Sys.time()
  out <- tryCatch(withTimeout(expr, timeout = timeout, onTimeout = "silent"), error = function(e) NULL)
  list(value = out, sec = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# Marginal variance path of a DCC/CCC fit: h_{t+1} = omega + alpha e_t^2 + beta h_t,
# updated with the realised data (parameters fixed).
dcc_onestep <- function(fit, Xfull, Tt, H) {
  cf <- coef(fit); sig <- sigma(fit); N <- nrow(Xfull)
  out <- matrix(NA_real_, H, N)
  for (i in seq_len(N)) {
    nm <- paste0("A", i)
    om <- cf[paste0("[", nm, "].omega")]; al <- cf[paste0("[", nm, "].alpha1")]; be <- cf[paste0("[", nm, "].beta1")]
    h <- as.numeric(sig[Tt, i])^2
    for (k in seq_len(H)) {
      h <- om + al * Xfull[i, Tt + k - 1]^2 + be * h
      out[k, i] <- h
    }
  }
  out
}

# GO-GARCH rolling one-step forecasts with fixed parameters, through rmgarch's
# own mechanism: the model is fitted on the full series with out.sample = H, so
# estimation (including the ICA rotation) uses the first Tt observations only,
# and gogarchforecast(n.ahead = 1, n.roll = H - 1) filters the hold-out period
# with those parameters.
go_rolling <- function(fit, H) {
  fc <- gogarchforecast(fit, n.ahead = 1, n.roll = H - 1)
  rc <- rcov(fc)
  t(sapply(seq_len(H), function(k) diag(as.matrix(rc[[k]][, , 1]))))
}

# Check of the first rolling forecast, computed from the fitted object:
# H_{T+1} = A diag(h_{T+1}) A', with factors f_T = A^{-1} (x_T - mu). Returns
# the relative differences with go_rolling() and with the fitted covariance.
go_check <- function(fit, X, Tt, pv1) {
  tryCatch({
    A   <- as.matrix(fit@mfit$A)
    uf  <- fit@mfit$ufit
    sig <- as.matrix(sigma(uf)); cf <- as.matrix(coef(uf))
    rcT <- rcov(fit)[, , Tt]
    c_in <- max(abs(A %*% diag(sig[Tt, ]^2) %*% t(A) - rcT)) / max(abs(rcT))
    mu  <- colMeans(t(X[, 1:Tt]) - as.matrix(residuals(fit)))
    fT  <- as.numeric(solve(A, X[, Tt] - mu))
    h1  <- cf["omega", ] + cf["alpha1", ] * fT^2 + cf["beta1", ] * sig[Tt, ]^2
    v1  <- diag(A %*% diag(h1) %*% t(A))
    c(in_sample = c_in, first_forecast = max(abs(v1 - pv1) / pv1))
  }, error = function(e) c(in_sample = NA, first_forecast = NA))
}

# One replication
spec_norm <- function(W) { diag(W) <- 0; W / max(Mod(eigen(W, only.values = TRUE)$values)) }
bounded_w <- function(R) spec_norm(1 / (2 * (1 - R) + 1))

run_rep <- function(dgp, N, r, Tt = 450, H = 52, whiten = TRUE,
                    timeout_dcc = 1800, timeout_go = 900) {
  p   <- if (dgp == "dcc") make_dcc_params(N) else make_net_params(N)
  sim <- if (dgp == "dcc") simulate_dcc(N, Tt + H + 1, p, seed = 1e5 * N + r)
         else              simulate_network(N, Tt + H + 1, p, seed = 2e5 * N + r)
  Yraw <- sim$Y                                             # N x (Tt + H + 1)

  # Modelling data: whitened or raw returns
  if (whiten) {
    wf <- fit_whitening(Yraw[, 1:(Tt + 1)], Yraw, mean_model = "AR1")
    X  <- wf$residuals[, -1, drop = FALSE]; Wh <- wf$whitening_matrix
  } else { X <- Yraw[, -1, drop = FALSE]; Wh <- diag(N) }
  Yw <- Yraw[, -1, drop = FALSE]                            # raw returns aligned with X (W source)

  # Oracle: true one-step conditional variance of the modelled series
  oracle <- NULL
  if (dgp == "dcc") {
    oracle <- t(sapply(seq_len(H), function(k) diag(Wh %*% sim$Sig[[Tt + 1 + k]] %*% t(Wh))))
  } else if (!whiten) {
    oracle <- t(sim$h[, Tt + 1 + seq_len(H)])
  }

  eps  <- as.numeric(quantile(X[, 1:Tt][X[, 1:Tt] != 0]^2, 0.01))
  act  <- t(X[, Tt + seq_len(H)]^2)                          # H x N realised
  alog <- log(pmax(act, eps)); av <- exp(alog)
  pv <- list(); tm <- list(); fail <- list()

  # --- Weights from raw returns of the training window
  W_corr <- bounded_w(cor(t(Yw[, 1:Tt])))
  ccc_w  <- timed(dccfit(dccspec(multispec(replicate(N, spec_u_mean)), dccOrder = c(0, 0),
                                 distribution = "mvt"), data = to_xts(Yw[, 1:Tt])), timeout_dcc)
  # dccfit() returns a uGARCHmultifit when a univariate fit fails, so the class is checked
  W_ccc  <- if (!inherits(ccc_w$value, "DCCfit")) NULL else
              tryCatch(bounded_w(rcor(ccc_w$value)[, , 1]), error = function(e) NULL)

  # --- Network models
  Ys <- log(pmax(X[, 1:Tt]^2, eps))
  for (nm in c("Net_Corr", "Net_CCC")) {
    W <- if (nm == "Net_Corr") W_corr else W_ccc
    if (is.null(W)) { fail[[nm]] <- TRUE; next }
    g <- timed(GMM_SDPD_2SLS_ARCH_ind_timelags(Ys, NULL, W, list(ksy = 10, ksx = 10, stl = 0, tl = 1, ted = 0)), 600)
    tm[[nm]] <- g$sec + if (nm == "Net_CCC") ccc_w$sec else 0      # Net_CCC pays for its CCC fit
    fail[[nm]] <- is.null(g$value)
    if (fail[[nm]]) next
    fc <- t(sapply(seq_len(H), function(k)
      as.numeric(forecast_gmm(g$value, W, log(pmax(X[, Tt + k - 1]^2, eps))))))
    fe <- alog - fc
    if (nm == "Net_Corr") pv$Net_Corr_rec <- predicted_var(fe, alog, nm)   # recursive smearing
    attr(fe, "smear") <- g$value$smear                       # in-sample smearing, Eq. (37)
    pv[[nm]] <- predicted_var(fe, alog, nm)
  }

  # --- DCC (marginal forecasts = CCC) and GO-GARCH on the modelling data
  d <- timed(dccfit(dccspec(multispec(replicate(N, spec_u)), dccOrder = c(1, 1), distribution = "mvt"),
                    data = to_xts(X[, 1:Tt])), timeout_dcc)
  tm$DCC <- d$sec; fail$DCC <- !inherits(d$value, "DCCfit")
  if (!fail$DCC) pv$DCC <- tryCatch(dcc_onestep(d$value, X, Tt, H), error = function(e) NULL)
  if (is.null(pv$DCC)) fail$DCC <- TRUE
  go <- timed(gogarchfit(gogarchspec(multispec(replicate(N, spec_u)), distribution.model = "mvnorm"),
                         data = to_xts(X[, 1:(Tt + H)]), out.sample = H), timeout_go)
  tm$GO <- go$sec; fail$GO <- !inherits(go$value, "goGARCHfit")
  if (!fail$GO) pv$GO <- tryCatch(go_rolling(go$value, H), error = function(e) NULL)
  go_chk <- if (!is.null(pv$GO)) go_check(go$value, X, Tt, pv$GO[1, ]) else c(in_sample = NA, first_forecast = NA)
  if (any(go_chk > 1e-4, na.rm = TRUE))
    cat(sprintf("    !! GO check: in-sample %.2e, first forecast %.2e\n", go_chk[1], go_chk[2]))
  if (is.null(pv$GO)) fail$GO <- TRUE
  if (!is.null(oracle)) pv$Oracle <- oracle

  pv <- pv[!sapply(pv, is.null)]
  ql <- lapply(pv, function(P) compute_qlike(av, P))
  list(dgp = dgp, N = N, rep = r, whiten = whiten,
       qlike = sapply(ql, mean, na.rm = TRUE),
       qlike_t = sapply(ql, rowMeans, na.rm = TRUE),           # per origin, for paired tests
       mse = sapply(pv, function(P) mean((av - P)^2, na.rm = TRUE)),
       time = unlist(tm), fail = unlist(fail), go_check = go_chk)
}

# Study driver (checkpointed) and summary
run_accuracy_study <- function(dgp, N_vals = c(20, 50, 100), R_vals = c(10, 10, 5), whiten = TRUE,
                               file = file.path(sim_dir, sprintf("sim_accuracy_v3_%s_%s.rds", dgp,
                                                                 if (whiten) "wh" else "raw"))) {
  res <- if (file.exists(file)) readRDS(file) else list()
  for (k in seq_along(N_vals)) for (r in seq_len(R_vals[k])) {
    key <- sprintf("N%d_r%d", N_vals[k], r)
    if (!is.null(res[[key]])) next
    t0 <- Sys.time()
    res[[key]] <- run_rep(dgp, N_vals[k], r, whiten = whiten)
    saveRDS(res, file)
    q <- res[[key]]$qlike
    cat(sprintf("[%s N=%d rep %d] %.0fs | QLIKE %s\n", dgp, N_vals[k], r,
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                paste(sprintf("%s %.3f", names(q), q), collapse = " | ")))
    flush.console()
  }
  res
}

summarize_accuracy <- function(res) {
  Ns <- sort(unique(sapply(res, `[[`, "N")))
  models <- c("Oracle", "Net_Corr", "Net_Corr_rec", "Net_CCC", "DCC", "GO")
  do.call(rbind, lapply(Ns, function(N) {
    sub <- res[sapply(res, function(x) x$N == N)]
    g <- function(field, m) mean(sapply(sub, function(x) if (m %in% names(x[[field]])) x[[field]][[m]] else NA), na.rm = TRUE)
    fr <- function(m) mean(sapply(sub, function(x) if (m %in% names(x$fail)) x$fail[[m]] else NA), na.rm = TRUE)
    # Paired comparison across all replications and origins: Net_Corr vs DCC (QLIKE)
    dq <- unlist(lapply(sub, function(x) if (all(c("Net_Corr", "DCC") %in% colnames(x$qlike_t)))
      x$qlike_t[, "Net_Corr"] - x$qlike_t[, "DCC"] else NULL))
    win <- mean(sapply(sub, function(x) if (all(c("Net_Corr", "DCC") %in% names(x$qlike)))
      x$qlike[["Net_Corr"]] < x$qlike[["DCC"]] else NA), na.rm = TRUE)
    data.frame(N = N, R = length(sub),
               setNames(as.list(sapply(models, function(m) g("qlike", m))), paste0("QLIKE_", models)),
               setNames(as.list(sapply(models[-1], function(m) g("mse", m))), paste0("MSE_", models[-1])),
               NetCorr_minus_DCC = if (length(dq)) mean(dq) else NA,
               p_paired = if (length(dq) > 2) t.test(dq)$p.value else NA,
               NetCorr_wins = win,
               Time_Net = g("time", "Net_Corr"), Time_NetCCC = g("time", "Net_CCC"),
               Time_DCC = g("time", "DCC"), Time_GO = g("time", "GO"),
               GO_check_max = suppressWarnings(max(unlist(lapply(sub, `[[`, "go_check")), na.rm = TRUE)),
               Fail_NetCCC = fr("Net_CCC"), Fail_DCC = fr("DCC"), Fail_GO = fr("GO"), check.names = FALSE)
  }))
}
res_dcc <- run_accuracy_study("dcc")
tab_dcc <- summarize_accuracy(res_dcc)
cat("\n=== Table D.3, Panel A: DCC-GARCH DGP ===\n")
print(tab_dcc, digits = 4, row.names = FALSE)
write.csv(tab_dcc, file.path(sim_dir, "Table_sim_accuracy_DCC_DGP.csv"), row.names = FALSE)

res_net <- run_accuracy_study("network", whiten = FALSE)   # oracle defined on the raw scale
tab_net <- summarize_accuracy(res_net)
cat("\n=== Table D.3, Panel B: network log-ARCH DGP ===\n")
print(tab_net, digits = 4, row.names = FALSE)
write.csv(tab_net, file.path(sim_dir, "Table_sim_accuracy_network_DGP.csv"), row.names = FALSE)

# Paired QLIKE comparisons with DCC- and GO-GARCH for both smearing factors
summarize_paired <- function(res, label) {
  Ns <- sort(unique(sapply(res, `[[`, "N")))
  do.call(rbind, lapply(Ns, function(N) {
    sub <- res[sapply(res, function(x) x$N == N)]
    one <- function(m, b) {
      d <- unlist(lapply(sub, function(x) if (all(c(m, b) %in% colnames(x$qlike_t)))
        x$qlike_t[, m] - x$qlike_t[, b] else NULL))
      w <- mean(sapply(sub, function(x) if (all(c(m, b) %in% names(x$qlike)))
        x$qlike[[m]] < x$qlike[[b]] else NA), na.rm = TRUE)
      c(diff = if (length(d)) mean(d) else NA, p = if (length(d) > 2) t.test(d)$p.value else NA, wins = w)
    }
    a <- one("Net_Corr", "DCC"); b <- one("Net_Corr_rec", "DCC")
    c2 <- one("Net_Corr", "GO"); d2 <- one("Net_Corr_rec", "GO")
    data.frame(DGP = label, N = N, R = length(sub),
               InSample_vs_DCC = a["diff"], p_InSample_DCC = a["p"], wins_InSample_DCC = a["wins"],
               Rec_vs_DCC = b["diff"], p_Rec_DCC = b["p"], wins_Rec_DCC = b["wins"],
               InSample_vs_GO = c2["diff"], Rec_vs_GO = d2["diff"], row.names = NULL)
  }))
}
tab_paired <- rbind(summarize_paired(res_dcc, "DCC-GARCH"), summarize_paired(res_net, "Network log-ARCH"))
print(tab_paired, digits = 3, row.names = FALSE)
write.csv(tab_paired, file.path(sim_dir, "Table_sim_paired_both_smearing.csv"), row.names = FALSE)


# ==============================================================================
# 15. Extended panel with Iraq and Kuwait (Section 6.3.7, Appendix G)
# ==============================================================================

# Eight producers from February 1998 (323 log-returns, 322 whitened months).
# T0 = 270 gives the same 52 evaluation months as the main panel. Run in its
# own environment so that no object of the main panel is overwritten. Outputs
# and stored fits are written to "Robustness N8".
results_N8 <- local({

  T0_N8 <- 270
  n8_dir <- "Robustness N8"
  if (!dir.exists(n8_dir)) dir.create(n8_dir, recursive = TRUE)

  lr8 <- read.csv("log_returns_data_8.csv", check.names = FALSE)
  Y_raw8 <- t(as.matrix(lr8[, 2:9]))
  rownames(Y_raw8) <- make.names(colnames(lr8)[2:9])
  Y_W8 <- Y_raw8[, -1, drop = FALSE]
  N <- nrow(Y_raw8)
  stopifnot(N == 8, ncol(Y_raw8) == 323, !anyNA(Y_raw8))
  start8 <- "1998-03-01"

  wh <- fit_whitening(Y_in = Y_raw8[, 1:(T0_N8 + 1), drop = FALSE], Y = Y_raw8, mean_model = "AR1")
  E  <- wh$residuals[, -1, drop = FALSE]
  dimnames(E) <- list(rownames(Y_raw8), colnames(Y_raw8)[-1])
  cat("Max |off-diagonal correlation| after whitening (training window):",
      round(max(abs(wh$correlation_after[upper.tri(wh$correlation_after)])), 4), "\n")
  eps <- as.numeric(quantile(E[E != 0]^2, 0.01))

  W_list <- generate_weights_for_T0(Y_W8, T0_N8, start_date = start8,
                                    go_cache_file = file.path(n8_dir, sprintf("N8_W_go_cached_T0_%d.rds", T0_N8)))
  W_all <- list(Euclidean = W_list$Euclidean, Correlation = W_list$Correlation,
                Piccolo = W_list$Piccolo, Net_CCC = W_list$CCC, Net_DCC = W_list$DCC,
                Net_GO = W_list$GO, Net_DY = W_list$DY)
  failed <- names(W_all)[sapply(W_all, is.null)]
  if (length(failed)) stop("Weight construction failed for: ", paste(failed, collapse = ", "))
  saveRDS(W_all, file.path(n8_dir, "N8_weight_matrices.rds"))

  # GMM estimates (Table G.1)
  Y_star <- log(pmax(E[, 1:T0_N8]^2, eps))
  fits <- lapply(W_all, function(W)
    GMM_SDPD_2SLS_ARCH_ind_timelags(Y = Y_star, X = NULL, W = W,
                                    info = list(ksy = 10, ksx = 10, stl = 0, tl = 1, ted = 0)))
  tabG1 <- t(sapply(names(fits), function(nm) {
    f <- fits[[nm]]; rho <- f$theta[1]
    c(rho = rho, t_rho = f$tstat[1],
      spec_radius_rhoW = max(Mod(eigen(rho * W_all[[nm]], only.values = TRUE)$values)))
  }))
  cat("\n=== Table G.1 ===\n"); print(round(tabG1, 4))

  # Network forecasts, 52 origins, W and parameters fixed at T0
  n_fc <- ncol(E) - T0_N8
  ferr_net <- lapply(names(fits), function(nm) {
    f <- fits[[nm]]; W <- W_all[[nm]]
    err <- matrix(NA_real_, n_fc, N, dimnames = list(NULL, rownames(E)))
    for (i in seq_len(n_fc)) {
      tt <- T0_N8 + i - 1
      pv <- log(pmax(E[, tt]^2, eps)); act <- log(pmax(E[, tt + 1]^2, eps))
      err[i, ] <- act - as.numeric(forecast_gmm(f, W, pv))
    }
    attr(err, "smear") <- f$smear
    err
  })
  names(ferr_net) <- names(fits)

  # MGARCH benchmarks on the whitened panel, re-estimated at each origin
  spec_univ_b <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                            mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                            distribution.model = "std")
  ccc_spec <- dccspec(uspec = multispec(replicate(N, spec_univ_b)), dccOrder = c(0, 0), distribution = "mvt")
  dcc_spec <- dccspec(uspec = multispec(replicate(N, spec_univ_b)), dccOrder = c(1, 1), distribution = "mvt")
  go_spec  <- gogarchspec(multispec(replicate(N, spec_univ_b)), distribution.model = "mvnorm")

  ccc_forecast_sigma <- function(fit, asset_names) {
    cf <- coef(fit); sig <- sigma(fit); res <- residuals(fit); Tt <- nrow(sig)
    vapply(seq_along(asset_names), function(k) {
      nm <- asset_names[k]
      sqrt(unname(cf[paste0("[", nm, "].omega")] +
                  cf[paste0("[", nm, "].alpha1")] * as.numeric(res[Tt, k])^2 +
                  cf[paste0("[", nm, "].beta1")]  * as.numeric(sig[Tt, k])^2))
    }, numeric(1))
  }
  # Each function returns the one-step forecast of the return volatilities and
  # the centring constant c = training-window mean of ln max(e^2, eps) - ln sigma^2,
  # as in the main-panel MGARCH script.
  # A failure of the centring step only affects the centred errors, not the
  # uncentred forecast (so the variance-scale results are unchanged).
  centring <- function(resid, sigma_in) tryCatch(
    colMeans(log(pmax(as.matrix(resid)^2, eps)) - log(as.matrix(sigma_in)^2), na.rm = TRUE),
    error = function(e) rep(NA_real_, N))
  fit_fc <- list(
    CCC = function(X) { f <- dccfit(ccc_spec, data = X)
                        list(sfc = ccc_forecast_sigma(f, colnames(X)), c = centring(residuals(f), sigma(f))) },
    DCC = function(X) { f <- dccfit(dcc_spec, data = X)
                        list(sfc = as.numeric(sigma(dccforecast(f, n.ahead = 1))[, , 1]),
                             c = centring(residuals(f), sigma(f))) },
    GO  = function(X) { f <- gogarchfit(go_spec, data = X)
                        cv <- rcov(gogarchforecast(f, n.ahead = 1)); if (is.list(cv)) cv <- cv[[1]]
                        sfc <- sqrt(diag(as.matrix(cv[, , 1])))    # return, not factor, volatilities
                        cc <- tryCatch({
                          ci  <- rcov(f)
                          res <- tryCatch(as.matrix(residuals(f)), error = function(e) as.matrix(X))
                          centring(res, sqrt(t(apply(ci, 3, diag))))
                        }, error = function(e) rep(NA_real_, N))
                        list(sfc = sfc, c = cc) })

  X_all <- xts(t(E), order.by = seq(as.Date(start8), by = "months", length.out = ncol(E)))
  ferr_mg <- list(); ferr_mg_c <- list()
  for (m in names(fit_fc)) {
    set.seed(123)    # the ICA step of GO-GARCH starts from random values
    f_csv  <- file.path(n8_dir, sprintf("N8_fERR_%s_T0_%d_eps_one_pct.csv", m, T0_N8))
    fc_csv <- file.path(n8_dir, sprintf("N8_fERR_%s_T0_%d_eps_one_pct_centred.csv", m, T0_N8))
    if (file.exists(f_csv) && file.exists(fc_csv)) {
      ferr_mg[[m]] <- as.matrix(read.csv(f_csv)); ferr_mg_c[[m]] <- as.matrix(read.csv(fc_csv)); next
    }
    err <- err_c <- matrix(NA_real_, n_fc, N, dimnames = list(NULL, rownames(E)))
    pb <- txtProgressBar(min = 0, max = n_fc, style = 3)
    for (i in seq_len(n_fc)) {
      tt  <- T0_N8 + i - 1
      out <- tryCatch(fit_fc[[m]](X_all[1:tt, ]), error = function(e) NULL)
      if (!is.null(out) && all(is.finite(out$sfc)) && length(out$sfc) == N) {
        err[i, ]   <- log(pmax(E[, tt + 1]^2, eps)) - log(pmax(out$sfc^2, eps))
        err_c[i, ] <- err[i, ] - out$c
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)
    cat(m, ": failed origins =", sum(!complete.cases(err)), "; without centring constant =",
        sum(!complete.cases(err_c)) - sum(!complete.cases(err)), "\n")
    write.csv(err, f_csv, row.names = FALSE); write.csv(err_c, fc_csv, row.names = FALSE)
    ferr_mg[[m]] <- err; ferr_mg_c[[m]] <- err_c
  }

  # Evaluation (Tables G.2 and G.3), benchmark Net_GO
  actual_log <- matrix(NA_real_, n_fc, N)
  for (i in seq_len(n_fc)) actual_log[i, ] <- log(pmax(E[, T0_N8 + i]^2, eps))
  av <- exp(actual_log)
  bench <- "Net_GO"
  mcs_names <- function(m) if (!is.null(m@Info$model.names)) m@Info$model.names else rownames(m@show)

  tabG2 <- t(sapply(names(ferr_net), function(nm) {
    e  <- ferr_net[[nm]]
    pr <- if (nm == bench) list(p.value = NA, mean_diff = NA) else panel_robust_dm_test(e, ferr_net[[bench]])
    c(RMSFE = sqrt(mean(e^2, na.rm = TRUE)), MAFE = mean(abs(e), na.rm = TRUE),
      Bias = mean(e, na.rm = TRUE), LossDiff_vs_NetGO = pr$mean_diff, DM_panel_p = pr$p.value)
  }))
  set.seed(123)
  mcs_log <- MCS::MCSprocedure(sapply(ferr_net, function(e) rowMeans(e^2, na.rm = TRUE)),
                               alpha = 0.05, B = 5000, statistic = "Tmax", verbose = FALSE)
  tabG2 <- cbind(tabG2, In_MCS = as.numeric(rownames(tabG2) %in% mcs_names(mcs_log)))

  # Centred MGARCH benchmarks on the log scale, tested against Net. GO-GARCH
  tabG2_mg <- t(sapply(c("DCC", "GO"), function(m) {
    e  <- ferr_mg_c[[m]]
    pr <- panel_robust_dm_test(e, ferr_net[[bench]])
    c(RMSFE = sqrt(mean(e^2, na.rm = TRUE)), MAFE = mean(abs(e), na.rm = TRUE),
      Bias = mean(e, na.rm = TRUE), LossDiff_vs_NetGO = pr$mean_diff, DM_panel_p = pr$p.value)
  }))
  rownames(tabG2_mg) <- c("Std_CCC_DCC_centred", "Std_GO_centred")
  set.seed(123)
  L_all <- sapply(c(ferr_net, list(DCC_c = ferr_mg_c$DCC, GO_c = ferr_mg_c$GO)),
                  function(e) rowMeans(e^2))
  ok <- complete.cases(L_all)    # origins at which GO-GARCH failed are dropped for all models
  cat("Log-scale MCS with the centred benchmarks uses", sum(ok), "of", nrow(L_all), "origins\n")
  mcs_log_all <- MCS::MCSprocedure(L_all[ok, ], alpha = 0.05, B = 5000, statistic = "Tmax", verbose = FALSE)
  cat("Log-scale MCS with the centred benchmarks retains:", paste(mcs_names(mcs_log_all), collapse = ", "), "\n")
  cat("CCC and DCC centred errors identical:",
      isTRUE(all.equal(ferr_mg_c$CCC, ferr_mg_c$DCC, check.attributes = FALSE)), "\n")

  all_e <- c(ferr_net, ferr_mg)
  Lq <- lapply(names(all_e), function(nm) compute_qlike(av, predicted_var(all_e[[nm]], actual_log, nm)))
  Lm <- lapply(names(all_e), function(nm) (av - predicted_var(all_e[[nm]], actual_log, nm))^2)
  names(Lq) <- names(Lm) <- names(all_e)
  tabG3 <- t(sapply(names(all_e), function(nm) {
    pr <- if (nm == bench) list(p.value = NA, mean_diff = NA) else panel_robust_dm_loss(Lq[[nm]], Lq[[bench]])
    c(QLIKE = mean(Lq[[nm]], na.rm = TRUE), MSE_var = mean(Lm[[nm]], na.rm = TRUE),
      QLIKE_diff_vs_NetGO = pr$mean_diff, DM_panel_p_QLIKE = pr$p.value)
  }))
  keep <- c("Net_CCC", "Net_DCC", "Net_GO", "Net_DY", "CCC", "GO")
  Lq_t <- sapply(Lq[keep], function(L) rowMeans(L, na.rm = TRUE))
  Lq_t <- Lq_t[, apply(Lq_t, 2, function(x) all(is.finite(x))), drop = FALSE]
  set.seed(123)
  mcs_q <- MCS::MCSprocedure(Lq_t, alpha = 0.05, B = 5000, statistic = "Tmax", verbose = FALSE)
  tabG3 <- cbind(tabG3, In_MCS_QLIKE = ifelse(rownames(tabG3) %in% colnames(Lq_t),
                                              as.numeric(rownames(tabG3) %in% mcs_names(mcs_q)), NA))

  cat("\n=== Table G.1: log scale ===\n");      print(round(tabG2, 4))
  cat("\n=== Table G.1: centred MGARCH benchmarks ===\n"); print(round(tabG2_mg, 4))
  write.csv(tabG2_mg, file.path(n8_dir, "TableG2_log_scale_MGARCH_centred.csv"))
  cat("\n=== Table G.1: variance scale ===\n"); print(round(tabG3, 4))
  write.csv(tabG1, file.path(n8_dir, "TableG1_rho.csv"))
  write.csv(tabG2, file.path(n8_dir, "TableG2_log_scale.csv"))
  write.csv(tabG3, file.path(n8_dir, "TableG3_variance_scale.csv"))
  list(rho = tabG1, log_scale = tabG2, log_scale_mgarch = tabG2_mg, variance_scale = tabG3)
})
