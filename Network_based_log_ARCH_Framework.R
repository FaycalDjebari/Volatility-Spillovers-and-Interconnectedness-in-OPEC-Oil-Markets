#======= Load the data ========

library(MASS)
library(Rsolnp)
library(forecast)
library(MTS)
library(TSclust)

# =============================
# Set working directory
# =============================

setwd("/Users/faycal/Library/CloudStorage/Dropbox/R code for Network log-ARCH framework/")

log_returns <- read.csv("log_returns_data.csv")



# ======== Load GMM + QML estimation function ============

GMM_SDPD_2SLS_ARCH_ind_timelags <- function(Y, X, W, info){
  
  ksy = info$ksy # spatial exapansion order of y
  ksx = info$ksx # spatial exapansion order of x
  
  if(is.null(Y)){
    stop("Y is missing")
  }
  # if(is.null(X)){
  # stop("X is missing")
  # }
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
  
  # cat("Number of cross-sectional units:", n,  "\n")
  # cat("Length of the time series:", t,  "\n")
  # cat("Number of spatial lags (weight matrices):", p,  "\n")
  
  
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
  # SIG <- 1/ns * solve(DSiD);
  
  std <- sqrt(abs(diag(SIG)));
  tstat <- theta/std;
  
  SIG_alt <- tryCatch(1/ns * solve(DSiD_alt), error = function(e){cat("DSiD_alt not invertible \n"); return(array(NA, dim = dim(DSiD_alt)))})
  # SIG_alt <- 1/ns * solve(DSiD_alt);
  std_alt <- sqrt(abs(diag(SIG_alt)));
  tstat_alt <- theta_alt/std_alt;
  
  results <- list(theta = theta, std = std, SIG = SIG, tstat = tstat, sigma2 = sigma2,
                  theta_alt = theta_alt, std_alt = std_alt, SIG_alt = SIG_alt, tstat_alt = tstat_alt, sigma2_alt = sigma2_alt,
                  e = array(e, dim = c(n, t)), e_alt = array(e_alt, dim = c(n, t)), hyt = array(hyt, dim = c(n, t)))
  
  return(results)
  
}


mldivide <- function(A, b){
  return(ginv(t(A) %*% A) %*% t(A) %*% b)
  # return(qr.solve(A, b))
}


# ========== Analysis of the oil price log-returns ==========

# ----- Descriptive analysis

plot(as.Date(log_returns$Date), log_returns$Algeria, type = "l")
lines(as.Date(log_returns$Date), log_returns$Iran, col = "red")
lines(as.Date(log_returns$Date), log_returns$Libya, col = "blue")
lines(as.Date(log_returns$Date), log_returns$Nigeria, col = "darkgreen")
lines(as.Date(log_returns$Date), log_returns$Saudi_arabia, col = "orange")
lines(as.Date(log_returns$Date), log_returns$UAE, col = "purple")

cor(log_returns[,2:7])


#  ========= Construction of the Weight matrices ================

Y <- t(as.matrix(log_returns[,2:7]))


# Alternative A: standard Euclidean distance
Wmat1 <- as.matrix(diss(Y, "EUCL"))
Wmat1 <- 1/Wmat1
diag(Wmat1) <- 0
Wmat1 <- Wmat1 / max(eigen(Wmat1, only.values = TRUE)$values)
print(Wmat1)



# Alternative B: correlation-based distance
Wmat2 <- as.matrix(diss(Y, "COR"))
Wmat2 <- 1/Wmat2
diag(Wmat2) <- 0
Wmat2 <- Wmat2 / max(eigen(Wmat2, only.values = TRUE)$values)
print(Wmat2)

# ------- Alternative C: log-ARCH approach -------- 
Wmat3 <- as.matrix(diss(Y, "AR.PIC"))
Wmat3 <- 1/Wmat3
diag(Wmat3) <- 0
Wmat3 <- Wmat3 / max(eigen(Wmat3, only.values = TRUE)$values)
print(Wmat3)


n <- dim(Y)[1]
t <- dim(Y)[2]


#=================Univariate and Multivariate GARCH Models=================

# ====== Load Required Libraries ========
library(rmgarch)
library(rugarch)
library(xts)
library(reshape2)
library(dplyr)



# ====== Data Preparation ========

start_date <- as.Date("1983-02-01")
date_index <- seq(from = start_date, by = "months", length.out = ncol(Y))
Y_xts <- xts(t(Y), order.by = date_index)
print(Y_xts)





# ====== Define Univariate GARCH(1,1) Model ========
spec_univ <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0)),
  distribution.model = "std"
)


# === Fit GARCH(1,1) to Each Series and Collect Results ===

fit_list <- list()
results_df <- data.frame(
  Variable = character(),
  Mu = numeric(),
  Omega = numeric(),
  Alpha1 = numeric(),
  Beta1 = numeric(),
  Shape = numeric(),        
  LogLik = numeric(),
  stringsAsFactors = FALSE
)


pval_df <- data.frame(
  Variable = character(),
  Mu = numeric(),
  Omega = numeric(),
  Alpha1 = numeric(),
  Beta1 = numeric(),
  Shape = numeric(),       
  stringsAsFactors = FALSE
)


for (i in 1:ncol(Y_xts)) {
  series_name <- colnames(Y_xts)[i]
  cat("Fitting series:", series_name, "\n")
  
  fit <- ugarchfit(spec = spec_univ, data = Y_xts[, i])
  fit_list[[series_name]] <- fit
  
  coef_vals <- coef(fit)
  coef_table <- as.data.frame(fit@fit$matcoef)
  
  
  results_df <- rbind(results_df, data.frame(
    Variable = series_name,
    Mu = coef_vals["mu"],
    Omega = coef_vals["omega"],
    Alpha1 = coef_vals["alpha1"],
    Beta1 = coef_vals["beta1"],
    Shape = coef_vals["shape"],      
    LogLik = likelihood(fit)
  ))
  
  
  pval_df <- rbind(pval_df, data.frame(
    Variable = series_name,
    Mu = coef_table["mu", "Pr(>|t|)"],
    Omega = coef_table["omega", "Pr(>|t|)"],
    Alpha1 = coef_table["alpha1", "Pr(>|t|)"],
    Beta1 = coef_table["beta1", "Pr(>|t|)"],
    Shape = coef_table["shape", "Pr(>|t|)"]   
  ))
}


# Format Final Table: Estimate (p-value) for Each Parameter 
combined_df <- data.frame(Variable = results_df$Variable)
for (param in c("Mu", "Omega", "Alpha1", "Beta1", "Shape")) {
  combined_df[[param]] <- sprintf("%.4f (%.4f)", results_df[[param]], pval_df[[param]])
}
combined_df$LogLik <- round(results_df$LogLik, 4)

# Print Final Table
print(combined_df)






# =========== CCC-GARCH Model using DCC framework ========

# ===== CCC-GARCH Model =====
spec_ccc <- dccspec(
  uspec = multispec(replicate(ncol(Y_xts), spec_univ)),  
  dccOrder = c(0, 0),
  distribution = "mvt"
)
ccc_fit <- dccfit(spec_ccc, data = Y_xts)
print(ccc_fit)
coef(ccc_fit)

# Extract CONSTANT correlation matrix (no time dimension)
ccc_corr <- rcor(ccc_fit)[,,1] 
print(ccc_corr)
# Calculate distances and weights
distance_matrix_ccc <- sqrt(2 * (1 - ccc_corr))
Wmat_ccc <- 1 / (distance_matrix_ccc^2 + 1)  

# Standardization 
diag(Wmat_ccc) <- 0  
Wmat_ccc <- Wmat_ccc / max(eigen(Wmat_ccc, symmetric=TRUE, only.values=TRUE)$values)

print(Wmat_ccc)



# ====== DCC-GARCH Model ========

spec_dcc <- dccspec(
  uspec = multispec(replicate(ncol(Y_xts), spec_univ)), 
  dccOrder = c(1, 1),
  distribution = "mvt"
)

# Fit the model
dcc_fit <- dccfit(spec_dcc, data = Y_xts)
print(dcc_fit)
coef(dcc_fit)

# Extract dynamic correlations
corr_array_dcc <- rcor(dcc_fit)  
corr_array_dcc <- aperm(corr_array_dcc, perm = c(3, 1, 2))


# Compute time-averaged correlations
avg_corr_dcc <- apply(corr_array_dcc, c(2,3), mean, na.rm = TRUE)

# Calculate distances and weights 
distance_matrix_dcc <- sqrt(2 * (1 - avg_corr_dcc))
Wmat_dcc <- 1 / (distance_matrix_dcc^2 + 1)  

# Cleanup 
diag(Wmat_dcc) <- 0  
Wmat_dcc <- (Wmat_dcc + t(Wmat_dcc)) / 2 
Wmat_dcc <- Wmat_dcc / max(eigen(Wmat_dcc, only.values=TRUE)$values) 

print(Wmat_dcc)



# ====== GO-GARCH Model ========
spec_go <- gogarchspec(
  multispec(replicate(nrow(Y), spec_univ)),
  distribution.model = "mvnorm"
)
go_fit <- gogarchfit(spec_go, data = t(Y))
print(go_fit)



likelihood(go_fit)
loglik <- likelihood(go_fit)
n <- ncol(Y)  
k <- 3 * (nrow(Y))
AIC <- (-2 * loglik + 2 * k) / n
BIC <- (-2 * loglik + k * log(n)) / n

# Print AIC and BIC
cat("AIC:", AIC, "\n")
cat("BIC:", BIC, "\n")




# Install rugarch version 1.5.3 if not already installed
# This specific version is required for reproducibility of the weight matrix
# The new version currently available on CRAN (1.5.4) does not not work properly with the GO-GARCH model
# Uncomment the following lines to install the specific version
#install.packages("remotes")
#remotes::install_version("rugarch", version = "1.5.3", repos = "http://cran.us.r-project.org", force = TRUE)



# Calculate time-averaged correlation matrix from GO-GARCH
T <- dim(rcov(go_fit))[3]  # Total time periods
R_go_list <- lapply(1:T, function(t) {
  V_t <- rcov(go_fit)[,,t]  # Extract time-varying covariance
  D_inv <- diag(1 / sqrt(diag(V_t)))  # Diagonal inverse standard deviations
  D_inv %*% V_t %*% D_inv  # Convert to correlation matrix
})
R_go_avg <- Reduce("+", R_go_list) / T  # Average correlations across time

# Compute distance matrix using averaged correlations
distance_go <- sqrt(2 * (1 - R_go_avg))

# Transform to weight matrix
W_go <- 1 / (distance_go^2 + 1)

# Standardization steps
diag(W_go) <- 0 
W_go <- (W_go + t(W_go)) / 2  
W_go <- W_go / max(eigen(W_go, only.values = TRUE)$values)

print(W_go)











#=============== Network Plots  =============================



# Required libraries
library(igraph)
library(wesanderson)
library(classInt)

# Country Names and Colors
country_names <- c("Algeria", "Iran", "Libya", "Nigeria", "Saudi Arabia", "UAE")
country_colors <- wes_palette("Zissou1", length(country_names), type = "continuous")

# Consistent Layout
get_network_layout <- function(Wmat) {
  Wmat <- as.matrix(Wmat)
  colnames(Wmat) <- country_names
  
  Wmat_sym <- (Wmat + t(Wmat)) / 2
  
  dist_matrix <- 1 / (Wmat_sym + 1e-6) 
  diag(dist_matrix) <- 0 
  
  # Create an undirected graph for layout calculation
  g_undir <- graph_from_adjacency_matrix(Wmat_sym, mode = "undirected", weighted = TRUE)
  
  layout_with_mds(g_undir, dist = dist_matrix)
}

plot_network <- function(Wmat) { 
  Wmat <- as.matrix(Wmat)
  colnames(Wmat) <- country_names
  
  # For plotting, ensure the matrix is symmetric to correctly represent undirected edges
  Wmat_plot <- (Wmat + t(Wmat)) / 2 
  
  # Create an undirected graph
  g <- graph_from_adjacency_matrix(Wmat_plot, mode = "undirected", weighted = TRUE)
  
  V(g)$color <- country_colors
  V(g)$size <- 35
  V(g)$label <- NA
  V(g)$frame.color <- "white"
  
  edge_weights <- E(g)$weight
  classes <- classIntervals(edge_weights, n = 8, style = "quantile")
  E(g)$color <- findColours(classes, pal = gray(seq(0.9, 0.1, length = 9))) 
  E(g)$width <- edge_weights * 6 
  
  plot(g,
       layout = get_network_layout(Wmat), 
       main = NULL, 
       vertex.label = NA)
}


plot_all_networks <- function(Wlist) { 
  # Define layout for multiple plots and a legend row
  layout_matrix <- matrix(c(1, 2, 3,
                            4, 5, 6,
                            7, 7, 7), nrow = 3, byrow = TRUE)
  
  layout(mat = layout_matrix, heights = c(4, 4, 1)) 
  par(mar = c(1, 1, 1, 1)) 
  
  for (i in seq_along(Wlist)) {
    plot_network(Wlist[[i]]) 
    label <- paste0("(", letters[i], ")") 
    mtext(label, side = 1, line = -1.5, adj = 0.5, cex = 1, font = 2)
  }
  
  par(mar = c(0, 0, 0, 0)) 
  plot.new()
  
  legend("center",
         legend = country_names,
         horiz = FALSE, 
         ncol = 3,    
         pch = 21,    
         pt.bg = country_colors, 
         pt.cex = 1.8, 
         cex = 0.95,   
         x.intersp = 0.8, 
         inset = c(0, 0.05), 
         bty = "n")
}


Wlist <- list(Wmat1, Wmat2, Wmat3, Wmat_ccc, Wmat_dcc, W_go)

# Plot all networks 
plot_all_networks(Wlist)









# ====== GMM-SDPD-2SLS-ARCH Model Estimation ======


# Named list of Network weight matrices
weight_matrices <- list(
  Net_CCC     = Wmat_ccc,
  Net_DCC     = Wmat_dcc,
  Net_GO      = W_go,
  Euclidean   = Wmat1,
  Correlation = Wmat2,
  Piccolo     = Wmat3
)

# Clear the environment of any previous model results
rm(list = names(weight_matrices)[names(weight_matrices) %in% ls()])

# -----Estimation and Assignment ------

for (name in names(weight_matrices)) {
  cat("Estimating model for:", name, "\n")
  
  model_output <- GMM_SDPD_2SLS_ARCH_ind_timelags(
    Y = log(Y^2 + min(Y[Y != 0]^2)),  # stabilizing log transformation
    X = NULL,
    W = weight_matrices[[name]],
    info = list(ksy = 10, ksx = 10, stl = 0, tl = 1, ted = 0)
  )
  
  # Assign result to an object named after the weight matrix label
  assign(name, model_output)
}

# Display Structured Results 

for (name in names(weight_matrices)) {
  res <- get(name)  # retrieve the model result (e.g., Net_GO)
  
  cat("\n=========== Results for", name, "===========\n")
  
  
  if (!is.null(res$theta)) {
    cat("Theta:\n"); print(res$theta)
  } else {
    cat("Theta: not found\n")
  }
  
  if (!is.null(res$std)) {
    cat("Standard Errors:\n"); print(res$std)
  } else {
    cat("Standard Errors: not found\n")
  }
  
  if (!is.null(res$tstat)) {
    cat("T-Statistics:\n"); print(round(res$tstat, 3))
  } else {
    cat("T-Statistics: not found\n")
  }
  
  if (!is.null(res$sigma2)) {
    cat("Sigma2:\n"); print(res$sigma2)
  } else {
    cat("Sigma2: not found\n")
  }
}












# ====== forecasting part of Network-based Log-ARCH Framework ======


# ====================================
# Define Relative Forecast Folder Path
# ====================================
forecast_subfolder <- "Standard MGARCH models"
# =============================
# Define models and available T0/epsilon combinations
# =============================
models <- c("ccc", "dcc", "go")
eps_filenames <- list(min = "min", one_pct = "one_pct.1%", small = "small")

# T0 where all three eps are available
T0_all_eps <- 300

# Other T0s where only "min" is available
T0_min_only <- c(200, 250, 350)

# =============================
# Helper function to read CSV safely
# =============================
read_ferr_matrix <- function(model, T0, eps) {
  folder <- paste0("fERR_", model, "_T0")
  filename <- sprintf("fERR_%s_T0_%d_eps_%s.csv", model, T0, eps_filenames[[eps]])
  path <- file.path(forecast_subfolder, folder, filename)
  if (file.exists(path)) {
    message("Reading: ", path)
    mat <- read.csv(path)
    return(as.matrix(sapply(mat, as.numeric)))
  } else {
    message("File not found: ", path)
    return(NULL)
  }
}

# =============================
# Load the available files only
# =============================
standard_fERR_list <- list()

for (model in models) {
  model_upper <- toupper(model)  
  standard_fERR_list[[model_upper]] <- list()
  
  # T0 = 300: all three eps values
  T0_chr <- as.character(T0_all_eps)
  standard_fERR_list[[model_upper]][[T0_chr]] <- list()
  for (eps in names(eps_filenames)) {
    mat <- read_ferr_matrix(model, T0_all_eps, eps)
    if (!is.null(mat)) {
      standard_fERR_list[[model_upper]][[T0_chr]][[eps]] <- mat
    }
  }
  
  # Other T0s: only eps = "min"
  for (T0 in T0_min_only) {
    T0_chr <- as.character(T0)
    standard_fERR_list[[model_upper]][[T0_chr]] <- list()
    mat <- read_ferr_matrix(model, T0, "min")
    if (!is.null(mat)) {
      standard_fERR_list[[model_upper]][[T0_chr]][["min"]] <- mat
    }
  }
}



# ===== 0. Setup and Dependencies =====
# Load required packages
dependencies <- c("forecast", "lmtest", "parallel", "future.apply", "MCS")
for (pkg in dependencies) { # Using a for loop to avoid parsing issues with lapply(..., character.only = TRUE)
  library(pkg, character.only = TRUE)
}

# reproducibility for all random processes
set.seed(123)


# ===== 1. Helper Functions =====

#' One‐step ahead forecast via GMM‐based spatial model
#' @param est_obj list returned by gmm estimation (theta vector)
#' @param W spatial weights matrix (n x n)
#' @param prev_vol numeric vector of log-squared volumes at time t
#' @return forecast for time t+1 (n-vector)
forecast_gmm <- function(est_obj, W, prev_vol) {
  rho <- est_obj$theta[1]
  gamma <- est_obj$theta[-1]
  n <- nrow(W)
  intercept_vec <- colMeans(prev_vol - rho * (W %*% prev_vol) - gamma * prev_vol,
                            na.rm = TRUE)
  solve(diag(n) - rho * W) %*% (gamma * prev_vol + intercept_vec)
}

#' Compute RMSFE and MAFE for forecasting errors matrix
#' @param err_mat matrix of errors (n_fore x n_series)
#' @return list with RMSFE (vector) and MAFE (vector)
calculate_metrics <- function(err_mat) {
  RMSFE <- apply(err_mat, 2, function(x) sqrt(mean(x^2, na.rm = TRUE)))
  MAFE <- apply(err_mat, 2, function(x) mean(abs(x), na.rm = TRUE))
  list(RMSFE = RMSFE, MAFE = MAFE)
}

#' Clark–West test for nested model comparison
#' @param e_small numeric vector of squared errors for small (nested) model
#' @param e_large numeric vector of squared errors for large model
#' @return list with statistic and p.value
cw_test <- function(e_small, e_large) {
  # Ensure inputs are numeric vectors and handle NAs consistently
  e_small <- as.vector(e_small)
  e_large <- as.vector(e_large)
  
  # Remove NA pairs
  common_nona <- !is.na(e_small) & !is.na(e_large)
  e_small <- e_small[common_nona]
  e_large <- e_large[common_nona]
  
  if (length(e_small) < 2) { # Need at least 2 observations for sd
    return(list(statistic = NA, p.value = NA))
  }
  
  d <- e_small^2 - e_large^2 + (e_large - e_small)^2
  mean_d <- mean(d, na.rm = TRUE)
  se_d <- sd(d, na.rm = TRUE) / sqrt(length(d))
  
  if (is.na(se_d) || se_d == 0) { # Handle cases where sd is NA or 0
    return(list(statistic = NA, p.value = NA))
  }
  
  cw_stat <- mean_d / se_d
  pval <- 2 * pnorm(-abs(cw_stat))
  list(statistic = cw_stat, p.value = pval)
}

#' Bootstrap confidence intervals for RMSFE & MAFE
#' @param err_mat matrix of forecast errors (n_fore x n_series)
#' @param R number of bootstrap replications (default 1000)
#' @param seed seed for reproducibility
#' @return list with RMSFE_CI and MAFE_CI
bootstrap_CI <- function(err_mat, R = 2000, seed = 123) {
  set.seed(seed)
  n_fore <- nrow(err_mat)
  rms_b <- numeric(R)
  maf_b <- numeric(R)
  for (b in seq_len(R)) {
    idx <- sample(n_fore, n_fore, replace = TRUE)
    errs <- as.vector(err_mat[idx, , drop = FALSE])
    rms_b[b] <- sqrt(mean(errs^2, na.rm = TRUE)) 
    maf_b[b] <- mean(abs(errs), na.rm = TRUE)    
  }
  list(
    RMSFE_CI = quantile(rms_b, c(0.025, 0.975), na.rm = TRUE),
    MAFE_CI = quantile(maf_b, c(0.025, 0.975), na.rm = TRUE)
  )
}


# ===== 2. Rolling‐Window Forecast Routine =====

models_list <- list(
  Euclidean = list(obj = Euclidean, W = Wmat1),
  Correlation = list(obj = Correlation, W = Wmat2),
  Piccolo = list(obj = Piccolo, W = Wmat3),
  Net_CCC = list(obj = Net_CCC, W = Wmat_ccc),
  Net_DCC = list(obj = Net_DCC, W = Wmat_dcc),
  Net_GO = list(obj = Net_GO, W = W_go)
)

#' Run rolling‐window forecasts and compute errors
#' @param Y data matrix (n_series x T_total)
#' @param models named list of model objects and weights
#' @param T0 initial window size
#' @param eps small constant for log-squared transform
#' @return list with fERR matrices, RMSFE and MAFE vectors
run_forecasts <- function(Y, models, T0, eps, standard_errors = NULL) {
  T_total <- ncol(Y)
  n_forecasts <- T_total - T0
  n_series <- nrow(Y)
  
  fERR <- lapply(models, function(x) matrix(NA, nrow = n_forecasts, ncol = n_series))
  
  pb <- txtProgressBar(min = 1, max = n_forecasts, style = 3)
  for (i in seq_len(n_forecasts)) {
    setTxtProgressBar(pb, i)
    t_train <- T0 + i - 1
    prev_vol <- log(Y[, t_train]^2 + eps)
    actual <- log(Y[, t_train + 1]^2 + eps)
    for (name in names(models)) {
      est <- models[[name]]$obj
      Wmat <- models[[name]]$W
      fcast <- forecast_gmm(est, Wmat, prev_vol)
      fERR[[name]][i, ] <- actual - fcast
    }
  }
  close(pb)
  
  # Add precomputed standard model fERRs if provided
  if (!is.null(standard_errors)) {
    for (model in names(standard_errors)) {
      fERR[[model]] <- standard_errors[[model]]
    }
  }
  
  metrics <- lapply(fERR, calculate_metrics)
  avg_RMSFE <- sapply(metrics, function(m) mean(m$RMSFE, na.rm = TRUE))
  avg_MAFE <- sapply(metrics, function(m) mean(m$MAFE, na.rm = TRUE))
  
  list(fERR = fERR, RMSFE = avg_RMSFE, MAFE = avg_MAFE)
}


# ===== 3. Experiments Over T0 and ε =====


eps_vals <- c(
  min = min(Y[Y != 0]^2),
  one_pct = as.numeric(quantile(Y[Y != 0]^2, 0.01)), 
  small = 1e-6
)
T0_vec <- c(200, 250, 300, 350)
results <- list()





results <- list()
rmsfe_values <- list()
mafe_baseline <- list()

T0 <- 300 

results[[as.character(T0)]]      <- list()
rmsfe_values[[as.character(T0)]] <- list()

for (eps_name in names(eps_vals)) {
  cat("Running forecasts with T0 =", T0, "and epsilon =", eps_name, "\n")
  
  # === Add only valid GARCH standard models ===
  garch_models <- list()
  for (model in c("CCC", "DCC", "GO")) {
    ferr_mat <- standard_fERR_list[[model]][[as.character(T0)]][[eps_name]]
    if (!is.null(ferr_mat)) {
      garch_models[[model]] <- ferr_mat
    }
  }
  
  # Run forecasts: with dynamic network + preloaded GARCH
  res <- run_forecasts(Y, models_list, T0, eps_vals[eps_name], standard_errors = garch_models)
  
  results[[as.character(T0)]][[eps_name]] <- res
  rmsfe_values[[as.character(T0)]][[eps_name]] <- round(res$RMSFE, 4)
  
  cat("  RMSFE:\n")
  print(rmsfe_values[[as.character(T0)]][[eps_name]])
  
  if (eps_name == "min") {
    mafe_baseline[["300"]] <- round(res$MAFE, 4)
    cat("\nBaseline MAFE (T0 = 300, eps = 'min'):\n")
    print(mafe_baseline[["300"]])
  }
}


# ===== 4. Additional T0 Values with eps = min =====

additional_T0 <- c(200, 250, 350)
eps_name <- "min"

for (T0 in additional_T0) {
  cat("Running forecasts with T0 =", T0, "and epsilon =", eps_name, "\n")
  garch_models <- list()
  for (model in c("CCC", "DCC", "GO")) {
    ferr_mat <- standard_fERR_list[[model]][[as.character(T0)]][[eps_name]]
    if (!is.null(ferr_mat)) {
      garch_models[[model]] <- ferr_mat
    }
  }
  
  res <- run_forecasts(Y, models_list, T0, eps_vals[[eps_name]], standard_errors = garch_models)
  results[[as.character(T0)]] <- list()
  results[[as.character(T0)]][[eps_name]] <- res
  
  rmsfe_values[[as.character(T0)]] <- list()
  rmsfe_values[[as.character(T0)]][[eps_name]] <- round(res$RMSFE, 4)
  
  cat("RMSFE for T0 =", T0, ", eps = min:\n")
  cat(paste0("  ", paste(names(res$RMSFE), collapse = " & "), " \\\\\n"))
  cat(paste0("  ", paste(sprintf("%.4f", res$RMSFE), collapse = " & "), " \\\\\n\n"))
}



# ===== Combine Standard and GMM-Based fERRs into a Unified List =====
all_ferr <- results[["300"]][["min"]]$fERR
model_names <- names(all_ferr)



# ===== 2. Bootstrap Confidence Intervals =====
ci_list <- lapply(all_ferr, function(mat) {
  bootstrap_CI(mat, R = 2000, seed = 123)
})

cat("\nBootstrap CIs (RMSFE, MAFE) for T0 = 300 (R=2000) including all models:\n")
print(ci_list)

# ===== 3. Diebold–Mariano (DM) Test =====

dm_p <- matrix(NA, length(model_names), length(model_names),
               dimnames = list(model_names, model_names))

for (i in seq_along(all_ferr)) {
  for (j in seq_along(all_ferr)) {
    if (i < j) {
      e1 <- as.vector(all_ferr[[i]])
      e2 <- as.vector(all_ferr[[j]])
      
      err_diff <- e1 - e2
      if (sd(err_diff, na.rm = TRUE) < 1e-10) {
        warning(sprintf("DM test skipped for models (%s, %s): near-zero variance in error differences.",
                        model_names[i], model_names[j]))
        next
      }
      
      dm_result <- tryCatch({
        forecast::dm.test(e1, e2, h = 1, power = 2)$p.value
      }, error = function(e) {
        warning(sprintf("DM test failed for models (%s, %s): %s",
                        model_names[i], model_names[j], e$message))
        return(NA)
      })
      
      dm_p[i, j] <- dm_result
    }
  }
}

cat("\nRobust Diebold–Mariano p-values (T0 = 300, ε = 'min'):\n")
print(round(dm_p, 4))


# ===== 4. Clark–West (CW) Test =====
cw_p <- matrix(NA, length(model_names), length(model_names),
               dimnames = list(model_names, model_names))

for (i in seq_along(all_ferr)) {
  for (j in seq_along(all_ferr)) {
    if (i < j) {
      e1 <- as.vector(all_ferr[[i]])
      e2 <- as.vector(all_ferr[[j]])
      
      # CW test adjustment
      fe1_sq <- e1^2
      fe2_sq <- e2^2
      adj_term <- fe1_sq - fe2_sq - (e1 - e2)^2
      
      # CW statistic and p-value
      cw_stat <- mean(adj_term, na.rm = TRUE) / sd(adj_term, na.rm = TRUE)
      cw_p[i, j] <- 2 * (1 - pnorm(abs(cw_stat)))
    }
  }
}

cat("\nClark–West p-values (GMM vs Standard and among all models):\n")
print(round(cw_p, 4))

# ===== 5. Model Confidence Set (MCS) Test =====
cat("\n--- Model Confidence Set (MCS) Results (alpha = 0.05) for All Models ---\n")

all_ferr <- results[["300"]][["min"]]$fERR

min_rows <- min(sapply(all_ferr, nrow))

all_err_list_aligned <- lapply(all_ferr, function(mat) {
  if (nrow(mat) > min_rows) {
    mat[1:min_rows, , drop = FALSE]
  } else {
    mat
  }
})

loss_sq_all_models <- sapply(all_err_list_aligned[model_names], function(err_mat) {
  rowMeans(err_mat^2, na.rm = TRUE)
})


set.seed(123)
mcs_result_all <- MCS::MCSprocedure(
  loss_sq_all_models,
  alpha = 0.05,
  B = 5000,
  statistic = "Tmax"
)
print(mcs_result_all)



# ===== 6. Visualization of Model Rankings =====

library(ggplot2)

model_ranks <- data.frame(
  Model = c("Net-CCC-GARCH", "Net-DCC-GARCH", "Net-GO-GARCH", 
            "Standard-DCC", "Standard-CCC"),
  Rank = c(3, 4, 2, 1, 5),
  Loss = c(5.8709, 5.8710, 5.8638, 5.7403, 6.8804),
  In_MCS = c(TRUE, TRUE, TRUE, TRUE, FALSE)
)

model_ranks <- model_ranks[order(model_ranks$Loss, decreasing = TRUE), ]
model_ranks$Model <- factor(model_ranks$Model, levels = model_ranks$Model)

ggplot(model_ranks, aes(x = Model, y = Loss, fill = In_MCS)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c("TRUE" = "#4CAF50", "FALSE" = "#B0BEC5"),
    name = expression("In MCS (" * alpha * " = 0.05)")
  ) +
  theme_minimal(base_size = 14) +
  coord_flip() +
  labs(
    y = "Average Forecast Loss",
    x = "Model"
  )




