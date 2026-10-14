# Reproducibility Guide – Network-Based Log-ARCH Framework

**Title:**  
_Forecasting Oil Volatility through Network Models with GARCH-Informed Correlation Weights_  
**Authors:** Faycal Djebari, Kahina Mehidi, Khelifa Mazouz, Philipp Otto

**Cite as:** 
Djebari, F., Mehidi, K., Mazouz, K., Otto, P. (2025). Forecasting Oil Volatility through Network Models with GARCH-Informed Correlation Weights. arXiv preprint arXiv:2507.15046. URL: https://arxiv.org/abs/2507.15046

---

## Overview

This repository contains the data and R code of the study, which forecasts the volatility of monthly OPEC oil prices with network log-ARCH models (Mattera and Otto, 2024) estimated by GMM, and compares them with CCC-, DCC- and GO-GARCH benchmarks.

The code in this version corresponds to the revised manuscript. The main changes relative to the first version are:

- All volatility models are estimated on pre-whitened returns (AR(1) filter and ZCA whitening, estimated within each training window). The weight matrices are still built from the log-returns.
- The forecast intercept of the network models follows Mattera and Otto (2024, Eq. 9), and the variance-scale forecasts use the smearing estimator.
- A seventh weight matrix, based on generalized forecast error variance decompositions (connectedness), is added.
- The MGARCH benchmarks are compared with the network models on the log-squared scale after centring, and on the variance scale with QLIKE and MSE.
- The forecast evaluation reports the standard, heteroskedasticity-robust (Harvey et al., 2024) and panel-robust Diebold-Mariano tests, and the Model Confidence Set.
- New analyses: forecast combinations, weight stability and structural breaks, economic value, simulations of estimation cost and forecast accuracy, and an extended panel with Iraq and Kuwait.

---

## Files

### Data

- `log_returns_data.csv`: monthly log-returns of the six producers, 1983 to 2024. Computed from the official OPEC prices (`t71.xlsx` in `/Official Oil Price Data from OPEC database/`).
- `log_returns_data_8.csv`: the same panel with Iraq and Kuwait, from 1998 (Appendix G).

### Scripts (run in this order)

**1. `MGARCH_benchmarks.R`**

Rolling one-step forecasts of the CCC-, DCC- and GO-GARCH benchmarks on the pre-whitened returns, for training windows T0 = 200, ..., 450. Writes to `Standard MGARCH models (whitened)/`:

- `fERR_<m>_T0/`: forecast errors of ln(sigma^2)
- `fERR_<m>_T0_centred/`: forecast errors of the centred forecasts (Section 6.3.1)
- `Computation_Metrics_<m>/`: computing time per origin (Table D.1)

This script is computationally intensive. Its output is read by the two other scripts.

**2. `Network_based_log_ARCH_Framework.R`**

Main script. Writes its tables to `Results/` and the extended panel to `Robustness N8/`.

| Section | Content | Paper |
|---|---|---|
| 1 | Packages and input files | |
| 2 | GMM estimator, forecast intercept | Eq. 17 |
| 3 | Distance measures and pre-whitening | Section 5.3 |
| 4 | Data and descriptive statistics | Tables 1-2 |
| 5 | Full-sample GARCH models and weight matrices | Tables 3, 4, C.1 |
| 6 | Weight matrices for each training window | |
| 7 | Network figures | Figures 2, 3, F.1, F.2 |
| 8 | GMM estimates at T0 = 450 | Table 5 |
| 9 | Rolling one-step forecasts | Tables 6, E.1 |
| 10 | Bootstrap intervals, DM tests, MCS | Tables 6, E.1, E.2 |
| 11 | Weight stability, subperiods, structural breaks | Table 8 |
| 12 | Economic value (volatility timing, VaR) | Table 9 |
| 13 | Simulation: estimation cost | Tables D.1, D.2 |
| 14 | Simulation: forecast accuracy | Table D.3 |
| 15 | Extended panel with Iraq and Kuwait | Table G.1 |

The seven weight matrices are Euclidean, correlation and Piccolo (on log-squared returns), time-averaged conditional correlations from CCC-, DCC- and GO-GARCH, and connectedness. All are normalized by their spectral radius.

**3. `Forecast_combination.R`**

Forecast combinations (equal weights, minimum variance, constrained OLS) on the log-squared scale (Table 7). Reads `Results/network_fERR_T0_450.rds`, saved by the main script, and the centred MGARCH errors.

### Notes

- The scripts save the forecast errors (`fERR_*.csv`) and the GO-GARCH fits (`W_go_cached_*.rds`) on the first run and reuse them afterwards. Delete them to re-estimate the models.
- The ICA step of `gogarchfit()` is not fully controlled by `set.seed()`, so GO-GARCH results from a new run can differ slightly from those in the paper.
- Random steps (GO-GARCH ICA, bootstrap, MCS, simulations) use `set.seed(123)`.
- Set the working directory at the top of each script (`setwd()`).

---

## Contact

For questions, feedback, or collaboration inquiries, please contact:

**Faycal Djebari**  
PhD Student in Quantitative Economics, University of Bejaia  
📧 Email: [faycal.djebari@univ-bejaia.dz](mailto:faycal.djebari@univ-bejaia.dz)
🔗 ORCID: [https://orcid.org/0009-0002-9265-9541](https://orcid.org/0009-0002-9265-9541)

or

**Prof. Philipp Otto**  
Professor of Statistics and Data Science, University of Glasgow  
📧 Email: [Philipp.Otto@glasgow.ac.uk](mailto:Philipp.Otto@glasgow.ac.uk)
🔗 ORCID: [https://orcid.org/0000-0002-9796-6682](https://orcid.org/0000-0002-9796-6682)

---

## License

This project is licensed under the Apache License, Version 2.0. See the [LICENSE](./LICENSE) file for details.

---

## ⚙️ How to Run the Code

### 1. Software Requirements

- **R version**: 4.2 or later
- **R packages:**

```r
install.packages(c("rugarch", "rmgarch", "xts", "vars", "MCS", "strucchange",
                   "igraph", "wesanderson", "classInt", "R.utils", "MASS",
                   "pbapply"))
```

> ⚠️ **Important:** Use `rugarch` version **1.5-3**.
> The newer version (1.5-4) introduces internal changes that may cause errors when extracting covariance matrices (`rcov`) from **GO-GARCH** models.
> To ensure full compatibility, please install the stable version:
>
> ```r
> if (!require("devtools")) install.packages("devtools")
> devtools::install_version("rugarch", version = "1.5-3")
> ```

### 2. Running

1. Run `MGARCH_benchmarks.R`.
2. Run `Network_based_log_ARCH_Framework.R`.
3. Run `Forecast_combination.R` in a fresh R session.
