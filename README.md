# Reproducibility Guide – Network-Based Log-ARCH Framework

**Preprint Title:**  
_Volatility Spillovers and Interconnectedness in OPEC Oil Markets: A Network-Based Log-ARCH Approach_  
**Authors:** Faycal Djebari, Kahina Mehidi, Khelifa Mazouz, Philipp Otto (2025)

**Cite as:** 
Djebari, F., Mehidi, K., Mazouz, K., Otto, P. (2025). Volatility Spillovers and Interconnectedness in OPEC Oil Markets: A Network-Based Log-ARCH Approach. arXiv preprint arXiv:2507.15046. URL: https://arxiv.org/abs/2507.15046

---

## Overview

This repository contains the complete R code and data preprocessing pipeline developed for the study, which investigates volatility spillovers among six selected OPEC countries using both standard Multivariate GARCH models and a novel **Network-Based Log-ARCH framework** estimated via GMM.

All results are derived from monthly official OPEC oil prices, ensuring full **reproducibility and transparency**.

---

## Folder and Script Structure

### 1. `Log Transformation and Descriptive Statistics/`

- **Objective:** Transform monthly OPEC crude oil prices into log returns.
- **Source Data:** Extracted from `T71.xlsx` in `/Official Oil Price Data from OPEC database/`.
- **Processed Data:** `Monthly_OPEC_OIL_PRICE_DATA.xlsx`

**Scripts:**

- `Calculate Log Returns.R` – Computes log returns from raw prices.
- `Plotting Log Returns and Squared Log Returns.R` – Plots time series of log and squared returns.
- `Histograms of Monthly Log Returns.R` – Generates histograms for each country.
- `Descriptive Statistics Summary.R` – Computes summary statistics (mean, sd, skewness, kurtosis, JB test).

**Outputs:**

- `log_returns_data.csv` – Primary dataset for all modeling stages.
- Figures and tables for descriptive analysis.

---

### 2. `Standard MGARCH models/`

**Scripts:**

- `CCC GARCH.R` – Forecasting error generation for CCC-GARCH
- `DCC GARCH.R` – Forecasting error generation for DCC-GARCH
- `GO GARCH.R` – Forecasting error generation for GO-GARCH

**Subfolders (auto-created):**

- `fERR_ccc_T0/` – Forecast error matrices from CCC-GARCH
- `fERR_dcc_T0/` – Forecast error matrices from DCC-GARCH
- `fERR_go_T0/` – Forecast error matrices from GO-GARCH

**Notes:**

- Forecast errors (`fERR_*.csv`) are saved and reused to save computational time.
- Rolling windows and multiple epsilon values are used for robustness.
- These scripts are computationally intensive and are separated for efficiency.

---

### 3. `Network_based_log_ARCH_Framework.R`

This is the **main script** implementing the proposed methodology.

**Key Components:**

- Loads `log_returns_data.csv`
- Constructs six types of weight matrices:
  - Euclidean, Correlation, Piccolo (AR-based)
  - Time-averaged correlations from CCC-, DCC-, and GO-GARCH
- Visualizes networks based on each matrix
- Estimates the GMM-based Network Log-ARCH models
- Performs out-of-sample forecasting with a rolling window
- Computes:
  - Bootstrap confidence intervals
  - Diebold-Mariano statistics
  - Model Confidence Set (MCS)

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

````r
install.packages(c("rugarch", "rmgarch", "tseries", "e1071", "igraph", "MCS",
                  "future", "parallel", "readxl", "dplyr", "ggplot2", "MASS",
                  "Rsolnp", "forecast", "MTS", "TSclust", "xts", "reshape2",
                  "wesanderson", "classInt", "lmtest", "future.apply"))

> ⚠️ **Important:** Use `rugarch` version **1.5-3**.
> The newer version (1.5-4) introduces internal changes that may cause errors when extracting covariance matrices (`rcov`) from **GO-GARCH** models.
> To ensure full compatibility, please install the stable version:
>
> ```r
> if (!require("devtools")) install.packages("devtools")
> devtools::install_version("rugarch", version = "1.5-3")
> ```


````
