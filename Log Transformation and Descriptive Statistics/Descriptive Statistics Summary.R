# Copyright 2025 Fayçal Djebari
# Licensed under the Apache License, Version 2.0
# See the LICENSE file in the repository root for full license information.

# ===============================
# Descriptive Statistics Summary
# ===============================

# Load necessary libraries
library(tidyverse)   # for data manipulation
library(tseries)     # for Jarque-Bera test
library(e1071)       # for skewness and kurtosis

# Set working directory (adjust if needed)
setwd("/Users/faycal/Library/CloudStorage/Dropbox/R code for Network log-ARCH framework/")

# Load your log returns data
df <- read.csv("log_returns_data.csv")

# Convert Date column to Date format
df$Date <- as.Date(df$Date)

# Extract the return series (excluding Date column)
log_data <- df[, -1]

# Compute summary statistics for each country
summary_stats <- data.frame(
  Country    = names(log_data),
  Mean       = sapply(log_data, mean, na.rm = TRUE),
  Median     = sapply(log_data, median, na.rm = TRUE),
  Std_Dev    = sapply(log_data, sd, na.rm = TRUE),
  Minimum    = sapply(log_data, min, na.rm = TRUE),
  Maximum    = sapply(log_data, max, na.rm = TRUE),
  Skewness   = sapply(log_data, skewness, na.rm = TRUE),
  Kurtosis   = sapply(log_data, kurtosis, na.rm = TRUE),
  JB_pvalue  = sapply(log_data, function(x) {
    pval <- jarque.bera.test(x)$p.value
    if (pval < 2.2e-16) {
      return("< 2.2e-16")
    } else {
      return(sprintf("%.4f", pval))
    }
  })
)

# Round numeric columns for display
summary_stats[, 2:8] <- round(summary_stats[, 2:8], 8)

# Print the result as a plain table
print(summary_stats, row.names = FALSE)
















