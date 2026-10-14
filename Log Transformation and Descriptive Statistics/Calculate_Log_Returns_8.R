# Copyright 2025 Fayçal Djebari
# Licensed under the Apache License, Version 2.0
# See the LICENSE file in the repository root for full license information.

# ==============================================
# Calculate Log Returns (OPEC, extended panel N = 8)
# Robustness panel for Referee 1, Comment 7: the six main countries plus Iraq
# and Kuwait. The Iraqi series has gaps up to December 1997 ("na" in the
# Excel file), so the common sample starts in January 1998 (first log-return:
# February 1998). The main file log_returns_data.csv is NOT modified.
# ==============================================

# Load Required Libraries
library(tidyverse)
library(quantmod)
library(readxl)

# Set your working directory
setwd("/Users/faycal/Library/CloudStorage/Dropbox/Code Paper 1/Volatility-Spillovers-and-Interconnectedness-in-OPEC-Oil-Markets-main")

# Load Excel Data ("na" marks the missing months of Iraq and Kuwait)
data <- read_excel("Log Transformation and descriptive statistics/Monthly_OPEC_OIL_PRICE_DATA_ROBUSTNESS.xlsx",
                   na = c("", "na", "NA"))

# Convert Date column to Date format
data$Date <- as.Date(data$Date)

# Same six countries in the same order as the main panel, then Iraq and Kuwait
country_columns <- c("Algeria", "Iran", "Libya", "Nigeria", "Saudi Arabia", "UAE",
                     "Iraq", "Kuwait")

# Common sample: from January 1998, the first month after the last Iraqi gap
data <- data %>%
  filter(Date >= as.Date("1998-01-01")) %>%
  select(Date, all_of(country_columns)) %>%
  mutate(across(all_of(country_columns), as.numeric))

stopifnot(!anyNA(data), all(data[country_columns] > 0))

# Function to compute log returns (same as the main panel)
calculate_log_returns <- function(x) {
  Delt(x, type = "log")
}

# Apply log returns calculation to country columns (column names preserved)
returns_only <- data %>%
  select(all_of(country_columns)) %>%
  mutate(across(everything(), ~ as.numeric(calculate_log_returns(.x))))

# Remove the first row (NA from differencing), keep original column names
returns_data <- cbind(Date = data$Date[-1], returns_only[-1, ])

# Checks: 323 monthly returns, February 1998 to December 2024, no missing values
stopifnot(!anyNA(returns_data))
cat("Rows:", nrow(returns_data), " from", format(min(returns_data$Date)),
    "to", format(max(returns_data$Date)), "\n")

# The six original columns must match the main file over the common months
main <- read.csv("log_returns_data.csv", check.names = FALSE)
main$Date <- as.Date(main$Date)
common <- merge(main, returns_data, by = "Date", suffixes = c("_main", "_ext"))
max_diff <- sapply(country_columns[1:6], function(c)
  max(abs(common[[paste0(c, "_main")]] - common[[paste0(c, "_ext")]])))
print(max_diff)   # should all be 0 (or ~1e-15)

# Save the output (column names unchanged)
write.csv(returns_data, "log_returns_data_8.csv", row.names = FALSE)

# View first few rows
head(returns_data)
