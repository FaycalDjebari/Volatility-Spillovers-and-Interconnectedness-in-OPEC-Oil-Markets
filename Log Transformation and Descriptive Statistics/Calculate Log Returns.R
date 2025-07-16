# ================================
# Calculate Log Returns (OPEC)
# ================================

# Load Required Libraries
library(tidyverse)
library(quantmod)
library(readxl)


# Set your working directory
setwd("/Users/faycal/Library/CloudStorage/Dropbox/R code for Network log-ARCH framework/")

# Load Excel Data 
data <- read_excel("Log Transformation and descriptive statistics/Monthly_OPEC_OIL_PRICE_DATA.xlsx")

# Convert Date column to Date format
data$Date <- as.Date(data$Date)

# List of country columns (do not modify names like "Saudi Arabia")
country_columns <- c("Algeria", "Iran", "Libya", "Nigeria", "Saudi Arabia", "UAE")

# Function to compute log returns
calculate_log_returns <- function(x) {
  Delt(x, type = "log")
}

# Apply log returns calculation to country columns (column names preserved)
returns_only <- data %>%
  select(all_of(country_columns)) %>%
  mutate(across(everything(), calculate_log_returns))

# Remove the first row (NA from differencing), keep original column names
returns_data <- cbind(Date = data$Date[-1], returns_only[-1, ])

# Save the output (column names unchanged)
write.csv(returns_data, "log_returns_data.csv", row.names = FALSE)

# View first few rows
head(returns_data)


