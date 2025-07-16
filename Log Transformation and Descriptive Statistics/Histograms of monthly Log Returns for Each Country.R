# Copyright 2025 Fayçal Djebari
# Licensed under the Apache License, Version 2.0
# See the LICENSE file in the repository root for full license information.

# Load necessary libraries for visualization
library(ggplot2)
library(reshape2)


setwd("/Users/faycal/Library/CloudStorage/Dropbox/R code for Network log-ARCH framework/")

# Read the CSV file
df <- read.csv("log_returns_data.csv")

# Convert Date to proper date format
df$Date <- as.Date(df$Date)

# Melt the dataframe for easier plotting
df_melted <- melt(df, id.vars = 'Date', variable.name = 'Country', value.name = 'LogReturn')

# Remove the "LogReturn_" prefix from the country names
df_melted$Country <- sub("LogReturn_", "", df_melted$Country)

# Plot the distribution of log returns for each country
p <- ggplot(df_melted, aes(x = LogReturn, fill = Country)) +
  geom_histogram(binwidth = 0.01, alpha = 0.6, position = 'identity') +
  facet_wrap(~ Country, scales = 'free') +
  theme_minimal() +
  labs( x = 'Log Return', y = 'Frequency') +
  coord_cartesian(xlim = c(-0.5, 0.5))  # Set x-axis limits

# Print the plot
print(p)
