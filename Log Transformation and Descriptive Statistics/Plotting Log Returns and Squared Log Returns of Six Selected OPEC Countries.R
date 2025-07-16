# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Set working directory
setwd("/Users/faycal/Library/CloudStorage/GoogleDrive-fayceldjebari@gmail.com/My Drive/Données/OPEC DATA/")

# Read the CSV file
df <- read.csv("log_returns_data.csv")

# Convert Date to proper date format
df$Date <- as.Date(df$Date)

# Rename "Saudi.Arabia" to "Saudi Arabia" for clarity in plots
colnames(df)[colnames(df) == "Saudi.Arabia"] <- "Saudi Arabia"

# Reshape the data for plotting (all columns except Date are countries)
returns_long <- df %>%
  pivot_longer(cols = -Date, names_to = "Country", values_to = "Return")

# Plot the returns
ggplot(returns_long, aes(x = Date, y = Return, color = Country)) +
  geom_line() +
  labs(x = "Date", y = "Log Return") +
  theme_minimal()

# Calculate squared returns
returns_long <- returns_long %>%
  mutate(Squared_Return = Return^2)

# Plot the squared returns
ggplot(returns_long, aes(x = Date, y = Squared_Return, color = Country)) +
  geom_line() +
  labs(x = "Date", y = "Squared Log Return") +
  theme_minimal()

