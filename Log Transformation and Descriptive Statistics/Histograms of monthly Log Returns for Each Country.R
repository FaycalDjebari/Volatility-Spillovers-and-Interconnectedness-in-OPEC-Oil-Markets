# Copyright 2026 Fayçal Djebari
# Licensed under the Apache License, Version 2.0
# See the LICENSE file in the repository root for full license information.

# Load necessary libraries for visualization
library(ggplot2)
library(reshape2)
library(qqplotr)
library(scales)



#setwd("/Users/faycal/Library/CloudStorage/Dropbox/R code for Network log-ARCH framework/")
setwd("/Users/faycal/Library/CloudStorage/Dropbox/Code Paper 1/Volatility-Spillovers-and-Interconnectedness-in-OPEC-Oil-Markets-main")

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


# -----------   Q-Q plot of log returns for each country --------



# Data prep
df_melted$ScaledReturn <- ave(df_melted$LogReturn, df_melted$Country, FUN = scale)
df_melted$Country <- gsub("\\.", " ", df_melted$Country)

# Publication-ready Q-Q plot
p_qq_paper <- ggplot(df_melted, aes(sample = ScaledReturn)) +
  
  # Confidence bands - slightly refined
  stat_qq_band(distribution = "t", dparams = list(df = 5), 
               alpha = 0.25, fill = "grey75", color = NA) +
  
  # Reference line - keep your vermillion but slightly thicker
  stat_qq_line(distribution = "t", dparams = list(df = 5), 
               color = "#D55E00", linewidth = 0.9) +
  
  # Points - keep your blue but adjust for print clarity
  stat_qq_point(distribution = "t", dparams = list(df = 5), 
                color = "#0072B2", alpha = 0.55, size = 0.7) +
  
  # FIXED scales (critical for academic comparison)
  facet_wrap(~ Country, scales = "fixed") +
  
  # Theme refinement
  theme_light(base_size = 9) +
  
  # NO title/subtitle (goes in LaTeX caption)
  labs(x = "Theoretical Quantiles (Student-t, df = 5)",
       y = "Standardized Sample Quantiles") +
  
  # Polished typography
  theme(
    strip.background = element_rect(fill = "black", color = "grey20"),
    strip.text = element_text(face = "bold", color = "white", size = 8.5),
    axis.title = element_text(face = "bold", size = 9.5),
    axis.text = element_text(size = 7.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
    panel.spacing = unit(0.35, "lines"),
    panel.border = element_rect(color = "grey70", linewidth = 0.5)
  )

print(p_qq_paper)





# ----- Density plot of log returns for each country with a Normal distribution benchmark  ----

# ==============================================================================
# 1. SETUP: "EXECUTIVE" COLOR PALETTE
# ==============================================================================
# Deep, saturated colors that look good in print and screen.
# Navy, Burgundy, Forest Green, Deep Orange, Slate Grey, Royal Purple
executive_palette <- c(
  "#003049", # Prussian Blue (Deep Navy)
  "#d62828", # Fire Engine Red (Deep Red)
  "#f77f00", # Orange (Saturated)
  "#2a9d8f", # Persian Green (Teal/Forest)
  "#8e44ad", # Wisteria (Deep Purple)
  "#505050"  # Dark Grey (Neutral but strong)
)

# Ensure country names are clean
df_melted$Country <- gsub("\\.", " ", df_melted$Country)

# Calculate Normal Benchmark Params
mean_ret <- mean(df_melted$LogReturn, na.rm = TRUE)
sd_ret   <- sd(df_melted$LogReturn, na.rm = TRUE)
peak_val <- dnorm(mean_ret, mean_ret, sd_ret) # Max height for annotation

# ==============================================================================
# 2. GENERATE THE "HIGH JOURNAL" PLOT
# ==============================================================================

p_executive <- ggplot(df_melted, aes(x = LogReturn)) +
  
  # A. The Benchmark (Subtle Background)
  # A thin, dotted grey line for the Normal distribution. 
  # It is there for reference but doesn't fight for attention.
  stat_function(fun = dnorm, args = list(mean = mean_ret, sd = sd_ret),
                linetype = "dotted", color = "black", linewidth = 0.6, alpha = 0.7) +
  
  # B. The Data (Lines ONLY - No Fill)
  # This prevents the "muddy" look of overlapping colors.
  geom_line(stat="density", aes(color = Country, linetype = Country), 
            linewidth = 0.8) +
  
  # C. Colors & Linetypes
  scale_color_manual(values = executive_palette) +
  scale_linetype_manual(values = rep("solid", 6)) + # Keep all solid for clarity
  
  # D. Scales
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  
  # E. Annotation (Elegant Serif Font)
  annotate("text", x = 0.16, y = peak_val * 0.4, 
           label = "Gaussian\nBenchmark", 
           family = "serif", fontface = "italic", color = "gray20", size = 3.5, hjust = 0) +
  
  # F. THEME: The "LaTeX" Look
  # 'base_family = "serif"' makes all text look like Times New Roman
  theme_classic(base_size = 12, base_family = "serif") +
  
  labs(x = "Monthly Log Returns", 
       y = "Density") +
  
  coord_cartesian(xlim = c(-0.25, 0.25)) +
  
  theme(
    # Legend: Clean and Minimal
    legend.position = c(0.85, 0.75),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = NA), # No box around legend
    legend.key.width = unit(1.5, "cm"), # Longer lines in legend look more elegant
    
    # Axes: Crisp black lines
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    axis.text = element_text(color = "black")
  )

print(p_executive)






















