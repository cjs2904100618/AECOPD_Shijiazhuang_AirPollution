# =============================================================================
# SCRIPT: sensitivity_analysis.R
# PURPOSE: Sensitivity analysis for GAM-DLNM models
# CORRESPONDS TO: Manuscript sensitivity analysis section
# =============================================================================

# Load required packages
library(dlnm)    # Distributed Lag Nonlinear Models
library(mgcv)    # Generalized Additive Models  
library(dplyr)   # Data manipulation
library(ggplot2) # Visualization
library(splines) # Spline functions

cat("=== Starting Sensitivity Analysis ===\n")

# -----------------------------------------------------------------------------
# Parameter settings
# -----------------------------------------------------------------------------
pollutants <- c("PM2.5", "PM10", "CO", "NO2", "SO2", "O3")
df_rh <- 3
df_temp <- 3
df_lag <- 4
max_lag <- 10
delta <- 10
df_times <- seq(3, 11, 2)  # Time spline degrees of freedom: 3, 5, 7, 9, 11

# -----------------------------------------------------------------------------
# Create day of year variable if not exists
# -----------------------------------------------------------------------------
if (!"doy" %in% colnames(data)) {
  data$doy <- as.numeric(format(data$Date, "%j"))
  cat("Created doy variable\n")
}

# Calculate number of years in data
years <- length(unique(format(data$Date, "%Y")))
cat("Number of years in data:", years, "\n")

# -----------------------------------------------------------------------------
# Sensitivity analysis: testing different time spline df
# -----------------------------------------------------------------------------
sensitivity_list <- list()

for (p in pollutants) {
  cat("Analyzing pollutant:", p, "\n")
  
  # Use median as centering value
  cen <- median(data[[p]], na.rm = TRUE)
  
  for (df_time_per_year in df_times) {
    # Calculate total degrees of freedom
    df_time_total <- df_time_per_year * years
    
    cat("  Testing df_time =", df_time_per_year, "per year (total =", df_time_total, ")\n")
    
    # Create crossbasis
    cb <- crossbasis(data[[p]], 
                     lag = max_lag,
                     argvar = list(fun = "lin"),
                     arglag = list(fun = "ns", df = df_lag))
    
    # Build model (FIXED: use df_time_total instead of df_time_per_year)
    model <- gam(Dailyhospitalizations ~ cb +
                   ns(AT, df_temp) + 
                   ns(RH, df_rh) +
                   ns(doy, df_time_total),  # FIXED HERE
                 family = quasipoisson(), 
                 data = data)
    
    # Predict
    pred <- crosspred(cb, model, 
                      cen = cen, 
                      at = c(cen, cen + delta), 
                      by = 1)
    
    # Extract lag curve
    rr_vec <- pred$matRRfit[2, ]
    rr_lo  <- pred$matRRlow[2, ]
    rr_hi  <- pred$matRRhigh[2, ]
    
    # Find lag with maximum effect
    max_lag_idx <- which.max(rr_vec)
    
    # Save RR information at maximum lag
    sensitivity_list[[paste(p, df_time_per_year, sep = "_")]] <- tibble(
      pollutant = p,
      df_time_per_year = df_time_per_year,
      df_time_total = df_time_total,
      lag = max_lag_idx - 1,
      RR = rr_vec[max_lag_idx],
      RR_low = rr_lo[max_lag_idx],
      RR_high = rr_hi[max_lag_idx]
    )
  }
}

# Combine results
sensitivity_df <- bind_rows(sensitivity_list)

# -----------------------------------------------------------------------------
# Plotting (using your original color scheme)
# -----------------------------------------------------------------------------
cat("Generating sensitivity analysis plot...\n")

# Create sensitivity plot
sensitivity_plot <- ggplot(sensitivity_df, aes(x = df_time_per_year, y = RR)) +
  # Error bars (steelblue color as in your code)
  geom_errorbar(aes(ymin = RR_low, ymax = RR_high), 
                width = 0.3,
                color = "steelblue", 
                alpha = 0.8) +
  # Line (steelblue color)
  geom_line(color = "steelblue", size = 1) +
  # Points (steelblue color)
  geom_point(color = "steelblue", size = 2) +
  # Reference line y=1
  geom_hline(yintercept = 1, 
             linetype = "dashed", 
             color = "red") +
  # Facet by pollutant
  facet_wrap(~ pollutant, 
             scales = "free_y", 
             ncol = 2) +
  # Labels and title
  labs(title = "Sensitivity Analysis by Time Spline Degrees of Freedom",
       x = "Time Spline Degrees of Freedom (per year)",
       y = "RR and 95% CI") +
  # Theme settings
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 13, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.spacing = unit(1.5, "lines"),
    axis.title = element_text(size = 12)
  )

# Display plot
print(sensitivity_plot)

# Save plot
ggsave("sensitivity_by_df_per_pollutant.png", 
       plot = sensitivity_plot, 
       width = 14, 
       height = 10, 
       dpi = 300)
cat("Plot saved: sensitivity_by_df_per_pollutant.png\n")

# -----------------------------------------------------------------------------
# Export results to CSV
# -----------------------------------------------------------------------------
write.csv(sensitivity_df, "sensitivity_analysis_results.csv", row.names = FALSE)
cat("Results exported: sensitivity_analysis_results.csv\n")

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n=== Sensitivity Analysis Completed ===\n")
cat("Files created:\n")
cat("1. sensitivity_by_df_per_pollutant.png\n")
cat("2. sensitivity_analysis_results.csv\n")

cat("\nSummary of maximum RR changes:\n")
for (p in pollutants) {
  subset_df <- sensitivity_df %>% filter(pollutant == p)
  cat(p, ": RR range =", round(min(subset_df$RR), 3), "-", 
      round(max(subset_df$RR), 3), 
      "(range =", round(max(subset_df$RR) - min(subset_df$RR), 3), ")\n")
}