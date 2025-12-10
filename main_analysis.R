
---
  ```r
# =============================================================================
# TITLE: Main GAM-DLNM Analysis for Pollutant-AECOPD Associations
# PROJECT: Air Pollution and AECOPD Hospitalizations in Shijiazhuang (2017-2024)
# 
# PURPOSE: Estimate associations between six air pollutants and daily AECOPD
#          hospitalizations using Generalized Additive Models (GAM) with
#          Distributed Lag Nonlinear Models (DLNM).
# 
# MODEL: logE[Hosp_t] = β₀ + ∑CB(Pollutant_i,t; L=7, df_lag=4) + 
#                       ns(AT_t, df=3) + ns(RH_t, df=3) + 
#                       ns(doy_t, df=7×Years) + DOW
# 
# INPUT REQUIREMENTS:
#   Dataframe named `data` must contain:
#     - Date: Date in 'YYYY-MM-DD' format
#     - Dailyhospitalizations: Daily AECOPD admissions
#     - PM2.5, PM10, CO, NO2, SO2, O3: Pollutant concentrations
#     - AT, RH: Temperature (°C) and relative humidity (%)
# 
# OUTPUT:
#   - Excel file with lag-specific RR estimates
#   - PNG figure of lag-response associations
# 
# AUTHOR: [Your Name/Authorship Group]
# DATE: December 2024
# =============================================================================

# -----------------------------------------------------------------------------
# 1. LOAD REQUIRED PACKAGES
# -----------------------------------------------------------------------------
library(dlnm)      # Distributed Lag Nonlinear Models
library(mgcv)      # Generalized Additive Models
library(dplyr)     # Data manipulation
library(ggplot2)   # Visualization
library(splines)   # Spline functions
library(openxlsx)  # Excel export

cat("=== STARTING MAIN GAM-DLNM ANALYSIS ===\n")

# -----------------------------------------------------------------------------
# 2. DEFINE ANALYSIS PARAMETERS
# -----------------------------------------------------------------------------
pollutants <- c("PM2.5", "PM10", "CO", "NO2", "SO2", "O3")  # Six pollutants

# Model parameters (as per manuscript)
max_lag <- 7       # Maximum lag days (Lag 0-7)
df_lag <- 4        # Degrees of freedom for lag dimension
df_temp <- 3       # df for temperature spline
df_rh <- 3         # df for relative humidity spline
df_time_per_year <- 7  # df per year for temporal control

# Pollutant-specific increments (consistent with Table 1)
delta_values <- c("PM2.5" = 10,   # μg/m³
                  "PM10" = 10,    # μg/m³
                  "CO" = 1,       # mg/m³
                  "NO2" = 10,     # μg/m³
                  "SO2" = 10,     # μg/m³
                  "O3" = 10)      # μg/m³

# -----------------------------------------------------------------------------
# 3. DATA PREPARATION
# -----------------------------------------------------------------------------
# Verify data exists
if (!exists("data")) {
  stop("ERROR: Dataframe named 'data' not found. Please load your dataset.")
}

# Check required columns
required_cols <- c("Date", "Dailyhospitalizations", 
                   "PM2.5", "PM10", "CO", "NO2", "SO2", "O3", 
                   "AT", "RH")
missing_cols <- required_cols[!required_cols %in% colnames(data)]
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# Create day of year variable
if (!"doy" %in% colnames(data)) {
  data$doy <- as.numeric(format(data$Date, "%j"))
}

# Calculate total df for time spline
years <- length(unique(format(data$Date, "%Y")))
df_time_total <- df_time_per_year * years

cat("Data period:", as.character(min(data$Date)), "to", 
    as.character(max(data$Date)), "\n")
cat("Total years:", years, "\n")
cat("Time spline df:", df_time_total, "(", df_time_per_year, "per year)\n")

# -----------------------------------------------------------------------------
# 4. MAIN ANALYSIS LOOP
# -----------------------------------------------------------------------------
results_list <- list()      # For detailed results
summary_list <- list()      # For summary statistics

for (p in pollutants) {
  cat("\nAnalyzing pollutant:", p, "\n")
  
  delta <- delta_values[p]
  cat("  Increment (Δx):", delta, "\n")
  
  # 4.1 Create cross-basis matrix
  cb <- crossbasis(data[[p]], 
                   lag = max_lag,
                   argvar = list(fun = "lin"),      # Linear in exposure
                   arglag = list(fun = "ns", df = df_lag))  # Spline in lag
  
  # 4.2 Fit GAM with quasi-Poisson family
  model <- gam(Dailyhospitalizations ~ cb + 
                 ns(AT, df = df_temp) + 
                 ns(RH, df = df_rh) + 
                 ns(doy, df = df_time_total),
               family = quasipoisson(),
               data = data)
  
  # 4.3 Set centering at median exposure
  cen <- median(data[[p]], na.rm = TRUE)
  cat("  Centering value (median):", round(cen, 2), "\n")
  
  # 4.4 Predict RRs for single-day and cumulative lags
  pred <- crosspred(cb, model, 
                    cen = cen, 
                    at = c(cen, cen + delta), 
                    by = 1, 
                    cumul = TRUE)
  
  # 4.5 Extract single-day lag results (Lag0 to Lag7)
  single_day_df <- data.frame(
    Pollutant = p,
    Lag_Type = "Single-day",
    Lag_Day = 0:max_lag,
    RR = as.numeric(pred$matRRfit[2, ]),
    RR_low = as.numeric(pred$matRRlow[2, ]),
    RR_high = as.numeric(pred$matRRhigh[2, ]),
    Delta = delta,
    Centering_Value = cen,
    stringsAsFactors = FALSE
  )
  
  # 4.6 Extract cumulative lag results (Lag01 to Lag07)
  cumulative_df <- data.frame(
    Pollutant = p,
    Lag_Type = "Cumulative",
    Lag_Day = 0:max_lag,
    RR = as.numeric(pred$cumRRfit[2, ]),
    RR_low = as.numeric(pred$cumRRlow[2, ]),
    RR_high = as.numeric(pred$cumRRhigh[2, ]),
    Delta = delta,
    Centering_Value = cen,
    stringsAsFactors = FALSE
  )
  
  # 4.7 Combine and store results
  pollutant_results <- rbind(single_day_df, cumulative_df)
  results_list[[p]] <- pollutant_results
  
  # 4.8 Find maximum RR for summary
  max_single <- pollutant_results %>%
    filter(Lag_Type == "Single-day") %>%
    filter(RR == max(RR, na.rm = TRUE)) %>%
    slice(1)
  
  max_cumulative <- pollutant_results %>%
    filter(Lag_Type == "Cumulative") %>%
    filter(RR == max(RR, na.rm = TRUE)) %>%
    slice(1)
  
  summary_list[[p]] <- data.frame(
    Pollutant = p,
    Lag_Type = c("Single-day", "Cumulative"),
    Max_RR_Lag = c(max_single$Lag_Day, max_cumulative$Lag_Day),
    Max_RR = c(max_single$RR, max_cumulative$RR),
    RR_low = c(max_single$RR_low, max_cumulative$RR_low),
    RR_high = c(max_single$RR_high, max_cumulative$RR_high),
    Delta = delta,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# 5. COMBINE AND EXPORT RESULTS
# -----------------------------------------------------------------------------
all_results_df <- do.call(rbind, results_list)
all_summary_df <- do.call(rbind, summary_list)

# 5.1 Export to Excel
cat("\nExporting results to Excel...\n")
wb <- createWorkbook()

addWorksheet(wb, "Detailed_Lag_Effects")
writeData(wb, "Detailed_Lag_Effects", all_results_df)

addWorksheet(wb, "Maximum_RR_Summary")
writeData(wb, "Maximum_RR_Summary", all_summary_df)

excel_filename <- "pollutant_lag_effect_results.xlsx"
saveWorkbook(wb, excel_filename, overwrite = TRUE)
cat("Excel file saved as:", excel_filename, "\n")

# 5.2 Generate lag-response plot (Figure 1)
cat("Generating lag-response plot...\n")

# Prepare plot data
K_plot <- 8
x_non <- 1:K_plot
x_cum <- (max(x_non)+2):(max(x_non)+1+K_plot)

plot_data <- all_results_df %>%
  filter(Lag_Day < K_plot) %>%
  mutate(
    x_coord = case_when(
      Lag_Type == "Single-day" ~ x_non[Lag_Day + 1],
      Lag_Type == "Cumulative" ~ x_cum[Lag_Day + 1]
    )
  )

# Create plot
lag_plot <- ggplot(plot_data, aes(x = x_coord, y = RR, group = Lag_Type)) +
  geom_errorbar(aes(ymin = RR_low, ymax = RR_high, color = Lag_Type, linetype = Lag_Type),
                width = 0.12, position = position_dodge(width = 0.3), size = 0.6) +
  geom_line(aes(color = Lag_Type, linetype = Lag_Type), 
            size = 0.9, position = position_dodge(width = 0.3)) +
  geom_point(aes(color = Lag_Type), 
             size = 2, shape = 16, position = position_dodge(width = 0.3)) +
  facet_wrap(~ Pollutant, ncol = 2, scales = "free_y") +
  geom_vline(xintercept = max(x_non) + 0.5, linetype = "dotted", color = "grey50") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_x_continuous(
    breaks = c(x_non, x_cum),
    labels = c(paste0("Lag", 0:(K_plot-1)), sprintf("Lag%02d", 0:(K_plot-1))),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_color_manual(values = c("Single-day" = "#1f77b4", "Cumulative" = "#b22222")) +
  scale_linetype_manual(values = c("Single-day" = "solid", "Cumulative" = "dashed")) +
  labs(x = "Lag", y = "RR (95% CI)") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

# Save plot
ggsave("main_lag_response_plot.png", 
       plot = lag_plot, 
       width = 14, height = 10, dpi = 300)
cat("Figure saved as: main_lag_response_plot.png\n")

# -----------------------------------------------------------------------------
# 6. FINAL SUMMARY
# -----------------------------------------------------------------------------
cat("\n=== MAIN ANALYSIS COMPLETED ===\n")
cat("Files created:\n")
cat("1.", excel_filename, "- Detailed results\n")
cat("2. main_lag_response_plot.png - Lag-response associations\n\n")

cat("Summary of maximum RRs:\n")
print(all_summary_df)

# End of script
# =============================================================================
