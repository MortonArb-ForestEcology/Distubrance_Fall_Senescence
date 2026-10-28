# =============================================================================
# Script 4: MidGreendown Phenology Analysis
# Purpose: Analyze how forest disturbance events affect autumn senescence
#          timing (MidGreendown DOY), comparing individual disturbed plots
#          against the multi-plot mean trend and examining deviation patterns
#          relative to disturbance year.
#
# Primary response variable: MidGreendown DOY
# Comparison baseline: mean MidGreendown DOY across all plots by year
#
# Hypothesis: Disturbance events will shift MidGreendown DOY relative to
#             pre-disturbance baseline, with the direction and magnitude of
#             the shift reflecting productivity changes induced by disturbance.
#
# Input:  modis_pheno_ndvi.csv (from Script 3)
# =============================================================================

library(ggplot2)
library(dplyr)
library(nlme)

# -----------------------------------------------------------------------------
# 1. File paths
# -----------------------------------------------------------------------------
path.google <- "/Users/brendonreidy/Library/CloudStorage/GoogleDrive-breidy@mortonarb.org/My Drive"
path.dat    <- file.path(path.google, "Reidy_research/Hansen exploration/2026")

# -----------------------------------------------------------------------------
# 2. Read and prepare data
# -----------------------------------------------------------------------------
flossdvi <- read.csv(file.path(path.dat, "modis_pheno_ndvi.csv"))

# Drop empty .geo column
flossdvi$.geo <- NULL

# Convert date columns
flossdvi$MidGreendown_Date <- as.Date(flossdvi$MidGreendown_Date)
flossdvi$Year              <- as.numeric(flossdvi$Year)

summary(flossdvi)
head(flossdvi)
dim(flossdvi)
length(unique(flossdvi$Label))

# -----------------------------------------------------------------------------
# 3. Extract disturbance year from Label
#    Label format: GDyyyy_n — uses regex to extract 4 consecutive digits
# -----------------------------------------------------------------------------
distyr <- function(label) {
  year_match <- regmatches(label, regexpr("\\d{4}", label))
  if (length(year_match) > 0) {
    return(as.numeric(year_match))
  } else {
    return(NA)
  }
}

flossdvi$dstrbyr <- sapply(flossdvi$Label, distyr)

# -----------------------------------------------------------------------------
# 4. NDVI disturbance verification
#    Purpose: confirm that NDVI drops at disturbance year for each plot,
#    verifying that Hansen loss polygons correspond to real canopy disturbance.
#    NDVI is NOT used as a productivity measure here.
# -----------------------------------------------------------------------------

# Calculate mean NDVI across all plots by year (comparison baseline)
mndviyr <- flossdvi %>%
  group_by(Year) %>%
  summarize(Mean.NDVI = mean(NDVI, na.rm = TRUE))

# Plot NDVI for individual plots vs mean trend — disturbance verification only
pltdistNDVI <- function(data, label, dstrbyr) {
  plot.dat <- data[data$Label == label, c("Label", "NDVI", "Year")]
  
  ggplot(plot.dat, aes(x = Year, y = NDVI)) +
    geom_line(color = "black", size = 1) +
    geom_vline(xintercept = dstrbyr, color = "red", linetype = "dashed", size = 1) +
    geom_line(data = mndviyr, aes(x = Year, y = Mean.NDVI),
              color = "blue", size = 1.2) +
    labs(title = paste("NDVI disturbance check — disturbance year:", dstrbyr),
         subtitle = "Verifying canopy disturbance signal in MODIS.",
         x = "Year", y = "NDVI") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    annotate("text", x = Inf, y = Inf,
             label = "Red: Disturbance year\nBlue: All-plot mean NDVI",
             hjust = 1.1, vjust = 1.1, size = 3)
}

# Run disturbance verification for a sample of plots
# Replace labels with actual labels from your dataset as needed
pltdistNDVI(flossdvi, unique(flossdvi$Label)[1], 
            unique(flossdvi$dstrbyr[flossdvi$Label == unique(flossdvi$Label)[1]]))

# -----------------------------------------------------------------------------
# 5. MidGreendown DOY time series per plot
#    Shows MidGreendown DOY by year for each disturbance plot with trend line.
#    Purpose: characterize baseline phenology trajectory and inter-annual
#    variability per site before centering on disturbance year.
# -----------------------------------------------------------------------------
unique.labels <- unique(flossdvi$Label)
nplots        <- length(unique.labels)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (i in 1:nplots) {
  plot.dat <- flossdvi[flossdvi$Label == unique.labels[i], ]
  
  plot(plot.dat$Year, plot.dat$MidGreendown_DOY,
       type = "n",
       ylim = range(plot.dat$MidGreendown_DOY, na.rm = TRUE),
       xlim = range(plot.dat$Year, na.rm = TRUE),
       xlab = "Year",
       ylab = "MidGreendown DOY",
       main = paste("MidGreendown DOY —", unique.labels[i]))
  
  lines(plot.dat$Year, plot.dat$MidGreendown_DOY, col = "goldenrod2", lwd = 2)
  points(plot.dat$Year, plot.dat$MidGreendown_DOY, col = "goldenrod2", pch = 16, cex = 1.2)
  
  # Add linear trend line if enough data
  if (sum(!is.na(plot.dat$MidGreendown_DOY)) > 2) {
    abline(lm(MidGreendown_DOY ~ Year, data = plot.dat, na.action = na.exclude),
           col = "goldenrod2", lty = 2, lwd = 1)
  }
  
  # Mark disturbance year
  abline(v = unique(plot.dat$dstrbyr), col = "red", lty = 2, lwd = 1.5)
  
  if (i == 1) {
    legend("topright",
           legend = c("MidGreendown DOY", "Trend", "Disturbance year"),
           col    = c("goldenrod2", "goldenrod2", "red"),
           lty    = c(1, 2, 2),
           lwd    = c(2, 1, 1.5),
           pch    = c(16, NA, NA),
           cex    = 0.8)
  }
}

par(mfrow = c(1, 1))

# -----------------------------------------------------------------------------
# 6. Pre-disturbance baseline and deviation calculation
#    Baseline: mean MidGreendown DOY across all plots in all years BEFORE
#    their respective disturbance year.
#    Deviation: how far each plot-year departs from that baseline.
# -----------------------------------------------------------------------------
predistub.dat <- flossdvi[flossdvi$Year < flossdvi$dstrbyr &
                            !is.na(flossdvi$MidGreendown_DOY), ]
pdistbmean <- mean(predistub.dat$MidGreendown_DOY, na.rm = TRUE)

cat("Pre-disturbance baseline MidGreendown DOY:", round(pdistbmean, 1), "\n")

flossdvi$yrfromdisturb        <- flossdvi$Year - flossdvi$dstrbyr
flossdvi$MidGreendown_Deviation <- flossdvi$MidGreendown_DOY - pdistbmean

# -----------------------------------------------------------------------------
# 7. MidGreendown deviation plot relative to disturbance year
#    Primary hypothesis visualization.
#    Shows whether senescence timing shifts after disturbance relative to
#    pre-disturbance baseline DOY.
#    Positive values = later than baseline | Negative = earlier than baseline
# -----------------------------------------------------------------------------
plot.dat <- flossdvi[!is.na(flossdvi$MidGreendown_Deviation) &
                       abs(flossdvi$yrfromdisturb) <= 5, ]

ggplot(plot.dat, aes(x = yrfromdisturb, y = MidGreendown_Deviation, color = Label)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, color = "red", linetype = "solid", size = 1.5, alpha = 0.7) +
  scale_x_continuous(breaks = seq(-5, 5, 1), minor_breaks = NULL) +
  labs(title    = "MidGreendown DOY Deviation Relative to Disturbance Year",
       subtitle = paste("Baseline: Pre-disturbance mean =", round(pdistbmean, 1), "DOY"),
       x        = "Years Relative to Disturbance (0 = Disturbance Year)",
       y        = "MidGreendown Deviation (Days)",
       color    = "Plot ID",
       caption  = "Positive = later than baseline | Negative = earlier than baseline") +
  theme_minimal() +
  theme(legend.position      = "right",
        panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5),
        panel.grid.minor.x   = element_blank(),
        axis.text.x          = element_text(size = 10),
        plot.title           = element_text(size = 14, face = "bold"),
        plot.subtitle        = element_text(size = 11)) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))


#coloring by disturbance year 
ggplot(plot.dat, aes(x = yrfromdisturb, y = MidGreendown_Deviation, 
                     color = as.factor(dstrbyr), group = Label)) +
  geom_line(size = 1, alpha = 0.6) +
  geom_point(size = 2, alpha = 0.8) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, color = "red", linetype = "solid", size = 1.5, alpha = 0.7) +
  scale_x_continuous(breaks = seq(-5, 5, 1), minor_breaks = NULL) +
  labs(title    = "MidGreendown DOY Deviation Relative to Disturbance Year",
       subtitle = paste("Baseline: Pre-disturbance mean =", round(pdistbmean, 1), "DOY"),
       x        = "Years Relative to Disturbance (0 = Disturbance Year)",
       y        = "MidGreendown Deviation (Days)",
       color    = "Disturbance Year",
       caption  = "Positive = later than baseline | Negative = earlier than baseline") +
  theme_minimal() +
  theme(legend.position      = "right",
        panel.grid.major.x   = element_line(color = "gray90", linewidth = 0.5),
        panel.grid.minor.x   = element_blank(),
        axis.text.x          = element_text(size = 10),
        plot.title           = element_text(size = 14, face = "bold"),
        plot.subtitle        = element_text(size = 11))
