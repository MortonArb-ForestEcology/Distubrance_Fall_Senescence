# =============================================================================
# Script 2: Polygon Processing and Filtering
# Purpose: Read forest loss polygons exported from Script 1, assign unique
#          labels, filter to polygons >= 500 m², and produce summary plots.
# Input:   forest_loss_polygons.csv (manually moved to path.dat after export)
# Output:  floss500_dat (filtered data frame, written out for Script 3)
# =============================================================================

library(ggplot2)
library(openxlsx)

# -----------------------------------------------------------------------------
# 1. File paths
# -----------------------------------------------------------------------------
path.google <- "~/Google Drive/My Drive"
path.dat    <- file.path(path.google, "Reidy_research/Hansen exploration/2026")

# -----------------------------------------------------------------------------
# 2. Read in CSV exported from Script 1
# -----------------------------------------------------------------------------
floss.dat <- read.csv(file.path(path.dat, "forest_loss_polygons.csv"))
head(floss.dat)

# -----------------------------------------------------------------------------
# 3. Create Year and unique Label columns
#    Label format: GDyyyy_n (e.g. GD2017_1)
#    'label' in Hansen = years since 2000, so add 2000 for full year
# -----------------------------------------------------------------------------
floss.dat$FullYear <- floss.dat$label + 2000

# Sort by area descending before assigning row numbers
# so Label suffix reflects size rank within the full dataset
floss.dat <- floss.dat[order(floss.dat$area_sqm, decreasing = TRUE), ]
floss.dat$RowNum <- 1:nrow(floss.dat)
floss.dat$Label  <- paste0("GD", floss.dat$FullYear, "_", floss.dat$RowNum)

# Select and rename columns
floss.dat <- floss.dat[, c("Label", "FullYear", ".geo", "latitude", "longitude",
                           "area_sqm", "system.index", "gain")]
names(floss.dat)[names(floss.dat) == "FullYear"] <- "Year"

# Set correct column types
floss.dat$Year     <- as.numeric(floss.dat$Year)
floss.dat$area_sqm <- as.numeric(floss.dat$area_sqm)

# -----------------------------------------------------------------------------
# 4. Filter to polygons >= 500 m²
#    Replaces the previous top-10 selection.
#    500 m² threshold retains all meaningful disturbance patches while
#    removing single-pixel noise artifacts from the Hansen vectorization.
# -----------------------------------------------------------------------------
floss500_dat <- floss.dat[floss.dat$area_sqm >= 500, ]

cat("Total polygons before filter:", nrow(floss.dat), "\n")
cat("Polygons >= 500 m²:", nrow(floss500_dat), "\n")
cat("Loss years represented:", paste(sort(unique(floss500_dat$Year)), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 5. Summary plots (all polygons, pre-filter, for QA)
# -----------------------------------------------------------------------------
# Area distribution histogram
ggplot(floss.dat, aes(x = area_sqm)) +
  geom_histogram(bins = 30, fill = "blue", color = "black") +
  labs(title = "Distribution of Area (sq m) — All Polygons",
       x = "Area (sq m)", y = "Frequency") +
  theme_minimal()

# Area by year boxplot
ggplot(floss.dat, aes(x = as.factor(Year), y = area_sqm)) +
  geom_boxplot() +
  labs(title = "Box Plot of Area by Year — All Polygons",
       x = "Year", y = "Area (sq m)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Area distribution histogram — filtered polygons only
ggplot(floss500_dat, aes(x = area_sqm)) +
  geom_histogram(bins = 30, fill = "darkgreen", color = "black") +
  labs(title = "Distribution of Area (sq m) — Polygons >= 500 m²",
       x = "Area (sq m)", y = "Frequency") +
  theme_minimal()

# -----------------------------------------------------------------------------
# 6. Write out filtered polygons for Script 3
# -----------------------------------------------------------------------------
write.csv(floss500_dat,
          file.path(path.dat, "floss500_dat.csv"),
          row.names = FALSE)


