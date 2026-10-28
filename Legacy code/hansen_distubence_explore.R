library(googlesheets4)
library(ggplot2)
library(terra)
library(MODISTools)
library(dplyr)
library(sf)
library(jsonlite)
library(readr)
library(tidyr)
library(readr)

# Read in the CSV exported
floss.dat <- read_csv("~/Google Drive/My Drive/Reidy_research/Hansen exploration/forest_loss_polygons.csv")
View(floss.dat)
head(floss.dat)

#  Create the Year column from the label, label is how hansen stores year of loss 
floss.dat$Year <- floss.dat$label


# Create a new   Label column in format GDyyyy_n this will be our unique polygon identifier
# Create a full year by adding 2000 to the year column
floss.dat$Year <- floss.dat$Year + 2000

# Add row number for the label suffix to account for polygons that may have the same year of loss 
floss.dat$RowNum <- 1:nrow(floss.dat)

# Create the Label: GD + year + _ + row number for a unique label
floss.dat$Label <- paste0("GD", floss.dat$Year, "_", floss.dat$RowNum)

# Select and reorder columns 
floss.dat <- floss.dat[, c("Label", "Year", ".geo", "latitude", "longitude", "area_sqm","system:index")]

summary(floss.dat)
head(floss.dat)
dim(floss.dat)

ggplot(data=floss.dat) +
  geom_point(aes(x=longitude, y=latitude))

head(data.frame(floss.dat))

#Setting year as a numeric
floss.dat$Year <- as.numeric(floss.dat$Year)


###### Examining polygons

# Convert year to factor for plotting
floss.dat$Year <- as.factor(floss.dat$Year)
floss.dat$area_sqm <- as.numeric(floss.dat$area_sqm)

# Area distribution histogram
ggplot(floss.dat, aes(x = area_sqm)) +
  geom_histogram(bins = 30, fill = "blue", color = "black") +
  labs(title = "Distribution of Area (sq m)", x = "Area (sq m)", y = "Frequency") +
  theme_minimal()

# Area by year boxplot
ggplot(floss.dat, aes(x = Year, y = area_sqm)) +
  geom_boxplot() +
  labs(title = "Box Plot of Area by Year", x = "Year", y = "Area (sqm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Get top 10 largest areas
floss10_dat <- floss.dat[order(floss.dat$area_sqm, decreasing = TRUE), ][1:10, ]

#write out top ten largest areas
write.xlsx(floss10_dat, "~/Google Drive/My Drive/Reidy_research/Hansen exploration/floss10_dat.xlsx", rowNames = FALSE)

###############
#analyzing NDVI
################
# Read NDVI data
flossdvi <- read.csv("~/Google Drive/My Drive/Reidy_research/Hansen exploration/modis_ndvi_mat_sen_gd50_data.csv")

# Fix label error
flossdvi$Label[flossdvi$Label == "GD202_8"] <- "GD2002_8"

# Convert date columns
flossdvi$Maturity_Date <- as.Date(flossdvi$Maturity_Date)
flossdvi$Senescence_Date <- as.Date(flossdvi$Senescence_Date)
flossdvi$MidGreendown_Date <- as.Date(flossdvi$MidGreendown_Date)
flossdvi$Year <- as.numeric(flossdvi$Year)

summary(flossdvi)
head(flossdvi)
dim(flossdvi)
length(unique(flossdvi$Label))


# NDVI by label(plot)
ggplot(flossdvi, aes(x = Year, y = NDVI, color = as.factor(Label))) +
  geom_line() +
  labs(title = "Mean NDVI by Year", x = "Year", y = "NDVI", color = "Plot ID") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# NDVI distribution across years
ggplot(flossdvi, aes(x = as.factor(Year), y = NDVI)) +
  geom_boxplot() +
  labs(title = "NDVI Distribution Across Years", x = "Year",y = "NDVI") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate mean NDVI by year
mndviyr <- flossdvi %>%
  group_by(Year) %>%
  summarize(Mean.NDVI = mean(NDVI, na.rm = TRUE))

summary(mndviyr)

# Individual disturbance plots with trend comparison
pltdistNDVI <- function(data, label, dstrbyr) {
  plot.dat <- data[data$Label == label, c("Label", "NDVI", "Year")]
  
  ggplot(plot.dat, aes(x = Year, y = NDVI)) +
    geom_line(color = "black", size = 1) +
    geom_vline(xintercept = dstrbyr, color = "red", linetype = "dashed", size = 1) +
    geom_line(data = mndviyr, aes(x = Year, y = Mean.NDVI), 
              color = "blue", size = 1.2) +
    labs(title = paste("NDVI series for", dstrbyr, "disturbance year"), x = "Year", y = "NDVI") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    annotate("text", x = Inf, y = Inf, 
             label = "Red: Disturbance\nBlue: Mean trend", 
             hjust = 1.1, vjust = 1.1, size = 3)
}

# Generating simpl plots for specific disturbance years
pltdistNDVI(flossdvi, "GD2017_11", 2017)
pltdistNDVI(flossdvi, "GD2015_4", 2015)
pltdistNDVI(flossdvi, "GD2002_8", 2002)

# Function to extract disturbance year from label(using regular expression of 4 consectutive digits to do so)
distyr <- function(label) {
  year_match <- regmatches(label, regexpr("\\d{4}", label))
  if(length(year_match) > 0) {
    return(as.numeric(year_match))
  } else {
    return(NA)
  }
}

# Adding disturbance year
flossdvi$dstrbyr <- sapply(flossdvi$Label, distyr)

# Calculate growing season metrics
flossdvi$growing <- flossdvi$Senescence_DOY - flossdvi$Maturity_DOY
flossdvi$early <- flossdvi$MidGreendown_DOY - flossdvi$Maturity_DOY
flossdvi$late <- flossdvi$Senescence_DOY - flossdvi$MidGreendown_DOY

# Dealing with negativity {in values not life} (cross-year boundaries)
flossdvi$growing[flossdvi$growing < 0] <- 
  flossdvi$growing[flossdvi$growing < 0] + 365
flossdvi$early[flossdvi$early < 0] <- 
  flossdvi$early[flossdvi$early < 0] + 365
flossdvi$late[flossdvi$late < 0] <- 
  flossdvi$late[flossdvi$late < 0] + 365

# Sens timing per plot - individual panels
unique.lables <- unique(flossdvi$Label)
nplots <- length(unique.lables)
pltcolor <- rainbow(nplots)

# Sens per timing plot
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (i in 1:nplots) {
  plot.dat <- flossdvi[flossdvi$Label == unique.lables[i], ]
  
  plot(plot.dat$Year, plot.dat$Maturity_DOY, 
       type = "n",
       ylim = range(c(plot.dat$Maturity_DOY, plot.dat$MidGreendown_DOY, plot.dat$Senescence_DOY), na.rm = TRUE),
       xlim = range(plot.dat$Year, na.rm = TRUE),
       xlab = "Year", 
       ylab = "Day of Year",
       main = paste("Phenology -", unique.lables[i]))
  
  # Add lines for each senesence metric
  lines(plot.dat$Year, plot.dat$Maturity_DOY, col = "green4", lwd = 2)
  points(plot.dat$Year, plot.dat$Maturity_DOY, col = "green4", pch = 16, cex = 1.2)
  
  lines(plot.dat$Year, plot.dat$MidGreendown_DOY, col = "goldenrod2", lwd = 2)
  points(plot.dat$Year, plot.dat$MidGreendown_DOY, col = "goldenrod2", pch = 16, cex = 1.2)
  
  lines(plot.dat$Year, plot.dat$Senescence_DOY, col = "red3", lwd = 2)
  points(plot.dat$Year, plot.dat$Senescence_DOY, col = "red3", pch = 16, cex = 1.2)
  
  # Add trend lines
  if (sum(!is.na(plot.dat$Maturity_DOY)) > 2) {
    abline(lm(Maturity_DOY ~ Year, data = plot.dat, na.action = na.exclude), 
           col = "green4", lty = 2, lwd = 1)
  }
  if (sum(!is.na(plot.dat$MidGreendown_DOY)) > 2) {
    abline(lm(MidGreendown_DOY ~ Year, data = plot.dat, na.action = na.exclude), 
           col = "goldenrod2", lty = 2, lwd = 1)
  }
  if (sum(!is.na(plot.dat$Senescence_DOY)) > 2) {
    abline(lm(Senescence_DOY ~ Year, data = plot.dat, na.action = na.exclude), 
           col = "red3", lty = 2, lwd = 1)
  }
  
  # Add legend to first plot
  if (i == 1) {
    legend("topright", 
           legend = c("Maturity", "MidGreendown", "Senescence"), 
           col = c("green4", "goldenrod2", "red3"),
           lwd = 2, pch = 16, cex = 0.8)
  }
}

par(mfrow = c(1, 1))

# Calculate pre-disturbance mean for MidGreendown analysis
predistub.dat <- flossdvi[flossdvi$Year < flossdvi$dstrbyr & !is.na(flossdvi$MidGreendown_DOY), ]
pdistbmean <- mean(predistub.dat$MidGreendown_DOY, na.rm = TRUE)

# Calculate deviations
flossdvi$yrfromdisturb <- flossdvi$Year - flossdvi$dstrbyr
flossdvi$MidGreendown_Deviation <- flossdvi$MidGreendown_DOY - pdistbmean

# MidGreendown deviation analysis relative to disturbance
plot.dat <- flossdvi[!is.na(flossdvi$MidGreendown_Deviation) & 
                        abs(flossdvi$yrfromdisturb) <= 10, ]

ggplot(plot.dat, aes(x = yrfromdisturb, y = MidGreendown_Deviation, color = Label)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, color = "red", linetype = "solid", size = 1.5, alpha = 0.7) +
  scale_x_continuous(breaks = seq(-10, 10, 1), 
                     minor_breaks = NULL,
                     limits = c(-3, 3)) +
  labs(title = "MidGreendown Deviation Relative to Disturbance Year",
       subtitle = paste("Baseline: Pre-disturbance mean =", round(pdistbmean, 1), "DOY"),
       x = "Years Relative to Disturbance (0 = Disturbance Year)", 
       y = "MidGreendown Deviation (Days)",
       color = "Plot ID",
       caption = "Positive values = Later than normal | Negative values = Earlier than normal") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.major.x = element_line(color = "gray90", size = 0.5),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11)) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

