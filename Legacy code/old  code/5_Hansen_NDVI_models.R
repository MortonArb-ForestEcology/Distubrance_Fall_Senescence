# =============================================================================
# Script 5: LME Models — MidGreendown Disturbance Effect
# Purpose: Test whether forest disturbance events significantly affect
#          MidGreendown DOY using linear mixed effects models.
#
# Three model variants test different random effect structures:
#   lmemod:    controls for both weather (Year) and site (Label) — primary model
#   lmemodwe:  controls for weather only — isolates site contribution
#   lmemodsit: controls for site only — isolates weather contribution
#
# Grouped model (timeGroupLME):
#   Collapses year-by-year factor into before/during/after periods.
#   Trades temporal resolution for statistical power — a more powerful
#   test of the overall disturbance effect with small sample sizes.
#
# NDVI verification models (NDVItest, NDVItest2):
#   Confirm that NDVI drops at disturbance year (yrfromdisturb == 0).
#   This is a data quality gate — verifying Hansen polygons correspond
#   to real canopy disturbance before interpreting phenology results.
#   These models appear BEFORE phenology models intentionally.

# =============================================================================

library(nlme)
library(MuMIn)
library(ggplot2)
library(readr)
# -----------------------------------------------------------------------------
# 1. File paths and data loading
# -----------------------------------------------------------------------------
path.google <- "~/Google Drive/My Drive"
path.dat    <- file.path(path.google, "Reidy_research/Hansen exploration/2026/")

flossdvi <- read.csv(file.path(path.dat, "flossdvi.csv"))
flossdvi$MidGreendown_Date <- as.Date(flossdvi$MidGreendown_Date)
flossdvi$Year              <- as.numeric(flossdvi$Year)

# -----------------------------------------------------------------------------
# 2. Filter to ±5 years around disturbance
# -----------------------------------------------------------------------------
dat.lme <- flossdvi[abs(flossdvi$yrfromdisturb) <= 5 &
                      !is.na(flossdvi$MidGreendown_DOY) &
                      !is.na(flossdvi$yrfromdisturb), ]

# Quick data checks
str(dat.lme[, c("Year", "Label", "MidGreendown_DOY", "dstrbyr", "yrfromdisturb")])
cat("Plots in model dataset:", length(unique(dat.lme$Label)), "\n")
cat("MidGreendown NAs:", sum(is.na(dat.lme$MidGreendown_DOY)), "\n")
cat("yrfromdisturb range:", range(dat.lme$yrfromdisturb, na.rm = TRUE), "\n")

# -----------------------------------------------------------------------------
# 3. NDVI disturbance verification — run BEFORE phenology models
#    Confirms Hansen loss polygons show detectable NDVI drop at disturbance year.
#    This is a quality gate, not a productivity analysis.
# -----------------------------------------------------------------------------
#removing NA values for NDVI 
dat.lme.ndvi <- dat.lme[!is.na(dat.lme$NDVI), ]
# Visual check — NDVI trajectory per plot centered on disturbance year
ggplot(data = dat.lme.ndvi) +
  geom_line(aes(x = yrfromdisturb, y = NDVI, color = Label), 
            alpha = 0.6, show.legend = FALSE) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(title   = "NDVI by Years Relative to Disturbance — Verification Only",
       subtitle = "Checking for NDVI drop at disturbance year (0)",
       x = "Years Relative to Disturbance",
       y = "NDVI") +
  theme_bw()

# LME: does NDVI differ significantly at disturbance year?
# Random effects: Year (weather) + Label (site)
NDVItest <- lme(NDVI ~ relevel(as.factor(yrfromdisturb), "0"),
                random = list(Year = ~1, Label = ~1),
                data   = dat.lme.ndvi)
summary(NDVItest)
anova(NDVItest)

# LME: site random effect only
NDVItest2 <- lme(NDVI ~ as.factor(yrfromdisturb),
                 random = list(Label = ~1),
                 data   = dat.lme.ndvi)
summary(NDVItest2)
anova(NDVItest2)

# -----------------------------------------------------------------------------
# 4. Primary LME models — MidGreendown DOY ~ disturbance year (factor)
#
#    lmemod:    Year + Label random effects (primary model)
#               Controls for both inter-annual weather variation and
#               site-level differences simultaneously.
#    lmemodwe:  Year random effect only
#               Tests whether disturbance effect holds when controlling
#               only for weather — isolates site contribution.
#    lmemodsit: Label random effect only
#               Tests whether disturbance effect holds when controlling
#               only for site — isolates weather contribution.
# -----------------------------------------------------------------------------

# Primary model: weather + site random effects
lmemod <- lme(MidGreendown_DOY ~ as.factor(yrfromdisturb),
              random = list(Year = ~1, Label = ~1),
              data   = dat.lme)
summary(lmemod)
r.squaredGLMM(lmemod)
anova(lmemod)

# Weather only
lmemodwe <- lme(MidGreendown_DOY ~ as.factor(yrfromdisturb),
                random = ~1 | Year,
                data   = dat.lme)
summary(lmemodwe)
r.squaredGLMM(lmemodwe)
anova(lmemodwe)

# Site only
lmemodsit <- lme(MidGreendown_DOY ~ as.factor(yrfromdisturb),
                 random = ~1 | Label,
                 data   = dat.lme)
summary(lmemodsit)
r.squaredGLMM(lmemodsit)
anova(lmemodsit)

# -----------------------------------------------------------------------------
# 5. Grouped disturbance period model — before / during / after
#    Collapses year-by-year factor into 3 periods for a more statistically
#    powerful test. Trades temporal resolution for sensitivity.
#    With n = 141 plots and ±5 years, this is a more appropriate test
#    than the 10-level factor model above.
# -----------------------------------------------------------------------------
dat.lme$yrGroup <- NA
dat.lme$yrGroup[dat.lme$yrfromdisturb < 0]  <- "before"
dat.lme$yrGroup[dat.lme$yrfromdisturb == 0] <- "during"
dat.lme$yrGroup[dat.lme$yrfromdisturb > 0]  <- "after"
dat.lme$yrGroup <- factor(dat.lme$yrGroup, levels = c("before", "during", "after"))

summary(dat.lme$yrGroup)

# MidGreendown grouped model
timeGroupLME <- lme(MidGreendown_DOY ~ yrGroup,
                    random = list(Year = ~1, Label = ~1),
                    data   = dat.lme)
summary(timeGroupLME)
r.squaredGLMM(timeGroupLME)
anova(timeGroupLME)

# NDVI grouped model — verification only
dat.lme.ndvi$yrGroup <- NA
dat.lme.ndvi$yrGroup[dat.lme.ndvi$yrfromdisturb < 0]  <- "before"
dat.lme.ndvi$yrGroup[dat.lme.ndvi$yrfromdisturb == 0] <- "during"
dat.lme.ndvi$yrGroup[dat.lme.ndvi$yrfromdisturb > 0]  <- "after"
dat.lme.ndvi$yrGroup <- factor(dat.lme.ndvi$yrGroup, levels = c("before", "during", "after"))

timeGroupLME2 <- lme(NDVI ~ yrGroup,
                     random = list(Year = ~1, Label = ~1),
                     data   = dat.lme.ndvi)
summary(timeGroupLME2)
r.squaredGLMM(timeGroupLME2)
anova(timeGroupLME2)
