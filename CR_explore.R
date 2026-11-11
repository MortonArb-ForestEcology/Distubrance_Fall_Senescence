library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(nlme)

PRIMARY  <- c("ndvi", "midgreendown")
REL_MIN  <- -5          # earliest year relative to the loss
REL_MAX  <- 10          # latest year relative to the loss
EARLY    <- 1:3         # "early" years after the loss
LATE     <- 4:10        # "late" years after the loss

outDir <- file.path("~/Google Drive/My Drive/Reidy_research", "scaleup_patches")
tsDir  <- file.path(outDir, "timeseries_full")
resDir <- file.path(outDir, "results_full")
# dir.create(resDir, recursive = TRUE, showWarnings = FALSE)

datRaw <- read.csv(file.path(tsDir, "patch_year_series.csv"))
datRaw <-datRaw[datRaw$in_analysis & datRaw$years_since_loss>=REL_MIN & datRaw$years_since_loss<=REL_MAX,]
datRaw$ndvi_diff <- datRaw$ndvi - datRaw$ctrl_ndvi
summary(datRaw)

pch_summary <- aggregate(cbind(ndvi, years_since_loss) ~ patch_id, data=datRaw, FUN=min)
names(pch_summary) <- c("patch_id", "ndvi.min", "years_since_loss.min")
pch_summary[,c("ndvi.max", "years_since_loss.max")] <- aggregate(cbind(ndvi, years_since_loss) ~ patch_id, data=datRaw, FUN=max)[,c("ndvi", "years_since_loss")]
summary(pch_summary)

# For now getting only the patches with full time series
pchUse <- pch_summary$patch_id[pch_summary$years_since_loss.min==-5 & pch_summary$years_since_loss.max==10]

datRaw <- datRaw[datRaw$patch_id %in% pchUse,]

# Creating my own summarized thing to graph
datAgg <- aggregate(cbind(ndvi, ctrl_ndvi, ndvi_diff, senescence_diff) ~ years_since_loss, data=datRaw, FUN=mean)
datAgg[,paste0(c("ndvi", "ctrl_ndvi", "ndvi_diff", "senescence_diff"), ".lwr")] <- aggregate(cbind(ndvi, ctrl_ndvi, ndvi_diff, senescence_diff) ~ years_since_loss, data=datRaw, FUN=quantile, 0.025)[,c("ndvi", "ctrl_ndvi", "ndvi_diff", "senescence_diff")]
datAgg[,paste0(c("ndvi", "ctrl_ndvi", "ndvi_diff", "senescence_diff"), ".upr")] <- aggregate(cbind(ndvi, ctrl_ndvi, ndvi_diff, senescence_diff) ~ years_since_loss, data=datRaw, FUN=quantile, 0.975)[,c("ndvi", "ctrl_ndvi", "ndvi_diff", "senescence_diff")]
summary(datAgg)

ggplot(data=datAgg) +
  geom_ribbon(aes(x=years_since_loss, ymin=ctrl_ndvi.lwr, ymax=ctrl_ndvi.upr, fill="reference"), alpha=0.5) +
  geom_line(aes(x=years_since_loss, y=ctrl_ndvi, color="reference")) +
  geom_point(aes(x=years_since_loss, y=ctrl_ndvi, color="reference")) +
  geom_ribbon(aes(x=years_since_loss, ymin=ndvi.lwr, ymax=ndvi.upr, fill="disturb"), alpha=0.5) +
  geom_line(aes(x=years_since_loss, y=ndvi, color="disturb")) +
  geom_point(aes(x=years_since_loss, y=ndvi, color="disturb")) +
  scale_color_manual(values=c("disturb"="darkorange", "reference"="darkgreen")) +
  scale_fill_manual(values=c("disturb"="darkorange", "reference"="darkgreen")) +
  theme_bw()

ggplot(data=datAgg) +
  geom_ribbon(aes(x=years_since_loss, ymin=ndvi_diff.lwr, ymax=ndvi_diff.upr, fill="disturb"), alpha=0.5) +
  geom_line(aes(x=years_since_loss, y=ndvi_diff, color="disturb")) +
  geom_point(aes(x=years_since_loss, y=ndvi_diff, color="disturb")) +
  scale_color_manual(values=c("disturb"="darkorange", "reference"="darkgreen")) +
  scale_fill_manual(values=c("disturb"="darkorange", "reference"="darkgreen")) +
  theme_bw()

# Is the difference significantly different from 0 for each year? 
# NOTE: Because we have years_since_loss here to identify different obs, we don't need a year random effect
ndviTime <- lme(ndvi_diff ~ as.factor(years_since_loss) - 1, random=c(us_l3code=~1, patch_id=~1), data=datRaw, na.action=na.omit)
summary(ndviTime)
# the "clean" subset shows NDVI has not recovered after 10 years



# Now doing the same thing with senescence

ggplot(data=datAgg) +
  geom_ribbon(aes(x=years_since_loss, ymin=senescence_diff.lwr, ymax=senescence_diff.upr, fill="disturb"), alpha=0.5) +
  geom_line(aes(x=years_since_loss, y=senescence_diff, color="disturb")) +
  geom_point(aes(x=years_since_loss, y=senescence_diff, color="disturb")) +
  scale_color_manual(values=c("disturb"="darkorange", "reference"="darkgreen")) +
  scale_fill_manual(values=c("disturb"="darkorange", "reference"="darkgreen")) +
  theme_bw()


senesceTime <- lme(senescence_diff ~ as.factor(years_since_loss) - 1, random=c(us_l3code=~1, patch_id=~1), data=datRaw, na.action=na.omit)
summary(senesceTime)
