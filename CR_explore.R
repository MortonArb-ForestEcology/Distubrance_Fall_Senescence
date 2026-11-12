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
datRaw$senescence_diff <- datRaw$senescence - datRaw$ctrl_senescence
datRaw$midgreendown_diff <- datRaw$midgreendown - datRaw$ctrl_midgreendown
summary(datRaw)

pch_summary <- aggregate(cbind(ndvi, years_since_loss) ~ patch_id, data=datRaw, FUN=min)
names(pch_summary) <- c("patch_id", "ndvi.min", "years_since_loss.min")
pch_summary[,c("ndvi.max", "years_since_loss.max")] <- aggregate(cbind(ndvi, years_since_loss) ~ patch_id, data=datRaw, FUN=max)[,c("ndvi", "years_since_loss")]
summary(pch_summary)

# For now getting only the patches with full time series
pchUse <- pch_summary$patch_id[pch_summary$years_since_loss.min==-5 & pch_summary$years_since_loss.max==10]

datRaw <- datRaw[datRaw$patch_id %in% pchUse,]

# Creating my own summarized thing to graph
# Changed to show the middle 50%
datAgg <- aggregate(cbind(ndvi, ctrl_ndvi, ndvi_diff, midgreendown_diff) ~ years_since_loss, data=datRaw, FUN=mean)
datAgg[,paste0(c("ndvi", "ctrl_ndvi", "ndvi_diff", "midgreendown_diff"), ".lwr")] <- aggregate(cbind(ndvi, ctrl_ndvi, ndvi_diff, midgreendown_diff) ~ years_since_loss, data=datRaw, FUN=quantile, 0.25)[,c("ndvi", "ctrl_ndvi", "ndvi_diff", "midgreendown_diff")]
datAgg[,paste0(c("ndvi", "ctrl_ndvi", "ndvi_diff", "midgreendown_diff"), ".upr")] <- aggregate(cbind(ndvi, ctrl_ndvi, ndvi_diff, midgreendown_diff) ~ years_since_loss, data=datRaw, FUN=quantile, 0.75)[,c("ndvi", "ctrl_ndvi", "ndvi_diff", "midgreendown_diff")]
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
  geom_ribbon(aes(x=years_since_loss, ymin=midgreendown_diff.lwr, ymax=midgreendown_diff.upr, fill="disturb"), alpha=0.5) +
  geom_line(aes(x=years_since_loss, y=midgreendown_diff, color="disturb")) +
  geom_point(aes(x=years_since_loss, y=midgreendown_diff, color="disturb")) +
  scale_color_manual(values=c("disturb"="darkorange", "reference"="darkgreen")) +
  scale_fill_manual(values=c("disturb"="darkorange", "reference"="darkgreen")) +
  theme_bw()

set.seed(354)
pchSamp <- sample(unique(datRaw$patch_id), 100)

ggplot(data=datRaw[datRaw$patch_id %in% pchSamp,]) +
  # geom_ribbon(aes(x=years_since_loss, ymin=midgreendown_diff.lwr, ymax=midgreendown_diff.upr, fill="disturb"), alpha=0.5) +
  geom_line(aes(x=years_since_loss, y=midgreendown_diff, color="disturb", group=patch_id), size=0.2, alpha=0.5) +
  # geom_point(aes(x=years_since_loss, y=midgreendown_diff, color="disturb", group=patch_id)) +
  scale_color_manual(values=c("disturb"="darkorange", "reference"="darkgreen")) +
  scale_fill_manual(values=c("disturb"="darkorange", "reference"="darkgreen")) +
  theme_bw()


midgreendownTime <- lme(midgreendown_diff ~ as.factor(years_since_loss) - 1, random=c(us_l3code=~1, patch_id=~1), data=datRaw[,], na.action=na.omit)
summary(midgreendownTime)

summary(datRaw[datRaw$patch_id %in% pchSamp,])
head(datRaw)

# Are there differences by ecoregion?
midgreendownTime2 <- lme(midgreendown_diff ~ as.factor(us_l3code)*as.factor(years_since_loss) - 1, random=c(us_l3code=~1, patch_id=~1), data=datRaw[,], na.action=na.omit)
anova(midgreendownTime2)


#brendon attempts looking at differences by ecoregion 

#checking to see if there are enough unique patches per ecoregion to run estiamtes
table(datRaw$na_l2code[!duplicated(datRaw$patch_id)])

#eco region 8.2 and 8.5  the central plains and the south east coastals plain (and mississippi alluvial plain) have very a very few number
#unique patches  so filter those out
datReg <- subset(datRaw, na_l2code != 8.2 & na_l2code != 8.5)
datReg$region <- as.factor(datReg$na_l2code)
#creating a column for the years since loss only 
datReg$after <- as.numeric(datReg$years_since_loss >= 0)

#model looking to see if there is the same yearly pattern in every region
mgreenregsamefx <- lme(midgreendown_diff~ region + as.factor(years_since_loss), random = ~after | patch_id, data = datReg, method = "ML", na.action = na.omit, control = lmeControl(opt = "optim"))

#model where ever region gets its own yearly pattern
mgreenregdiffx <- lme(midgreendown_diff~ region * as.factor(years_since_loss), random = ~after | patch_id, data = datReg, method = "ML", na.action = na.omit, control = lmeControl(opt = "optim"))

# does the yearly pattern  of midgreenn down timing differ by region?
anova(mgreenregsamefx, mgreenregdiffx)

#does the yearly pattern  of mid green down timing differ by region when it comes to early and late timing after loss?
#setting early and late time, and removing year of disturbance 
datReg <- datReg[datReg$years_since_loss != 0, ]
datReg$early <- as.numeric(datReg$years_since_loss >= min(EARLY) & datReg$years_since_loss <= max(EARLY))
datReg$late  <- as.numeric(datReg$years_since_loss >= min(LATE)  & datReg$years_since_loss <= max(LATE))

# model where the early and late change is the same in every region
mgreenpersamefx <- lme(midgreendown_diff ~ region + early + late, random = list(patch_id = pdDiag(~ early + late)), data = datReg, method = "ML", na.action = na.omit, control = lmeControl(opt = "optim"))

# model where every region gets its own early and late change
mgreenperdiffx <- lme(midgreendown_diff ~ region * (early + late), random = list(patch_id = pdDiag(~ early + late)), data = datReg, method = "ML", na.action = na.omit, control = lmeControl(opt = "optim"))

# does the change after the loss differ by region?
anova(mgreenpersamefx, mgreenperdiffx)
