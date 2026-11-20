# Figures for the poster
# -----------------------------------------------------------------------------
library(ggplot2)
library(sf)

outDir <- file.path("~/Google Drive/My Drive/Reidy_research", "scaleup_patches")
tsDir  <- file.path(outDir, "timeseries_full")
resDir <- file.path(outDir, "results_full")
figDir <- file.path(outDir, "poster_figures")
dir.create(figDir, showWarnings=FALSE)

REL_MIN <- -5      # earliest year relative to the loss
REL_MAX <- 10      # latest year relative to the loss
LATE    <- 4:10    # "late" years after the loss

# Poster colors
COL_TEXT   <- "white"      
COL_GRID   <- "#3f8a52"   
COL_GREEN  <- "#9ed39a"    
COL_ORANGE <- "#f08a4b"    
COL_LAND  <- "#e5f2df"

# Setting a theme so these will all show up on the poster in the same way and large enough size
themePoster <- theme_bw() +
  theme(plot.background=element_blank(), panel.background=element_blank(), panel.border=element_blank(),
        panel.grid.minor=element_blank(), panel.grid.major=element_line(color=COL_GRID, linewidth=0.4),
        axis.line=element_line(color=COL_TEXT), axis.ticks=element_line(color=COL_TEXT),
        text=element_text(color=COL_TEXT, size=18), axis.text=element_text(color=COL_TEXT, size=16),
        plot.title=element_text(face="bold", size=22), plot.subtitle=element_text(size=15),
        plot.caption=element_text(size=13, hjust=0),
        legend.background=element_blank(), legend.key=element_blank())

# resulrts to draw upong
evt     <- read.csv(file.path(resDir, "event_study_coefficients.csv"))
mainRes <- read.csv(file.path(resDir, "main_results.csv"))
summary(evt)
mainRes[,c("outcome", "period", "estimate", "ci_lo", "ci_hi")]


# "How we read fall from space"
datCurve <- data.frame(doy=1:365)
datCurve$uncut  <- 1/(1+exp(-(datCurve$doy-125)/9)) - 1/(1+exp(-(datCurve$doy-285)/11))
datCurve$regrow <- 1/(1+exp(-(datCurve$doy-125)/9)) - 1/(1+exp(-(datCurve$doy-277)/5.1))

datPhase <- data.frame(phase=c("Greenup", "Senescence", "Mid-greendown", "Dormancy"),
                       doy=c(109, 261, 285, 304),
                       hjust=c(-0.2, -0.1, 1.05, -0.15),  
                       vjust=c(0.5, -0.8, 1.8, -0.8))
datPhase$evi <- datCurve$uncut[datPhase$doy]
datPhase

ggplot(data=datCurve) +
  geom_ribbon(data=datCurve[datCurve$doy>=109 & datCurve$doy<=261,], aes(x=doy, ymin=0, ymax=uncut), fill=COL_GREEN, alpha=0.2) +
  geom_ribbon(data=datCurve[datCurve$doy>=261 & datCurve$doy<=285,], aes(x=doy, ymin=0, ymax=uncut), fill=COL_ORANGE, alpha=0.35) +
  geom_hline(yintercept=c(0.9, 0.5), linetype="dashed", color="gray70") +
  annotate("text", x=365, y=c(0.94, 0.54), label=c("90%", "50%"), hjust=1, color=COL_TEXT, size=5) +
  geom_line(aes(x=doy, y=uncut, color="uncut forest"), linewidth=1.5) +
  geom_line(data=datCurve[datCurve$doy>=236,], aes(x=doy, y=regrow, color="regrowing stand"), linewidth=1.2, linetype="dashed") +
  geom_point(data=datPhase, aes(x=doy, y=evi), color=COL_GREEN, size=4) +
  geom_text(data=datPhase, aes(x=doy, y=evi, label=phase, hjust=hjust, vjust=vjust), color=COL_TEXT, size=6, fontface="bold") +
  annotate("text", x=273, y=0.07, label="fall\ndecline", color=COL_ORANGE, size=5, fontface="bold") +
  scale_color_manual(values=c("uncut forest"=COL_GREEN, "regrowing stand"=COL_ORANGE)) +
  scale_x_continuous(breaks=c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349), labels=month.abb) +
  themePoster +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid=element_blank(),
        legend.position="top", legend.justification="left") +
  labs(title="How we read fall from space", subtitle="Schematic, MODIS MCD12Q2 (EVI2)",
       x=NULL, y="Greenness", color=NULL)

ggsave(file.path(figDir, "fig0_how_we_read_fall.png"), width=14, height=5.5, dpi=300, bg="transparent")


# MidGreendown by year since loss
evtMGD <- evt[evt$outcome=="midgreendown",]
evtMGD$side <- "before"
evtMGD$side[evtMGD$rel_year>=0] <- "after"   # so the line breaks at the loss
evtMGD[,c("rel_year", "estimate", "ci_lo", "ci_hi")]

# The late estimate and its CI for the label years 4-10 vs the 5 years before
lateMGD <- mainRes[mainRes$outcome=="midgreendown" & mainRes$period=="late",]
lateMGD
labMGD <- paste0("Years 4-10 vs. the 5 years before loss:\n",
                 round(-lateMGD$estimate, 1), " days earlier (95% CI ", round(-lateMGD$ci_hi, 1), " to ", round(-lateMGD$ci_lo, 1), ")")
labMGD

ggplot(data=evtMGD) +
  annotate("rect", xmin=-5.5, xmax=-0.5, ymin=-Inf, ymax=Inf, fill="white", alpha=0.07) +
  annotate("text", x=-3, y=1.4, label="before loss", color=COL_TEXT, size=5) +
  geom_hline(yintercept=0, color=COL_TEXT) +
  geom_vline(xintercept=0, linetype="dashed", color=COL_TEXT) +
  geom_errorbar(aes(x=rel_year, ymin=ci_lo, ymax=ci_hi), width=0, color=COL_ORANGE, linewidth=1) +
  geom_line(aes(x=rel_year, y=estimate, group=side), color=COL_ORANGE, linewidth=1) +
  geom_point(aes(x=rel_year, y=estimate), color=COL_ORANGE, size=4) +
  annotate("text", x=10, y=-0.7, label=labMGD, hjust=1, color=COL_TEXT, size=6, fontface="bold") +
  scale_x_continuous(breaks=REL_MIN:REL_MAX) +
  themePoster +
  labs(title="Mid-greendown comes earlier after canopy loss, and stays earlier",
       subtitle="Mid-greendown date, loss patch minus nearby uncut forest. MODIS MCD12Q2 (500 m, daily)",
       x="Years since forest loss", y="Change in mid-greendown (days)\nnegative = earlier, vs. year before loss",)
ggsave(file.path(figDir, "fig2_midgreendown.png"), width=14, height=7, dpi=300, bg="transparent")

#July NDVI 
evtNDVI <- evt[evt$outcome=="ndvi",]
evtNDVI$side <- "before"
evtNDVI$side[evtNDVI$rel_year>=0] <- "after"
evtNDVI[,c("rel_year", "estimate", "ci_lo", "ci_hi")]

# Getting early and  early and late estimates for the label
earlyNDVI <- mainRes[mainRes$outcome=="ndvi" & mainRes$period=="early",]
lateNDVI  <- mainRes[mainRes$outcome=="ndvi" & mainRes$period=="late",]
labNDVI <- paste0("Model estimates vs. the 5 years before loss:\n",
                  "years 1-3: ", round(earlyNDVI$estimate, 3), " (95% CI ", round(earlyNDVI$ci_lo, 3), " to ", round(earlyNDVI$ci_hi, 3), ")\n",
                  "years 4-10: ", round(lateNDVI$estimate, 3), " (95% CI ", round(lateNDVI$ci_lo, 3), " to ", round(lateNDVI$ci_hi, 3), ")")
labNDVI

ggplot(data=evtNDVI) +
  annotate("rect", xmin=-5.5, xmax=-0.5, ymin=-Inf, ymax=Inf, fill="white", alpha=0.07) +
  annotate("text", x=-3, y=-0.11, label="before loss", color=COL_TEXT, size=5) +
  geom_hline(yintercept=0, color=COL_TEXT) +
  geom_vline(xintercept=0, linetype="dashed", color=COL_TEXT) +
  geom_line(aes(x=rel_year, y=estimate, group=side), color=COL_GREEN, linewidth=1) +
  geom_point(aes(x=rel_year, y=estimate), color=COL_GREEN, size=4) +
  annotate("text", x=10, y=-0.09, label=labNDVI, hjust=1, color=COL_TEXT, size=5) +
  scale_x_continuous(breaks=REL_MIN:REL_MAX) +
  themePoster +
  labs(title="Summer greenness drops after canopy loss, then recovers",
       subtitle="July NDVI, loss patch minus nearby uncut forest. MODIS MOD13A1 v6.1 (500 m, 16-day composites)",
       x="Years since forest loss", y="Change in July NDVI\nvs. year before loss",)

ggsave(file.path(figDir, "fig1_ndvi.png"), width=14, height=6.5, dpi=300, bg="transparent")


#  MidGreendown for  ecoregions
regRes <- read.csv(file.path(resDir, "checks", "region_effects.csv"))

# Late years only droppin 8.2, 8.5, 9.2
regMGD <- regRes[regRes$outcome=="midgreendown" & regRes$period=="late" & regRes$n>=100,]
regMGD$regLab <- paste0(regMGD$region, "\n", prettyNum(regMGD$n, big.mark=","), " patches")
regMGD$regLab <- factor(regMGD$regLab, levels=regMGD$regLab[order(regMGD$estimate, decreasing=TRUE)])   # biggest shift on top
regMGD$valLab <- paste0(round(regMGD$estimate, 2), " [", round(regMGD$ci_lo, 2), ", ", round(regMGD$ci_hi, 2), "]")
regMGD[,c("region", "n", "estimate", "ci_lo", "ci_hi")]

ggplot(data=regMGD) +
  geom_vline(xintercept=0, color=COL_TEXT) +
  geom_segment(aes(x=ci_lo, xend=ci_hi, y=regLab, yend=regLab), color=COL_ORANGE, linewidth=1.5) +
  geom_point(aes(x=estimate, y=regLab), color=COL_ORANGE, size=5) +
  geom_text(aes(x=ci_hi, y=regLab, label=valLab), hjust=-0.15, color=COL_TEXT, size=5) +
  scale_x_continuous(limits=c(-5.5, 2), breaks=-5:0) +
  themePoster +
  theme(panel.grid.major.y=element_blank()) +
  labs(title="Earlier in all five main ecoregions",
       subtitle="Mid-greendown shift in years 4-10 after loss, by EPA Level II ecoregion",
       x="Change in mid-greendown, years 4-10 (days; negative = earlier)", y=NULL,)

ggsave(file.path(figDir, "fig3_regions.png"), width=14, height=7.5, dpi=300, bg="transparent")



# 5. Study area map -
eco <- st_read(file.path(outDir, "boundary_ecoregions.gpkg"), quiet=TRUE)
eco <- st_transform(eco, 5070)
eco$l2 <- as.character(eco$na_l2code)
table(eco$l2)

# The patches in the analysis
pchInfo <- read.csv(file.path(tsDir, "patch_summary.csv"))
pchInfo <- pchInfo[pchInfo$in_analysis,]
nrow(pchInfo)   # should be 9,097
pchSF <- st_transform(st_as_sf(pchInfo, coords=c("lon", "lat"), crs=4326), 5070)

ECO_COLORS <- c("5.2"="#b3eced", "5.3"="#fdd6d3", "8.1"="#b3ecd2", "8.2"="#f0e0b3",
                "8.3"="#eed8ff", "8.4"="#b3e5ff", "8.5"="#d8e7b3", "9.2"="#ffd0f0")
ECO_NAMES  <- c("5.2"="Mixed Wood Shield", "5.3"="Atlantic Highlands", "8.1"="Mixed Wood Plains",
                "8.2"="Central USA Plains", "8.3"="Southeastern USA Plains", "8.4"="Ozark-Appalachian Forests",
                "8.5"="Miss. Alluvial & SE Coastal Plains", "9.2"="Temperate Prairies")

ggplot() +
  geom_sf(data=eco, aes(fill=l2), color="gray50", linewidth=0.1) +   # Level III outlines, colored by Level II
  geom_sf(data=pchSF, color="darkgreen", size=0.3, alpha=0.6) +      # one dot per patch
  scale_fill_manual(values=ECO_COLORS, breaks=names(ECO_NAMES), labels=ECO_NAMES, na.value="gray80") +
  theme_void() +
  theme(text=element_text(color=COL_TEXT, size=16)) +
  labs(fill="Level II ecoregion")

ggsave(file.path(figDir, "map_study_area.png"), width=11, height=9, dpi=300, bg="transparent")


# 6. Map of each patchs change in MidGreendown for every patch mean difference = patch - reference in the late years minus
# the mean in the years before the loss. This is just  raw change not  model.
# Single patches are wayyy too noisy so average them into hexagons.
datMap <- read.csv(file.path(tsDir, "patch_year_series.csv"))
datMap <- datMap[datMap$in_analysis & datMap$years_since_loss>=REL_MIN & datMap$years_since_loss<=REL_MAX,]

# Mean before the loss for each patch
pchChange <- aggregate(midgreendown_diff ~ patch_id, data=datMap[datMap$years_since_loss<0,], FUN=mean)
names(pchChange) <- c("patch_id", "before")

# Mean in the late years for each patch
pchLate <- aggregate(midgreendown_diff ~ patch_id, data=datMap[datMap$years_since_loss>=min(LATE) & datMap$years_since_loss<=max(LATE),], FUN=mean)
names(pchLate) <- c("patch_id", "late")

# Putting them together patches with no late years yet were dropped
pchChange <- merge(pchChange, pchLate, by="patch_id")
pchChange$change <- pchChange$late - pchChange$before
pchChange <- merge(pchChange, pchInfo[,c("patch_id", "lon", "lat")], by="patch_id")
summary(pchChange)
# should be in the ballpark of the -3.3 days from the model
mean(pchChange$change)   

# Checking how many patches land in each hexagon 
ggplot(data=pchChange) +
  geom_hex(aes(x=lon, y=lat), binwidth=c(1, 1)) +
  borders("state", colour="gray40", fill=COL_LAND)  +
  coord_quickmap(xlim=range(pchChange$lon), ylim=range(pchChange$lat)) +
  theme_bw() +
  labs(title="Patches per hexagon", x="Longitude", y="Latitude")

# Average change in each hexagon, only hexagons with at least 10 patches 
ggplot(data=pchChange) +
  stat_summary_hex(aes(x=lon, y=lat, z=change), fun=function(z) ifelse(length(z)>=10, mean(z), NA), binwidth=c(1, 1)) +
  borders("state", colour="gray40", fill=COL_LAND) +
  scale_fill_gradient2(low="darkorange", mid="white", high="steelblue", midpoint=0, limits=c(-8, 8), oob=scales::squish, na.value=NA) +
  coord_quickmap(xlim=range(pchChange$lon), ylim=range(pchChange$lat)) +
  theme_void() +
  theme(text=element_text(color=COL_TEXT, size=16)) +
  labs(fill="Change in\nmid-greendown\n(days)")

ggsave(file.path(figDir, "map_midgreendown_change.png"), width=11, height=9, dpi=300, bg="transparent")

