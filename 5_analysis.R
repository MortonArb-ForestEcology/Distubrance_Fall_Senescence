# Patch Step 5: 
# -----------------------------------------------------------------------------
# Runs entirely in R. Run line by lineor use Source with Echo
# Needs: dplyr, tidyr, ggplot2, lme4, lmerTest
# Reads:  scaleup_patches/timeseries/patch_year_series.csv   (Patch Step 4)
# Saves (in scaleup_patches/results/):
#   event_study_main.png, event_study_all.png, trajectories_patch_vs_control.png,
#   region_effects_exploratory.png
#   period_effects.csv, event_study_coefficients.csv, pretrend_tests.csv,
#   region_effects_exploratory.csv, analysis_set_patches.csv
#
# THE COMPARISON
#   For every patch and year, Step 4 already subtracted the patch's own nearby
#   uncut forest ("controls"). Difference = patch value - control value.
#   That removes weather and other year toyear swings shared by the patch and
#   its surroundings. Each patch is compared with its OWN never cut controls,
#   not with other patches cut in other years
#
# What it does (for the 100 patches in the analysis set, years -5 to +10):
#   1. Describes the analysis set (regions, loss years, sizes).
#   2. One model per variable:
#        difference ~ year-since-loss + (1 | ecoregion) + (1 | patch)
#      Each year is compared with the year before the loss (-1). The pre-loss
#      years (-5 to -2) should sit near 0: .
#   3. Summary effects: early (years 1-3) and late (4-10) vs before the loss,
#      for each variable, with p-values adjusted for testing 5 variables x 2
#      periods (Holm). NDVI and MidGreendown are the main outcomes.
#   4. I do have an exploratory by region section looking  early and late effects within each Level II
#      region that has at least 5 patches. But needs work as of 10/27
#   5. Figures: Model plots, raw patch vs control trajectories, region
#      effects(if it works)
#
#  NDVI differences are in NDVI units. Phenology
# differences are in days; positive = LATER than the controls, negative = earlier.
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)

# ---- Settings ---------------------------------------------------------------
OUTCOMES <- c(ndvi         = "July NDVI",
              greenup      = "Greenup (day of year)",
              senescence   = "Senescence (day of year)",
              midgreendown = "MidGreendown (day of year)",
              dormancy     = "Dormancy (day of year)")
PRIMARY  <- c("ndvi", "midgreendown")
REL_MIN  <- -5          # earliest year relative to the loss
REL_MAX  <- 10          # latest year relative to the loss
EARLY    <- 1:3         # "early" years after the loss
LATE     <- 4:10        # "late" years after the loss
MIN_PATCHES_REGION <- 5 # region results only where at least this many patches

outDir <- file.path("~/Google Drive/My Drive/Reidy_research", "scaleup_patches")
tsDir  <- file.path(outDir, "timeseries")
resDir <- file.path(outDir, "results")
dir.create(resDir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. The analysis set ----------------------------------------------------
ser <- read.csv(file.path(tsDir, "patch_year_series.csv"))
d <- ser %>%
  filter(as.logical(in_analysis), years_since_loss >= REL_MIN, years_since_loss <= REL_MAX) %>%
  mutate(patch_id  = factor(patch_id),
         us_l3code = factor(us_l3code),
         na_l2code = factor(na_l2code))

patches <- d %>% distinct(patch_id, na_l2code, us_l3code, loss_year, area_ha)
nrow(patches)                       # patches in the analysis set
nrow(d)                             # patch-years
table(patches$na_l2code)            # patches by Level II region
table(patches$loss_year)            # patches by year of loss
round(summary(patches$area_ha), 1)  # patch size (ha)
write.csv(patches, file.path(resDir, "analysis_set_patches.csv"), row.names = FALSE)

prep <- function(var) {
  d %>%
    mutate(y = .data[[paste0(var, "_diff")]]) %>%
    filter(!is.na(y))
}

# ---- 2. Event study ---------------------------------------------------------
event_rows <- list()
pretrend   <- list()
for (var in names(OUTCOMES)) {
  dd <- prep(var) %>% mutate(rel = relevel(factor(years_since_loss), ref = "-1"))
  m  <- lmer(y ~ rel + (1 | us_l3code) + (1 | patch_id), data = dd)
  co <- summary(m)$coefficients
  k  <- grep("^rel", rownames(co))
  event_rows[[var]] <- bind_rows(
    tibble(outcome = var, rel_year = as.integer(sub("^rel", "", rownames(co)[k])),
           estimate = co[k, "Estimate"], se = co[k, "Std. Error"], p = co[k, "Pr(>|t|)"]),
    tibble(outcome = var, rel_year = -1L, estimate = 0, se = 0, p = NA_real_))
  
  # Joint test: are the pre-loss years (-5 to -2) all ~0?
  pre_terms <- intersect(paste0("rel", -5:-2), names(fixef(m)))
  L <- matrix(0, nrow = length(pre_terms), ncol = length(fixef(m)),
              dimnames = list(pre_terms, names(fixef(m))))
  for (j in seq_along(pre_terms)) L[j, pre_terms[j]] <- 1
  ct <- lmerTest::contest(m, L, joint = TRUE)
  pretrend[[var]] <- tibble(outcome = var, pre_years_tested = length(pre_terms),
                            F_value = ct$`F value`, p = ct$`Pr(>F)`,
                            patches = n_distinct(dd$patch_id), patch_years = nrow(dd))
}
event    <- bind_rows(event_rows) %>%
  mutate(ci_lo = estimate - 1.96 * se, ci_hi = estimate + 1.96 * se,
         label = OUTCOMES[outcome]) %>%
  arrange(outcome, rel_year)
pretrend <- bind_rows(pretrend)
write.csv(event,    file.path(resDir, "event_study_coefficients.csv"), row.names = FALSE)
write.csv(pretrend, file.path(resDir, "pretrend_tests.csv"),           row.names = FALSE)

# Pre-loss check (p > 0.05 = no sign the patches were already drifting before the loss)
pretrend %>% mutate(across(c(F_value, p), ~ round(.x, 3))) %>% as.data.frame()

# ---- 3. Summary effects: early and late vs before ---------------------------
period_of <- function(r) case_when(r < 0 ~ "pre", r %in% EARLY ~ "early", r %in% LATE ~ "late")

fit_period <- function(dd, re_ecoregion = TRUE) {
  dd <- dd %>%
    mutate(period = factor(period_of(years_since_loss), levels = c("pre", "early", "late"))) %>%
    filter(!is.na(period))
  f <- if (re_ecoregion) y ~ period + (1 | us_l3code) + (1 | patch_id) else y ~ period + (1 | patch_id)
  m  <- lmer(f, data = dd)
  co <- summary(m)$coefficients
  tibble(period = c("early", "late"),
         estimate = co[c("periodearly", "periodlate"), "Estimate"],
         se       = co[c("periodearly", "periodlate"), "Std. Error"],
         p        = co[c("periodearly", "periodlate"), "Pr(>|t|)"],
         patches  = n_distinct(dd$patch_id))
}

period_eff <- bind_rows(lapply(names(OUTCOMES), function(var) {
  fit_period(prep(var)) %>% mutate(outcome = var, .before = 1)
})) %>%
  mutate(ci_lo   = estimate - 1.96 * se,
         ci_hi   = estimate + 1.96 * se,
         p_holm  = p.adjust(p, method = "holm"),
         primary = outcome %in% PRIMARY,
         label   = OUTCOMES[outcome])
write.csv(period_eff, file.path(resDir, "period_effects.csv"), row.names = FALSE)

# Effect of forest loss (patch minus controls, compared with before the loss)
period_eff %>%
  transmute(outcome, period, estimate = signif(estimate, 3),
            ci = sprintf("%.3g to %.3g", ci_lo, ci_hi),
            p = signif(p, 3), p_holm = signif(p_holm, 3), primary) %>%
  as.data.frame()

# ---- 4. By region (exploratory) ---------------------------------------------
regions <- names(which(table(patches$na_l2code) >= MIN_PATCHES_REGION))
region_eff <- bind_rows(lapply(PRIMARY, function(var) {
  bind_rows(lapply(regions, function(r) {
    dd <- prep(var) %>% filter(na_l2code == r)
    out <- tryCatch(fit_period(dd, re_ecoregion = FALSE), error = function(e) NULL)
    if (is.null(out)) return(NULL)
    out %>% mutate(outcome = var, na_l2code = r, .before = 1)
  }))
}))
if (nrow(region_eff) > 0) {
  region_eff <- region_eff %>%
    mutate(ci_lo = estimate - 1.96 * se, ci_hi = estimate + 1.96 * se, label = OUTCOMES[outcome])
}
write.csv(region_eff, file.path(resDir, "region_effects_exploratory.csv"), row.names = FALSE)

# By Level II region (exploratory; regions with >= MIN_PATCHES_REGION patches)
if (nrow(region_eff) > 0) {
  region_eff %>%
    transmute(outcome, region = na_l2code, period, patches, estimate = signif(estimate, 3),
              ci = sprintf("%.3g to %.3g", ci_lo, ci_hi)) %>%
    as.data.frame()
}

# ---- 5. Figures -------------------------------------------------------------
# Event study, main outcomes
ggplot(data=filter(event, outcome %in% PRIMARY)) +
  facet_wrap(~ label, ncol=1, scales="free_y") + # one panel per variable
  annotate("rect", xmin=REL_MIN - 0.5, xmax=-0.5, ymin=-Inf, ymax=Inf, fill="grey92") + # grey = before the loss
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dashed", color="red") + # year of loss
  geom_ribbon(aes(x=rel_year, ymin=ci_lo, ymax=ci_hi), fill="darkgreen", alpha=0.2) + # 95% CI
  geom_line(aes(x=rel_year, y=estimate), color="darkgreen") +
  geom_point(aes(x=rel_year, y=estimate), color="darkgreen", size=1.3) +
  scale_x_continuous(breaks=REL_MIN:REL_MAX) +
  theme_bw() +
  labs(title=sprintf("Ndvi/ Midgreendown vs year (n = %d patches)", nrow(patches)),
       x="Years since loss", y="Change relative to year -1")

ggsave(file.path(resDir, "event_study_main.png"), width=8, height=7, dpi=200)

# Event study, all five outcomes
ggplot(data=event) +
  facet_wrap(~ label, ncol=1, scales="free_y") + # one panel per variable
  annotate("rect", xmin=REL_MIN - 0.5, xmax=-0.5, ymin=-Inf, ymax=Inf, fill="grey92") + # grey = before the loss
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dashed", color="red") + # year of loss
  geom_ribbon(aes(x=rel_year, ymin=ci_lo, ymax=ci_hi), fill="darkgreen", alpha=0.2) + # 95% CI
  geom_line(aes(x=rel_year, y=estimate), color="darkgreen") +
  geom_point(aes(x=rel_year, y=estimate), color="darkgreen", size=1.3) +
  scale_x_continuous(breaks=REL_MIN:REL_MAX) +
  theme_bw() +
  labs(title=sprintf("All outcomes vs year (n = %d patches)", nrow(patches)),
       x="Years since loss", y="Change relative to year -1")

ggsave(file.path(resDir, "event_study_all.png"), width=8, height=12, dpi=200)

# Raw values: patches vs their controls
traj <- bind_rows(lapply(PRIMARY, function(var) {
  d %>%
    select(patch_id, years_since_loss, Patch = all_of(var), Controls = all_of(paste0("ctrl_", var))) %>%
    pivot_longer(c(Patch, Controls), names_to = "series", values_to = "value") %>%
    filter(!is.na(value)) %>%
    group_by(series, years_since_loss) %>%
    summarise(mean = mean(value), se = sd(value) / sqrt(n()), .groups = "drop") %>%
    mutate(label = OUTCOMES[var])
}))

ggplot(data=traj) +
  facet_wrap(~ label, ncol=1, scales="free_y") + # one panel per variable
  geom_vline(xintercept=0, linetype="dashed", color="red") + # year of loss
  geom_ribbon(aes(x=years_since_loss, ymin=mean - 1.96*se, ymax=mean + 1.96*se, fill=series), alpha=0.15) + # 95% CI
  geom_line(aes(x=years_since_loss, y=mean, color=series)) +
  scale_fill_manual(name="", values=c("Patch"="darkorange", "Controls"="darkgreen")) +
  scale_color_manual(name="", values=c("Patch"="darkorange", "Controls"="darkgreen")) +
  scale_x_continuous(breaks=REL_MIN:REL_MAX) +
  theme_bw() +
  labs(title="Raw values: disturbed patches vs their nearby uncut forest", x="Years since loss", y=NULL)

ggsave(file.path(resDir, "trajectories_patch_vs_control.png"), width=8, height=7, dpi=200)

# Effect by region (IDK abut this yet)
ggplot(data=region_eff) +
  facet_wrap(~ label, ncol=1, scales="free_y") + # one panel per variable
  geom_hline(yintercept=0, linetype="dotted") +
  geom_pointrange(aes(x=na_l2code, y=estimate, ymin=ci_lo, ymax=ci_hi, color=period),
                  position=position_dodge(width=0.5)) + # bars = 95% CI
  scale_color_manual(name="After loss", values=c("early"="darkorange", "late"="darkred")) +
  theme_bw() +
  labs(title="Effect by Level II region (exploratory)",
       subtitle="Early = years 1-3, late = years 4-10, each vs before the loss. Bars = 95% CI.",
       x="Level II ecoregion", y="Patch minus controls, vs before")

ggsave(file.path(resDir, "region_effects_exploratory.png"), width=8, height=7, dpi=200)