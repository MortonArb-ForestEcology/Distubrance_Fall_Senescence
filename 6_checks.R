# Patch Step 6: checks from the review
# -----------------------------------------------------------------------------
# Run line by line  or use Source with Echo
#  so the tabl and plots show. It fits about 60 mixed
# models on every patch-year, so expect it to take a while.
# Needs: dplyr, tidyr, ggplot2, lme4, lmerTest
# Reads (scaleup_patches/timeseries_full/, from Patch Step 4):
#   patch_summary.csv, patch_year_series.csv
# Saves: tables (.csv) and figures (.png) in scaleup_patches/results_full/checks/
#   (each section lists its files)
#
# Which check each section answers:
#   1. Comparison pixels: how many each patch has, and which rule removes the
#      pixels around patches in each region (pine plantations?)
#   2. Year of loss by region: is there a bias in when the losses happened?
#   3. Pseudoreplication: the old Step 5 model next to models that test the
#      loss effect across patches instead of across patch-years, and a simple
#      patch-level t-test. Model B is the main model (used in Step 5); C
#      (adds calendar year) is kept to show why calendar year is left out.
#   4. Spatial non-stationarity: Level II ecoregion as a fixed effect with an
#      interaction, and latitude (south to north)
#   5. Covariates: patch size, and how much of the core MODIS pixels was lost
#   6. One model for all five variables (in standard-deviation units)
# The effect by years since loss (every year with SE, CI and df, and the
# recovery rate) is in Patch Step 5.
#
# All shaded bands and error bars are 95% confidence intervals.
# NDVI effects are in NDVI units. Phenology effects are in days; negative =
# earlier than the comparison pixels.
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
PRIMARY      <- c("ndvi", "midgreendown")
REL_MIN      <- -5     # earliest year relative to the loss
REL_MAX      <- 10     # latest year relative to the loss
EARLY        <- 1:3    # "early" years after the loss
LATE         <- 4:10   # "late" years after the loss
MIN_CONTROLS <- 3      # must match Patch Step 4

outDir <- file.path("~/Google Drive/My Drive/Reidy_research", "scaleup_patches")
tsDir  <- file.path(outDir, "timeseries_full")
chkDir <- file.path(outDir, "results_full", "checks")
dir.create(chkDir, recursive = TRUE, showWarnings = FALSE)

# ---- 0. Read the inputs -----------------------------------------------------
# One row per patch (every patch from Step 3), with its comparison pixel counts
psum <- read.csv(file.path(tsDir, "patch_summary.csv")) %>%
  mutate(na_l2code   = as.character(na_l2code),
         in_analysis = as.logical(in_analysis),
         downloaded  = as.logical(downloaded))

# The analysis set, with the covariates centred on their averages
a_patches <- psum %>%
  filter(in_analysis) %>%
  mutate(lat_c       = lat - mean(lat),                                # degrees north of the average
         log2_area_c = log2(area_ha) - mean(log2(area_ha)),            # doublings of patch size
         core_loss_c = 10 * (core_loss_frac - mean(core_loss_frac)))   # per 10% more of the core pixels lost
nrow(a_patches)   # patches in the analysis set

region_n <- a_patches %>% count(na_l2code, region)
region_n          # analysis set by Level II region

# Patch-years for the models
d <- read.csv(file.path(tsDir, "patch_year_series.csv")) %>%
  filter(as.logical(in_analysis), years_since_loss >= REL_MIN, years_since_loss <= REL_MAX) %>%
  left_join(a_patches %>% select(patch_id, region, eco = na_l2code, lat_c, log2_area_c, core_loss_c),
            by = "patch_id") %>%
  mutate(patch_id  = factor(patch_id),
         us_l3code = factor(us_l3code),
         eco       = factor(eco),
         year_f    = factor(year))

# One variable's differences (patch minus comparison pixels)
prep <- function(var) {
  d %>%
    mutate(y = .data[[paste0(var, "_diff")]]) %>%
    filter(!is.na(y))
}
# ... with the before / early / late periods (the loss year itself is left out)
prep_period <- function(var) {
  prep(var) %>%
    mutate(period = case_when(years_since_loss < 0       ~ "pre",
                              years_since_loss %in% EARLY ~ "early",
                              years_since_loss %in% LATE  ~ "late")) %>%
    filter(!is.na(period)) %>%
    mutate(early = as.numeric(period == "early"),
           late  = as.numeric(period == "late"))
}

# ---- 1. Comparison pixels ---------------------------------------------------
# How the comparison works: for each loss patch, EVERY MODIS pixel 1-10 km
# away that looks like uncut hardwood forest is used: >= 75% deciduous or
# mixed forest in 2001, <= 5% of the pixel ever lost, <= 10% lost in it and
# its 8 neighbours, and elevation within 100 m of the patch. Each year, patch
# value minus the mean of those pixels is the difference used in the models.
# A patch needs at least 3 comparison pixels to be used.
# Saves: comparison_pixels_per_patch.csv, comparison_pixels_by_region.csv,
#        comparison_rules_by_region.csv, comparison_pixels_per_patch.png,
#        comparison_rules_by_region.png

# 1a. Comparison pixels for every patch, and how many pixels in its 1-10 km
#     ring passed each rule
comp_per_patch <- psum %>%
  filter(downloaded) %>%
  select(patch_id, region, loss_year, area_ha, n_ring, n_hardwood, n_pixel_loss, n_near_loss,
         n_elevation, n_controls, status, in_analysis) %>%
  mutate(ctrl_class = cut(n_controls, c(-Inf, 0, 2, 9, 99, Inf),
                          labels = c("0", "1-2", "3-9", "10-99", "100+"))) %>%
  arrange(region, patch_id) %>%
  as.data.frame()
head(comp_per_patch, 20)   # first 20; every patch is in the csv
write.csv(comp_per_patch, file.path(chkDir, "comparison_pixels_per_patch.csv"), row.names = FALSE)

# 1b. Summary by Level II region
comp_by_region <- comp_per_patch %>%
  group_by(region) %>%
  summarise(patches         = n(),
            with_3_or_more  = sum(n_controls >= MIN_CONTROLS),
            pct_3_or_more   = round(100 * mean(n_controls >= MIN_CONTROLS), 1),
            median_controls = median(n_controls),
            `0 pixels`      = sum(ctrl_class == "0"),
            `1-2 pixels`    = sum(ctrl_class == "1-2"),
            `3-9 pixels`    = sum(ctrl_class == "3-9"),
            `10-99 pixels`  = sum(ctrl_class == "10-99"),
            `100+ pixels`   = sum(ctrl_class == "100+"),
            .groups = "drop") %>%
  as.data.frame()
comp_by_region
write.csv(comp_by_region, file.path(chkDir, "comparison_pixels_by_region.csv"), row.names = FALSE)

ggplot(data=comp_per_patch) +
  geom_bar(aes(x=ctrl_class, fill=region)) +
  geom_vline(xintercept=2.5, linetype="dashed", color="red") + # left of the line = too few to use
  scale_x_discrete(drop=FALSE) +
  theme_bw() +
  labs(title="Comparison pixels per patch", x="Comparison pixels", y="Number of patches",
       fill="Level II ecoregion")

ggsave(file.path(chkDir, "comparison_pixels_per_patch.png"), width=9, height=5, dpi=200)

# 1c. Which rule removes the pixels? % of the pixels 1-10 km from the patches
#     that pass each rule on its own, by region. If "hardwood" is low in the
#     southern plains, that supports the pine plantation explanation.
rules_by_region <- comp_per_patch %>%
  group_by(region) %>%
  summarise(ring_pixels = sum(n_ring),
            hardwood    = round(100 * sum(n_hardwood)   / sum(n_ring), 1),
            pixel_loss  = round(100 * sum(n_pixel_loss) / sum(n_ring), 1),
            near_loss   = round(100 * sum(n_near_loss)  / sum(n_ring), 1),
            elevation   = round(100 * sum(n_elevation)  / sum(n_ring), 1),
            all_rules   = round(100 * sum(n_controls)   / sum(n_ring), 1),
            .groups = "drop") %>%
  as.data.frame()
rules_by_region
write.csv(rules_by_region, file.path(chkDir, "comparison_rules_by_region.csv"), row.names = FALSE)

rules_long <- rules_by_region %>%
  select(-ring_pixels) %>%
  pivot_longer(-region, names_to = "rule", values_to = "pct_pass") %>%
  mutate(rule = factor(rule, levels = c("hardwood", "pixel_loss", "near_loss", "elevation", "all_rules")))

ggplot(data=rules_long) +
  geom_col(aes(x=rule, y=pct_pass, fill=region), position="dodge") +
  theme_bw() +
  labs(title="Pixels 1-10 km from the patches passing each rule, by region",
       x="Rule", y="% of pixels passing", fill="Level II ecoregion")

ggsave(file.path(chkDir, "comparison_rules_by_region.png"), width=10, height=5, dpi=200)

# ---- 2. Year of loss by region ----------------------------------------------
# Late-loss patches can't reach years 6-10, so if one region's losses are
# mostly recent, the late-period effect leans on the other regions.
# Saves: loss_year_by_region.csv, loss_year_by_region.png, patches_per_year_since_loss.png
loss_year_tbl <- a_patches %>%
  count(loss_year, region) %>%
  pivot_wider(names_from = region, values_from = n, values_fill = 0) %>%
  arrange(loss_year) %>%
  as.data.frame()
loss_year_tbl
write.csv(loss_year_tbl, file.path(chkDir, "loss_year_by_region.csv"), row.names = FALSE)

# Typical year of loss in each region
a_patches %>%
  group_by(region) %>%
  summarise(patches = n(), median_loss_year = median(loss_year),
            earliest = min(loss_year), latest = max(loss_year), .groups = "drop") %>%
  as.data.frame()

kruskal.test(loss_year ~ factor(region), data = a_patches)   # p < 0.05 = regions differ in when their losses happened

ggplot(data=a_patches) +
  geom_bar(aes(x=factor(loss_year), fill=region)) +
  theme_bw() +
  labs(title="Year of loss in the analysis set, by region", x="Year of loss",
       y="Number of patches", fill="Level II ecoregion")

ggsave(file.path(chkDir, "loss_year_by_region.png"), width=9, height=5, dpi=200)

# How many patches have data at each year since the loss, by region
contrib <- d %>% filter(!is.na(ndvi_diff)) %>% distinct(patch_id, years_since_loss, region)

ggplot(data=contrib) +
  geom_bar(aes(x=factor(years_since_loss), fill=region)) +
  geom_vline(xintercept=which(sort(unique(contrib$years_since_loss)) == 0), linetype="dashed", color="red") + # year of loss
  theme_bw() +
  labs(title="Patches with July NDVI data at each year since the loss", x="Years since loss",
       y="Number of patches", fill="Level II ecoregion")

ggsave(file.path(chkDir, "patches_per_year_since_loss.png"), width=9, height=5, dpi=200)

# ---- 3. Pseudoreplication ---------------------------------------------------
# Four ways of testing the early (years 1-3) and late (4-10) effects, each
# against all years before the loss (-5 to -1):
#   A  Old Step 5 model: random intercepts only. Treats every patch-year as a
#      replicate, so the p-values are too small.
#   B  + random slopes: each patch gets its own loss effect, so the effect is
#      tested across patches (df near the number of patches). THIS IS THE
#      MAIN MODEL (Step 5).
#   C  + calendar year: patches measured in the same year share its weather,
#      so year is a random effect too. With every patch this moved the late
#      effects away from the raw data and the t-test: calendar year = year of
#      loss + years since loss, so the model can trade the two off. Kept as a
#      check only.
#   T  Simple cross-check: one number per patch (mean after minus mean
#      before), then a t-test across patches.
# If a model with correlated random slopes is "singular" (can't estimate all
# the correlations), it is refitted with uncorrelated slopes; the
# random_effects column says which was used, and note lists any warnings.
# Saves: model_comparison.csv, main_results.csv, model_comparison.png

MODEL_LABELS <- c(A = "A: old Step 5 model",
                  B = "B: + random slopes by patch (main)",
                  C = "C: + calendar year (check)",
                  T = "T: patch-level t-test")

ctrl_lmer <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

# Fit one model; keep any warnings in $note instead of printing them
fit <- function(f, data) {
  notes <- character(0)
  m <- withCallingHandlers(
    lmer(as.formula(f), data = data, control = ctrl_lmer),
    warning = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(x) { notes <<- c(notes, trimws(conditionMessage(x))); invokeRestart("muffleMessage") })
  list(m = m, note = paste(unique(notes), collapse = "; "))
}

# Fit with random slopes by patch: correlated if possible, otherwise uncorrelated
fit_slopes <- function(fixed, rest, data) {
  r <- fit(paste(fixed, "+ (1 + early + late | patch_id)", rest), data)
  r$re <- "correlated random slopes"
  if (isSingular(r$m)) {
    r <- fit(paste(fixed, "+ (1 + early + late || patch_id)", rest), data)
    r$re <- "uncorrelated random slopes"
  }
  r
}

# Likelihood-ratio test: does adding the interaction terms improve the model?
lrt <- function(f_add, f_int, rest, data) {
  int <- fit_slopes(f_int, rest, data)
  bar <- if (int$re == "correlated random slopes") "|" else "||"
  add <- fit(paste(f_add, "+ (1 + early + late", bar, "patch_id)", rest), data)
  a   <- suppressMessages(anova(add$m, int$m))   # refits both with maximum likelihood
  list(int  = int,
       test = data.frame(chisq = a$Chisq[2], df = a$Df[2], p = a$`Pr(>Chisq)`[2]))
}

# Estimate, SE, df (Satterthwaite), t and p for some terms of a model
coefs <- function(m, terms = NULL) {
  co <- summary(m)$coefficients
  if (is.null(terms)) terms <- rownames(co)
  data.frame(term = terms, estimate = co[terms, "Estimate"], se = co[terms, "Std. Error"],
             df = co[terms, "df"], t = co[terms, "t value"], p = co[terms, "Pr(>|t|)"],
             row.names = NULL)
}

cmp_rows <- list()
sd_pre   <- c()
for (var in names(OUTCOMES)) {
  dd <- prep_period(var)
  sd_pre[var] <- sd(dd$y[dd$period == "pre"])   # SD of the difference before the loss
  
  fits <- list(
    A = c(fit("y ~ early + late + (1 | us_l3code) + (1 | patch_id)", dd), re = "random intercepts only"),
    B = fit_slopes("y ~ early + late", "+ (1 | us_l3code)", dd),
    C = fit_slopes("y ~ early + late", "+ (1 | us_l3code) + (1 | year_f)", dd))
  for (nm in names(fits)) {
    cmp_rows[[paste(var, nm)]] <- coefs(fits[[nm]]$m, c("early", "late")) %>%
      mutate(outcome = var, period = term, model = nm,
             random_effects = fits[[nm]]$re, note = fits[[nm]]$note, .before = 1)
  }
  
  # T: one number per patch and period, then a t-test across patches
  pc <- dd %>%
    group_by(patch_id) %>%
    summarise(pre   = mean(y[period == "pre"]),
              early = mean(y[period == "early"]) - pre,
              late  = mean(y[period == "late"])  - pre,
              .groups = "drop")
  for (per in c("early", "late")) {
    x  <- pc[[per]][is.finite(pc[[per]])]
    tt <- t.test(x)
    cmp_rows[[paste(var, "T", per)]] <- data.frame(
      outcome = var, period = per, model = "T", term = per,
      estimate = unname(tt$estimate), se = tt$stderr, df = unname(tt$parameter),
      t = unname(tt$statistic), p = tt$p.value,
      random_effects = "none (one value per patch)", note = paste(length(x), "patches"))
  }
}

cmp <- bind_rows(cmp_rows) %>%
  group_by(model) %>%
  mutate(p_holm = p.adjust(p, method = "holm")) %>%   # Holm across the 10 tests in each model
  ungroup() %>%
  mutate(ci_lo       = estimate - qt(0.975, df) * se,
         ci_hi       = estimate + qt(0.975, df) * se,
         label       = factor(OUTCOMES[outcome], levels = OUTCOMES),
         period      = factor(period, levels = c("early", "late")),
         model_label = factor(MODEL_LABELS[model], levels = rev(MODEL_LABELS)))
write.csv(cmp, file.path(chkDir, "model_comparison.csv"), row.names = FALSE)

# Effect of loss under each model (look at how df, the CI and p change from A to C and T)
cmp %>%
  arrange(outcome, period, model) %>%
  transmute(outcome, period, model, estimate = signif(estimate, 3), se = signif(se, 2),
            df = round(df), ci = sprintf("%.3g to %.3g", ci_lo, ci_hi),
            p = signif(p, 2), p_holm = signif(p_holm, 2)) %>%
  as.data.frame()

# Main results (model B): estimate, SE, df, 95% CI, p, and the effect in
# pre-loss SD units (std_effect: effect divided by the SD of the difference
# before the loss, so NDVI and days can be compared)
main <- cmp %>%
  filter(model == "B") %>%
  mutate(sd_pre     = unname(sd_pre[outcome]),
         std_effect = estimate / sd_pre,
         std_ci_lo  = ci_lo / sd_pre,
         std_ci_hi  = ci_hi / sd_pre,
         primary    = outcome %in% PRIMARY) %>%
  select(outcome, label, period, estimate, se, df, t, ci_lo, ci_hi, p, p_holm,
         sd_pre, std_effect, std_ci_lo, std_ci_hi, primary, random_effects, note)
write.csv(main, file.path(chkDir, "main_results.csv"), row.names = FALSE)

main %>%
  transmute(outcome, period, estimate = signif(estimate, 3), se = signif(se, 2), df = round(df),
            ci = sprintf("%.3g to %.3g", ci_lo, ci_hi), p = signif(p, 2), p_holm = signif(p_holm, 2),
            sd_pre = signif(sd_pre, 2), std_effect = round(std_effect, 2), random_effects) %>%
  as.data.frame()

ggplot(data=cmp) +
  facet_wrap(~ label + period, ncol=2, scales="free_x") + # one panel per variable and period, each with its own effect axis
  geom_vline(xintercept=0, linetype="dotted") +
  geom_pointrange(aes(y=model_label, x=estimate, xmin=ci_lo, xmax=ci_hi)) + # bars = 95% CI
  theme_bw() +
  labs(title="Effect of forest loss under each model (vs before the loss)",
       y=NULL, x="Patch minus comparison pixels",
       caption="Bars = 95% confidence intervals. Early = years 1-3, late = years 4-10 after the loss.")

ggsave(file.path(chkDir, "model_comparison.png"), width=10, height=12, dpi=200)

# ---- 4. Spatial non-stationarity --------------------------------------------
# Does the loss effect differ from place to place?
#   - Level II ecoregion as a FIXED effect, with a period x ecoregion
#     interaction: one early and one late effect per region. The test asks if
#     the regions differ more than chance would give.
#   - Latitude: does the effect change from south to north? (per degree north)
# Both keep the main model's random effects (slopes by patch, Level III ecoregion).
# Saves: region_effects.csv, latitude_effects.csv, region_effects.png
#        (tests go in interaction_tests.csv, section 5)
rest       <- "+ (1 | us_l3code)"
test_rows  <- list()
region_rows <- list()
lat_rows   <- list()
for (var in names(OUTCOMES)) {
  dd <- prep_period(var)
  
  R <- lrt("y ~ eco + early + late", "y ~ eco + eco:early + eco:late", rest, dd)
  co <- coefs(R$int$m)
  region_rows[[var]] <- co[grepl(":(early|late)$", co$term), ] %>%
    mutate(outcome   = var,
           na_l2code = sub("^eco(.*):(early|late)$", "\\1", term),
           period    = sub(".*:", "", term), .before = 1)
  test_rows[[paste(var, "eco")]] <- R$test %>%
    mutate(outcome = var, test = "period x Level II ecoregion", random_effects = R$int$re, .before = 1)
  
  L <- lrt("y ~ early + late + lat_c", "y ~ early + late + lat_c + early:lat_c + late:lat_c", rest, dd)
  co <- coefs(L$int$m)
  lat_rows[[var]] <- co[grepl(":", co$term), ] %>%
    mutate(outcome = var, period = ifelse(grepl("early", term), "early", "late"),
           meaning = "change in the effect per degree north", .before = 1)
  test_rows[[paste(var, "lat")]] <- L$test %>%
    mutate(outcome = var, test = "period x latitude", random_effects = L$int$re, .before = 1)
}

region_eff <- bind_rows(region_rows) %>%
  left_join(region_n, by = "na_l2code") %>%
  mutate(ci_lo     = estimate - qt(0.975, df) * se,
         ci_hi     = estimate + qt(0.975, df) * se,
         label     = factor(OUTCOMES[outcome], levels = OUTCOMES),
         region_lab = paste0(na_l2code, "\n(", n, " patches)"))
write.csv(region_eff, file.path(chkDir, "region_effects.csv"), row.names = FALSE)

# Early and late effects in each Level II region
region_eff %>%
  transmute(outcome, region, period, patches = n, estimate = signif(estimate, 3),
            se = signif(se, 2), df = round(df), ci = sprintf("%.3g to %.3g", ci_lo, ci_hi)) %>%
  as.data.frame()

lat_eff <- bind_rows(lat_rows) %>%
  mutate(ci_lo = estimate - qt(0.975, df) * se, ci_hi = estimate + qt(0.975, df) * se)
write.csv(lat_eff, file.path(chkDir, "latitude_effects.csv"), row.names = FALSE)

# Change in the early and late effects per degree north
lat_eff %>%
  transmute(outcome, period, estimate = signif(estimate, 3), se = signif(se, 2), df = round(df),
            ci = sprintf("%.3g to %.3g", ci_lo, ci_hi), p = signif(p, 2)) %>%
  as.data.frame()

ggplot(data=filter(region_eff, outcome %in% PRIMARY)) +
  facet_wrap(~ label, ncol=1, scales="free_y") + # one panel per variable
  geom_hline(yintercept=0, linetype="dotted") +
  geom_pointrange(aes(x=region_lab, y=estimate, ymin=ci_lo, ymax=ci_hi, color=period),
                  position=position_dodge(width=0.5)) + # bars = 95% CI
  scale_color_manual(name="After loss", values=c("early"="darkorange", "late"="darkred")) +
  theme_bw() +
  labs(title="Effect of forest loss by Level II ecoregion (fixed effect)",
       x="Level II ecoregion", y="Patch minus comparison pixels, vs before",
       caption="Bars = 95% confidence intervals. Early = years 1-3, late = years 4-10 after the loss.")

ggsave(file.path(chkDir, "region_effects.png"), width=8, height=7, dpi=200)

# ---- 5. Covariates ----------------------------------------------------------
# Does the effect depend on patch size or on how much of the core MODIS
# pixels was lost? (A patch whose core pixels are only half lost mixes in
# uncut forest, which should weaken the signal.)
# Saves: covariate_effects.csv, interaction_tests.csv
cov_rows <- list()
for (var in names(OUTCOMES)) {
  dd <- prep_period(var)
  K <- lrt("y ~ early + late + log2_area_c + core_loss_c",
           paste("y ~ early + late + log2_area_c + core_loss_c",
                 "+ early:log2_area_c + late:log2_area_c + early:core_loss_c + late:core_loss_c"),
           rest, dd)
  co <- coefs(K$int$m)
  cov_rows[[var]] <- co[grepl(":", co$term), ] %>%
    mutate(outcome = var, period = ifelse(grepl("early", term), "early", "late"),
           meaning = ifelse(grepl("area", term), "change in the effect per doubling of patch size",
                            "change in the effect per 10% more of the core pixels lost"), .before = 1)
  test_rows[[paste(var, "cov")]] <- K$test %>%
    mutate(outcome = var, test = "period x (patch size + core loss)", random_effects = K$int$re, .before = 1)
}

cov_eff <- bind_rows(cov_rows) %>%
  mutate(ci_lo = estimate - qt(0.975, df) * se, ci_hi = estimate + qt(0.975, df) * se)
write.csv(cov_eff, file.path(chkDir, "covariate_effects.csv"), row.names = FALSE)

cov_eff %>%
  transmute(outcome, period, meaning, estimate = signif(estimate, 3), se = signif(se, 2),
            df = round(df), ci = sprintf("%.3g to %.3g", ci_lo, ci_hi), p = signif(p, 2)) %>%
  as.data.frame()

# Likelihood-ratio tests for sections 4 and 5 (p < 0.05 = the effect differs
# by region / latitude / covariates)
tests <- bind_rows(test_rows) %>% mutate(chisq = round(chisq, 2), p = signif(p, 2))
write.csv(tests, file.path(chkDir, "interaction_tests.csv"), row.names = FALSE)
tests %>% arrange(test, outcome) %>% as.data.frame()


# ---- 6. One model for all five variables ------------------------------------
# Each variable is divided by its SD before the loss, so all five are in the
# same units, then stacked and fitted in one model with an early and a late
# effect per variable. Compared with the standardized effects from the
# separate main models (section 3).
# Saves: one_model_all_outcomes.csv, standardized_effects.png
dl <- bind_rows(lapply(names(OUTCOMES), function(var) {
  prep_period(var) %>%
    transmute(patch_id, us_l3code, year_f, early, late, outcome = var, z = y / sd_pre[var])
})) %>%
  mutate(outcome = factor(outcome, levels = names(OUTCOMES)))

G <- tryCatch(
  fit(paste("z ~ 0 + outcome + outcome:early + outcome:late",
            "+ (1 + early + late || patch_id) + (1 | patch_id:outcome)"), dl),
  error = function(e) { message("One-model fit failed: ", conditionMessage(e)); NULL })

std_rows <- main %>%
  transmute(outcome, period = as.character(period), source = "Separate models (model B)",
            estimate = std_effect, ci_lo = std_ci_lo, ci_hi = std_ci_hi, p = p, p_holm = p_holm)
if (!is.null(G)) {
  co <- coefs(G$m)
  co <- co[grepl(":(early|late)$", co$term), ]
  std_rows <- bind_rows(std_rows, co %>%
                          transmute(outcome = sub("^outcome(.*):(early|late)$", "\\1", term),
                                    period  = sub(".*:", "", term),
                                    source  = "One model, all variables",
                                    estimate, ci_lo = estimate - qt(0.975, df) * se, ci_hi = estimate + qt(0.975, df) * se,
                                    p, p_holm = p.adjust(p, method = "holm")))
}
std_all <- std_rows %>%
  mutate(label = factor(OUTCOMES[outcome], levels = rev(OUTCOMES)))
write.csv(std_all, file.path(chkDir, "one_model_all_outcomes.csv"), row.names = FALSE)
G$note   # any warnings from the one-model fit (empty = none)

# Effects in pre-loss SD units, separate models vs one model
std_all %>%
  transmute(outcome, period, source, estimate = round(estimate, 2),
            ci = sprintf("%.2f to %.2f", ci_lo, ci_hi), p_holm = signif(p_holm, 2)) %>%
  arrange(outcome, period, source) %>%
  as.data.frame()

ggplot(data=std_all) +
  facet_wrap(~ period) + # early and late
  geom_hline(yintercept=0, linetype="dotted") +
  geom_pointrange(aes(x=label, y=estimate, ymin=ci_lo, ymax=ci_hi, color=source),
                  position=position_dodge(width=0.5)) + # bars = 95% CI
  scale_color_manual(name="", values=c("Separate models (model B)"="darkgreen", "One model, all variables"="purple")) +
  coord_flip() +
  theme_bw() +
  theme(legend.position="bottom") +
  labs(title="Effect of forest loss in standard-deviation units", x=NULL,
       y="Effect / SD of the difference before the loss",
       caption="Bars = 95% confidence intervals. Early = years 1-3, late = years 4-10 after the loss.")

ggsave(file.path(chkDir, "standardized_effects.png"), width=9, height=5, dpi=200)