# Patch Step 7: robustness checks for the main results
# -----------------------------------------------------------------------------
# Runs entirely in R. Run line by line (Ctrl+Enter) or use Source with Echo
# (Ctrl+Shift+Enter) so the tables and plots show. It fits about 35 mixed
# models on every patch-year, so expect it to take a while.
# Needs: dplyr, tidyr, ggplot2, lme4, lmerTest
# Reads (scaleup_patches/timeseries_full/, from Patch Step 4):
#   patch_year_series.csv, patch_summary.csv
# Saves (scaleup_patches/results_full/robustness/):
#   robustness_results.csv, robustness_primary.png,
#   fall_decline_by_year.png
#
# Each check reruns the main model (model B, Step 5) on a different version
# of the data. If the main results are sound, the estimates should stay close.
#   main             as in Step 5
#   before -5 to -2  "before" = years -5 to -2, year -1 left out. NDVI and
#                    MidGreendown already move a little in year -1, maybe
#                    because the disturbance starts before Hansen's loss year.
#   main regions     only Level II regions with at least MIN_REGION_PATCHES
#                    patches (drops the 3-12 patch regions 8.2, 8.5, 9.2)
#   loss to 2016     only losses up to 2016, whose regrowth is visible in NLCD
#                    2021 (catches pine replanting the evergreen rule missed)
#   strict hardwood  only patches < 25% evergreen and >= 50% deciduous or
#                    mixed forest in NLCD 2021
#   balanced -3..5   years -3 to +5 only, and only patches with data in every
#                    one of those years (late = years 4-5), so every year
#                    since loss has the same patches
# Plus one extra outcome, "fall decline" = MidGreendown - Senescence (days
# from 90% to 50% of the seasonal EVI range). If Senescence moves later while
# MidGreendown moves earlier, this gets shorter: the autumn drop got steeper,
# or the curve's shape changed (Patch Step 8 checks the shape directly).
#
# All error bars are 95% confidence intervals.
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
              dormancy     = "Dormancy (day of year)",
              decline      = "Fall decline (MidGreendown - Senescence, days)")
PRIMARY            <- c("ndvi", "midgreendown")
MIN_REGION_PATCHES <- 30     # "main regions" check: regions with at least this many patches
STRICT_EVERGREEN   <- 0.25   # "strict hardwood": evergreen share below this in 2021 ...
STRICT_DECID       <- 0.50   # ... and deciduous/mixed share at least this

outDir <- file.path("~/Google Drive/My Drive/Reidy_research", "scaleup_patches")
tsDir  <- file.path(outDir, "timeseries_full")
robDir <- file.path(outDir, "results_full", "robustness")
dir.create(robDir, recursive = TRUE, showWarnings = FALSE)

# ---- Model helpers ----------------------------------------------------------
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

# Model B: random slopes by patch (correlated if possible), Level III ecoregion
model_b <- function(dd) {
  r  <- fit("y ~ early + late + (1 + early + late | patch_id) + (1 | us_l3code)", dd)
  re <- "correlated random slopes"
  if (isSingular(r$m)) {
    r  <- fit("y ~ early + late + (1 + early + late || patch_id) + (1 | us_l3code)", dd)
    re <- "uncorrelated random slopes"
  }
  co <- summary(r$m)$coefficients[c("early", "late"), , drop = FALSE]
  data.frame(period = c("early", "late"), estimate = co[, "Estimate"], se = co[, "Std. Error"],
             df = co[, "df"], p = co[, "Pr(>|t|)"], patches = n_distinct(dd$patch_id),
             patch_years = nrow(dd), random_effects = re, note = r$note, row.names = NULL)
}

# ---- 1. Data ----------------------------------------------------------------
psum <- read.csv(file.path(tsDir, "patch_summary.csv")) %>%
  mutate(in_analysis = as.logical(in_analysis), na_l2code = as.character(na_l2code))

d <- read.csv(file.path(tsDir, "patch_year_series.csv")) %>%
  filter(as.logical(in_analysis)) %>%
  select(-na_l2code) %>%
  left_join(psum %>% select(patch_id, na_l2code, region, after_decid, after_evergreen), by = "patch_id") %>%
  mutate(decline_diff = (midgreendown - senescence) - (ctrl_midgreendown - ctrl_senescence),
         patch_id     = factor(patch_id),
         us_l3code    = factor(us_l3code))
n_distinct(d$patch_id)   # patches in the analysis set

region_counts <- psum %>% filter(in_analysis) %>% count(na_l2code, region) %>% as.data.frame()
region_counts            # patches per Level II region
main_regions <- region_counts$na_l2code[region_counts$n >= MIN_REGION_PATCHES]

# One variable, labelled into before / early / late (other years dropped)
prep_period <- function(var, data, pre = -5:-1, early = 1:3, late = 4:10) {
  data %>%
    mutate(y = .data[[paste0(var, "_diff")]]) %>%
    filter(!is.na(y)) %>%
    mutate(period = case_when(years_since_loss %in% pre   ~ "pre",
                              years_since_loss %in% early ~ "early",
                              years_since_loss %in% late  ~ "late")) %>%
    filter(!is.na(period)) %>%
    mutate(early = as.numeric(period == "early"),
           late  = as.numeric(period == "late"))
}

# Balanced: years -3 to +5, only patches with a value in every one of those
# years except the loss year (8 years)
prep_balanced <- function(var, data) {
  dd <- prep_period(var, data, pre = -3:-1, early = 1:3, late = 4:5)
  full <- dd %>% count(patch_id) %>% filter(n == 8) %>% pull(patch_id)
  dd %>% filter(patch_id %in% full)
}

# ---- 2. The checks ----------------------------------------------------------
CHECKS <- c(main            = "main (Step 5)",
            before_m5_m2    = "before = years -5 to -2",
            main_regions    = "main regions only",
            loss_to_2016    = "losses up to 2016",
            strict_hardwood = "strict hardwood regrowth",
            balanced        = "balanced, years -3 to 5")

check_data <- function(check, var) {
  switch(check,
         main            = prep_period(var, d),
         before_m5_m2    = prep_period(var, d, pre = -5:-2),
         main_regions    = prep_period(var, filter(d, na_l2code %in% main_regions)),
         loss_to_2016    = prep_period(var, filter(d, loss_year <= 2016)),
         strict_hardwood = prep_period(var, filter(d, after_evergreen < STRICT_EVERGREEN,
                                                   after_decid >= STRICT_DECID)),
         balanced        = prep_balanced(var, d))
}

rows <- list()
for (ck in names(CHECKS)) {
  for (var in names(OUTCOMES)) {
    message(format(Sys.time(), "%H:%M"), "  ", CHECKS[ck], ": ", OUTCOMES[var])   # progress
    rows[[paste(ck, var)]] <- model_b(check_data(ck, var)) %>%
      mutate(check = ck, outcome = var, .before = 1)
  }
}

rob <- bind_rows(rows) %>%
  group_by(check) %>%
  mutate(p_holm = p.adjust(p, method = "holm")) %>%   # Holm across the 12 tests in each check
  ungroup() %>%
  mutate(ci_lo       = estimate - qt(0.975, df) * se,
         ci_hi       = estimate + qt(0.975, df) * se,
         label       = factor(OUTCOMES[outcome], levels = OUTCOMES),
         period      = factor(period, levels = c("early", "late")),
         check_label = factor(CHECKS[check], levels = rev(CHECKS)))
write.csv(rob, file.path(robDir, "robustness_results.csv"), row.names = FALSE)

# Every check, every outcome (estimate and 95% CI)
rob %>%
  arrange(outcome, period, match(check, names(CHECKS))) %>%
  transmute(outcome, period, check, estimate = signif(estimate, 3),
            ci = sprintf("%.3g to %.3g", ci_lo, ci_hi), patches, p_holm = signif(p_holm, 2)) %>%
  as.data.frame()

# Main outcomes side by side: one column per check
rob %>%
  filter(outcome %in% PRIMARY) %>%
  transmute(outcome, period, check, estimate = signif(estimate, 3)) %>%
  pivot_wider(names_from = check, values_from = estimate) %>%
  as.data.frame()

# ---- 3. Figures -------------------------------------------------------------
ggplot(data=filter(rob, outcome %in% c(PRIMARY, "decline"))) +
  facet_wrap(~ label + period, ncol=2, scales="free_x") + # one panel per variable and period
  geom_vline(xintercept=0, linetype="dotted") +
  geom_pointrange(aes(y=check_label, x=estimate, xmin=ci_lo, xmax=ci_hi), color="darkgreen") + # bars = 95% CI
  theme_bw() +
  labs(title="Main results under each robustness check (model B)", y=NULL,
       x="Patch minus comparison, vs before the loss",
       caption="Bars = 95% confidence intervals. Early = years 1-3, late = years 4-10 (4-5 for the balanced check).")

ggsave(file.path(robDir, "robustness_primary.png"), width=10, height=9, dpi=200)

# Fall decline (MidGreendown - Senescence) by year since loss: patch minus comparison
decline_by_year <- d %>%
  filter(years_since_loss >= -5, years_since_loss <= 10, !is.na(decline_diff)) %>%
  group_by(years_since_loss) %>%
  summarise(n = n(), mean = mean(decline_diff), se = sd(decline_diff) / sqrt(n()), .groups = "drop")
decline_by_year

ggplot(data=decline_by_year) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dashed", color="red") + # year of loss
  geom_ribbon(aes(x=years_since_loss, ymin=mean - 1.96*se, ymax=mean + 1.96*se), fill="darkgreen", alpha=0.2) + # 95% CI
  geom_line(aes(x=years_since_loss, y=mean), color="darkgreen") +
  geom_point(aes(x=years_since_loss, y=mean), color="darkgreen", size=1.3) +
  scale_x_continuous(breaks=-5:10) +
  theme_bw() +
  labs(title="Fall decline (Senescence to MidGreendown): patch minus comparison forest",
       x="Years since loss", y="Days (negative = shorter decline in the patch)",
       caption="Shaded band = 95% confidence interval of the mean across patches (not model-adjusted).")

ggsave(file.path(robDir, "fall_decline_by_year.png"), width=8, height=5, dpi=200)