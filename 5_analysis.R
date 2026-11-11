# Patch Step 5: does forest loss change greenness and fall timing?
# -----------------------------------------------------------------------------
# Runs entirely in R. Run line by line (Ctrl+Enter) or use Source with Echo
# (Ctrl+Shift+Enter) so the tables and plots show. A plain "Source" still
# saves every file but won't show them. With every patch this can take a
# while (the models have many thousands of patch-years).
# Needs: dplyr, tidyr, ggplot2, lme4, lmerTest
# Reads:  scaleup_patches/timeseries_full/patch_year_series.csv   (Patch Step 4)
# Saves (in scaleup_patches/results_full/):
#   main_results.csv            early and late effects: estimate, SE, df, 95% CI, p, Holm p,
#                               SD before the loss and the effect in SD units
#   event_study_coefficients.csv, pretrend_tests.csv, recovery_rate.csv, analysis_set_patches.csv
#   event_study_main.png, event_study_all.png, main_effects.png,
#   trajectories_patch_vs_control.png
#
# THE COMPARISON
#   For every patch and year, Step 4 already subtracted the patch's nearby
#   uncut forest (every qualifying MODIS pixel 1-10 km away):
#   difference = patch value - comparison value. That removes weather and
#   other year-to-year swings shared by the patch and its surroundings.
#
# What it does (every patch in the analysis set, years -5 to +10):
#   1. Describes the analysis set (regions, loss years, sizes).
#   2. EVENT STUDY, one model per variable: every year compared with the year
#      before the loss (-1), with random effects for patch, Level III
#      ecoregion and calendar year. The pre-loss years (-5 to -2) should sit
#      near 0: that's the check that patches and comparison forest were
#      behaving alike before the loss. This is the "effect by years since loss"
#      picture.
#   3. MAIN RESULTS (model C from the review checks): early (years 1-3) and
#      late (4-10) vs all years before the loss (-5 to -1), the loss year left
#      out. Each patch gets its own loss effect (random slopes), so the effect
#      is tested across patches, not patch-years (no pseudoreplication).
#      Calendar year and Level III ecoregion are random effects too. Reported
#      with SE, df, 95% CI, p and Holm-adjusted p (5 variables x 2 periods).
#      NDVI and MidGreendown are the main outcomes, chosen beforehand.
#   4. RECOVERY RATE: from year 1 to year 10, how much the effect changes per
#      year (each patch compared with its own years before the loss).
#   5. Figures.
#
# Reading the numbers: NDVI differences are in NDVI units. Phenology
# differences are in days; positive = LATER than the comparison forest.
# The region, latitude, covariate and model checks are in Patch Step 6.
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

outDir <- file.path("~/Google Drive/My Drive/Reidy_research", "scaleup_patches")
tsDir  <- file.path(outDir, "timeseries_full")
resDir <- file.path(outDir, "results_full")
dir.create(resDir, recursive = TRUE, showWarnings = FALSE)

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

# Estimate, SE, df (Satterthwaite), t and p for some terms of a model
coefs <- function(m, terms = NULL) {
  co <- summary(m)$coefficients
  if (is.null(terms)) terms <- rownames(co)
  data.frame(term = terms, estimate = co[terms, "Estimate"], se = co[terms, "Std. Error"],
             df = co[terms, "df"], t = co[terms, "t value"], p = co[terms, "Pr(>|t|)"],
             row.names = NULL)
}

# ---- 1. The analysis set ----------------------------------------------------
ser <- read.csv(file.path(tsDir, "patch_year_series.csv"))
d <- ser %>%
  filter(as.logical(in_analysis), years_since_loss >= REL_MIN, years_since_loss <= REL_MAX) %>%
  mutate(patch_id  = factor(patch_id),
         us_l3code = factor(us_l3code),
         na_l2code = factor(na_l2code),
         year_f    = factor(year))

patches <- d %>% distinct(patch_id, na_l2code, us_l3code, loss_year, area_ha, n_controls)
nrow(patches)                          # patches in the analysis set
nrow(d)                                # patch-years
table(patches$na_l2code)               # patches by Level II region
table(patches$loss_year)               # patches by year of loss
round(summary(patches$area_ha), 1)     # patch size (ha)
summary(patches$n_controls)            # comparison pixels per patch
write.csv(patches, file.path(resDir, "analysis_set_patches.csv"), row.names = FALSE)

# One variable's differences (patch minus comparison)
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

# ---- 2. Event study: every year vs the year before the loss -----------------
event_rows <- list()
pretrend   <- list()
for (var in names(OUTCOMES)) {
  dd <- prep(var) %>% mutate(rel = relevel(factor(years_since_loss), ref = "-1"))
  E  <- fit("y ~ rel + (1 | us_l3code) + (1 | patch_id) + (1 | year_f)", dd)
  co <- coefs(E$m)
  co <- co[grepl("^rel", co$term), ] %>% mutate(rel_year = as.integer(sub("^rel", "", term)))
  event_rows[[var]] <- bind_rows(co, data.frame(term = "rel-1", estimate = 0, se = 0, df = NA,
                                                t = NA, p = NA, rel_year = -1L)) %>%
    mutate(outcome = var, note = E$note, .before = 1)
  
  # Joint test: are the pre-loss years (-5 to -2) all ~0?
  fe <- names(fixef(E$m))
  pre_terms <- intersect(paste0("rel", -5:-2), fe)
  L <- matrix(0, nrow = length(pre_terms), ncol = length(fe), dimnames = list(pre_terms, fe))
  for (j in seq_along(pre_terms)) L[j, pre_terms[j]] <- 1
  ct <- lmerTest::contest(E$m, L, joint = TRUE)
  pretrend[[var]] <- data.frame(outcome = var, pre_years_tested = length(pre_terms),
                                F_value = ct$`F value`, p = ct$`Pr(>F)`,
                                patches = n_distinct(dd$patch_id), patch_years = nrow(dd))
}
event <- bind_rows(event_rows) %>%
  mutate(ci_lo = ifelse(is.na(df), estimate, estimate - qt(0.975, df) * se),
         ci_hi = ifelse(is.na(df), estimate, estimate + qt(0.975, df) * se),
         label = factor(OUTCOMES[outcome], levels = OUTCOMES)) %>%
  arrange(outcome, rel_year)
pretrend <- bind_rows(pretrend)
write.csv(event,    file.path(resDir, "event_study_coefficients.csv"), row.names = FALSE)
write.csv(pretrend, file.path(resDir, "pretrend_tests.csv"),           row.names = FALSE)

# Pre-loss check (p > 0.05 = no sign the patches were already drifting before the loss)
pretrend %>% mutate(across(c(F_value, p), ~ round(.x, 3))) %>% as.data.frame()

# Main outcomes, every year vs year -1
event %>%
  filter(outcome %in% PRIMARY) %>%
  transmute(outcome, rel_year, estimate = signif(estimate, 3), se = signif(se, 2),
            ci = sprintf("%.3g to %.3g", ci_lo, ci_hi), p = signif(p, 2)) %>%
  as.data.frame()

# ---- 3. Main results: early and late vs before (model C) --------------------
main_rows <- list()
sd_pre    <- c()
for (var in names(OUTCOMES)) {
  dd <- prep_period(var)
  sd_pre[var] <- sd(dd$y[dd$period == "pre"])   # SD of the difference before the loss
  
  # Random slopes by patch: correlated if the model can estimate them, otherwise uncorrelated
  C  <- fit("y ~ early + late + (1 + early + late | patch_id) + (1 | us_l3code) + (1 | year_f)", dd)
  re <- "correlated random slopes"
  if (isSingular(C$m)) {
    C  <- fit("y ~ early + late + (1 + early + late || patch_id) + (1 | us_l3code) + (1 | year_f)", dd)
    re <- "uncorrelated random slopes"
  }
  main_rows[[var]] <- coefs(C$m, c("early", "late")) %>%
    mutate(outcome = var, period = term, patches = n_distinct(dd$patch_id),
           patch_years = nrow(dd), random_effects = re, note = C$note, .before = 1)
}

main <- bind_rows(main_rows) %>%
  mutate(p_holm     = p.adjust(p, method = "holm"),   # Holm across the 10 tests
         ci_lo      = estimate - qt(0.975, df) * se,
         ci_hi      = estimate + qt(0.975, df) * se,
         sd_pre     = unname(sd_pre[outcome]),
         std_effect = estimate / sd_pre,               # effect in SDs of the pre-loss difference
         primary    = outcome %in% PRIMARY,
         label      = factor(OUTCOMES[outcome], levels = OUTCOMES))
write.csv(main, file.path(resDir, "main_results.csv"), row.names = FALSE)

# Effect of forest loss (patch minus comparison, vs before the loss)
main %>%
  transmute(outcome, period, estimate = signif(estimate, 3), se = signif(se, 2), df = round(df),
            ci = sprintf("%.3g to %.3g", ci_lo, ci_hi), p = signif(p, 2), p_holm = signif(p_holm, 2),
            sd_pre = signif(sd_pre, 2), std_effect = round(std_effect, 2), patches, primary) %>%
  as.data.frame()

# ---- 4. Recovery rate -------------------------------------------------------
recov_rows <- list()
for (var in names(OUTCOMES)) {
  dd <- prep(var)
  pre_mean <- dd %>%
    filter(years_since_loss < 0) %>%
    group_by(patch_id) %>%
    summarise(pre_mean = mean(y), .groups = "drop")
  dp <- dd %>%
    filter(years_since_loss >= 1) %>%
    inner_join(pre_mean, by = "patch_id") %>%
    mutate(y_rel = y - pre_mean,              # change from the patch's own years before the loss
           yr    = years_since_loss - 1)      # 0 = year 1
  Rv <- fit("y_rel ~ yr + (1 + yr | patch_id) + (1 | us_l3code) + (1 | year_f)", dp)
  re <- "correlated random slopes"
  if (isSingular(Rv$m)) {
    Rv <- fit("y_rel ~ yr + (1 + yr || patch_id) + (1 | us_l3code) + (1 | year_f)", dp)
    re <- "uncorrelated random slopes"
  }
  recov_rows[[var]] <- coefs(Rv$m, c("(Intercept)", "yr")) %>%
    mutate(outcome = var, quantity = c("effect in year 1", "change per year after that"),
           patches = n_distinct(dp$patch_id), random_effects = re, note = Rv$note, .before = 1)
}
recov <- bind_rows(recov_rows) %>%
  mutate(ci_lo = estimate - qt(0.975, df) * se, ci_hi = estimate + qt(0.975, df) * se)
write.csv(recov, file.path(resDir, "recovery_rate.csv"), row.names = FALSE)

# Effect in year 1 and how much it changes each year after (NDVI: positive = recovering)
recov %>%
  transmute(outcome, quantity, estimate = signif(estimate, 3), se = signif(se, 2), df = round(df),
            ci = sprintf("%.3g to %.3g", ci_lo, ci_hi), p = signif(p, 2), patches) %>%
  as.data.frame()

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
       x="Years since loss", y="Change relative to year -1",
       caption="Shaded band = 95% confidence interval from the model. Grey = years before the loss.")

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
       x="Years since loss", y="Change relative to year -1",
       caption="Shaded band = 95% confidence interval from the model. Grey = years before the loss.")

ggsave(file.path(resDir, "event_study_all.png"), width=8, height=12, dpi=200)

# Main results, early and late
ggplot(data=main) +
  facet_wrap(~ label, ncol=1, scales="free_x") + # one panel per variable, each with its own effect axis
  geom_vline(xintercept=0, linetype="dotted") +
  geom_pointrange(aes(y=factor(period, levels=c("late", "early")), x=estimate, xmin=ci_lo, xmax=ci_hi),
                  color="darkgreen") + # bars = 95% CI
  theme_bw() +
  labs(title="Effect of forest loss (model C)", y=NULL, x="Patch minus comparison, vs before the loss",
       caption="Bars = 95% confidence intervals. Early = years 1-3, late = years 4-10 after the loss.")

ggsave(file.path(resDir, "main_effects.png"), width=7, height=9, dpi=200)

# Raw values: patches vs their comparison forest
traj <- bind_rows(lapply(PRIMARY, function(var) {
  d %>%
    select(patch_id, years_since_loss, Patch = all_of(var), Comparison = all_of(paste0("ctrl_", var))) %>%
    pivot_longer(c(Patch, Comparison), names_to = "series", values_to = "value") %>%
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
  scale_fill_manual(name="", values=c("Patch"="darkorange", "Comparison"="darkgreen")) +
  scale_color_manual(name="", values=c("Patch"="darkorange", "Comparison"="darkgreen")) +
  scale_x_continuous(breaks=REL_MIN:REL_MAX) +
  theme_bw() +
  labs(title="Raw values: disturbed patches vs their nearby uncut forest", x="Years since loss", y=NULL,
       caption="Shaded band = 95% confidence interval of the mean across patches (not model-adjusted).")

ggsave(file.path(resDir, "trajectories_patch_vs_control.png"), width=8, height=7, dpi=200)