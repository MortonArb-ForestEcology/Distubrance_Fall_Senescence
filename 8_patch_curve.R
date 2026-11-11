# Patch Step 8: is the timing shift real, or a change in the greenness curve's shape?
# -----------------------------------------------------------------------------
# Normal run or  Source with Echo 
# Needs: rgee, googledrive, dplyr, tidyr, ggplot2, lme4, lmerTest
# Reads:  scaleup_patches/all_patches.csv                     (Patch Step 3)
#         scaleup_patches/timeseries_full/patch_year_series.csv (Patch Step 4)
#         the patch and core-pixel assets Step 4 uploaded
# Saves (scaleup_patches/results_full/curve_shape/):
#   curve_values_by_patch_year.csv, curve_shape_effects.csv,
#   timing_with_curve_checks.csv, curve_shape_by_year.png, timing_with_curve_checks.png
#
# WHY THIS CHECK
#   MCD12Q2 dates are where the EVI curve crosses fractions of its OWN seasonal
#   range (amplitude): Senescence = last day at or above 90% of the range,
#   MidGreendown = last day at or above 50%. After a clearcut the curve is
#   lower and flatter, and that alone can move those dates without any change
#   in when leaves actually turn. The patches show Senescence LATER but
#   MidGreendown EARLIER, which is what a change in curve shape could do.
#
# What it does:
#   1. Earth Engine, same patches, core pixels and comparison rules as Step 4:
#      MCD12Q2's own curve-shape numbers each year -- EVI minimum, EVI
#      amplitude (seasonal range) and EVI area (greenness summed from
#      greenup to dormancy) -- for each patch's core pixels and the mean of
#      its comparison pixels.
#   2. Does the curve change after the loss? Model B (Step 5) on patch minus
#      comparison for minimum, amplitude and area.
#   3. Do the timing effects hold up?
#        a. the main timing model again, with this year's amplitude and minimum
#           differences added as covariates
#        b. only patch-years whose amplitude is within AMP_TOL of the
#           comparison forest's (curves of similar shape)
#      If the timing effects mostly vanish in (a) and (b), they're largely a
#      curve-shape effect. If they stay, the timing shift is more believable.
#      (b) picks years using something the loss itself changes, so treat it as
#      a check, not the main result.
# =============================================================================

library(rgee)
library(googledrive)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)

ee_Initialize(drive = TRUE)

# ---- Settings ---------------------------------------------------------------
RUN_TAG            <- "v3"       # same as Patch Steps 1 and 4 (finds Step 4's assets)
CURVE_TAG          <- "v3curve"  # names the Earth Engine jobs for this step
CONTROL_MIN_M      <- 1000       # comparison rules: same as Step 4
CONTROL_MAX_M      <- 10000
CTRL_DECID_MIN     <- 0.75
CTRL_LOSS_MAX      <- 0.05
CTRL_NEAR_LOSS_MAX <- 0.10
ELEV_TOL_M         <- 100
MCD_LAST_YEAR      <- 2025
PART_SIZE          <- 1000       # must match Step 4 (the uploaded parts)
TILE_SCALE         <- 16
AMP_TOL            <- 0.15       # check 3b: patch amplitude within 15% of the comparison's
PRIMARY            <- c("midgreendown", "senescence")
TIMING             <- c(greenup      = "Greenup (day of year)",
                        senescence   = "Senescence (day of year)",
                        midgreendown = "MidGreendown (day of year)",
                        dormancy     = "Dormancy (day of year)")
CURVE              <- c(evimin  = "EVI minimum",
                        eviamp  = "EVI amplitude (seasonal range)",
                        eviarea = "EVI area (greenup to dormancy)")

outDir      <- file.path("~/Google Drive/My Drive/Reidy_research", "scaleup_patches")
tsDir       <- file.path(outDir, "timeseries_full")
curveDir    <- file.path(outDir, "results_full", "curve_shape")
partDir     <- file.path(curveDir, paste0("by_part_", CURVE_TAG))
dir.create(partDir, recursive = TRUE, showWarnings = FALSE)
driveFolder <- "Reidy_research"
assetFolder <- "projects/breidyee/assets"
patch_asset <- function(k) paste0(assetFolder, "/full_", RUN_TAG, "_patches_", k)
core_asset  <- function(k) paste0(assetFolder, "/full_", RUN_TAG, "_core_", k)

# Earth Engine helpers (same as Step 4)
gi <- function(x, tries = 4) {   # getInfo, retrying when the servers are busy
  for (k in seq_len(tries)) {
    out <- tryCatch(x$getInfo(), error = function(e) e)
    if (!inherits(out, "error")) return(out)
    if (k == tries || !grepl("timed out|Too many|429|503|Internal", conditionMessage(out))) stop(out)
    message("  Earth Engine is busy; trying again in 60 seconds...")
    Sys.sleep(60)
  }
}
asset_exists <- function(id) !is.null(tryCatch(ee$data$getAsset(id), error = function(e) NULL))
ee_states <- function() {   # newest state of each task, if Earth Engine will tell us
  tl <- tryCatch(ee$data$getTaskList(), error = function(e) NULL)
  if (is.null(tl) || length(tl) == 0) {
    return(data.frame(description = character(), state = character(), error = character()))
  }
  d <- data.frame(
    description = vapply(tl, function(x) if (is.null(x$description)) "" else x$description, ""),
    state       = vapply(tl, function(x) if (is.null(x$state)) "" else x$state, ""),
    error       = vapply(tl, function(x) if (is.null(x$error_message)) "" else x$error_message, ""),
    stringsAsFactors = FALSE)
  d[!duplicated(d$description), ]   # the list is newest first
}

# ---- 1. Patches and Step 4's uploaded parts ---------------------------------
allp <- read.csv(file.path(outDir, "all_patches.csv")) %>%
  arrange(patch_n) %>%
  mutate(part = ceiling(patch_n / PART_SIZE))
all_parts <- sort(unique(allp$part))
length(all_parts)   # Earth Engine jobs (same parts as Step 4)

no_asset <- all_parts[!vapply(all_parts, function(k) asset_exists(patch_asset(k)) &&
                                asset_exists(core_asset(k)), TRUE)]
if (length(no_asset) > 0) {
  stop("Step 4's uploaded parts are missing for part(s) ", paste(no_asset, collapse = ", "),
       ". Run Patch Step 4 first.", call. = FALSE)
}

# ---- 2. Comparison rules (same as Step 4) and the curve-shape image ---------
modisProj <- ee$ImageCollection("MODIS/061/MOD13A1")$first()$select("NDVI")$projection()
to_modis  <- function(img) img$reduceResolution(reducer = ee$Reducer$mean(), maxPixels = 65535L)

hansen  <- ee$Image("UMD/hansen/global_forest_change_2025_v1_13")
lossAny <- hansen$select("loss")$unmask(0)
decid01 <- ee$ImageCollection("USGS/NLCD_RELEASES/2019_REL/NLCD")$
  filter(ee$Filter$eq("system:index", "2001"))$first()$
  select("landcover")$remap(list(41L, 43L), list(1L, 1L), 0)

anyLossM  <- to_modis(lossAny)
nearLossM <- anyLossM$reduceNeighborhood(reducer = ee$Reducer$mean(),
                                         kernel  = ee$Kernel$square(radius = 1, units = "pixels"))
decidM    <- to_modis(decid01)
elevM     <- to_modis(ee$Image("USGS/SRTMGL1_003")$select("elevation"))$rename("elev")
validCtrl <- decidM$gte(CTRL_DECID_MIN)$And(anyLossM$lte(CTRL_LOSS_MAX))$And(nearLossM$lte(CTRL_NEAR_LOSS_MAX))

# MCD12Q2 curve shape each year (cycle 1). Stored values are scaled:
# minimum and amplitude x 0.0001, area x 0.1. Values outside the valid range are masked.
mcd <- ee$ImageCollection("MODIS/061/MCD12Q2")
try_years <- 2001:MCD_LAST_YEAR
n_imgs <- unlist(gi(ee$List(lapply(try_years, function(y) {
  mcd$filter(ee$Filter$calendarRange(y, y, "year"))$size()
}))))
mcd_years <- try_years[as.numeric(n_imgs) > 0]
range(mcd_years)   # MCD12Q2 years Earth Engine has

bands <- list()
for (yr in mcd_years) {
  img <- mcd$filter(ee$Filter$calendarRange(yr, yr, "year"))$first()
  mn  <- img$select("EVI_Minimum_1");   mn <- mn$updateMask(mn$lte(10000))$multiply(0.0001)
  am  <- img$select("EVI_Amplitude_1"); am <- am$updateMask(am$lte(10000))$multiply(0.0001)
  ar  <- img$select("EVI_Area_1");      ar <- ar$updateMask(ar$lte(3700))$multiply(0.1)
  bands[[length(bands) + 1]] <- mn$rename(paste0("evimin_", yr))
  bands[[length(bands) + 1]] <- am$rename(paste0("eviamp_", yr))
  bands[[length(bands) + 1]] <- ar$rename(paste0("eviarea_", yr))
}
curveImg <- Reduce(function(a, b) a$addBands(b), bands)

curve_names <- paste0(rep(names(CURVE), times = length(mcd_years)), "_", rep(mcd_years, each = 3))
ctrl_names  <- paste0("ctrl_", curve_names)

per_patch <- ee_utils_pyfunc(function(f) {
  core_pts <- ee$FeatureCollection(ee$List(f$get("core")))$geometry()   # the patch's core pixel centres
  ring     <- f$geometry()$buffer(CONTROL_MAX_M, 50)$
    difference(f$geometry()$buffer(CONTROL_MIN_M, 50), 50)               # 1-10 km from the patch
  pElev    <- ee$Number(elevM$reduceRegion(
    reducer = ee$Reducer$mean()$unweighted(), geometry = core_pts, crs = modisProj,
    tileScale = TILE_SCALE)$get("elev"))
  ctrlMask <- validCtrl$And(elevM$subtract(pElev)$abs()$lte(ELEV_TOL_M))
  
  pv <- curveImg$reduceRegion(reducer = ee$Reducer$mean()$unweighted(), geometry = core_pts,
                              crs = modisProj, maxPixels = 1e9)
  cv <- curveImg$updateMask(ctrlMask)$rename(as.list(ctrl_names))$
    reduceRegion(reducer = ee$Reducer$mean()$unweighted(), geometry = ring, crs = modisProj,
                 maxPixels = 1e9, tileScale = TILE_SCALE)
  ee$Feature(NULL, pv)$set(cv)$set("patch_id", f$get("patch_id"))
})

part_job <- function(k) {
  joined <- ee$Join$saveAll(matchesKey = "core")$apply(
    primary   = ee$FeatureCollection(patch_asset(k)),
    secondary = ee$FeatureCollection(core_asset(k)),
    condition = ee$Filter$equals(leftField = "patch_id", rightField = "patch_id"))
  ee$FeatureCollection(joined)$map(per_patch)
}

# ---- 3. Run the jobs and download (finished parts are skipped) --------------
export_name <- function(k) paste0("full_", CURVE_TAG, "_values_", k)
part_folder <- function(k) file.path(partDir, paste0("part_", k))
is_done     <- function(k) file.exists(file.path(part_folder(k), "done.txt"))
selectors   <- c("patch_id", curve_names, ctrl_names)

download_ready <- function(parts) {
  f <- drive_ls(path = driveFolder, pattern = paste0("full_", CURVE_TAG, "_values_"))
  for (k in parts) {
    nm  <- export_name(k)
    hit <- f[f$name == paste0(nm, ".csv") | startsWith(f$name, paste0(nm, "-")), ]
    if (nrow(hit) == 0) next
    dir.create(part_folder(k), recursive = TRUE, showWarnings = FALSE)
    for (j in seq_len(nrow(hit))) {
      drive_download(as_id(hit$id[j]), path = file.path(part_folder(k), hit$name[j]), overwrite = TRUE)
    }
    writeLines(hit$name, file.path(part_folder(k), "done.txt"))
  }
}

todo <- all_parts[!vapply(all_parts, is_done, TRUE)]
length(todo)   # parts not downloaded yet (0 = all done, skips straight to section 4)
failed <- data.frame()

if (length(todo) > 0) {
  download_ready(todo)                       # anything an earlier run finished
  todo <- todo[!vapply(todo, is_done, TRUE)]
  
  st <- ee_states()
  waiting <- st$description[st$state %in% c("READY", "RUNNING")]
  n_started <- 0
  for (k in todo) {
    if (export_name(k) %in% waiting) next   # still running from an earlier session
    t <- ee_table_to_drive(collection = part_job(k), description = export_name(k),
                           folder = driveFolder, fileFormat = "CSV",
                           selectors = selectors, timePrefix = FALSE)
    t$start()
    n_started <- n_started + 1
  }
  message(sprintf("Started %d Earth Engine jobs (%d were already running). Checking every 3 minutes...",
                  n_started, length(todo) - n_started))
  
  t0 <- Sys.time()
  repeat {
    Sys.sleep(180)
    download_ready(todo)
    todo <- todo[!vapply(todo, is_done, TRUE)]
    st <- ee_states()
    st <- st[st$description %in% export_name(todo), ]
    failed <- st[st$state %in% c("FAILED", "CANCELLED"), ]
    message(sprintf("%s  parts left: %d  (running: %d, waiting: %d, failed: %d)",
                    format(Sys.time(), "%H:%M"), length(todo), sum(st$state == "RUNNING"),
                    sum(st$state == "READY"), nrow(failed)))
    if (length(todo) == 0 || length(todo) == nrow(failed)) break
    if (difftime(Sys.time(), t0, units = "hours") > 48) break
  }
}

failed   # parts that failed (empty = none). Source again to retry them.

have <- all_parts[vapply(all_parts, is_done, TRUE)]
if (length(have) < length(all_parts)) {
  warning(length(all_parts) - length(have), " part(s) not downloaded yet. Source again to retry them; ",
          "the analysis below uses only the parts that are done.", call. = FALSE)
}

# ---- 4. Join the curve values to the patch-year table -----------------------
read_part <- function(k) {
  files <- list.files(part_folder(k), pattern = "\\.csv$", full.names = TRUE)
  bind_rows(lapply(files, function(p) read.csv(p, colClasses = c(patch_id = "character"))))
}
w <- bind_rows(lapply(have, read_part)) %>% distinct(patch_id, .keep_all = TRUE)
nrow(w)   # patches downloaded

curve_long <- w %>%
  select(patch_id, matches("_\\d{4}$")) %>%
  pivot_longer(-patch_id, names_to = c("metric", "year"), names_pattern = "^(.*)_(\\d{4})$") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(year = as.integer(year))

d <- read.csv(file.path(tsDir, "patch_year_series.csv")) %>%
  filter(as.logical(in_analysis), years_since_loss >= -5, years_since_loss <= 10) %>%
  inner_join(curve_long, by = c("patch_id", "year")) %>%
  mutate(evimin_diff  = evimin  - ctrl_evimin,
         eviamp_diff  = eviamp  - ctrl_eviamp,
         eviarea_diff = eviarea - ctrl_eviarea,
         amp_ratio    = eviamp / ctrl_eviamp,       # 1 = same seasonal range as the comparison forest
         patch_id     = factor(patch_id),
         us_l3code    = factor(us_l3code))
write.csv(d, file.path(curveDir, "curve_values_by_patch_year.csv"), row.names = FALSE)
n_distinct(d$patch_id)   # patches in the analysis set with curve values

# ---- 5. Models --------------------------------------------------------------
ctrl_lmer <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
fit <- function(f, data) {   # keep warnings in $note instead of printing them
  notes <- character(0)
  m <- withCallingHandlers(
    lmer(as.formula(f), data = data, control = ctrl_lmer),
    warning = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(x) { notes <<- c(notes, trimws(conditionMessage(x))); invokeRestart("muffleMessage") })
  list(m = m, note = paste(unique(notes), collapse = "; "))
}
# Model B (Step 5), optionally with extra covariates
model_b <- function(dd, covars = "") {
  r  <- fit(paste("y ~ early + late", covars, "+ (1 + early + late | patch_id) + (1 | us_l3code)"), dd)
  if (isSingular(r$m)) {
    r <- fit(paste("y ~ early + late", covars, "+ (1 + early + late || patch_id) + (1 | us_l3code)"), dd)
  }
  co <- summary(r$m)$coefficients[c("early", "late"), , drop = FALSE]
  data.frame(period = c("early", "late"), estimate = co[, "Estimate"], se = co[, "Std. Error"],
             df = co[, "df"], p = co[, "Pr(>|t|)"], patches = n_distinct(dd$patch_id),
             patch_years = nrow(dd), note = r$note, row.names = NULL)
}
prep_period <- function(var, data) {   # before = -5..-1, early = 1-3, late = 4-10
  data %>%
    mutate(y = .data[[paste0(var, "_diff")]]) %>%
    filter(!is.na(y)) %>%
    mutate(period = case_when(years_since_loss < 0       ~ "pre",
                              years_since_loss %in% 1:3  ~ "early",
                              years_since_loss %in% 4:10 ~ "late")) %>%
    filter(!is.na(period)) %>%
    mutate(early = as.numeric(period == "early"), late = as.numeric(period == "late"))
}
add_ci <- function(x) x %>% mutate(ci_lo = estimate - qt(0.975, df) * se, ci_hi = estimate + qt(0.975, df) * se)

# 5a. Does the curve change after the loss? (patch minus comparison)
curve_eff <- bind_rows(lapply(names(CURVE), function(var) {
  model_b(prep_period(var, d)) %>% mutate(outcome = var, .before = 1)
})) %>% add_ci() %>% mutate(label = CURVE[outcome])
write.csv(curve_eff, file.path(curveDir, "curve_shape_effects.csv"), row.names = FALSE)

curve_eff %>%
  transmute(outcome, period, estimate = signif(estimate, 3), ci = sprintf("%.3g to %.3g", ci_lo, ci_hi),
            p = signif(p, 2), patches) %>%
  as.data.frame()

# How similar are the curves? Share of patch-years with amplitude within AMP_TOL of the comparison's
d %>%
  mutate(when = case_when(years_since_loss < 0 ~ "before", years_since_loss <= 3 ~ "years 0-3",
                          TRUE ~ "years 4-10")) %>%
  group_by(when) %>%
  summarise(patch_years = sum(!is.na(amp_ratio)),
            median_amp_ratio = round(median(amp_ratio, na.rm = TRUE), 3),
            pct_similar = round(100 * mean(abs(amp_ratio - 1) <= AMP_TOL, na.rm = TRUE), 1),
            .groups = "drop") %>%
  as.data.frame()

# 5b. Timing effects three ways: main model, with curve covariates, similar-curve years only
d_cov <- d %>%
  filter(!is.na(eviamp_diff), !is.na(evimin_diff)) %>%
  mutate(amp_c = eviamp_diff - mean(eviamp_diff), min_c = evimin_diff - mean(evimin_diff))
d_sim <- d %>% filter(abs(amp_ratio - 1) <= AMP_TOL)

timing_rows <- list()
for (var in names(TIMING)) {
  message(format(Sys.time(), "%H:%M"), "  timing models: ", TIMING[var])   # progress
  timing_rows[[var]] <- bind_rows(
    model_b(prep_period(var, d_cov))                         %>% mutate(version = "main model"),
    model_b(prep_period(var, d_cov), "+ amp_c + min_c")      %>% mutate(version = "+ amplitude and minimum"),
    model_b(prep_period(var, d_sim))                         %>% mutate(version = "similar-curve years only")) %>%
    mutate(outcome = var, .before = 1)
}
timing <- bind_rows(timing_rows) %>%
  add_ci() %>%
  mutate(label   = factor(TIMING[outcome], levels = TIMING),
         period  = factor(period, levels = c("early", "late")),
         version = factor(version, levels = rev(c("main model", "+ amplitude and minimum",
                                                  "similar-curve years only"))))
write.csv(timing, file.path(curveDir, "timing_with_curve_checks.csv"), row.names = FALSE)

# Timing effects (days, patch minus comparison vs before), three ways
timing %>%
  arrange(outcome, period, desc(version)) %>%
  transmute(outcome, period, version, estimate = signif(estimate, 3),
            ci = sprintf("%.3g to %.3g", ci_lo, ci_hi), patches, patch_years) %>%
  as.data.frame()

# ---- 6. Figures -------------------------------------------------------------
curve_by_year <- d %>%
  select(patch_id, years_since_loss, all_of(paste0(names(CURVE), "_diff"))) %>%
  pivot_longer(-c(patch_id, years_since_loss), names_to = "outcome", values_to = "diff") %>%
  filter(!is.na(diff)) %>%
  mutate(outcome = sub("_diff$", "", outcome), label = factor(CURVE[outcome], levels = CURVE)) %>%
  group_by(label, years_since_loss) %>%
  summarise(mean = mean(diff), se = sd(diff) / sqrt(n()), .groups = "drop")

ggplot(data=curve_by_year) +
  facet_wrap(~ label, ncol=1, scales="free_y") + # one panel per curve measure
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dashed", color="red") + # year of loss
  geom_ribbon(aes(x=years_since_loss, ymin=mean - 1.96*se, ymax=mean + 1.96*se), fill="darkgreen", alpha=0.2) + # 95% CI
  geom_line(aes(x=years_since_loss, y=mean), color="darkgreen") +
  geom_point(aes(x=years_since_loss, y=mean), color="darkgreen", size=1.3) +
  scale_x_continuous(breaks=-5:10) +
  theme_bw() +
  labs(title="Shape of the seasonal EVI curve: patch minus comparison forest",
       x="Years since loss", y="Patch minus comparison",
       caption="Shaded band = 95% confidence interval of the mean across patches (not model-adjusted).")

ggsave(file.path(curveDir, "curve_shape_by_year.png"), width=8, height=8, dpi=200)

ggplot(data=timing) +
  facet_wrap(~ label + period, ncol=2, scales="free_x") + # one panel per timing variable and period
  geom_vline(xintercept=0, linetype="dotted") +
  geom_pointrange(aes(y=version, x=estimate, xmin=ci_lo, xmax=ci_hi), color="darkgreen") + # bars = 95% CI
  theme_bw() +
  labs(title="Timing effects with and without accounting for curve shape (model B)", y=NULL,
       x="Days, patch minus comparison, vs before the loss",
       caption="Bars = 95% confidence intervals. Early = years 1-3, late = years 4-10 after the loss.")

ggsave(file.path(curveDir, "timing_with_curve_checks.png"), width=10, height=9, dpi=200)