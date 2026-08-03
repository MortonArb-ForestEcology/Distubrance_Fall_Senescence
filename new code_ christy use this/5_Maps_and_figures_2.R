# Step 5: Before/After Extent Maps + NDVI Trends
# -----------------------------------------------------------------------------
# Part A — quick GEE-side preview (interactive, stays in rgee, good for
#          exploring the layers yourself)
# Part B — polished ggplot maps built from downloaded layers (good for
#          figures — full control over styling, legends, panel layout)
# Part C — NDVI yearly trends from Step 4's output, brought in two ways:
#          (1) an aggregated recovery trajectory (all patches pooled,
#              relative to each patch's own disturbance year), and
#          (2) a spatial map of per-patch NDVI change (post vs pre).
#          Relativizing to loss_year happens HERE, not in Step 4 — Step 4
#          deliberately stored absolute calendar years so that decision
#          could be made at analysis time.
# Part D — MODIS phenology (MidGreendown) trends from Step 4b's output,
#          mirroring Part C's structure but for fall-color transition
#          timing rather than NDVI. See Part D's header below for details.
# Part E — Pre/post histograms (NDVI + MidGreendown) and baseline-normalized
#          trajectories, where each patch's own pre-disturbance window is
#          rescaled to a common reference point (1 for NDVI, 0 for
#          MidGreendown — see Part E's header for why these differ).
# Part F — Deeper senescence-timing analyses: a regional (climate-controlled)
#          reference comparison (needs Step 4c's output), a per-patch
#          pre-disturbance TREND baseline instead of a flat mean, senescence
#          duration/rate, interannual + spatial variability, and a QA-
#          stratified check for data-quality confounds. See Part F's header.
# Part G — Mixed-effects model (patch + year crossed random effects) that
#          consolidates Part E's baseline normalization and Part F #2's
#          regional-climate comparison into one statistically principled
#          framework -- see Part G's header for why these turn out to be
#          the same idea. Also includes disturbance-year and patch-area
#          histograms, and runs the model on several subsets (all patches,
#          a random 10%, the 10 largest patches, and the patches with the
#          biggest NDVI change), surfacing only the statistically
#          significant results.
# Part H — Does patch size moderate the disturbance effect? A continuous
#          post*log(area_ha) interaction model on the full population, for
#          both NDVI and MidGreendown, with full residual and influence
#          diagnostics -- includes the robustness check that found NDVI's
#          size-effect holds up but MidGreendown's does not (hinges on a
#          single ordinary patch). See Part H's header for the full story.
# Part I — Disaggregates Part G's binary pre/post into individual years
#          (-5 to +5) as a factor, so a delayed or strengthening effect
#          isn't hidden by averaging all post-years together. Includes a
#          simple one-way ANOVA + Tukey HSD, the equivalent mixed-model
#          ANOVA, per-year effect-size plots, and a superimposed epoch
#          analysis (SEA) for both NDVI and MidGreendown.
# Part J — Does the MAGNITUDE of NDVI loss (rather than patch area, per
#          Part H) moderate the effect? Built as a genuine, non-circular
#          continuous interaction for MidGreendown; built as a clearly-
#          flagged DESCRIPTIVE-ONLY severity-tercile check for NDVI itself,
#          since using NDVI's own change to moderate NDVI would be
#          circular. See Part J's header for why these are treated
#          differently.
# Part K — Resolves the two questions Parts F/I/J left open. K1 re-runs
#          Part I's year-by-year period model on the CLIMATE-CONTROLLED
#          regional anomaly (patch value minus that year's undisturbed-
#          forest reference from Step 4c), which supersedes Part F #1's
#          coarser two-window version and settles whether the delayed
#          MidGreendown shift is real or regional drift. K2 runs the formal
#          per-patch severity regression that Part J only plotted.
#          REQUIRES Step 4c to have produced a non-empty
#          reference_phenology_by_year.csv.
#
# "Before" = NLCD forest cover extent / Hansen treecover2000, i.e. what was
#            there before disturbance.
# "After"  = current persistent-loss extent, and which patches/pixels meet
#            the >=0.75 forest-cover criterion (from Step 3's attribute
#            tables).
# -----------------------------------------------------------------------------
# Requires: rgee, sf, terra, tidyterra, ggplot2, patchwork, dplyr, tidyr,
#           lme4, lmerTest, purrr, scales, influence.ME
# =============================================================================

library(rgee)
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

ee_Initialize()

roi <- ee$Geometry$Point(c(-79.862539, 37.829550))$buffer(100000)

modisProjection <- ee$ImageCollection('MODIS/061/MOD13A1')$
  first()$select('NDVI')$projection()

# =============================================================================
# PART A — Quick GEE-side interactive preview
# =============================================================================
hansen <- ee$Image('UMD/hansen/global_forest_change_2024_v1_12')$clip(roi)
forestFractional <- ee$Image('projects/breidyee/assets/nlcd_forest_fraction_modis')

# .toInt32() cast: images built via reduceToImage (as hansen_loss_forest_masked
# was, in Step 3) often come out as 64-bit Long pixels, which GEE's export
# pipeline refuses outright ("Pixel type not supported: Type<Long>"). The
# actual values here are just lossyear codes (1-24), so int32 loses nothing.
hansenMasked <- ee$Image('projects/breidyee/assets/hansen_loss_forest_masked')$toInt32()

Map$setCenter(-79.862539, 37.829550, 10)

Map$addLayer(
  hansen$select('treecover2000')$updateMask(hansen$select('treecover2000')$gt(0)),
  list(min = 0, max = 100, palette = c('black', 'green')),
  'BEFORE: Tree Cover 2000'
) +
  Map$addLayer(
    forestFractional$select('forest_fraction_mean')$selfMask(),
    list(min = 0, max = 1, palette = c('white', 'lightgreen', 'darkgreen')),
    'BEFORE: Forest Fraction (NLCD, MODIS 500m)'
  ) +
  Map$addLayer(
    hansenMasked,
    list(min = 1, max = 24, palette = c('yellow', 'orange', 'red')),
    'AFTER: Persistent Loss, Forest-Masked (30m)'
  ) +
  Map$addLayer(
    forestFractional$select('forest_fraction_mean')$gte(0.75)$selfMask(),
    list(palette = c('blue')),
    'AFTER: Pixels Meeting >=0.75 Cover Threshold'
  )

# =============================================================================
# PART B — Polished ggplot maps from downloaded layers
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Download rasters locally
#    Using ee_as_rast (terra-native) rather than the deprecated ee_as_raster,
#    which requires the 'stars' package. ee_as_rast returns a SpatRaster
#    directly, so no rast() wrapper is needed either.
#    Note: hansenMasked at 30m over a 100km buffer is a large download —
#    if it's too slow, drop the scale or shrink the ROI for the figure.
# -----------------------------------------------------------------------------
forestFrac_r <- ee_as_rast(
  forestFractional, region = roi, scale = 500, crs = 'EPSG:4326', via = 'drive'
)
# <-- FIX: reproject via GEE at export time (crs argument above) rather than
# locally afterward with terra::project(). A local terra::project() call on
# the downloaded sinusoidal raster was introducing a real ~18km northward
# shift in the data itself (confirmed by comparing forestFractional's own
# GEE-side bounds against forestFrac_r's post-reprojection bounds) — not
# just an extent/display issue. Letting GEE do the reprojection during
# export avoids that shift entirely.

hansenMasked_r <- ee_as_rast(
  hansenMasked, region = roi, scale = 30, via = 'drive'  # native Hansen resolution
)

treecover2000_r <- ee_as_rast(
  hansen$select('treecover2000'), region = roi, scale = 30, via = 'drive'  # native resolution
)

# -----------------------------------------------------------------------------
# 2. Load patches + join the patch attribute table (from Step 3) for the
#    meets_forest_threshold flag
# -----------------------------------------------------------------------------
# Shared output directory — the SAME local folder Steps 1, 3, 4, 4b and 4c
# write to. Defined directly rather than sourced from a config.R: the config
# file earlier drafts referenced doesn't exist at that path, so every
# source() call failed, left outputDir undefined, and broke every read/write
# below it. Set once, here.
outputDir <- file.path("~/Google Drive/My Drive", "Reidy_research")
if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)

patches_sf <- st_read(file.path(outputDir, 'hansen_persistent_loss_patches.gpkg'))
patch_attrs <- read.csv(file.path(outputDir, 'patch_attribute_table.csv'))
patch_attrs$meets_forest_threshold <- as.logical(patch_attrs$meets_forest_threshold) 

patches_sf <- merge(
  patches_sf, patch_attrs[, c('patch_uuid', 'meets_forest_threshold', 'forest_cover_mean')],
  by = 'patch_uuid', all.x = TRUE
)

# forestFractional's underlying GEE asset is wider than roi (this is real,
# confirmed via forestFractional$geometry()$bounds()$getInfo() vs.
# roi$bounds()$getInfo() — not an artifact of the crs fix above), so crop
# forestFrac_r down to patches_sf's own verified-correct extent before any
# plotting uses it.
bbox_patches <- st_bbox(patches_sf)
crop_ext <- terra::ext(bbox_patches["xmin"], bbox_patches["xmax"],
                       bbox_patches["ymin"], bbox_patches["ymax"])
forestFrac_r <- terra::crop(forestFrac_r, crop_ext)

# -----------------------------------------------------------------------------
# 3. BEFORE / AFTER panel
# -----------------------------------------------------------------------------
p_before <- ggplot() +
  geom_spatraster(data = treecover2000_r) +
  scale_fill_gradient(low = 'grey95', high = 'darkgreen', na.value = NA,
                      name = '% tree\ncover 2000') +
  labs(title = 'Before: Tree Cover (2000)') +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank())

p_after <- ggplot() +
  geom_spatraster(data = treecover2000_r, fill = 'grey90') +
  geom_spatraster(data = hansenMasked_r) +
  scale_fill_gradient(low = 'yellow', high = 'red', na.value = NA,
                      name = 'loss\nyear (code)') +
  labs(title = 'After: Persistent Forest Loss (forest-masked)') +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank())

before_after_panel <- p_before + p_after
ggsave(file.path(outputDir, 'map_before_after.png'), before_after_panel,
       width = 12, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 4. % Cover choropleth
# -----------------------------------------------------------------------------
p_cover <- ggplot() +
  geom_spatraster(data = forestFrac_r, aes(fill = forest_fraction_mean)) +
  scale_fill_viridis_c(name = 'Forest\nfraction', na.value = NA, limits = c(0, 1)) +
  labs(title = 'Fractional Forest Cover (MODIS 500m)') +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank())

ggsave(file.path(outputDir, 'map_forest_cover_fraction.png'), p_cover,
       width = 8, height = 7, dpi = 300)

# -----------------------------------------------------------------------------
# 5. Patches meeting criteria
# -----------------------------------------------------------------------------
p_criteria <- ggplot() +
  geom_spatraster(data = forestFrac_r, aes(fill = forest_fraction_mean), alpha = 0.5) +
  scale_fill_viridis_c(name = 'Forest\nfraction', na.value = NA, limits = c(0, 1)) +
  ggnewscale::new_scale_fill() +
  geom_sf(data = patches_sf, aes(fill = meets_forest_threshold), color = NA) +
  scale_fill_manual(values = c(`TRUE` = 'red', `FALSE` = 'grey40'),
                    name = 'Meets >=0.75\ncover threshold') +
  labs(title = 'Persistent-Loss Patches Meeting Forest-Cover Criterion') +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank())

ggsave(file.path(outputDir, 'map_patches_meeting_criteria.png'), p_criteria,
       width = 8, height = 7, dpi = 300)

cat(sprintf('Maps saved to %s/:\n', outputDir),
    '  map_before_after.png\n',
    '  map_forest_cover_fraction.png\n',
    '  map_patches_meeting_criteria.png\n')

# =============================================================================
# PART C — NDVI Yearly Trends (from Step 4's ndvi_yearly_july_by_patch.csv)
# =============================================================================

ndvi_long <- read.csv(file.path(outputDir, 'ndvi_yearly_july_by_patch.csv'))

# <-- FIX: every figure in this project that touches NDVI has captioned
# itself "Patches >= 25 ha only (MODIS-resolvable filter, see Step 4)" --
# but that filter was never actually enforced anywhere. Confirmed directly:
# ndvi_long contains 114,095 distinct patches with area down to 0.12 ha,
# not the 615 patches that actually meet area_ha >= 25 in patch_attrs. Step
# 4's CSV export was apparently never restricted by area at the source, so
# every NDVI figure built so far (Part C's trajectory, Part E's histograms
# and normalized trajectory, Part G's models) was silently drawing on the
# full, unfiltered catalog -- including many sub-MODIS-pixel patches whose
# July NDVI is heavily diluted by whatever's in the surrounding
# neighborhood, not a clean signal from the patch itself. Enforcing the
# filter HERE, right at load, means the fix propagates to every part of
# the script that uses ndvi_long downstream, rather than needing to patch
# each usage site separately.
ndvi_long <- ndvi_long %>%
  select(-any_of('area_ha')) %>%
  left_join(patch_attrs %>% select(patch_uuid, area_ha), by = 'patch_uuid') %>%
  filter(area_ha >= 25)

# Drop rows with no valid NDVI (ndvi_july_count == 0 -> ndvi_july_mean is NA,
# per the fix in Step 4 that guarantees the column exists even when blank).
ndvi_long <- ndvi_long %>% filter(!is.na(ndvi_july_mean))

# Relativize to each patch's own loss_year HERE — this is the deliberate
# downstream step Step 4 was designed to leave open, rather than baking a
# fixed pre/post window into the extraction itself.
ndvi_long <- ndvi_long %>%
  mutate(years_since_loss = year - loss_year)

# -----------------------------------------------------------------------------
# 6. Aggregate NDVI trajectory relative to disturbance (all patches pooled)
#    Mean line + IQR ribbon by years_since_loss. Years with very few patches
#    contributing (e.g. the earliest/latest edges of the record) are dropped
#    so the line isn't driven by a handful of patches.
# -----------------------------------------------------------------------------
ndvi_trajectory <- ndvi_long %>%
  group_by(years_since_loss) %>%
  summarise(
    n_patches = n(),
    ndvi_mean = mean(ndvi_july_mean, na.rm = TRUE),
    ndvi_q25  = quantile(ndvi_july_mean, 0.25, na.rm = TRUE),
    ndvi_q75  = quantile(ndvi_july_mean, 0.75, na.rm = TRUE),
    .groups   = 'drop'
  ) %>%
  filter(n_patches >= 5)

p_trajectory <- ggplot(ndvi_trajectory, aes(x = years_since_loss)) +
  geom_ribbon(aes(ymin = ndvi_q25, ymax = ndvi_q75), fill = 'darkgreen', alpha = 0.2) +
  geom_line(aes(y = ndvi_mean), color = 'darkgreen', linewidth = 1) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  labs(
    title    = 'July NDVI Relative to Year of Disturbance (all patches pooled)',
    subtitle = 'Line = mean across patches, band = IQR. Dashed line = disturbance year.\nPatches ≥ 25 ha only (MODIS-resolvable filter, see Step 4).',
    x        = 'Years since disturbance (0 = loss_year)',
    y        = 'July NDVI'
  ) +
  theme_minimal()

ggsave(file.path(outputDir, 'ndvi_trajectory_relative_to_disturbance.png'), p_trajectory,
       width = 9, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 6b. Sample-size panel + cohort-faceted trajectory
#     The pooled trajectory above hides a real confound: only early-loss-year
#     patches (2001-2003) can contribute data out at years_since_loss = +20,
#     and only late-loss-year patches (2020+) contribute at the negative end
#     -- the tails are different subsets of patches, not the same patches
#     followed longer, so apparent multi-decade "recovery" could just be a
#     cohort artifact (space-for-time substitution bias). These two figures
#     make that visible instead of leaving it implicit in one pooled line.
# -----------------------------------------------------------------------------

# Sample-size sub-panel, stacked under the trajectory so thin-data tails are
# obvious at a glance rather than requiring a separate look at n_patches.
p_n <- ggplot(ndvi_trajectory, aes(x = years_since_loss, y = n_patches)) +
  geom_col(fill = 'grey60') +
  labs(x = 'Years since disturbance (0 = loss_year)', y = 'n patches') +
  theme_minimal()

p_trajectory_with_n <- (p_trajectory + labs(x = NULL) +
                          theme(axis.text.x = element_blank())) /
  p_n +
  plot_layout(heights = c(3, 1))

ggsave(file.path(outputDir, 'ndvi_trajectory_relative_to_disturbance.png'),
       p_trajectory_with_n, width = 9, height = 7.5, dpi = 300)

# Cohort-faceted trajectory: split patches into early/mid/late loss-year
# terciles and plot each cohort's trajectory separately. A real recovery
# signal should look broadly similar across cohorts; sharp divergence
# between them flags the pooled line above as a pooling artifact rather
# than genuine recovery.
#
# Breaks are computed from the loss_year distribution of the patches that
# actually appear in ndvi_long (one row per distinct patch_uuid), NOT from
# the full patch_attrs catalog — ndvi_long already reflects Step 4's >= 25 ha
# filter, and that filtered population's loss_year distribution need not
# match the full catalog's, so terciles drawn from the wrong population
# could produce badly unbalanced (rather than even) cohort groups here.
cohort_breaks <- quantile(
  distinct(ndvi_long, patch_uuid, loss_year)$loss_year,
  probs = c(0, 1/3, 2/3, 1), na.rm = TRUE
)
cohort_breaks <- unique(cohort_breaks)

if (length(cohort_breaks) == 4) {
  ndvi_long <- ndvi_long %>%
    mutate(loss_cohort = cut(loss_year, breaks = cohort_breaks,
                             include.lowest = TRUE,
                             labels = c('early', 'mid', 'late')))
  
  ndvi_trajectory_cohort <- ndvi_long %>%
    group_by(loss_cohort, years_since_loss) %>%
    summarise(
      n_patches = n(),
      ndvi_mean = mean(ndvi_july_mean, na.rm = TRUE),
      .groups   = 'drop'
    ) %>%
    filter(n_patches >= 5)
  
  p_trajectory_cohort <- ggplot(ndvi_trajectory_cohort,
                                aes(x = years_since_loss, y = ndvi_mean, color = loss_cohort)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
    labs(
      title    = 'July NDVI Relative to Disturbance, by Loss-Year Cohort',
      subtitle = 'Cohorts are loss_year terciles -- divergence flags pooling bias.\nPatches ≥ 25 ha only (MODIS-resolvable filter, see Step 4).',
      x        = 'Years since disturbance (0 = loss_year)',
      y        = 'July NDVI',
      color    = 'Loss cohort'
    ) +
    theme_minimal()
  
  ggsave(file.path(outputDir, 'ndvi_trajectory_by_cohort.png'), p_trajectory_cohort,
         width = 9, height = 6, dpi = 300)
} else {
  warning('loss_year quantile breaks were not unique (only ', length(cohort_breaks),
          ' distinct value(s)) -- skipping the cohort-faceted trajectory figure.')
}

# -----------------------------------------------------------------------------
# 7. Per-patch NDVI change (post vs pre), mapped spatially
#    pre  = mean NDVI across all available years BEFORE loss_year
#    post = mean NDVI across years 1-3 AFTER loss_year (early recovery window)
#    These windows are a starting point — adjust once you've looked at the
#    trajectory plot above (e.g. if recovery clearly isn't visible until
#    later years, widen the post window).
#    Patches are required to have at least 3 pre-years and 1 post-year of
#    valid data — anything thinner than that isn't a reliable baseline.
# -----------------------------------------------------------------------------
ndvi_change <- ndvi_long %>%
  group_by(patch_uuid) %>%
  summarise(
    pre_ndvi_mean  = mean(ndvi_july_mean[years_since_loss < 0], na.rm = TRUE),
    post_ndvi_mean = mean(ndvi_july_mean[years_since_loss >= 1 & years_since_loss <= 3], na.rm = TRUE),
    n_pre_years    = sum(years_since_loss < 0),
    n_post_years   = sum(years_since_loss >= 1 & years_since_loss <= 3),
    .groups        = 'drop'
  ) %>%
  mutate(ndvi_change = post_ndvi_mean - pre_ndvi_mean) %>%
  filter(n_pre_years >= 3, n_post_years >= 1)

patches_ndvi_sf <- merge(patches_sf, ndvi_change, by = 'patch_uuid', all.x = FALSE)

p_ndvi_change <- ggplot(patches_ndvi_sf) +
  geom_sf(aes(fill = ndvi_change), color = NA) +
  scale_fill_gradient2(
    low = 'darkred', mid = 'blue', high = 'darkgreen', midpoint = 0,
    # Fixed, tighter limits with out-of-bounds squishing: the raw data range
    # (~ -0.65 to +0.25) is dominated by a handful of outlier patches, which
    # compressed the vast majority of near-zero values into one indistinct
    # color and hid the real spatial signal.
    limits = c(-0.15, 0.15), oob = scales::squish,
    name = 'NDVI change\n(post - pre)'
  ) +
  labs(title    = 'NDVI Change Pre- vs Post-Disturbance, by Patch',
       subtitle = 'Patches ≥ 25 ha only (MODIS-resolvable filter, see Step 4)') +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank())

ggsave(file.path(outputDir, 'map_ndvi_change_by_patch.png'), p_ndvi_change,
       width = 8, height = 7, dpi = 300)

cat(sprintf('NDVI figures saved to %s/:\n', outputDir),
    '  ndvi_trajectory_relative_to_disturbance.png\n',
    '  ndvi_trajectory_by_cohort.png\n',
    '  map_ndvi_change_by_patch.png\n',
    sprintf('(%d of %d patches had enough pre/post NDVI coverage to include in the change map)\n',
            nrow(ndvi_change), n_distinct(ndvi_long$patch_uuid)))

# =============================================================================
# PART D — MODIS Phenology Trends (from Step 4b's phenology_yearly_by_patch.csv)
# =============================================================================
# Mirrors Part C's structure (pooled trajectory / cohort split / per-patch
# change map), but for phenology-transition TIMING rather than NDVI level.
#
# MidGreendown_1 is used as the primary "fall color" metric here, per Step
# 4b's header comment: it's the User Guide's recommended, more stable
# phenophase estimate (50% amplitude crossing, vs. Greenup/Dormancy's 15%
# crossing), and the most representative single peak-transition day for
# this project specifically. Senescence_1 and Dormancy_1 are already in the
# CSV if you want to add matching panels for the onset or near-total-leaf-off
# ends of the transition -- same pattern as below, just swap the column.
#
# Step 4b deliberately ran over the FULL, unfiltered patch catalog (no
# area_ha cutoff), unlike Step 4's NDVI extraction -- but both share the
# same 500m MODIS grid, so the same neighborhood-dilution caveat applies to
# small patches here too (per Step 4b's own header comment). Applying the
# same >=25 ha filter below (via patch_attrs, already loaded in Part B)
# keeps this a patch-level rather than neighborhood-level analysis, for a
# fair comparison against Part C's NDVI figures.

phenology_long <- read.csv(file.path(outputDir, 'phenology_yearly_by_patch.csv'))

phenology_long <- phenology_long %>%
  left_join(patch_attrs[, c('patch_uuid', 'area_ha')], by = 'patch_uuid') %>%
  filter(area_ha >= 25, !is.na(midgreendown_doy_mean))

# Relativize to each patch's own loss_year, same as Part C.
phenology_long <- phenology_long %>%
  mutate(years_since_loss = year - loss_year)

# -----------------------------------------------------------------------------
# 8. Aggregate MidGreendown-DOY trajectory relative to disturbance (pooled)
# -----------------------------------------------------------------------------
phenology_trajectory <- phenology_long %>%
  group_by(years_since_loss) %>%
  summarise(
    n_patches         = n(),
    midgreendown_mean = mean(midgreendown_doy_mean, na.rm = TRUE),
    midgreendown_q25  = quantile(midgreendown_doy_mean, 0.25, na.rm = TRUE),
    midgreendown_q75  = quantile(midgreendown_doy_mean, 0.75, na.rm = TRUE),
    .groups           = 'drop'
  ) %>%
  filter(n_patches >= 5)

p_pheno_trajectory <- ggplot(phenology_trajectory, aes(x = years_since_loss)) +
  geom_ribbon(aes(ymin = midgreendown_q25, ymax = midgreendown_q75),
              fill = 'darkorange', alpha = 0.2) +
  geom_line(aes(y = midgreendown_mean), color = 'darkorange', linewidth = 1) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  labs(
    title    = 'MidGreendown Timing Relative to Year of Disturbance (all patches pooled)',
    subtitle = 'Line = mean across patches, band = IQR. Dashed line = disturbance year.\nPatches >= 25 ha only (MODIS-resolvable filter, matching Part C).',
    x        = 'Years since disturbance (0 = loss_year)',
    y        = 'MidGreendown day-of-year'
  ) +
  theme_minimal()

# Sample-size sub-panel, same rationale as Part C 6b: thin-data tails at the
# edges of the record shouldn't be read as a real signal on their own.
p_pheno_n <- ggplot(phenology_trajectory, aes(x = years_since_loss, y = n_patches)) +
  geom_col(fill = 'grey60') +
  labs(x = 'Years since disturbance (0 = loss_year)', y = 'n patches') +
  theme_minimal()

p_pheno_trajectory_with_n <- (p_pheno_trajectory + labs(x = NULL) +
                                theme(axis.text.x = element_blank())) /
  p_pheno_n +
  plot_layout(heights = c(3, 1))

ggsave(file.path(outputDir, 'phenology_trajectory_relative_to_disturbance.png'),
       p_pheno_trajectory_with_n, width = 9, height = 7.5, dpi = 300)

# -----------------------------------------------------------------------------
# 9. Cohort-faceted phenology trajectory — same space-for-time substitution
#    caveat as Part C 6b applies here (see that section's comment): only
#    early-loss-year patches can contribute at large positive
#    years_since_loss, and only late-loss-year patches at the negative end,
#    so a real fall-color-timing shift should look broadly similar across
#    cohorts; sharp divergence flags the pooled line as a pooling artifact.
# -----------------------------------------------------------------------------
pheno_cohort_breaks <- quantile(
  distinct(phenology_long, patch_uuid, loss_year)$loss_year,
  probs = c(0, 1/3, 2/3, 1), na.rm = TRUE
)
pheno_cohort_breaks <- unique(pheno_cohort_breaks)

if (length(pheno_cohort_breaks) == 4) {
  phenology_long <- phenology_long %>%
    mutate(loss_cohort = cut(loss_year, breaks = pheno_cohort_breaks,
                             include.lowest = TRUE,
                             labels = c('early', 'mid', 'late')))
  
  phenology_trajectory_cohort <- phenology_long %>%
    group_by(loss_cohort, years_since_loss) %>%
    summarise(
      n_patches         = n(),
      midgreendown_mean = mean(midgreendown_doy_mean, na.rm = TRUE),
      .groups           = 'drop'
    ) %>%
    filter(n_patches >= 5)
  
  p_pheno_trajectory_cohort <- ggplot(phenology_trajectory_cohort,
                                      aes(x = years_since_loss, y = midgreendown_mean,
                                          color = loss_cohort)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
    labs(
      title    = 'MidGreendown Timing Relative to Disturbance, by Loss-Year Cohort',
      subtitle = 'Cohorts are loss_year terciles -- divergence flags pooling bias.\nPatches >= 25 ha only (MODIS-resolvable filter, matching Part C).',
      x        = 'Years since disturbance (0 = loss_year)',
      y        = 'MidGreendown day-of-year',
      color    = 'Loss cohort'
    ) +
    theme_minimal()
  
  ggsave(file.path(outputDir, 'phenology_trajectory_by_cohort.png'),
         p_pheno_trajectory_cohort, width = 9, height = 6, dpi = 300)
} else {
  warning('loss_year quantile breaks were not unique (only ', length(pheno_cohort_breaks),
          ' distinct value(s)) -- skipping the cohort-faceted phenology trajectory figure.')
}

# -----------------------------------------------------------------------------
# 10. Per-patch MidGreendown shift (post vs pre), mapped spatially
#     Same pre/post window logic and minimum-coverage thresholds as Part C's
#     NDVI change map (section 7), applied to MidGreendown timing instead.
# -----------------------------------------------------------------------------
phenology_change <- phenology_long %>%
  group_by(patch_uuid) %>%
  summarise(
    pre_midgreendown_mean  = mean(midgreendown_doy_mean[years_since_loss < 0], na.rm = TRUE),
    post_midgreendown_mean = mean(midgreendown_doy_mean[years_since_loss >= 1 & years_since_loss <= 3], na.rm = TRUE),
    n_pre_years            = sum(years_since_loss < 0),
    n_post_years           = sum(years_since_loss >= 1 & years_since_loss <= 3),
    .groups                = 'drop'
  ) %>%
  mutate(midgreendown_shift = post_midgreendown_mean - pre_midgreendown_mean) %>%
  filter(n_pre_years >= 3, n_post_years >= 1)

patches_pheno_sf <- merge(patches_sf, phenology_change, by = 'patch_uuid', all.x = FALSE)

p_pheno_change <- ggplot(patches_pheno_sf) +
  geom_sf(aes(fill = midgreendown_shift), color = NA) +
  scale_fill_gradient2(
    low = 'purple', mid = 'white', high = 'darkorange', midpoint = 0,
    name = 'MidGreendown shift\n(post - pre, days)'
  ) +
  labs(title    = 'MidGreendown Timing Shift Pre- vs Post-Disturbance, by Patch',
       subtitle = 'Patches >= 25 ha only (MODIS-resolvable filter, matching Part C)') +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank())

ggsave(file.path(outputDir, 'map_phenology_change_by_patch.png'), p_pheno_change,
       width = 8, height = 7, dpi = 300)

cat(sprintf('Phenology figures saved to %s/:\n', outputDir),
    '  phenology_trajectory_relative_to_disturbance.png\n',
    '  phenology_trajectory_by_cohort.png\n',
    '  map_phenology_change_by_patch.png\n',
    sprintf('(%d of %d patches had enough pre/post MidGreendown coverage to include in the change map)\n',
            nrow(phenology_change), n_distinct(phenology_long$patch_uuid)))

# =============================================================================
# PART E — Pre/Post Histograms + Baseline-Normalized Trajectories
# =============================================================================
# Two things bundled together, both built from ndvi_change/phenology_change
# (per-patch pre/post summaries, Parts C/D) and ndvi_long/phenology_long
# (per-patch-per-year series, Parts C/D):
#
#   1. Histograms comparing the PRE vs POST distribution across patches, for
#      both NDVI and MidGreendown timing.
#   2. Baseline-normalized trajectories: every patch's own pre-disturbance
#      mean is rescaled to a common reference point (1 for NDVI, 0 for
#      MidGreendown -- see note below on why these differ), so patches with
#      very different absolute starting values become comparable on one
#      relative scale.
#
# NOTE on the NDVI-vs-MidGreendown normalization difference: NDVI is a
# ratio-scale quantity (true zero = no vegetation), so dividing by the
# pre-disturbance mean is meaningful -- a patch at 0.5x its baseline really
# did lose half its relative greenness. MidGreendown is a day-of-year value
# on an INTERVAL scale (Jan 1 is an arbitrary reference point, not a
# meaningful zero), so dividing one DOY by another doesn't have a clean
# physical interpretation. The equivalent, well-defined operation for an
# interval-scale variable is SUBTRACTION, not division: expressing each
# year's MidGreendown as an ANOMALY IN DAYS relative to the patch's own
# pre-disturbance mean, with the pre-disturbance baseline set to 0 rather
# than 1. Both are provided below so you can compare, but the anomaly-in-days
# version is the one we'd recommend leaning on.

# -----------------------------------------------------------------------------
# 11. Pre/Post histograms — overlaid distributions across patches
# -----------------------------------------------------------------------------
ndvi_prepost_long <- ndvi_change %>%
  select(patch_uuid, pre_ndvi_mean, post_ndvi_mean) %>%
  pivot_longer(cols = c(pre_ndvi_mean, post_ndvi_mean),
               names_to = 'period', values_to = 'ndvi') %>%
  mutate(period = recode(period, pre_ndvi_mean = 'Pre-disturbance',
                         post_ndvi_mean = 'Post-disturbance (yr 1-3)'))

p_ndvi_hist <- ggplot(ndvi_prepost_long, aes(x = ndvi, fill = period)) +
  geom_histogram(position = 'identity', alpha = 0.5, bins = 30, color = 'white') +
  scale_fill_manual(values = c('Pre-disturbance' = 'grey40',
                               'Post-disturbance (yr 1-3)' = 'darkgreen'), name = NULL) +
  labs(title = 'Distribution of Per-Patch July NDVI, Pre vs Post Disturbance',
       subtitle = 'Patches >= 25 ha only (MODIS-resolvable filter, see Step 4).',
       x = 'July NDVI (patch mean)', y = 'Number of patches') +
  theme_minimal()

ggsave(file.path(outputDir, 'hist_ndvi_pre_post.png'), p_ndvi_hist,
       width = 8, height = 5, dpi = 300)

pheno_prepost_long <- phenology_change %>%
  select(patch_uuid, pre_midgreendown_mean, post_midgreendown_mean) %>%
  pivot_longer(cols = c(pre_midgreendown_mean, post_midgreendown_mean),
               names_to = 'period', values_to = 'midgreendown_doy') %>%
  mutate(period = recode(period, pre_midgreendown_mean = 'Pre-disturbance',
                         post_midgreendown_mean = 'Post-disturbance (yr 1-3)'))

p_pheno_hist <- ggplot(pheno_prepost_long, aes(x = midgreendown_doy, fill = period)) +
  geom_histogram(position = 'identity', alpha = 0.5, bins = 30, color = 'white') +
  scale_fill_manual(values = c('Pre-disturbance' = 'grey40',
                               'Post-disturbance (yr 1-3)' = 'darkorange'), name = NULL) +
  labs(title = 'Distribution of Per-Patch MidGreendown Timing, Pre vs Post Disturbance',
       subtitle = 'Patches >= 25 ha only (MODIS-resolvable filter, matching Part C).',
       x = 'MidGreendown day-of-year (patch mean)', y = 'Number of patches') +
  theme_minimal()

ggsave(file.path(outputDir, 'hist_midgreendown_pre_post.png'), p_pheno_hist,
       width = 8, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 12. Baseline-normalized trajectories
# -----------------------------------------------------------------------------
ndvi_long_norm <- ndvi_long %>%
  left_join(ndvi_change %>% select(patch_uuid, pre_ndvi_mean), by = 'patch_uuid') %>%
  filter(!is.na(pre_ndvi_mean), pre_ndvi_mean != 0) %>%
  mutate(ndvi_relative = ndvi_july_mean / pre_ndvi_mean)

ndvi_trajectory_norm <- ndvi_long_norm %>%
  group_by(years_since_loss) %>%
  summarise(
    n_patches     = n(),
    ndvi_rel_mean = mean(ndvi_relative, na.rm = TRUE),
    ndvi_rel_q25  = quantile(ndvi_relative, 0.25, na.rm = TRUE),
    ndvi_rel_q75  = quantile(ndvi_relative, 0.75, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(n_patches >= 5)

p_ndvi_trajectory_norm <- ggplot(ndvi_trajectory_norm, aes(x = years_since_loss)) +
  geom_ribbon(aes(ymin = ndvi_rel_q25, ymax = ndvi_rel_q75), fill = 'darkgreen', alpha = 0.2) +
  geom_line(aes(y = ndvi_rel_mean), color = 'darkgreen', linewidth = 1) +
  geom_hline(yintercept = 1, linetype = 'dotted', color = 'grey30') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  labs(
    title    = "NDVI Relative to Each Patch's Own Pre-Disturbance Baseline",
    subtitle = "Baseline (dotted line at 1.0) = each patch's own pre-loss mean NDVI.\nPatches >= 25 ha only.",
    x = 'Years since disturbance (0 = loss_year)',
    y = 'NDVI relative to pre-disturbance baseline'
  ) +
  theme_minimal()

ggsave(file.path(outputDir, 'ndvi_trajectory_normalized.png'), p_ndvi_trajectory_norm,
       width = 9, height = 6, dpi = 300)

phenology_long_norm <- phenology_long %>%
  left_join(phenology_change %>% select(patch_uuid, pre_midgreendown_mean), by = 'patch_uuid') %>%
  filter(!is.na(pre_midgreendown_mean)) %>%
  mutate(
    midgreendown_anomaly_days = midgreendown_doy_mean - pre_midgreendown_mean,  # recommended
    midgreendown_relative     = midgreendown_doy_mean / pre_midgreendown_mean   # literal ratio -- see caveat above
  )

phenology_trajectory_norm <- phenology_long_norm %>%
  group_by(years_since_loss) %>%
  summarise(
    n_patches    = n(),
    anomaly_mean = mean(midgreendown_anomaly_days, na.rm = TRUE),
    anomaly_q25  = quantile(midgreendown_anomaly_days, 0.25, na.rm = TRUE),
    anomaly_q75  = quantile(midgreendown_anomaly_days, 0.75, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(n_patches >= 5)

p_pheno_trajectory_norm <- ggplot(phenology_trajectory_norm, aes(x = years_since_loss)) +
  geom_ribbon(aes(ymin = anomaly_q25, ymax = anomaly_q75), fill = 'darkorange', alpha = 0.2) +
  geom_line(aes(y = anomaly_mean), color = 'darkorange', linewidth = 1) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'grey30') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  labs(
    title    = "MidGreendown Timing Relative to Each Patch's Own Pre-Disturbance Baseline",
    subtitle = "Baseline (dotted line at 0) = each patch's own pre-loss mean MidGreendown DOY.\nPatches >= 25 ha only.",
    x = 'Years since disturbance (0 = loss_year)',
    y = 'MidGreendown anomaly (days relative to pre-disturbance baseline)'
  ) +
  theme_minimal()

ggsave(file.path(outputDir, 'phenology_trajectory_normalized.png'), p_pheno_trajectory_norm,
       width = 9, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 13. Single-distribution histograms of the change itself
# -----------------------------------------------------------------------------
ndvi_change_norm <- ndvi_change %>%
  mutate(ndvi_ratio = post_ndvi_mean / pre_ndvi_mean)

p_ndvi_ratio_hist <- ggplot(ndvi_change_norm, aes(x = ndvi_ratio)) +
  geom_histogram(bins = 30, fill = 'darkgreen', alpha = 0.7, color = 'white') +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'grey30') +
  labs(title = 'Per-Patch NDVI Change, Relative to Own Pre-Disturbance Baseline',
       subtitle = 'Dashed line at 1.0 = no change from baseline.\nPatches >= 25 ha only.',
       x = 'Post / Pre NDVI ratio', y = 'Number of patches') +
  theme_minimal()

ggsave(file.path(outputDir, 'hist_ndvi_ratio.png'), p_ndvi_ratio_hist,
       width = 8, height = 5, dpi = 300)

p_pheno_anomaly_hist <- ggplot(phenology_change, aes(x = midgreendown_shift)) +
  geom_histogram(bins = 30, fill = 'darkorange', alpha = 0.7, color = 'white') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey30') +
  labs(title = 'Per-Patch MidGreendown Shift, Relative to Own Pre-Disturbance Baseline',
       subtitle = 'Dashed line at 0 = no change from baseline.\nPatches >= 25 ha only.',
       x = 'MidGreendown shift (days, post - pre)', y = 'Number of patches') +
  theme_minimal()

ggsave(file.path(outputDir, 'hist_midgreendown_shift.png'), p_pheno_anomaly_hist,
       width = 8, height = 5, dpi = 300)

cat(sprintf('Part E figures saved to %s/:\n', outputDir),
    '  hist_ndvi_pre_post.png\n',
    '  hist_midgreendown_pre_post.png\n',
    '  ndvi_trajectory_normalized.png\n',
    '  phenology_trajectory_normalized.png\n',
    '  hist_ndvi_ratio.png\n',
    '  hist_midgreendown_shift.png\n')

# =============================================================================
# PART F — Deeper Senescence-Timing Analyses
# =============================================================================
# F1: regional (climate) reference comparison -- needs Step 4c's output
# F2: patch's own pre-disturbance TREND (not just flat mean) as baseline
# F3: senescence duration/rate (Senescence_1 -> Dormancy_1), not just a date
# F4: interannual variability, pre vs post (does timing get less predictable?)
# F5: QA-stratified check (is the shift a data-quality artifact?)
# =============================================================================

# -----------------------------------------------------------------------------
# F1. Regional reference comparison (climate-controlled)
#     Requires reference_phenology_by_year.csv from Step 4c. Computes each
#     patch-year's anomaly relative to that YEAR's regional reference mean,
#     then compares each patch's own pre- vs post-disturbance anomaly level
#     (a difference-in-differences estimate): this nets out BOTH (a) the
#     region's common year-to-year climate variation, since we subtract the
#     matching year's regional mean first, and (b) each patch's own
#     idiosyncratic offset from the regional average (e.g. due to elevation,
#     aspect, forest type), since we compare the patch's own anomaly level
#     before vs after -- not the anomaly's absolute value.
# -----------------------------------------------------------------------------
reference_phenology <- read.csv(file.path(outputDir, 'reference_phenology_by_year.csv'))

phenology_long_ref <- phenology_long %>%
  left_join(reference_phenology %>% select(year, ref_midgreendown_mean), by = 'year') %>%
  filter(!is.na(ref_midgreendown_mean)) %>%
  mutate(regional_anomaly = midgreendown_doy_mean - ref_midgreendown_mean)

regional_anomaly_trajectory <- phenology_long_ref %>%
  group_by(years_since_loss) %>%
  summarise(
    n_patches   = n(),
    anomaly_mean = mean(regional_anomaly, na.rm = TRUE),
    anomaly_q25  = quantile(regional_anomaly, 0.25, na.rm = TRUE),
    anomaly_q75  = quantile(regional_anomaly, 0.75, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(n_patches >= 5)

p_regional_anomaly <- ggplot(regional_anomaly_trajectory, aes(x = years_since_loss)) +
  geom_ribbon(aes(ymin = anomaly_q25, ymax = anomaly_q75), fill = 'steelblue', alpha = 0.2) +
  geom_line(aes(y = anomaly_mean), color = 'steelblue', linewidth = 1) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  labs(
    title    = 'MidGreendown Timing Relative to the Regional (Undisturbed-Forest) Baseline',
    subtitle = paste0('Each point = patch MidGreendown DOY minus that YEAR\'s regional reference mean\n',
                      '(Step 4c, ', reference_phenology$year[1], '-', tail(reference_phenology$year, 1),
                      ', n ref cells varies by year). This nets out common year-to-year climate\n',
                      'variability; a level SHIFT around loss_year (not the absolute level, which\n',
                      'reflects each patch\'s own baseline offset from the region) indicates a real\n',
                      'disturbance effect. Patches >= 25 ha only.'),
    x = 'Years since disturbance (0 = loss_year)',
    y = 'MidGreendown timing relative to regional reference (days)'
  ) +
  theme_minimal()

ggsave(file.path(outputDir, 'phenology_regional_anomaly_trajectory.png'), p_regional_anomaly,
       width = 9, height = 6.5, dpi = 300)

# Difference-in-differences per patch: change in the patch's OWN anomaly
# level, pre vs post -- this is the actual disturbance-attributable shift,
# controlling for both regional climate and the patch's own baseline offset.
regional_anomaly_change <- phenology_long_ref %>%
  group_by(patch_uuid) %>%
  summarise(
    pre_anomaly_mean  = mean(regional_anomaly[years_since_loss < 0], na.rm = TRUE),
    post_anomaly_mean = mean(regional_anomaly[years_since_loss >= 1 & years_since_loss <= 3], na.rm = TRUE),
    n_pre_years  = sum(years_since_loss < 0),
    n_post_years = sum(years_since_loss >= 1 & years_since_loss <= 3),
    .groups = 'drop'
  ) %>%
  mutate(regional_anomaly_shift = post_anomaly_mean - pre_anomaly_mean) %>%
  filter(n_pre_years >= 3, n_post_years >= 1)

p_regional_anomaly_hist <- ggplot(regional_anomaly_change, aes(x = regional_anomaly_shift)) +
  geom_histogram(bins = 30, fill = 'steelblue', alpha = 0.7, color = 'white') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey30') +
  labs(title = 'Climate-Controlled MidGreendown Shift, Per Patch',
       subtitle = paste0('Change in each patch\'s own gap from the regional reference, post vs pre\n',
                         'disturbance -- controls for both regional climate variability and each\n',
                         'patch\'s own baseline offset from the region. Dashed line at 0 = no shift\n',
                         'beyond regional climate. Patches >= 25 ha only.'),
       x = 'Regional-anomaly shift (days, post - pre)', y = 'Number of patches') +
  theme_minimal()

ggsave(file.path(outputDir, 'hist_phenology_regional_anomaly_shift.png'), p_regional_anomaly_hist,
       width = 8, height = 5.5, dpi = 300)

# -----------------------------------------------------------------------------
# F2. Patch's own pre-disturbance TREND as baseline (not just a flat mean)
#     If senescence timing has been drifting across the record for climate
#     reasons, a flat pre-mean baseline misattributes that drift to
#     disturbance. Fit year vs DOY on each patch's PRE-loss years only, then
#     compare post-loss values to the trend's own extrapolation.
#     Requires >=5 pre-loss years (a flat-mean baseline elsewhere in this
#     project only requires >=3 -- a slope needs more points to be stable).
# -----------------------------------------------------------------------------
phenology_trend_baseline <- phenology_long %>%
  filter(years_since_loss < 0, !is.na(midgreendown_doy_mean)) %>%
  group_by(patch_uuid) %>%
  filter(n() >= 5) %>%
  group_modify(~ {
    fit <- lm(midgreendown_doy_mean ~ year, data = .x)
    tibble(intercept = unname(coef(fit)[1]), slope = unname(coef(fit)[2]))
  }) %>%
  ungroup()

phenology_long_trend <- phenology_long %>%
  inner_join(phenology_trend_baseline, by = 'patch_uuid') %>%
  filter(!is.na(midgreendown_doy_mean)) %>%
  mutate(
    expected_midgreendown = intercept + slope * year,
    trend_residual        = midgreendown_doy_mean - expected_midgreendown
  )

trend_residual_trajectory <- phenology_long_trend %>%
  group_by(years_since_loss) %>%
  summarise(
    n_patches     = n(),
    residual_mean = mean(trend_residual, na.rm = TRUE),
    residual_q25  = quantile(trend_residual, 0.25, na.rm = TRUE),
    residual_q75  = quantile(trend_residual, 0.75, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(n_patches >= 5)

p_trend_residual <- ggplot(trend_residual_trajectory, aes(x = years_since_loss)) +
  geom_ribbon(aes(ymin = residual_q25, ymax = residual_q75), fill = 'darkred', alpha = 0.2) +
  geom_line(aes(y = residual_mean), color = 'darkred', linewidth = 1) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'grey30') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  labs(
    title    = "MidGreendown Deviation from Each Patch's Own Pre-Disturbance TREND",
    subtitle = paste0('0 = matches the linear trend extrapolated from that patch\'s pre-loss years\n',
                      '(not just its flat pre-loss mean) -- controls for any pre-existing drift.\n',
                      'Requires >=5 pre-loss years per patch to fit a trend; patches with fewer\n',
                      'are excluded, so this uses a smaller patch set than Part D\'s figures.'),
    x = 'Years since disturbance (0 = loss_year)',
    y = 'Residual (days, observed - trend-expected)'
  ) +
  theme_minimal()

ggsave(file.path(outputDir, 'phenology_trend_residual_trajectory.png'), p_trend_residual,
       width = 9, height = 6.5, dpi = 300)

# -----------------------------------------------------------------------------
# F3. Senescence duration/rate: does disturbance change HOW FAST a patch
#     senesces, not just when it starts? Senescence_1 (onset, 90% crossing)
#     to Dormancy_1 (near-total leaf-off, 15% crossing) already exist in
#     phenology_long from Step 4b.
#     Sanity bounds (0 < duration < 200 days) guard against edge cases where
#     a transition date is missing or wraps oddly across the year boundary.
# -----------------------------------------------------------------------------
phenology_long <- phenology_long %>%
  mutate(senescence_duration_days = dormancy_doy_mean - senescence_doy_mean)

duration_valid <- phenology_long %>%
  filter(!is.na(senescence_duration_days),
         senescence_duration_days > 0, senescence_duration_days < 200)

duration_trajectory <- duration_valid %>%
  group_by(years_since_loss) %>%
  summarise(
    n_patches     = n(),
    duration_mean = mean(senescence_duration_days, na.rm = TRUE),
    duration_q25  = quantile(senescence_duration_days, 0.25, na.rm = TRUE),
    duration_q75  = quantile(senescence_duration_days, 0.75, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(n_patches >= 5)

p_duration_trajectory <- ggplot(duration_trajectory, aes(x = years_since_loss)) +
  geom_ribbon(aes(ymin = duration_q25, ymax = duration_q75), fill = 'chocolate4', alpha = 0.2) +
  geom_line(aes(y = duration_mean), color = 'chocolate4', linewidth = 1) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  labs(
    title    = 'Senescence Duration (Onset to Dormancy) Relative to Disturbance',
    subtitle = 'Senescence_1 to Dormancy_1, in days -- shorter = faster senescence.\nPatches >= 25 ha only; duration filtered to a plausible 0-200 day range.',
    x = 'Years since disturbance (0 = loss_year)',
    y = 'Senescence duration (days)'
  ) +
  theme_minimal()

ggsave(file.path(outputDir, 'senescence_duration_trajectory.png'), p_duration_trajectory,
       width = 9, height = 6, dpi = 300)

duration_change <- duration_valid %>%
  group_by(patch_uuid) %>%
  summarise(
    pre_duration_mean  = mean(senescence_duration_days[years_since_loss < 0], na.rm = TRUE),
    post_duration_mean = mean(senescence_duration_days[years_since_loss >= 1 & years_since_loss <= 3], na.rm = TRUE),
    n_pre_years  = sum(years_since_loss < 0),
    n_post_years = sum(years_since_loss >= 1 & years_since_loss <= 3),
    .groups = 'drop'
  ) %>%
  mutate(duration_change = post_duration_mean - pre_duration_mean) %>%
  filter(n_pre_years >= 3, n_post_years >= 1)

p_duration_hist <- ggplot(duration_change, aes(x = duration_change)) +
  geom_histogram(bins = 30, fill = 'chocolate4', alpha = 0.7, color = 'white') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey30') +
  labs(title = 'Change in Senescence Duration, Per Patch',
       subtitle = 'Negative = senescence became FASTER after disturbance.\nPatches >= 25 ha only.',
       x = 'Duration change (days, post - pre)', y = 'Number of patches') +
  theme_minimal()

ggsave(file.path(outputDir, 'hist_senescence_duration_change.png'), p_duration_hist,
       width = 8, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# F4. Variability, pre vs post
#     (a) INTERANNUAL variability: how much does the patch's own YEARLY mean
#         MidGreendown bounce around from year to year, pre vs post? This is
#         a genuinely different question from (b) below.
#     (b) WITHIN-YEAR spatial variability: midgreendown_doy_stdDev already
#         exists per patch-year from Step 4b's reducer (pixel-to-pixel spread
#         inside that year's 500m cell) -- a bonus comparison since it's
#         already in hand, distinct from (a).
#     Both are stdDev ratios (a ratio-scale quantity: true zero = no
#     variability), so post/pre RATIO normalization is appropriate here,
#     unlike the raw MidGreendown DOY itself (see Part E's note).
# -----------------------------------------------------------------------------
# MATCHED WINDOWS (-3:-1 vs 1:3), both required complete.
#
# The earlier version compared an UNBOUNDED pre-window (all pre-loss years,
# often 10-14 of them) against a 3-year post-window, and required only >=4
# pre / >=2 post. That guarantees a downward-biased ratio even when nothing
# real changes, because the sample SD is systematically low at small n:
# E[s] ~= 0.886*sigma at n=3 versus ~0.981*sigma at n=14. Expected ratio
# with ZERO true change ~= 0.90 -- essentially exactly where the old
# histogram peaked. The apparent "timing became more predictable after
# disturbance" was that artifact, not a finding.
#
# Equal window lengths make the bias identical in numerator and denominator
# so it cancels. Requiring BOTH windows complete (3 of 3 each) means every
# patch contributes n=3 vs n=3 -- no differential bias by construction.
interannual_variability <- phenology_long %>%
  group_by(patch_uuid) %>%
  summarise(
    pre_sd  = sd(midgreendown_doy_mean[years_since_loss %in% -3:-1], na.rm = TRUE),
    post_sd = sd(midgreendown_doy_mean[years_since_loss %in% 1:3],   na.rm = TRUE),
    n_pre_years  = sum(years_since_loss %in% -3:-1 & !is.na(midgreendown_doy_mean)),
    n_post_years = sum(years_since_loss %in% 1:3   & !is.na(midgreendown_doy_mean)),
    .groups = 'drop'
  ) %>%
  filter(n_pre_years == 3, n_post_years == 3, pre_sd > 0) %>%
  mutate(sd_ratio     = post_sd / pre_sd,
         log2_ratio   = log2(sd_ratio))

# Summarise on the LOG2 ratio and use the MEDIAN, not the mean.
# With n=3 per window (df=2), the variance ratio follows F(2,2) under the
# null, whose mean is undefined and whose right tail is extremely heavy --
# a mean ratio would be dominated by a few patches with a near-zero pre_sd
# denominator. log2 makes the null distribution symmetric about 0, and the
# median is well-defined (median of F(2,2) = 1, i.e. log2 = 0).
median_ratio <- median(interannual_variability$sd_ratio, na.rm = TRUE)

# Formal test: is the log-ratio distribution centred away from 0? Wilcoxon
# signed-rank rather than a t-test, since even in log space the tails stay
# heavy at this n.
wilcox_var <- wilcox.test(interannual_variability$log2_ratio, mu = 0)

cat(sprintf('\nInterannual variability (matched 3yr windows, n = %d patches):\n',
            nrow(interannual_variability)))
cat(sprintf('  Median post/pre SD ratio: %.3f\n', median_ratio))
cat(sprintf('  Wilcoxon signed-rank on log2 ratio vs 0: p = %.4g\n', wilcox_var$p.value))
if (wilcox_var$p.value >= 0.05) {
  cat('  -> No detectable change in year-to-year predictability of timing.\n')
} else if (median_ratio > 1) {
  cat('  -> Timing became LESS predictable year-to-year after disturbance.\n')
} else {
  cat('  -> Timing became MORE predictable year-to-year after disturbance.\n')
}

p_interannual_sd_hist <- ggplot(interannual_variability, aes(x = sd_ratio)) +
  geom_histogram(bins = 30, fill = 'purple', alpha = 0.7, color = 'white') +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'grey30') +
  geom_vline(xintercept = median_ratio, linetype = 'solid', color = 'darkred', linewidth = 0.8) +
  scale_x_continuous(trans = 'log2',
                     breaks = c(0.25, 0.5, 1, 2, 4),
                     labels = c('0.25', '0.5', '1', '2', '4')) +
  labs(title = 'Interannual Variability of MidGreendown Timing, Post vs Pre Disturbance',
       subtitle = paste0('Ratio of post-loss to pre-loss year-to-year stdDev in each patch\'s own\n',
                         'yearly mean MidGreendown DOY, using MATCHED 3-year windows (-3:-1 vs\n',
                         '1:3, both complete) so small-sample SD bias cancels rather than\n',
                         'manufacturing a shift. Log2 x-axis: the null distribution is symmetric\n',
                         'in log space, not linear. Dashed = no change; solid red = observed median.\n',
                         sprintf('Median = %.3f, Wilcoxon p = %.3g, n = %d patches.',
                                 median_ratio, wilcox_var$p.value,
                                 nrow(interannual_variability))),
       x = 'Post-loss SD / Pre-loss SD (interannual, log2 scale)',
       y = 'Number of patches') +
  theme_minimal()

ggsave(file.path(outputDir, 'hist_midgreendown_interannual_variability.png'), p_interannual_sd_hist,
       width = 8, height = 6, dpi = 300)

# Bonus (b): within-year spatial heterogeneity, pre vs post -- reuses the
# midgreendown_doy_stdDev column already computed by Step 4b's reducer.
spatial_heterogeneity_change <- phenology_long %>%
  filter(!is.na(midgreendown_doy_stdDev), midgreendown_doy_stdDev > 0) %>%
  group_by(patch_uuid) %>%
  summarise(
    pre_spatial_sd  = mean(midgreendown_doy_stdDev[years_since_loss < 0], na.rm = TRUE),
    post_spatial_sd = mean(midgreendown_doy_stdDev[years_since_loss >= 1 & years_since_loss <= 3], na.rm = TRUE),
    n_pre_years  = sum(years_since_loss < 0),
    n_post_years = sum(years_since_loss >= 1 & years_since_loss <= 3),
    .groups = 'drop'
  ) %>%
  filter(n_pre_years >= 3, n_post_years >= 1, pre_spatial_sd > 0) %>%
  mutate(spatial_sd_ratio = post_spatial_sd / pre_spatial_sd)

p_spatial_sd_hist <- ggplot(spatial_heterogeneity_change, aes(x = spatial_sd_ratio)) +
  geom_histogram(bins = 30, fill = 'orchid4', alpha = 0.7, color = 'white') +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'grey30') +
  labs(title = 'Within-Patch Spatial Heterogeneity of MidGreendown Timing, Post vs Pre',
       subtitle = paste0('Ratio of post-loss to pre-loss mean pixel-to-pixel stdDev WITHIN each\n',
                         'patch\'s own 500m cell(s) in a given year (not year-to-year -- that\'s the\n',
                         'panel above). Values > 1 = the patch became more spatially non-uniform\n',
                         'in timing after disturbance. Patches >= 25 ha only.'),
       x = 'Post-loss spatial SD / Pre-loss spatial SD', y = 'Number of patches') +
  theme_minimal()

ggsave(file.path(outputDir, 'hist_midgreendown_spatial_heterogeneity.png'), p_spatial_sd_hist,
       width = 8, height = 5.5, dpi = 300)

# -----------------------------------------------------------------------------
# F5. QA-stratified check: is the apparent shift a data-quality artifact?
#     qa_overall_mean already exists per patch-year from Step 4b
#     (QA_Overall_1: 0=best ... 3=poor). Two checks:
#       (i)  does retrieval quality itself change pre vs post disturbance?
#       (ii) does a bigger apparent MidGreendown shift co-occur with a
#            bigger QA degradation across patches?
# -----------------------------------------------------------------------------
qa_change <- phenology_long %>%
  filter(!is.na(qa_overall_mean)) %>%
  group_by(patch_uuid) %>%
  summarise(
    pre_qa_mean  = mean(qa_overall_mean[years_since_loss < 0], na.rm = TRUE),
    post_qa_mean = mean(qa_overall_mean[years_since_loss >= 1 & years_since_loss <= 3], na.rm = TRUE),
    n_pre_years  = sum(years_since_loss < 0),
    n_post_years = sum(years_since_loss >= 1 & years_since_loss <= 3),
    .groups = 'drop'
  ) %>%
  filter(n_pre_years >= 3, n_post_years >= 1) %>%
  mutate(qa_change = post_qa_mean - pre_qa_mean)

p_qa_change_hist <- ggplot(qa_change, aes(x = qa_change)) +
  geom_histogram(bins = 30, fill = 'grey40', alpha = 0.7, color = 'white') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  labs(title = 'Change in MCD12Q2 Retrieval Quality, Per Patch',
       subtitle = paste0('QA_Overall_1: 0=best ... 3=poor. Positive = quality got WORSE after\n',
                         'disturbance. Systematic degradation here would mean some of the\n',
                         'apparent phenology shift could be a data-quality artifact.\n',
                         'Patches >= 25 ha only.'),
       x = 'QA change (post - pre; higher = worse quality)', y = 'Number of patches') +
  theme_minimal()

ggsave(file.path(outputDir, 'hist_phenology_qa_change.png'), p_qa_change_hist,
       width = 8, height = 5.5, dpi = 300)

qa_vs_shift <- phenology_change %>%
  select(patch_uuid, midgreendown_shift) %>%
  inner_join(qa_change %>% select(patch_uuid, qa_change), by = 'patch_uuid')

p_qa_vs_shift <- ggplot(qa_vs_shift, aes(x = qa_change, y = midgreendown_shift)) +
  geom_point(alpha = 0.4, color = 'darkorange') +
  geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'grey50') +
  geom_vline(xintercept = 0, linetype = 'dotted', color = 'grey50') +
  labs(title = 'Does Retrieval Quality Explain the Apparent MidGreendown Shift?',
       subtitle = paste0('Each point = one patch. A positive relationship here would suggest some of\n',
                         'the apparent timing shift tracks data-quality degradation rather than real\n',
                         'phenology change. Patches >= 25 ha only.'),
       x = 'Change in QA_Overall_1 (post - pre; higher = worse quality)',
       y = 'MidGreendown shift (days, post - pre)') +
  theme_minimal()

ggsave(file.path(outputDir, 'scatter_qa_vs_phenology_shift.png'), p_qa_vs_shift,
       width = 8, height = 6, dpi = 300)

cat(sprintf('Part F figures saved to %s/:\n', outputDir),
    '  phenology_regional_anomaly_trajectory.png\n',
    '  hist_phenology_regional_anomaly_shift.png\n',
    '  phenology_trend_residual_trajectory.png\n',
    '  senescence_duration_trajectory.png\n',
    '  hist_senescence_duration_change.png\n',
    '  hist_midgreendown_interannual_variability.png\n',
    '  hist_midgreendown_spatial_heterogeneity.png\n',
    '  hist_phenology_qa_change.png\n',
    '  scatter_qa_vs_phenology_shift.png\n')

# =============================================================================
# PART G — Mixed-Effects Model: Patch + Year Random Effects
# =============================================================================
# CONCEPTUAL NOTE: this consolidates two things that were previously built
# separately (Part E's baseline-normalized trajectories, and Part F #2's
# regional-reference comparison against undisturbed forest pixels pulled via
# Step 4c) into a single mixed-effects model:
#
#   response ~ post + (1 | patch_uuid) + (1 | year)
#
# A per-patch random intercept absorbs each patch's own baseline level --
# the same job Part E's manual rescaling (NDVI ratio -> 1, MidGreendown
# anomaly -> 0) was doing by hand. A per-year random intercept absorbs
# common year-to-year climate variation shared by every patch that year --
# the same job Step 4c's synthetic undisturbed-reference-pixel pull was
# doing. Because patches were disturbed in staggered years (not all at
# once), any given calendar year contains a mix of patches that are
# currently "pre" and others that are currently "post" -- this is exactly
# the kind of staggered-treatment-timing structure that lets a year random
# effect legitimately separate "this year ran warm/cool for everyone" from
# "this patch was disturbed." That means Step 4c's regional-reference pull
# is no longer strictly necessary for this purpose -- it's left in the
# project as an independent cross-check if you want one, but this model is
# the primary, more rigorous route to the same question.
#
# ONE TERMINOLOGY CORRECTION vs. how this was originally described: patch
# and year are CROSSED random effects here, (1|patch_uuid) + (1|year), NOT
# nested (1|year/patch_uuid). Nesting is for group IDs that are only
# meaningful within their parent (e.g. "Plot 3" means something different
# in every different year). patch_uuid is the opposite -- it's the SAME
# physical patch tracked across many years, which is the textbook case for
# crossed, not nested, random effects.
#
# Requires: lme4, lmerTest (for Satterthwaite-approximation p-values on the
# fixed effect, which base lme4 does not provide on its own).
# =============================================================================
library(lme4)
library(lmerTest)

# -----------------------------------------------------------------------------
# G0a. How many patches were disturbed each year?
#      Uses the FULL patch catalog (patch_attrs), not just the >=25ha
#      MODIS-resolvable subset used in the NDVI/phenology models below --
#      this is about the disturbance record itself, not the modeling subset.
# -----------------------------------------------------------------------------
p_loss_year_hist <- ggplot(patch_attrs, aes(x = loss_year)) +
  geom_histogram(binwidth = 1, boundary = 0.5, fill = 'firebrick', alpha = 0.8, color = 'white') +
  labs(title = 'Number of Persistent-Loss Patches by Disturbance Year',
       subtitle = 'Full patch catalog (all sizes, not just the >=25 ha modeling subset below).',
       x = 'Loss year', y = 'Number of patches') +
  theme_minimal()

ggsave(file.path(outputDir, 'hist_patches_per_loss_year.png'), p_loss_year_hist,
       width = 9, height = 5.5, dpi = 300)

# -----------------------------------------------------------------------------
# G0b. Patch area distribution
#      Log10 x-axis: patch area is typically heavily right-skewed (many
#      small patches, a few very large ones), so a linear axis would
#      compress almost everything into the first bar.
# -----------------------------------------------------------------------------
p_area_hist <- ggplot(patch_attrs, aes(x = area_ha)) +
  geom_histogram(bins = 30, fill = 'darkslategray', alpha = 0.8, color = 'white') +
  scale_x_log10(labels = scales::label_number()) +
  labs(title = 'Distribution of Patch Area',
       subtitle = 'Full patch catalog. Log10 x-axis -- patch area is heavily right-skewed.',
       x = 'Patch area (ha, log scale)', y = 'Number of patches') +
  theme_minimal()

ggsave(file.path(outputDir, 'hist_patch_area.png'), p_area_hist,
       width = 8, height = 5.5, dpi = 300)

# -----------------------------------------------------------------------------
# G1. Prepare model data
#     "post" = years_since_loss >= 1. The transition year itself
#     (years_since_loss == 0, i.e. loss_year) is excluded from both pre and
#     post -- it's a partial-year, ambiguous observation (disturbance could
#     have happened before or after that year's MODIS composite date),
#     matching how pre/post windows were already handled elsewhere in this
#     project. Unlike Part C/D/E/F's pre/post SUMMARY tables (which used a
#     fixed years_since_loss 1-3 "post" window), this model uses ALL
#     available post-loss years, not just 1-3 -- the point of the model is
#     to use the full panel and let the patch/year random effects handle
#     the heterogeneity, rather than pre-restricting the window by hand.
# -----------------------------------------------------------------------------
model_data_ndvi <- ndvi_long %>%
  filter(years_since_loss != 0) %>%
  mutate(post = years_since_loss >= 1)

model_data_pheno <- phenology_long %>%
  filter(years_since_loss != 0, !is.na(midgreendown_doy_mean)) %>%
  mutate(post = years_since_loss >= 1)

# -----------------------------------------------------------------------------
# G2. Generic fit-and-summarize helper
#     Returns one row per subset: the "post" fixed-effect estimate (the
#     disturbance effect, controlling for patch baseline + year climate),
#     its standard error, p-value, and a flag for statistical significance.
#     Subsets too thin to support both random effects (fewer than 2 distinct
#     patches or fewer than 2 distinct years) are reported as such rather
#     than silently attempted.
# -----------------------------------------------------------------------------
fit_mixed_effect <- function(data, response_var, label) {
  n_patches <- n_distinct(data$patch_uuid)
  n_years   <- n_distinct(data$year)
  
  if (n_patches < 2 || n_years < 2 || nrow(data) < 10) {
    return(tibble(subset = label, response = response_var,
                  n_patches = n_patches, n_obs = nrow(data),
                  estimate = NA_real_, std_error = NA_real_, p_value = NA_real_,
                  significant = NA, note = 'Too few patches/years/observations to fit'))
  }
  
  form <- as.formula(paste0(response_var, ' ~ post + (1 | patch_uuid) + (1 | year)'))
  fit <- tryCatch(lmerTest::lmer(form, data = data), error = function(e) NULL,
                  warning = function(w) { message(label, ' (', response_var, '): ', conditionMessage(w)); NULL })
  
  if (is.null(fit)) {
    return(tibble(subset = label, response = response_var,
                  n_patches = n_patches, n_obs = nrow(data),
                  estimate = NA_real_, std_error = NA_real_, p_value = NA_real_,
                  significant = NA, note = 'Model failed to fit/converge'))
  }
  
  co <- summary(fit)$coefficients
  if (!'postTRUE' %in% rownames(co)) {
    return(tibble(subset = label, response = response_var,
                  n_patches = n_patches, n_obs = nrow(data),
                  estimate = NA_real_, std_error = NA_real_, p_value = NA_real_,
                  significant = NA, note = 'No post-disturbance contrast estimable (e.g. all pre or all post)'))
  }
  
  tibble(subset = label, response = response_var,
         n_patches = n_patches, n_obs = nrow(data),
         estimate = co['postTRUE', 'Estimate'],
         std_error = co['postTRUE', 'Std. Error'],
         p_value = co['postTRUE', 'Pr(>|t|)'],
         significant = co['postTRUE', 'Pr(>|t|)'] < 0.05,
         note = NA_character_)
}

# -----------------------------------------------------------------------------
# G3. Build the four subsets
# -----------------------------------------------------------------------------
set.seed(42)  # reproducible random 10% subset

# (a) Full dataset -- nothing to build, use model_data_ndvi / model_data_pheno directly.

# (b) Random 10% of patches
patches_all <- unique(model_data_ndvi$patch_uuid)
patches_10pct <- sample(patches_all, size = ceiling(0.10 * length(patches_all)))

# (c) 10 largest patches by area
patches_top10_area <- patch_attrs %>%
  arrange(desc(area_ha)) %>%
  slice_head(n = 10) %>%
  pull(patch_uuid)
# Caveat: a patch random effect with only 10 groups is thin (rule-of-thumb
# guidance generally wants >=5-10 groups minimum) -- treat this subset's
# variance-component estimates as indicative, not precise, even where the
# fixed-effect p-value comes out significant.

# (d) Top 20 patches by absolute NDVI change (post - pre) -- "greatest
#     change" is read here as largest magnitude in either direction, not
#     just largest decline. 20 rather than 10 for a bit more stability in
#     the random-effects estimates; change to 10 if you want strict parity
#     with the "10 largest patches" subset above.
patches_top20_change <- ndvi_change %>%
  mutate(abs_change = abs(ndvi_change)) %>%
  arrange(desc(abs_change)) %>%
  slice_head(n = 20) %>%
  pull(patch_uuid)

subset_defs <- list(
  'All patches'                       = patches_all,
  'Random 10% of patches'             = patches_10pct,
  '10 largest patches (by area)'      = patches_top10_area,
  'Top 20 patches (|NDVI change|)'    = patches_top20_change
)

# -----------------------------------------------------------------------------
# G4. Fit both response variables (NDVI, MidGreendown) across all four subsets
# -----------------------------------------------------------------------------
model_results <- purrr::imap_dfr(subset_defs, function(patch_ids, label) {
  bind_rows(
    fit_mixed_effect(filter(model_data_ndvi, patch_uuid %in% patch_ids), 'ndvi_july_mean', label),
    fit_mixed_effect(filter(model_data_pheno, patch_uuid %in% patch_ids), 'midgreendown_doy_mean', label)
  )
})

# Full results table (ALL subsets/responses, significant or not) -- saved so
# nothing is silently hidden, even though the figures below only surface the
# significant ones per your request.
write.csv(model_results, file.path(outputDir, 'mixed_model_results_all.csv'), row.names = FALSE)

# -----------------------------------------------------------------------------
# G5. Significant-only view, per your request: "only show things where we
#     can show a statistical impact on NDVI." Applied to both response
#     variables for consistency, but NDVI is the one explicitly asked for.
# -----------------------------------------------------------------------------
model_results_significant <- model_results %>%
  filter(significant == TRUE)

cat('Mixed-effects model results (all subsets, both responses) saved to:\n',
    file.path(outputDir, 'mixed_model_results_all.csv'), '\n\n')
cat('Subsets with a statistically significant (p < 0.05) disturbance effect:\n')
print(model_results_significant)

if (nrow(model_results_significant %>% filter(response == 'ndvi_july_mean')) == 0) {
  cat('\nNo subset showed a statistically significant disturbance effect on NDVI.\n',
      'This itself is informative -- check mixed_model_results_all.csv for the\n',
      'full estimates/p-values/notes (e.g. convergence issues, thin subsets)\n',
      'rather than assuming a null result means nothing happened.\n')
}

# -----------------------------------------------------------------------------
# G6. Effect-size comparison plots (one per response variable), showing the
#     "post" estimate +/- 95% CI for each subset. Only significant subsets
#     are plotted, per your request -- the full table above is the place to
#     check what was excluded and why.
# -----------------------------------------------------------------------------
plot_effect_sizes <- function(results_df, response_label, y_lab, fill_color, filename) {
  df <- results_df %>%
    filter(response == response_label, significant == TRUE) %>%
    mutate(ci_lo = estimate - 1.96 * std_error,
           ci_hi = estimate + 1.96 * std_error,
           subset = factor(subset, levels = names(subset_defs)))
  
  if (nrow(df) == 0) {
    cat(sprintf('\nSkipping effect-size plot for %s -- no significant subsets to show.\n',
                response_label))
    return(invisible(NULL))
  }
  
  p <- ggplot(df, aes(x = subset, y = estimate)) +
    geom_hline(yintercept = 0, linetype = 'dotted', color = 'grey40') +
    geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), color = fill_color, linewidth = 0.8, size = 0.7) +
    coord_flip() +
    labs(title = paste0('Disturbance Effect on ', y_lab, ' (Mixed-Effects Model)'),
         subtitle = 'Only subsets with a statistically significant (p < 0.05) effect are shown.\nPoint = estimate, bars = 95% CI. See mixed_model_results_all.csv for every subset.',
         x = NULL, y = paste0(y_lab, ' change (post vs pre, controlling for patch + year)')) +
    theme_minimal()
  
  ggsave(file.path(outputDir, filename), p, width = 9, height = 5, dpi = 300)
}

plot_effect_sizes(model_results, 'ndvi_july_mean', 'NDVI', 'darkgreen', 'effect_size_ndvi.png')
plot_effect_sizes(model_results, 'midgreendown_doy_mean', 'MidGreendown (days)', 'darkorange', 'effect_size_midgreendown.png')

cat(sprintf('\nPart G figures saved to %s/:\n', outputDir),
    '  hist_patches_per_loss_year.png\n',
    '  hist_patch_area.png\n',
    '  mixed_model_results_all.csv\n',
    '  effect_size_ndvi.png (only if any subset was significant)\n',
    '  effect_size_midgreendown.png (only if any subset was significant)\n')

# =============================================================================
# PART H — Does Patch Size Moderate the Disturbance Effect? (NDVI + MidGreendown)
# =============================================================================
# Tests whether the disturbance effect scales with patch area, using a
# continuous interaction model on the FULL qualifying population (>=25 ha)
# rather than the ad hoc subset comparisons in Part G:
#
#   response ~ post * log(area_ha) + (1 | patch_uuid) + (1 | year)
#
# The postTRUE:log(area_ha) coefficient is the effect of interest: does the
# disturbance effect (postTRUE) get bigger or smaller as patch area grows?
# log(area_ha) is used because patch area is heavily right-skewed (see the
# area histogram in Part G).
#
# IMPORTANT: (Intercept) and postTRUE alone are evaluated at log(area_ha)=0,
# i.e. a hypothetical 1 ha patch -- outside this project's actual data range
# (patches run 25-~280 ha). Only the COMBINED effect at a realistic area is
# meaningful: effect_at_area = postTRUE + postTRUE:log(area_ha) * log(area)
#
# For EACH response variable, this section:
#   1. Fits the full interaction model on all qualifying patches.
#   2. Saves a residual-vs-fitted diagnostic plot.
#   3. Identifies high-leverage patches -- the METHOD DIFFERS between the two
#      responses, deliberately, because their residual plots showed
#      different problems (see each subsection).
#   4. Refits EXCLUDING those patches, as a robustness check.
#   5. Reports both fits side by side in one CSV, so a fragile result
#      (significance that depends on a handful of patches) is visible
#      rather than hidden.
#
# HEADLINE FINDING: NDVI's size-interaction is robust -- it survives
# excluding the 16 most volatile patches (2.6% of the data), shrinking from
# -0.0297 to -0.0155 but staying significant at p<2e-16 either way.
# MidGreendown's is NOT robust -- excluding just ONE ordinary, non-QA-
# flagged patch (out of 615) is enough to flip p from 0.042 to non-
# significant. That single patch showed no data-quality problem of any
# kind; it's simply an large, complete-record patch whose ordinary
# statistical leverage happens to be large enough to decide the result.
# The honest conclusion is that patch size moderates the disturbance
# effect for NDVI, but this is NOT demonstrated for MidGreendown timing.
#
# Requires: influence.ME (for Cook's distance on a patch-level grouping
# factor in a mixed model -- base R's cooks.distance() doesn't support
# grouped/random-effect structures).
# =============================================================================
library(influence.ME)

# -----------------------------------------------------------------------------
# H1. NDVI: fit full model, diagnose residuals
# -----------------------------------------------------------------------------
model_data_ndvi_area <- model_data_ndvi %>%
  filter(!is.na(area_ha))

fit_ndvi_area <- lmerTest::lmer(
  ndvi_july_mean ~ post * log(area_ha) + (1 | patch_uuid) + (1 | year),
  data = model_data_ndvi_area
)

ndvi_diag <- model_data_ndvi_area %>%
  mutate(fit_val = fitted(fit_ndvi_area),
         resid_val = resid(fit_ndvi_area, type = 'pearson'))

p_resid_ndvi <- ggplot(ndvi_diag, aes(x = fit_val, y = resid_val)) +
  geom_point(alpha = 0.15, color = 'steelblue', shape = 1) +
  geom_hline(yintercept = 0, color = 'black') +
  labs(title = 'Residuals vs Fitted: NDVI Size-Interaction Model',
       subtitle = paste0('A visually distinct low-fitted-value cluster (fitted < 0.75) shows much wider\n',
                         'spread (+/- 0.3-0.4) than the main cloud -- flagged and excluded below.'),
       x = 'Fitted value', y = 'Pearson residual') +
  theme_minimal()

ggsave(file.path(outputDir, 'residuals_ndvi_size_interaction.png'), p_resid_ndvi,
       width = 7.5, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# H2. NDVI: identify the high-variance cluster, refit without it
#     Method: fitted < 0.75. This is a visual/exploratory threshold (there
#     is a genuinely distinct second population below it in the residual
#     plot), not a formal statistical rule -- Cook's distance is used for
#     MidGreendown instead below, since that residual plot didn't show an
#     equivalent clean visual break, just individual high-leverage points.
# -----------------------------------------------------------------------------
ndvi_low_fit_patches <- ndvi_diag %>%
  filter(fit_val < 0.75) %>%
  distinct(patch_uuid) %>%
  pull(patch_uuid)

model_data_ndvi_clean <- model_data_ndvi_area %>%
  filter(!patch_uuid %in% ndvi_low_fit_patches)

fit_ndvi_area_clean <- lmerTest::lmer(
  ndvi_july_mean ~ post * log(area_ha) + (1 | patch_uuid) + (1 | year),
  data = model_data_ndvi_clean
)

# -----------------------------------------------------------------------------
# H3. MidGreendown: fit full model, diagnose residuals
# -----------------------------------------------------------------------------
model_data_pheno_area <- model_data_pheno %>%
  filter(!is.na(area_ha))

fit_pheno_area <- lmerTest::lmer(
  midgreendown_doy_mean ~ post * log(area_ha) + (1 | patch_uuid) + (1 | year),
  data = model_data_pheno_area
)

pheno_diag <- model_data_pheno_area %>%
  mutate(fit_val = fitted(fit_pheno_area),
         resid_val = resid(fit_pheno_area, type = 'pearson'))

p_resid_pheno <- ggplot(pheno_diag, aes(x = fit_val, y = resid_val)) +
  geom_point(alpha = 0.15, color = 'darkorange', shape = 1) +
  geom_hline(yintercept = 0, color = 'black') +
  labs(title = 'Residuals vs Fitted: MidGreendown Size-Interaction Model',
       subtitle = paste0('No visually distinct second cluster (unlike NDVI) -- a single continuous\n',
                         'cloud with a mild funnel. Concern here is individual high-leverage points,\n',
                         'identified via Cook\'s distance below, not a subgroup.'),
       x = 'Fitted value', y = 'Pearson residual') +
  theme_minimal()

ggsave(file.path(outputDir, 'residuals_midgreendown_size_interaction.png'), p_resid_pheno,
       width = 7.5, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# H4. MidGreendown: Cook's distance per patch, refit without top-influence
#     Method: Cook's distance via influence.ME. Top 2 patches are excluded
#     here -- this specific number was validated directly: excluding just
#     these 2 (0.3% of patches) is enough to flip the postTRUE:log(area_ha)
#     p-value from 0.042 to 0.16. Both flagged patches were checked
#     individually: the #1 patch (highest Cook's distance) turned out to
#     have a perfect qa_overall_mean (0, best possible) throughout its
#     record -- no data-quality justification for exclusion on its own,
#     and removing ONLY that patch barely moves the result (p: 0.042 ->
#     0.0485). It's specifically the #2 patch -- an entirely ordinary,
#     complete-record, non-anomalous large patch -- whose presence/absence
#     is what actually decides significance. That is itself the finding:
#     this result is fragile, not contaminated by a data error.
# -----------------------------------------------------------------------------
infl_pheno <- influence(fit_pheno_area, group = 'patch_uuid')
cooks_d_pheno <- cooks.distance(infl_pheno)
cooks_df_pheno <- data.frame(patch_uuid = rownames(cooks_d_pheno),
                             cooks_d = as.numeric(cooks_d_pheno)) %>%
  arrange(desc(cooks_d))

write.csv(cooks_df_pheno, file.path(outputDir, 'cooks_distance_midgreendown_by_patch.csv'),
          row.names = FALSE)

pheno_top_influence_patches <- head(cooks_df_pheno$patch_uuid, 2)

model_data_pheno_clean <- model_data_pheno_area %>%
  filter(!patch_uuid %in% pheno_top_influence_patches)

fit_pheno_area_clean <- lmerTest::lmer(
  midgreendown_doy_mean ~ post * log(area_ha) + (1 | patch_uuid) + (1 | year),
  data = model_data_pheno_clean
)

# -----------------------------------------------------------------------------
# H5. Consolidated comparison table: full vs clean, both responses
# -----------------------------------------------------------------------------
extract_interaction_row <- function(fit, response_label, model_version,
                                    n_patches, n_excluded, excluded_patches) {
  if (is.null(fit)) {
    return(tibble(response = response_label, model_version = model_version,
                  n_patches = n_patches, n_patches_excluded = n_excluded,
                  excluded_patch_uuids = paste(excluded_patches, collapse = '; '),
                  interaction_estimate = NA_real_, interaction_se = NA_real_,
                  interaction_p_value = NA_real_, interaction_significant = NA))
  }
  co <- summary(fit)$coefficients
  tibble(response = response_label, model_version = model_version,
         n_patches = n_patches, n_patches_excluded = n_excluded,
         excluded_patch_uuids = paste(excluded_patches, collapse = '; '),
         interaction_estimate = co['postTRUE:log(area_ha)', 'Estimate'],
         interaction_se = co['postTRUE:log(area_ha)', 'Std. Error'],
         interaction_p_value = co['postTRUE:log(area_ha)', 'Pr(>|t|)'],
         interaction_significant = co['postTRUE:log(area_ha)', 'Pr(>|t|)'] < 0.05)
}

size_interaction_results <- bind_rows(
  extract_interaction_row(fit_ndvi_area, 'ndvi_july_mean', 'Full (all qualifying patches)',
                          n_distinct(model_data_ndvi_area$patch_uuid), 0, character(0)),
  extract_interaction_row(fit_ndvi_area_clean, 'ndvi_july_mean',
                          'Excl. low-fitted-value cluster (fitted < 0.75)',
                          n_distinct(model_data_ndvi_clean$patch_uuid),
                          length(ndvi_low_fit_patches), ndvi_low_fit_patches),
  extract_interaction_row(fit_pheno_area, 'midgreendown_doy_mean', 'Full (all qualifying patches)',
                          n_distinct(model_data_pheno_area$patch_uuid), 0, character(0)),
  extract_interaction_row(fit_pheno_area_clean, 'midgreendown_doy_mean',
                          "Excl. top-2 Cook's distance patches",
                          n_distinct(model_data_pheno_clean$patch_uuid),
                          length(pheno_top_influence_patches), pheno_top_influence_patches)
)

write.csv(size_interaction_results, file.path(outputDir, 'size_interaction_robustness_results.csv'),
          row.names = FALSE)

cat('Part H results saved to', outputDir, ':\n',
    '  residuals_ndvi_size_interaction.png\n',
    '  residuals_midgreendown_size_interaction.png\n',
    '  cooks_distance_midgreendown_by_patch.csv\n',
    '  size_interaction_robustness_results.csv\n\n')
print(size_interaction_results)

# =============================================================================
# PART I — Year-by-Year Disturbance Effect (-5 to +5), Not Just Pre/Post
# =============================================================================
# Part G's "post" variable was a single TRUE/FALSE indicator, pooling ALL
# post-loss years together. That can hide a real effect that only emerges
# a few years out -- if, say, MidGreendown doesn't shift in year+1 but does
# by year+3, lumping them into one flat "post" average dilutes the year+3
# signal with the null year+1 signal, and could be part of why some of
# Part G/H's results were weaker or more fragile than expected.
#
# This section replaces the binary "post" with a categorical FACTOR:
#   - "pre"    = years_since_loss -5 through -1, pooled into one reference
#                level (matching the project's existing convention of a
#                multi-year pre-disturbance baseline elsewhere)
#   - "post_1" through "post_5" = each individual year after disturbance,
#                kept SEPARATE rather than pooled, so a delayed or
#                strengthening effect is visible year by year instead of
#                averaged away.
# Years beyond +/-5 and the transition year itself (years_since_loss == 0)
# are excluded, matching this specific -5..+5 window.
#
# Three complementary views of the same restructured data, per your
# request:
#   1. A "straight ANOVA" -- both a simple one-way aov() ignoring the
#      patch/year structure (a classic means comparison across periods,
#      with Tukey HSD for pairwise comparisons), AND the equivalent overall
#      F-test from the proper mixed model (anova() on the lmer fit), so you
#      can see whether accounting for patch/year random effects changes the
#      simple picture.
#   2. The mixed model itself, with "period" as a factor -- gives one
#      coefficient (and p-value) per individual post-year, relative to the
#      pooled pre-baseline.
#   3. A superimposed epoch analysis (SEA) plot -- the classic event-study
#      technique of aligning every patch's series on relative time and
#      plotting the mean +/- 95% CI at each relative-time step. This is
#      conceptually the same idea as Parts C/D's pooled trajectories, but
#      restricted to this -5..+5 window and paired directly with the
#      factor model's significance results (points are marked when that
#      year is statistically different from pooled pre-disturbance).
# =============================================================================

# -----------------------------------------------------------------------------
# I1. Build the period factor for both response variables
# -----------------------------------------------------------------------------
build_period_factor <- function(data) {
  data %>%
    filter(years_since_loss %in% c(-5:-1, 1:5)) %>%
    mutate(period = case_when(
      years_since_loss %in% -5:-1 ~ 'pre',
      years_since_loss == 1        ~ 'post_1',
      years_since_loss == 2        ~ 'post_2',
      years_since_loss == 3        ~ 'post_3',
      years_since_loss == 4        ~ 'post_4',
      years_since_loss == 5        ~ 'post_5'
    )) %>%
    mutate(period = factor(period,
                           levels = c('pre', 'post_1', 'post_2', 'post_3', 'post_4', 'post_5')))
}

model_data_ndvi_period  <- build_period_factor(model_data_ndvi)
model_data_pheno_period <- build_period_factor(model_data_pheno)

# -----------------------------------------------------------------------------
# I2. Mixed model with "period" as a factor (one coefficient per post-year,
#     vs. pooled pre-disturbance) + the overall ANOVA F-test on that factor
# -----------------------------------------------------------------------------
fit_ndvi_period <- lmerTest::lmer(
  ndvi_july_mean ~ period + (1 | patch_uuid) + (1 | year),
  data = model_data_ndvi_period
)
fit_pheno_period <- lmerTest::lmer(
  midgreendown_doy_mean ~ period + (1 | patch_uuid) + (1 | year),
  data = model_data_pheno_period
)

anova_ndvi_period_mixed  <- anova(fit_ndvi_period)
anova_pheno_period_mixed <- anova(fit_pheno_period)

cat('Mixed-model overall F-test for "period" factor -- NDVI:\n')
print(anova_ndvi_period_mixed)
cat('\nMixed-model overall F-test for "period" factor -- MidGreendown:\n')
print(anova_pheno_period_mixed)

# -----------------------------------------------------------------------------
# I3. "Straight" one-way ANOVA (ignoring patch/year random effects) + Tukey
#     HSD pairwise comparisons -- the simpler, classical means-test version
# -----------------------------------------------------------------------------
aov_ndvi_period <- aov(ndvi_july_mean ~ period, data = model_data_ndvi_period)
aov_pheno_period <- aov(midgreendown_doy_mean ~ period, data = model_data_pheno_period)

cat('\n\nSimple one-way ANOVA (no random effects) -- NDVI:\n')
print(summary(aov_ndvi_period))
cat('\nTukey HSD pairwise comparisons -- NDVI:\n')
print(TukeyHSD(aov_ndvi_period))

cat('\n\nSimple one-way ANOVA (no random effects) -- MidGreendown:\n')
print(summary(aov_pheno_period))
cat('\nTukey HSD pairwise comparisons -- MidGreendown:\n')
print(TukeyHSD(aov_pheno_period))

# Save the Tukey tables to CSV for reference
write.csv(as.data.frame(TukeyHSD(aov_ndvi_period)$period),
          file.path(outputDir, 'tukey_hsd_ndvi_period.csv'))
write.csv(as.data.frame(TukeyHSD(aov_pheno_period)$period),
          file.path(outputDir, 'tukey_hsd_midgreendown_period.csv'))

# -----------------------------------------------------------------------------
# I4. Tidy per-period coefficients from the MIXED model (the version that
#     controls for patch + year), for both plotting and CSV export
# -----------------------------------------------------------------------------
extract_period_coefs <- function(fit, response_label) {
  co <- summary(fit)$coefficients
  rows <- rownames(co)[rownames(co) != '(Intercept)']
  tibble(
    response    = response_label,
    period      = gsub('^period', '', rows),
    estimate    = co[rows, 'Estimate'],
    std_error   = co[rows, 'Std. Error'],
    p_value     = co[rows, 'Pr(>|t|)'],
    significant = co[rows, 'Pr(>|t|)'] < 0.05
  )
}

period_results <- bind_rows(
  extract_period_coefs(fit_ndvi_period, 'ndvi_july_mean'),
  extract_period_coefs(fit_pheno_period, 'midgreendown_doy_mean')
)

write.csv(period_results, file.path(outputDir, 'period_factor_model_results.csv'), row.names = FALSE)
cat('\n\nPer-period estimates (vs. pooled pre-disturbance baseline), mixed model:\n')
print(period_results)

# -----------------------------------------------------------------------------
# I5. Forest-plot of per-year estimates -- directly shows whether the effect
#     is present immediately, or builds/strengthens over the post-loss years
# -----------------------------------------------------------------------------
plot_period_effects <- function(results_df, response_label, y_lab, filename) {
  df <- results_df %>%
    filter(response == response_label) %>%
    mutate(period = factor(period, levels = c('post_1', 'post_2', 'post_3', 'post_4', 'post_5')))
  
  p <- ggplot(df, aes(x = period, y = estimate)) +
    geom_hline(yintercept = 0, linetype = 'dotted', color = 'grey40') +
    geom_pointrange(aes(ymin = estimate - 1.96 * std_error, ymax = estimate + 1.96 * std_error,
                        color = significant), linewidth = 0.8, size = 0.7) +
    scale_color_manual(values = c(`TRUE` = 'darkred', `FALSE` = 'grey60'), name = 'p < 0.05') +
    labs(title = paste0(y_lab, ' Effect by Year Since Disturbance'),
         subtitle = 'Each post-year is its own factor level, compared against a pooled pre-disturbance\nbaseline (years -5 to -1) -- not lumped into one flat "post" average.',
         x = 'Year since disturbance', y = paste0(y_lab, ' difference from pre-disturbance baseline')) +
    theme_minimal()
  
  ggsave(file.path(outputDir, filename), p, width = 8, height = 5.5, dpi = 300)
}

plot_period_effects(period_results, 'ndvi_july_mean', 'NDVI', 'period_effects_ndvi.png')
plot_period_effects(period_results, 'midgreendown_doy_mean', 'MidGreendown', 'period_effects_midgreendown.png')

# -----------------------------------------------------------------------------
# I6. Superimposed epoch analysis (SEA): mean +/- 95% CI at each relative
#     year, -5 to +5, with points marked where that year is significantly
#     different from the pooled pre-baseline (per the factor model above).
# -----------------------------------------------------------------------------
compute_sea <- function(data, response_col) {
  data %>%
    group_by(years_since_loss) %>%
    summarise(
      n        = n(),
      mean_val = mean(.data[[response_col]], na.rm = TRUE),
      se       = sd(.data[[response_col]], na.rm = TRUE) / sqrt(n()),
      ci_lo    = mean_val - 1.96 * se,
      ci_hi    = mean_val + 1.96 * se,
      .groups  = 'drop'
    )
}

sea_ndvi  <- compute_sea(model_data_ndvi_period, 'ndvi_july_mean')
sea_pheno <- compute_sea(model_data_pheno_period, 'midgreendown_doy_mean')

sig_lookup <- period_results %>%
  mutate(years_since_loss = case_when(
    period == 'post_1' ~ 1, period == 'post_2' ~ 2, period == 'post_3' ~ 3,
    period == 'post_4' ~ 4, period == 'post_5' ~ 5, TRUE ~ NA_real_
  )) %>%
  filter(!is.na(years_since_loss))

sea_ndvi <- sea_ndvi %>%
  left_join(sig_lookup %>% filter(response == 'ndvi_july_mean') %>% select(years_since_loss, significant),
            by = 'years_since_loss') %>%
  mutate(significant = ifelse(is.na(significant), FALSE, significant))  # pre-years: not tested against themselves

sea_pheno <- sea_pheno %>%
  left_join(sig_lookup %>% filter(response == 'midgreendown_doy_mean') %>% select(years_since_loss, significant),
            by = 'years_since_loss') %>%
  mutate(significant = ifelse(is.na(significant), FALSE, significant))

plot_sea <- function(sea_df, y_lab, fill_color, filename) {
  p <- ggplot(sea_df, aes(x = years_since_loss, y = mean_val)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), fill = fill_color, alpha = 0.2) +
    geom_line(color = fill_color, linewidth = 1) +
    geom_point(aes(color = significant), size = 3) +
    scale_color_manual(values = c(`TRUE` = 'red', `FALSE` = 'grey30'), name = 'Sig. vs pre (p<0.05)') +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
    labs(title = paste0('Superimposed Epoch Analysis: ', y_lab, ' by Year Since Disturbance'),
         subtitle = 'Points = yearly mean +/- 95% CI. Pre-disturbance years (-5 to -1) are the pooled\nreference period; red points are significantly different from that pooled baseline.',
         x = 'Years since disturbance (0 excluded; negative = pre, positive = post)',
         y = y_lab) +
    theme_minimal()
  
  ggsave(file.path(outputDir, filename), p, width = 9, height = 6, dpi = 300)
}

plot_sea(sea_ndvi, 'July NDVI', 'darkgreen', 'sea_ndvi.png')
plot_sea(sea_pheno, 'MidGreendown day-of-year', 'darkorange', 'sea_midgreendown.png')

cat('\n\nPart I figures/tables saved to', outputDir, ':\n',
    '  period_effects_ndvi.png\n',
    '  period_effects_midgreendown.png\n',
    '  sea_ndvi.png\n',
    '  sea_midgreendown.png\n',
    '  period_factor_model_results.csv\n',
    '  tukey_hsd_ndvi_period.csv\n',
    '  tukey_hsd_midgreendown_period.csv\n')

# =============================================================================
# PART J — Does Severity of NDVI Loss Moderate the Effect? (vs. Part H's Area)
# =============================================================================
# Part H asked whether patch AREA moderates the disturbance effect. This
# section asks the same style of question using a different moderator:
# the MAGNITUDE OF NDVI LOSS each patch actually experienced (ndvi_change,
# from Part C section 7 -- post_ndvi_mean minus pre_ndvi_mean per patch).
#
# IMPORTANT ASYMMETRY, handled deliberately differently for the two
# responses:
#
#   MidGreendown: this is a clean, non-circular test. ndvi_change and
#   midgreendown_doy_mean are two independently measured quantities, so
#   asking "do patches with bigger NDVI declines also show bigger
#   senescence-timing shifts" is a genuine severity-moderation question --
#   built below as a continuous interaction model, the same style as
#   Part H's post*log(area_ha).
#
#   NDVI: using ndvi_change to moderate ndvi_july_mean itself would use the
#   outcome to explain the outcome -- a patch's post/pre NDVI difference
#   is arithmetically DERIVED from the same ndvi_july_mean values the model
#   is trying to explain, so grouping patches by their own NDVI change and
#   then asking "do groups differ in NDVI" is close to guaranteed to look
#   large by construction, not evidence of anything new. This is built
#   below as a labeled, clearly-flagged DESCRIPTIVE severity-tercile check
#   (mild/moderate/severe groups), not presented as equivalent evidence to
#   the MidGreendown version.
# =============================================================================

# -----------------------------------------------------------------------------
# J1. MidGreendown ~ post * NDVI-loss-severity (continuous, non-circular)
# -----------------------------------------------------------------------------
model_data_pheno_severity <- model_data_pheno %>%
  left_join(ndvi_change %>% select(patch_uuid, ndvi_change), by = 'patch_uuid') %>%
  filter(!is.na(ndvi_change)) %>%
  mutate(abs_ndvi_change = abs(ndvi_change))
# abs_ndvi_change used as the primary moderator (rather than signed
# ndvi_change) purely for interpretability: a NEGATIVE interaction
# coefficient on abs_ndvi_change then directly reads as "bigger NDVI loss
# -> bigger (earlier) MidGreendown shift," without having to mentally
# flip a sign. The signed version is left in the data if you want it.

fit_pheno_severity <- lmerTest::lmer(
  midgreendown_doy_mean ~ post * abs_ndvi_change + (1 | patch_uuid) + (1 | year),
  data = model_data_pheno_severity
)

cat('MidGreendown ~ post * |NDVI change| (severity moderator, continuous):\n')
print(summary(fit_pheno_severity))

# Direct patch-level scatter as a companion to the model -- same visual
# style as Part F's QA-vs-shift check, using the pre-computed per-patch
# summary tables (phenology_change, ndvi_change) rather than the full panel.
severity_vs_shift <- phenology_change %>%
  select(patch_uuid, midgreendown_shift) %>%
  inner_join(ndvi_change %>% select(patch_uuid, ndvi_change), by = 'patch_uuid')

p_severity_vs_shift <- ggplot(severity_vs_shift, aes(x = ndvi_change, y = midgreendown_shift)) +
  geom_point(alpha = 0.4, color = 'darkgreen') +
  geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = 'dotted', color = 'grey50') +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'grey50') +
  labs(title = 'MidGreendown Shift vs. NDVI Change, Per Patch',
       subtitle = paste0('Do patches with bigger NDVI declines (more negative, further left) also show\n',
                         'bigger senescence-timing shifts (more negative, further down)? Independent\n',
                         'measurements on each axis -- not a circular comparison.'),
       x = 'NDVI change (post - pre)', y = 'MidGreendown shift (days, post - pre)') +
  theme_minimal()

ggsave(file.path(outputDir, 'scatter_ndvi_severity_vs_midgreendown_shift.png'), p_severity_vs_shift,
       width = 8, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# J2. NDVI ~ post * severity-tercile -- DESCRIPTIVE ONLY, see header caveat.
#     Terciles are built from the same ndvi_change used to define severity,
#     which is itself derived from the response variable being modeled --
#     expect this to show large, close-to-definitional differences between
#     groups. Useful as a sanity check that the grouping behaves as
#     expected, not as independent confirmation of a severity effect.
# -----------------------------------------------------------------------------
ndvi_severity_terciles <- ndvi_change %>%
  mutate(abs_ndvi_change = abs(ndvi_change),
         severity_tercile = ntile(abs_ndvi_change, 3)) %>%
  mutate(severity_tercile = factor(severity_tercile, labels = c('mild', 'moderate', 'severe')))

model_data_ndvi_severity <- model_data_ndvi %>%
  left_join(ndvi_severity_terciles %>% select(patch_uuid, severity_tercile), by = 'patch_uuid') %>%
  filter(!is.na(severity_tercile))

fit_ndvi_severity <- lmerTest::lmer(
  ndvi_july_mean ~ post * severity_tercile + (1 | patch_uuid) + (1 | year),
  data = model_data_ndvi_severity
)

cat('\n\nNDVI ~ post * severity tercile (DESCRIPTIVE ONLY -- see Part J header, this is\n',
    'expected to look large by construction since terciles are built from the\n',
    'same variable being modeled):\n')
print(summary(fit_ndvi_severity))
print(anova(fit_ndvi_severity))

# -----------------------------------------------------------------------------
# J3. Consolidated results CSV -- both checks, clearly labeled by whether
#     each is a genuine independent test or a descriptive/circular check
# -----------------------------------------------------------------------------
extract_key_rows <- function(fit, response_label, check_type, term_pattern) {
  co <- summary(fit)$coefficients
  rows <- rownames(co)[grepl(term_pattern, rownames(co))]
  tibble(
    response = response_label,
    check_type = check_type,
    term = rows,
    estimate = co[rows, 'Estimate'],
    std_error = co[rows, 'Std. Error'],
    p_value = co[rows, 'Pr(>|t|)'],
    significant = co[rows, 'Pr(>|t|)'] < 0.05
  )
}

severity_results <- bind_rows(
  extract_key_rows(fit_pheno_severity, 'midgreendown_doy_mean',
                   'Independent test (non-circular)', 'post.*abs_ndvi_change'),
  extract_key_rows(fit_ndvi_severity, 'ndvi_july_mean',
                   'DESCRIPTIVE ONLY (circular -- severity defined from same response)',
                   'post.*severity_tercile')
)

write.csv(severity_results, file.path(outputDir, 'ndvi_severity_moderation_results.csv'), row.names = FALSE)

cat('\n\nPart J results saved to', outputDir, ':\n',
    '  scatter_ndvi_severity_vs_midgreendown_shift.png\n',
    '  ndvi_severity_moderation_results.csv\n\n')
print(severity_results)

# =============================================================================
# PART K — Climate-Controlled Period Model + Formal Severity Test
# =============================================================================
# Two tests that resolve the open questions left hanging by Parts F, I and J.
# Both are short; both change how the earlier results should be read.
#
# K1 supersedes Part F #1. Part F #1 compared each patch's pre- vs post-
# disturbance gap from the regional reference using a crude two-window
# (pre / post yr 1-3) summary. K1 does the same climate correction but with
# Part I's year-by-year period factor, which is the structure that revealed
# the delayed MidGreendown response in the first place. Run K1; treat Part
# F #1's two figures as the earlier, coarser version of the same idea.
#
# WHY THE YEAR-SPECIFIC REFERENCE, NOT A LINEAR YEAR TERM: an earlier plan
# was to add a fixed linear `year` covariate to Part I's model. Step 4c's
# output shows why that would be inadequate. The regional reference series
# ranges from 270.1 (2012) to 283.9 (2007) -- a 14-day peak-to-peak swing,
# between-year SD ~4.2 days -- while its LINEAR trend is only
# -0.12 +/- 0.14 days/yr (p = 0.41, R^2 = 0.03). A linear term would remove
# the (negligible, non-significant) trend and leave the large year-to-year
# swing untouched. Subtracting the year-SPECIFIC reference removes both.
#
# Note also that the linear-trend test alone could NOT settle the drift
# question: its 95% CI runs -0.40 to +0.16 days/yr, and the lower bound
# over the ~7 years from the pre-window centroid (yr -3) to post_4 implies
# -2.8 days -- numerically identical to the observed post_4 effect. Hence
# this stronger test.
#
# Requires: lme4, lmerTest, dplyr; phenology_long (Part D),
# severity_vs_shift (Part J), and reference_phenology_by_year.csv (Step 4c).
# =============================================================================

# -----------------------------------------------------------------------------
# K1. Period model on the REGIONAL ANOMALY (climate-controlled)
# -----------------------------------------------------------------------------
reference_phenology <- read.csv(file.path(outputDir, 'reference_phenology_by_year.csv'))

# Guard: Step 4c silently exported an all-NA file on its first run, which fed
# blank figures into Part F #1 with no error anywhere. Fail loudly instead.
if (sum(!is.na(reference_phenology$ref_midgreendown_mean)) == 0) {
  stop('reference_phenology_by_year.csv is empty (all years NA). Re-run Step 4c ',
       'and confirm its section-4 cell counts before running Part K.')
}

phenology_long_ref <- phenology_long %>%
  left_join(reference_phenology %>% select(year, ref_midgreendown_mean), by = 'year') %>%
  filter(!is.na(ref_midgreendown_mean)) %>%
  mutate(regional_anomaly = midgreendown_doy_mean - ref_midgreendown_mean)

model_data_anom_period <- phenology_long_ref %>%
  filter(years_since_loss %in% c(-5:-1, 1:5)) %>%
  mutate(period = factor(
    case_when(
      years_since_loss %in% -5:-1 ~ 'pre',
      years_since_loss == 1 ~ 'post_1', years_since_loss == 2 ~ 'post_2',
      years_since_loss == 3 ~ 'post_3', years_since_loss == 4 ~ 'post_4',
      years_since_loss == 5 ~ 'post_5'),
    levels = c('pre', 'post_1', 'post_2', 'post_3', 'post_4', 'post_5')))

# No (1|year) here, deliberately: the anomaly has already had the year effect
# subtracted out by the reference series, so a year random effect would be
# redundant with it. Leaving it out also makes the reference series visibly
# responsible for the climate control, rather than splitting that job
# ambiguously between two mechanisms.
fit_anom_period <- lmerTest::lmer(
  regional_anomaly ~ period + (1 | patch_uuid),
  data = model_data_anom_period
)

cat('\n=== K1: MidGreendown period effects, CLIMATE-CONTROLLED ===\n')
print(summary(fit_anom_period))

# Side-by-side against the raw-DOY period model from Part I. This comparison
# IS the result: if the anomaly estimates hold near the raw ones, the delayed
# senescence advance survives a real climate control; if they collapse toward
# zero, Part F #2's trend-residual figure was right and the effect was
# regional climate all along.
anom_coefs <- summary(fit_anom_period)$coefficients
anom_rows  <- rownames(anom_coefs)[rownames(anom_coefs) != '(Intercept)']

period_comparison <- tibble(
  period            = gsub('^period', '', anom_rows),
  anomaly_estimate  = anom_coefs[anom_rows, 'Estimate'],
  anomaly_se        = anom_coefs[anom_rows, 'Std. Error'],
  anomaly_p         = anom_coefs[anom_rows, 'Pr(>|t|)']
) %>%
  left_join(
    period_results %>%
      filter(response == 'midgreendown_doy_mean') %>%
      select(period, raw_estimate = estimate, raw_p = p_value),
    by = 'period'
  ) %>%
  select(period, raw_estimate, raw_p, anomaly_estimate, anomaly_se, anomaly_p) %>%
  mutate(anomaly_significant = anomaly_p < 0.05)

write.csv(period_comparison,
          file.path(outputDir, 'period_raw_vs_climate_controlled.csv'), row.names = FALSE)

cat('\nRaw-DOY vs climate-controlled period effects (MidGreendown):\n')
print(period_comparison)

p_anom_period <- ggplot(period_comparison,
                        aes(x = factor(period, levels = c('post_1','post_2','post_3','post_4','post_5')))) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'grey40') +
  geom_pointrange(aes(y = anomaly_estimate,
                      ymin = anomaly_estimate - 1.96 * anomaly_se,
                      ymax = anomaly_estimate + 1.96 * anomaly_se,
                      color = anomaly_significant),
                  linewidth = 0.8, size = 0.7) +
  geom_point(aes(y = raw_estimate), shape = 4, size = 3, color = 'grey45') +
  scale_color_manual(values = c(`TRUE` = 'steelblue', `FALSE` = 'grey60'), name = 'p < 0.05') +
  labs(title = 'MidGreendown Effect by Year, Climate-Controlled vs Raw',
       subtitle = paste0('Points + CI = effect on the regional anomaly (each patch-year minus that\n',
                         "YEAR's undisturbed-forest reference, Step 4c). Grey x = the raw-DOY\n",
                         'estimate from Part I, for comparison. Divergence between them is the\n',
                         'portion of the raw effect attributable to regional climate.'),
       x = 'Year since disturbance',
       y = 'MidGreendown anomaly difference from pre-disturbance baseline (days)') +
  theme_minimal()

ggsave(file.path(outputDir, 'period_effects_climate_controlled.png'), p_anom_period,
       width = 8.5, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# K2. Formal per-patch severity test
#     Part J's mixed model found no post:|NDVI change| interaction (p = 0.50),
#     but Part J's own scatter shows a clearly sloped fit. These aren't
#     contradictory so much as different questions: the panel model asks the
#     question with a patch random intercept absorbing exactly the
#     between-patch variation the moderator lives in, which is where its
#     power goes. This lm asks it directly -- one observation per patch, no
#     nesting, no random effect competing for the same variance.
# -----------------------------------------------------------------------------
fit_severity_lm <- lm(midgreendown_shift ~ ndvi_change, data = severity_vs_shift)

cat('\n\n=== K2: Per-patch severity test (one row per patch, no nesting) ===\n')
print(summary(fit_severity_lm))

sev_co <- summary(fit_severity_lm)$coefficients
severity_lm_results <- tibble(
  test        = 'lm(midgreendown_shift ~ ndvi_change), one row per patch',
  n_patches   = nrow(severity_vs_shift),
  estimate    = sev_co['ndvi_change', 'Estimate'],
  std_error   = sev_co['ndvi_change', 'Std. Error'],
  p_value     = sev_co['ndvi_change', 'Pr(>|t|)'],
  r_squared   = summary(fit_severity_lm)$r.squared,
  significant = sev_co['ndvi_change', 'Pr(>|t|)'] < 0.05
)

write.csv(severity_lm_results,
          file.path(outputDir, 'severity_lm_results.csv'), row.names = FALSE)
print(severity_lm_results)

if (severity_lm_results$significant) {
  cat('\nSignificant at the patch level. The correct framing is then: severity\n',
      'DOES predict timing shift patch-to-patch, but the Part J panel model is\n',
      'underpowered to detect it as an interaction -- NOT "no relationship".\n')
} else {
  cat('\nNot significant here either. Part J\'s null interaction stands, and the\n',
      'apparent slope in the scatter is within noise.\n')
}

cat(sprintf('\nPart K outputs saved to %s/:\n', outputDir),
    '  period_effects_climate_controlled.png\n',
    '  period_raw_vs_climate_controlled.csv\n',
    '  severity_lm_results.csv\n')