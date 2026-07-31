# Step 4b: MODIS Phenology Metrics by Patch — Greenup / Senescence / Dormancy
# -----------------------------------------------------------------------------
# Forked from Step 4 per Critic review: the "fall color project" name implies
# autumn senescence phenology, but Step 4 only ever samples July (peak
# growing-season) NDVI, which cannot speak to fall coloration at all. Rather
# than derive senescence timing ourselves from raw NDVI, this script pulls it
# directly from MODIS's purpose-built land-surface-phenology product,
# MCD12Q2, which already computes per-pixel, per-year Greenup/Senescence/
# Dormancy transition dates.
#
# Deliberately runs over the FULL, unfiltered patch catalog (patch_uuid +
# loss_year only, no area_ha filter) -- unlike Step 4's NDVI extraction,
# which restricts to patches >= 25 ha (one MODIS pixel) to avoid neighborhood
# dilution. Both scripts share the same 500m MODIS grid, so the same
# dilution caveat applies here for small patches; this script does not
# filter on it so the phenology dataset stays maximally complete. Apply an
# area filter downstream at analysis time if patch-level (rather than
# neighborhood-level) phenology is what's needed.
#
# Output: one long-format CSV -- one row per (patch_uuid, year).
#   patch_uuid, loss_year, year,
#   greenup_doy_mean, greenup_doy_stdDev,
#   senescence_doy_mean, senescence_doy_stdDev,
#   midgreendown_doy_mean, midgreendown_doy_stdDev,
#   dormancy_doy_mean, dormancy_doy_stdDev,
#   qa_overall_mean, qa_overall_stdDev (QA_Overall_1 diagnostic, masked to
#   its documented 0-3 valid range -- 0=best, 1=good, 2=fair, 3=poor;
#   not used to filter the phenology dates themselves, see convertToDoy())
# -----------------------------------------------------------------------------
# Requires: rgee installed and authenticated (ee_Initialize() must succeed)
# =============================================================================

library(rgee)

ee_Initialize(drive = TRUE)

# -----------------------------------------------------------------------------
# 1. Establish ROI (same as Steps 1-4)
# -----------------------------------------------------------------------------
roi <- ee$Geometry$Point(c(-79.862539, 37.829550))$buffer(100000)

# -----------------------------------------------------------------------------
# 2. Load Step 1 asset — full, unfiltered patch catalog
# -----------------------------------------------------------------------------
hansenVectors <- ee$FeatureCollection(
  'projects/breidyee/assets/hansen_persistent_loss_vectors_2'
)$select(list('patch_uuid', 'loss_year'))

# -----------------------------------------------------------------------------
# 3. MCD12Q2 Land Cover Dynamics (phenology), cycle 1 only
#    Cycle 2 exists for regions with two growing seasons per year
#    (double-cropping, some evergreen/grassland dynamics) -- not relevant to
#    Virginia deciduous forest, so only the "_1" bands are used.
#
#    MidGreendown_1 added per Critic review of the official MCD12Q2 v6.1 User
#    Guide: the product's own "Known Issues" section recommends MidGreenup/
#    MidGreendown (50% amplitude crossing) over Greenup/Dormancy (15%
#    crossing) as "the more realistic and stable" phenophase estimates. For
#    this project specifically, MidGreendown is also the most representative
#    single "peak fall-color transition" day-of-year available in this
#    product -- more central than Senescence (90% crossing, onset) or
#    Dormancy (15% crossing, near-total leaf-off).
# -----------------------------------------------------------------------------
phenoBands <- c('Greenup_1', 'Senescence_1', 'MidGreendown_1', 'Dormancy_1')

# MCD12Q2's date bands are stored as days-since-1970-01-01 (epoch days), NOT
# day-of-year -- treating the raw value as a DOY directly would silently
# produce nonsense dates. Convert explicitly, per image, using that image's
# own year.
convertToDoy <- function(image) {
  imgYear   <- ee$Date(image$get('system:time_start'))$get('year')
  yearStart <- ee$Date$fromYMD(imgYear, 1, 1)
  epochDaysAtYearStart <- yearStart$difference(ee$Date('1970-01-01'), 'day')

  doyImg <- image$select(phenoBands)$
    subtract(epochDaysAtYearStart)$
    rename(paste0(phenoBands, '_doy'))

  # Pixels with no detected phenological cycle that year (e.g. evergreen
  # conifer, persistent snow/urban, missing data) carry a large fill value
  # in the raw band. After the epoch-to-DOY conversion those fill values
  # land far outside a real day-of-year, so mask on the plausible [1, 366]
  # range rather than relying on the exact raw fill code.
  validMask <- doyImg$gte(1)$And(doyImg$lte(366))
  doyImg <- doyImg$updateMask(validMask)

  # QA_Overall_1: per-pixel confidence flag for the cycle-1 phenology
  # retrieval. Confirmed via the official MCD12Q2 v6.1 User Guide (Table 1):
  # valid range 0-3 (0=best, 1=good, 2=fair, 3=poor), fill value 32767 for
  # pixels with no valid vegetation cycle that year. Masked to the documented
  # valid range here -- addBands() does NOT propagate doyImg's mask onto this
  # band, so without an explicit mask a fill-valued pixel (32767) would
  # silently corrupt the qa_overall_mean/stdDev reduceRegions calculation
  # below. Carried through as its own diagnostic band (not used to filter
  # the phenology dates themselves) so downstream analysis can decide
  # whether/how to use it -- the same philosophy as Step 4's ndvi_july_count.
  qaImg <- image$select('QA_Overall_1')
  qaImg <- qaImg$updateMask(qaImg$gte(0)$And(qaImg$lte(3)))

  doyImg$addBands(qaImg)$copyProperties(image, list('system:time_start'))
}

phenologyCollection <- ee$ImageCollection('MODIS/061/MCD12Q2')$
  filterBounds(roi)$
  map(ee_utils_pyfunc(convertToDoy))

# -----------------------------------------------------------------------------
# 4. Build one year's phenology means, reduced over every patch
#    MCD12Q2 is already an annual product (one image per calendar year), so
#    -- unlike Step 4's 16-day MOD13A1 composites -- no within-year
#    filtering is needed beyond selecting that year's single image.
# -----------------------------------------------------------------------------
startYear <- 2001  # MCD12Q2 begins with the 2001 growing season
endYear   <- 2022  # most recent finalized vintage at time of writing --
                    # check the MCD12Q2 catalog page before assuming later
                    # years are available (dormancy detection needs data
                    # from into the following year to finalize)

yearlyPhenology <- function(yr) {
  yrImage <- phenologyCollection$
    filter(ee$Filter$calendarRange(yr, yr, 'year'))$
    first()

  # mean + stdDev combined into one reducer/one reduceRegions pass (same
  # pattern as Steps 2, 3, and the revised Step 4). With a multi-band input
  # image (4 phenology DOY bands + QA_Overall_1), EE's default per-band
  # output naming (e.g. 'Greenup_1_doy_mean') is already what we want, so no
  # setOutputs() renaming is needed here -- that was only required in Step 4
  # because its input was a single band.
  statsReducer <- ee$Reducer$mean()$combine(
    reducer2     = ee$Reducer$stdDev(),
    sharedInputs = TRUE
  )

  fc <- yrImage$reduceRegions(
    collection = hansenVectors, reducer = statsReducer, scale = 500, tileScale = 4
  )

  # Same guard as Step 4: explicitly re-set every expected property (even to
  # a null get()) so no feature can cause GEE's CSV exporter to drop a
  # column just because one feature happened to lack that key.
  fc$map(ee_utils_pyfunc(function(feature) {
    feature$set(
      'year',                     yr,
      'greenup_doy_mean',         feature$get('Greenup_1_doy_mean'),
      'greenup_doy_stdDev',       feature$get('Greenup_1_doy_stdDev'),
      'senescence_doy_mean',      feature$get('Senescence_1_doy_mean'),
      'senescence_doy_stdDev',    feature$get('Senescence_1_doy_stdDev'),
      'midgreendown_doy_mean',    feature$get('MidGreendown_1_doy_mean'),
      'midgreendown_doy_stdDev',  feature$get('MidGreendown_1_doy_stdDev'),
      'dormancy_doy_mean',        feature$get('Dormancy_1_doy_mean'),
      'dormancy_doy_stdDev',      feature$get('Dormancy_1_doy_stdDev'),
      'qa_overall_mean',          feature$get('QA_Overall_1_mean'),
      'qa_overall_stdDev',        feature$get('QA_Overall_1_stdDev')
    )
  }))
}

allYearsPhenology <- ee$FeatureCollection(
  lapply(startYear:endYear, yearlyPhenology)
)$flatten()

# -----------------------------------------------------------------------------
# 5. Export as one long CSV (single batch task)
#    Columns specified via `selectors` on the export call itself, for the
#    same reason as Step 4: relying on select() + the exporter's automatic
#    header inference is what silently dropped columns in Step 4's earlier
#    attempts.
# -----------------------------------------------------------------------------
phenologyColumns <- c('patch_uuid', 'loss_year', 'year',
                       'greenup_doy_mean', 'greenup_doy_stdDev',
                       'senescence_doy_mean', 'senescence_doy_stdDev',
                       'midgreendown_doy_mean', 'midgreendown_doy_stdDev',
                       'dormancy_doy_mean', 'dormancy_doy_stdDev',
                       'qa_overall_mean', 'qa_overall_stdDev')

task <- ee_table_to_drive(
  collection  = allYearsPhenology,
  description = 'phenology_yearly_by_patch',
  folder      = 'Reidy_research',
  fileFormat  = 'CSV',
  selectors   = phenologyColumns,
  timePrefix  = FALSE
)
task$start()
cat('MODIS phenology (MCD12Q2) yearly export started: phenology_yearly_by_patch\n')
ee_monitoring(task)

# -----------------------------------------------------------------------------
# 6. Pull it down locally
# -----------------------------------------------------------------------------
library(googledrive)

# Shared output directory (see config.R) — keeps every pipeline script
# writing to/reading from the same place. Sourced via a fixed absolute path
# (not a bare relative "config.R") so this doesn't depend on R's working
# directory at run time.
source(file.path("~/Google Drive/My Drive/Reidy_research/fall color code", "config.R"))

Sys.sleep(15)
phenoFile <- drive_ls(path = 'Reidy_research', pattern = 'phenology_yearly_by_patch')
if (nrow(phenoFile) >= 1) {
  drive_download(
    file      = phenoFile[1, ],
    path      = file.path(outputDir, 'phenology_yearly_by_patch.csv'),
    overwrite = TRUE
  )
  cat('MODIS phenology yearly series saved locally.\n')
}

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
4