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
#
# "Before" = NLCD forest cover extent / Hansen treecover2000, i.e. what was
#            there before disturbance.
# "After"  = current persistent-loss extent, and which patches/pixels meet
#            the >=0.75 forest-cover criterion (from Step 3's attribute
#            tables).
# -----------------------------------------------------------------------------
# Requires: rgee, sf, terra, tidyterra, ggplot2, patchwork, dplyr
# =============================================================================

library(rgee)
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)
library(patchwork)
library(dplyr)

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
  forestFractional, region = roi, scale = 500, via = 'drive'
)

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
# Set this to the SAME local folder Steps 1 and 3 saved their files to
path.google <- "~/Google Drive/My Drive"
path.dat    <- file.path(path.google, "Reidy_research")
outputDir   <- path.dat
if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)

patches_sf <- st_read(file.path(outputDir, 'hansen_persistent_loss_patches.gpkg'))
patch_attrs <- read.csv(file.path(outputDir, 'patch_attribute_table.csv'))
patch_attrs$meets_forest_threshold <- as.logical(patch_attrs$meets_forest_threshold) 

patches_sf <- merge(
  patches_sf, patch_attrs[, c('patch_uuid', 'meets_forest_threshold', 'forest_cover_mean')],
  by = 'patch_uuid', all.x = TRUE
)

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
    subtitle = 'Line = mean across patches, band = IQR. Dashed line = disturbance year.',
    x        = 'Years since disturbance (0 = loss_year)',
    y        = 'July NDVI'
  ) +
  theme_minimal()

ggsave(file.path(outputDir, 'ndvi_trajectory_relative_to_disturbance.png'), p_trajectory,
       width = 9, height = 6, dpi = 300)

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
    name = 'NDVI change\n(post - pre)'
  ) +
  labs(title = 'NDVI Change Pre- vs Post-Disturbance, by Patch') +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank())

ggsave(file.path(outputDir, 'map_ndvi_change_by_patch.png'), p_ndvi_change,
       width = 8, height = 7, dpi = 300)

cat(sprintf('NDVI figures saved to %s/:\n', outputDir),
    '  ndvi_trajectory_relative_to_disturbance.png\n',
    '  map_ndvi_change_by_patch.png\n',
    sprintf('(%d of %d patches had enough pre/post NDVI coverage to include in the change map)\n',
            nrow(ndvi_change), n_distinct(ndvi_long$patch_uuid)))

