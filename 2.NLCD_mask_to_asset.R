# Step 2: NLCD Forest Mask — Fractional Cover at MODIS Scale
# -----------------------------------------------------------------------------
# Order of operations:
#   1. Establish ROI
#   2. Load NLCD 2021, clip to ROI, mask to forest types 41 and 43
#   3. Reproject to MODIS 500m sinusoidal grid
#   4. Compute fractional forest cover (0-1) per MODIS pixel
#   5. Save as GEE asset
#   6. Visualize fractional cover
# -----------------------------------------------------------------------------
# Requires: rgee installed and authenticated (ee_Initialize() must succeed)
# =============================================================================

library(rgee)

ee_Initialize()

# -----------------------------------------------------------------------------
# 1. Establish ROI
#    Same ROI as Step 1 — 100km buffer around central Virginia point.
# -----------------------------------------------------------------------------
roi <- ee$Geometry$Point(c(-79.862539, 37.829550))$buffer(100000)

# -----------------------------------------------------------------------------
# 2. Load NLCD 2021, clip to ROI, mask to forest types 41 and 43
#    NLCD class 41 = deciduous forest
#    NLCD class 43 = mixed forest
#    Binary mask: 1 where pixel is forest type of interest, 0 elsewhere.
# -----------------------------------------------------------------------------
nlcd <- ee$Image("USGS/NLCD_RELEASES/2021_REL/NLCD/2021")$
  select('landcover')$
  clip(roi)

forestClasses <- c(41L, 43L)

# Remap forest classes to 1, everything else to 0
forestMask <- nlcd$remap(
  from         = as.list(forestClasses),
  to           = list(1L, 1L),
  defaultValue = 0
)$rename('forest_mask')

# -----------------------------------------------------------------------------
# 3. Get MODIS 500m sinusoidal projection
#    Used to anchor the reprojection and aggregation to the native MODIS grid.
# -----------------------------------------------------------------------------
modisProjection <- ee$ImageCollection('MODIS/061/MOD13A1')$
  first()$
  select('NDVI')$
  projection()

# -----------------------------------------------------------------------------
# 4. Compute fractional forest cover per MODIS 500m pixel
#    mean() reducer over the binary 30m forest mask gives the fraction of
#    30m pixels within each 500m cell that are forest type 41 or 43.
#    Output is a 0-1 continuous value per MODIS pixel where:
#      0   = no deciduous/mixed forest in this 500m cell
#      0.5 = half of the 30m pixels in this cell are forest
#      1   = entire 500m cell is deciduous/mixed forest
#    reproject() anchors the output to the MODIS grid server-side.
# -----------------------------------------------------------------------------
forestFractional <- forestMask$
  reduceResolution(
    reducer    = ee$Reducer$mean(),
    bestEffort = TRUE,
    maxPixels  = 1024
  )$
  reproject(
    crs   = modisProjection,
    scale = 500
  )$
  rename('forest_fraction')

# -----------------------------------------------------------------------------
# 5. Save as GEE asset
# -----------------------------------------------------------------------------
assetId <- 'projects/breidyee/assets/nlcd_forest_fraction_modis'

task <- ee$batch$Export$image$toAsset(
  image       = forestFractional,
  description = 'nlcd_forest_fraction_modis',
  assetId     = assetId,
  region      = roi,
  crs         = modisProjection,
  scale       = 500,
  maxPixels   = 1e10
)

task$start()
cat(sprintf('Export task started. Asset will be saved to:\n  %s\n', assetId))
ee_monitoring(task)

# -----------------------------------------------------------------------------
# 6. Visualize fractional cover
#    Color ramp from white (0 = no forest) to dark green (1 = full forest).
#    Also display the raw NLCD forest mask for comparison.
# -----------------------------------------------------------------------------
Map$setCenter(-79.862539, 37.829550, 10)

Map$addLayer(
  forestMask$selfMask(),
  list(min = 0, max = 1, palette = c('lightgreen')),
  'NLCD Forest Mask 30m (41 + 43)'
)

Map$addLayer(
  forestFractional$selfMask(),
  list(min = 0, max = 1, palette = c('white', 'lightgreen', 'darkgreen')),
  'Forest Fractional Cover at MODIS 500m'
)
