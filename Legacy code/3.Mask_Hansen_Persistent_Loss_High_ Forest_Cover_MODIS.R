# Step 3: Mask Hansen Persistent Loss to High Forest Cover MODIS Pixels
# -----------------------------------------------------------------------------
# Order of operations:
#   1. Establish ROI
#   2. Load Step 1 asset (Hansen persistent loss vectors)
#   3. Load Step 2 asset (NLCD forest fraction raster)
#   4. Create binary mask where forest_fraction >= 0.75
#   5. Rasterize Hansen loss vectors back to 30m
#   6. Mask Hansen loss raster to high forest cover pixels
#   7. Save masked Hansen loss as GEE asset
#   8. Visualize
# -----------------------------------------------------------------------------
# Threshold: only MODIS 500m pixels where >= 75% of 30m pixels are
# deciduous (41) or mixed (43) forest are retained for analysis.
# Requires: rgee installed and authenticated (ee_Initialize() must succeed)
# =============================================================================

library(rgee)

ee_Initialize()

# -----------------------------------------------------------------------------
# 1. Establish ROI
# -----------------------------------------------------------------------------
roi <- ee$Geometry$Point(c(-79.862539, 37.829550))$buffer(100000)

# -----------------------------------------------------------------------------
# 2. Load Step 1 asset — Hansen persistent loss vectors
#    One polygon per contiguous disturbance patch, labeled with lossyear
#    (Hansen encoding 1-24) and loss_year (calendar year 2001-2024).
# -----------------------------------------------------------------------------
hansenVectors <- ee$FeatureCollection(
  'projects/breidyee/assets/hansen_persistent_loss_vectors'
)

# -----------------------------------------------------------------------------
# 3. Load Step 2 asset — NLCD forest fraction raster
#    500m MODIS pixels with 0-1 fractional cover of classes 41 and 43.
# -----------------------------------------------------------------------------
forestFractional <- ee$Image(
  'projects/breidyee/assets/nlcd_forest_fraction_modis'
)

# -----------------------------------------------------------------------------
# 4. Get MODIS projection
# -----------------------------------------------------------------------------
modisProjection <- ee$ImageCollection('MODIS/061/MOD13A1')$
  first()$
  select('NDVI')$
  projection()

# -----------------------------------------------------------------------------
# 5. Create binary forest mask at MODIS 500m
#    Pixels where forest_fraction >= 0.75 are retained (value = 1).
#    All other pixels are masked out.
# -----------------------------------------------------------------------------
forestThresholdMask <- forestFractional$gte(0.75)$selfMask()

# -----------------------------------------------------------------------------
# 6. Rasterize Hansen loss vectors to 30m
#    Converts the Step 1 vector asset back to a raster so it can be
#    spatially masked by the 500m forest threshold mask.
#    Each 30m pixel gets the lossyear value of the polygon it falls within.
# -----------------------------------------------------------------------------
hansenLossRaster <- hansenVectors$
  filter(ee$Filter$notNull(list('lossyear')))$
  reduceToImage(
    properties = list('lossyear'),
    reducer    = ee$Reducer$first()
  )$
  rename('lossyear')$
  clip(roi)

# -----------------------------------------------------------------------------
# 7. Mask Hansen loss raster to high forest cover MODIS pixels
#    forestThresholdMask is at 500m — updateMask propagates it down to
#    the 30m Hansen pixels that fall within qualifying MODIS cells.
#    Only 30m loss pixels inside MODIS cells with >= 75% forest cover
#    are retained.
# -----------------------------------------------------------------------------
hansenMasked <- hansenLossRaster$updateMask(forestThresholdMask)

# -----------------------------------------------------------------------------
# 8. Save masked Hansen loss raster as GEE asset
# -----------------------------------------------------------------------------
assetId <- 'projects/breidyee/assets/hansen_loss_forest_masked'

task <- ee$batch$Export$image$toAsset(
  image       = hansenMasked,
  description = 'hansen_loss_forest_masked',
  assetId     = assetId,
  region      = roi,
  scale       = 30,
  maxPixels   = 1e10
)

task$start()
cat(sprintf('Export task started. Asset will be saved to:\n  %s\n', assetId))
ee_monitoring(task)

# -----------------------------------------------------------------------------
# 9. Visualize
# -----------------------------------------------------------------------------
Map$setCenter(-79.862539, 37.829550, 10)

Map$addLayer(
  forestThresholdMask,
  list(min = 0, max = 1, palette = c('darkgreen')),
  'Forest Fraction >= 0.75 (MODIS 500m)'
)

Map$addLayer(
  hansenLossRaster,
  list(min = 1, max = 24, palette = c('yellow', 'orange', 'red')),
  'Hansen Persistent Loss — Unmasked (30m)'
)

Map$addLayer(
  hansenMasked,
  list(min = 1, max = 24, palette = c('yellow', 'orange', 'red')),
  'Hansen Persistent Loss — Forest Masked (30m)'
)
