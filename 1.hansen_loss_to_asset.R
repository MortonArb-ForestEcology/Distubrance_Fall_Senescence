# Step 1: Hansen Persistent Loss Asset
# -----------------------------------------------------------------------------
# Order of operations:
#   1. Establish ROI
#   2. Load Hansen 2024, clip to ROI
#   3. Apply persistent loss filter (loss == 1 AND gain == 0)
#   4. Vectorize at native 30m resolution using labelProperty
#   5. Add area property and filter small polygons
#   6. Save as GEE asset
# -----------------------------------------------------------------------------
# Approach follows EEFA Book Chapter F5.1 (Nomura & Bowers) which uses
# labelProperty in reduceToVectors to directly label polygons with lossyear,
# and works at native scale via projection().nominalScale(). No dissolve step
# is needed — individual contiguous patch polygons are preserved, one per
# spatially connected loss event per year.
# Requires: rgee installed and authenticated (ee_Initialize() must succeed)
# =============================================================================

library(rgee)

ee_Initialize()

# -----------------------------------------------------------------------------
# 1. Establish ROI
# -----------------------------------------------------------------------------
roi <- ee$Geometry$Point(c(-79.862539, 37.829550))$buffer(100000)

# -----------------------------------------------------------------------------
# 2. Load Hansen Global Forest Change 2024, clip to ROI
#    Clipping is safe here — we are working at native 30m resolution and
#    not feeding into reduceResolution.
# -----------------------------------------------------------------------------
hansen <- ee$Image('UMD/hansen/global_forest_change_2024_v1_12')$clip(roi)

loss     <- hansen$select('loss')
gain     <- hansen$select('gain')
lossyear <- hansen$select('lossyear')

# -----------------------------------------------------------------------------
# 3. Apply persistent loss filter
#    Persistent loss: canopy removed (loss == 1) with no regrowth (gain == 0).
#    Note: Hansen gain band covers 2000-2012 only — pixels with loss after
#    ~2010 will pass the gain == 0 filter by default since gain is not
#    recorded for that period.
# -----------------------------------------------------------------------------
persistentLoss <- loss$eq(1)$And(gain$eq(0))

# Mask lossyear to persistent loss pixels only
lossYearMasked <- lossyear$updateMask(persistentLoss)

# -----------------------------------------------------------------------------
# 4. Vectorize at native 30m resolution
#    labelProperty passes the lossyear value directly to each output polygon
#    as a named property — no 2-band workaround needed. Each output polygon
#    represents one spatially contiguous patch of persistent loss, labeled
#    with the year that loss occurred (Hansen encoding: 1=2001 ... 24=2024).
#    tileScale = 4 splits computation internally to avoid memory errors.
# -----------------------------------------------------------------------------
lossVectors <- lossYearMasked$reduceToVectors(
  scale          = lossYearMasked$projection()$nominalScale(),
  geometry       = roi,
  geometryType   = 'polygon',
  eightConnected = FALSE,
  labelProperty  = 'lossyear',
  maxPixels      = 1e13,
  tileScale      = 4
)

# -----------------------------------------------------------------------------
# 5. Add area (ha) and calendar year properties, filter out slivers
#    area_ha: polygon area in hectares.
#    loss_year: Hansen lossyear (1-24) converted to calendar year (2001-2024).
#    Filter removes polygons smaller than 0.09 ha (~1 x 30m pixel) which
#    can appear at patch edges due to partial pixel overlap.
# -----------------------------------------------------------------------------
lossVectors <- lossVectors$map(ee_utils_pyfunc(function(feature) {
  feature$set(
    'area_ha',   feature$geometry()$area(maxError = 1)$divide(10000),
    'loss_year', ee$Number(feature$get('lossyear'))$add(2000)
  )
}))

lossVectors <- lossVectors$filter(ee$Filter$gt('area_ha', 0.09))

# -----------------------------------------------------------------------------
# 6. Save as GEE asset
# -----------------------------------------------------------------------------
assetId <- 'projects/breidyee/assets/hansen_persistent_loss_vectors'

task <- ee$batch$Export$table$toAsset(
  collection  = lossVectors,
  description = 'hansen_persistent_loss_vectors',
  assetId     = assetId
)

task$start()
cat(sprintf('Export task started. Asset will be saved to:\n  %s\n', assetId))
ee_monitoring(task)

# -----------------------------------------------------------------------------
# 7. Optional: Visualize in interactive map
# -----------------------------------------------------------------------------
Map$setCenter(-79.862539, 37.829550, 10)

Map$addLayer(
  hansen$select('treecover2000')$updateMask(hansen$select('treecover2000')$gt(0)),
  list(min = 0, max = 100, palette = c('black', 'green')),
  'Tree Cover 2000'
)

Map$addLayer(
  lossYearMasked,
  list(min = 1, max = 24, palette = c('yellow', 'orange', 'red')),
  'Persistent Loss Year (30m)'
)
