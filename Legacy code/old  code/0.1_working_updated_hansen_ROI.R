# =============================================================================
# Script 1: Forest Loss Identification via rgee
# Purpose: Identify persistent forest loss polygons (loss without recovery)
#          within the region of interest using Hansen GFC and NLCD data.
#          Each output polygon = one full 500 m MODIS pixel (non-overlapping).
#          Only pixels that were forested in both 2001 and 2021 (NLCD) are
#          included, ensuring disturbance occurred within a persistent forest
#          matrix rather than representing land use conversion.
# Output:  CSV of loss polygons with area, centroid coordinates, year of loss,
#          and gain flag — written to Google Drive.
# Requires: rgee installed and authenticated (ee_Initialize() must succeed)
# =============================================================================
### Start with a code that creates a clipped region to run this  

library(rgee)

ee_Initialize(drive = TRUE)

# -----------------------------------------------------------------------------
# 1. Define region of interest
# -----------------------------------------------------------------------------
roi <- ee$Geometry$Point(c(-79.862539, 37.829550))$buffer(100000)

# -----------------------------------------------------------------------------
# 2. Load Hansen Global Forest Change 2024
# -----------------------------------------------------------------------------
hansen <- ee$Image('UMD/hansen/global_forest_change_2024_v1_12')$clip(roi)

# -----------------------------------------------------------------------------
# 3. Get MODIS projection
#    Passed into reduceResolution and reduceToVectors so all aggregation
#    and vectorization occurs on the native MODIS sinusoidal grid.
#    No reproject() calls are used — reproject() triggers eager local
#    computation in rgee before tasks reach GEE servers.
# -----------------------------------------------------------------------------
modisProjection <- ee$ImageCollection('MODIS/061/MOD13A1')$
  first()$
  select('NDVI')$
  projection()

# -----------------------------------------------------------------------------
# 4. Build stable forest mask
#    Pixel must be deciduous (41) or mixed (43) forest in BOTH 2001 and
#    2021 NLCD. reduceResolution with min() requires every 30 m pixel
#    within a 500 m cell to be forest — ensures passing cells are
#    unambiguously forest-dominated. crs = modisProjection passed into
#    reduceResolution anchors aggregation to the MODIS sinusoidal grid
#    server-side without triggering local computation.
# -----------------------------------------------------------------------------
nlcd2001 <- ee$Image("USGS/NLCD_RELEASES/2019_REL/NLCD/2001")$select('landcover')
nlcd2021 <- ee$Image("USGS/NLCD_RELEASES/2021_REL/NLCD/2021")$select('landcover')

forestClasses <- c(41L, 43L)  # 41 = deciduous forest, 43 = mixed forest

forestMask2001 <- nlcd2001$remap(
  from         = as.list(forestClasses),
  to           = as.list(forestClasses),
  defaultValue = 0
)$neq(0)$
  reduceResolution(
    reducer    = ee$Reducer$min(),
    bestEffort = TRUE,
    maxPixels  = 1024
  )

forestMask2021 <- nlcd2021$remap(
  from         = as.list(forestClasses),
  to           = as.list(forestClasses),
  defaultValue = 0
)$neq(0)$
  reduceResolution(
    reducer    = ee$Reducer$min(),
    bestEffort = TRUE,
    maxPixels  = 1024
  )

stableForestMask <- forestMask2001$And(forestMask2021)

# -----------------------------------------------------------------------------
# 5. Identify persistent loss within stable forest mask
#    loss == 1 AND gain == 0: canopy removed and no regrowth detected.
#    Hansen gain band covers 2000-2012 only — the 2021 NLCD forest
#    requirement serves as the primary post-disturbance filter for loss
#    events after ~2010 where the gain band is uninformative.
# -----------------------------------------------------------------------------
loss <- hansen$select('loss')
gain <- hansen$select('gain')

persistentLoss <- loss$gt(0)$And(gain$eq(0))
persistentLossMasked <- persistentLoss$updateMask(stableForestMask)

# -----------------------------------------------------------------------------
# 6. Aggregate to 500 m MODIS pixel grid
#    persistentLoss500m: max() — any 30 m loss pixel marks the whole
#    500 m cell as loss.
#    lossYear500m: min() — earliest loss year within the cell, so the
#    pre-disturbance NDVI baseline window in downstream scripts is
#    anchored to the first detectable disturbance event.
# -----------------------------------------------------------------------------
persistentLoss500m <- persistentLossMasked$
  reduceResolution(
    reducer    = ee$Reducer$max(),
    bestEffort = TRUE,
    maxPixels  = 1024
  )

lossYear500m <- hansen$select('lossyear')$
  updateMask(persistentLossMasked)$
  reduceResolution(
    reducer    = ee$Reducer$min(),
    bestEffort = TRUE,
    maxPixels  = 1024
  )

# -----------------------------------------------------------------------------
# 7. Vectorize 500 m loss pixels to polygons
#    lossyear is the first band so reduceToVectors first() operates on it.
#    GEE always names the output property 'first' regardless of band name —
#    renamed to 'lossyear' in section 8 via feature$get('first').
#    crs = modisProjection pins output polygons to the MODIS sinusoidal
#    grid and cascades back through the computation graph to ensure all
#    reduceResolution calls evaluate in the same pixel space.
# -----------------------------------------------------------------------------
multiImage <- lossYear500m$rename('lossyear')$
  addBands(ee$Image$constant(1L)$updateMask(persistentLoss500m)$rename('loss_flag'))

lossPolygons <- multiImage$reduceToVectors(
  reducer      = ee$Reducer$first(),
  geometry     = roi,
  scale        = 500,
  crs          = modisProjection,
  geometryType = 'polygon',
  maxPixels    = 1e9
)

# -----------------------------------------------------------------------------
# 8. Add metadata to each polygon
# -----------------------------------------------------------------------------
lossWithMetadata <- lossPolygons$map(ee_utils_pyfunc(function(feature) {
  feature$set(
    'area_sqm',  feature$geometry()$area(maxError = 1),
    'longitude', feature$geometry()$centroid(maxError = 1)$coordinates()$get(0),
    'latitude',  feature$geometry()$centroid(maxError = 1)$coordinates()$get(1),
    'lossyear',  feature$get('first'),
    # gain_30m_presence: sampled from raw 30 m Hansen gain band.
    # NOT the same as the gain$eq(0) filter in section 5 — that operated
    # at 500 m. This column is for QA reference only, do not use as a
    # filter criterion downstream.
    'gain_30m_presence', gain$reduceRegion(
      reducer   = ee$Reducer$max(),
      geometry  = feature$geometry(),
      scale     = 30,
      maxPixels = 1e6
    )$get('gain'),
    'nlcd_2001_class', nlcd2001$reduceRegion(
      reducer   = ee$Reducer$mode(),
      geometry  = feature$geometry(),
      scale     = 500,
      maxPixels = 1e6
    )$get('landcover'),
    'nlcd_2021_class', nlcd2021$reduceRegion(
      reducer   = ee$Reducer$mode(),
      geometry  = feature$geometry(),
      scale     = 500,
      maxPixels = 1e6
    )$get('landcover')
  )
}))

# -----------------------------------------------------------------------------
# 9. QA filters
# -----------------------------------------------------------------------------
# Area filter: ~250,000 m² expected per MODIS pixel, allow 10% tolerance
lossWithMetadata <- lossWithMetadata$filter(
  ee$Filter$And(
    ee$Filter$gt('area_sqm', 225000),
    ee$Filter$lt('area_sqm', 275000)
  )
)

# NLCD filter: hard enforce stable forest mask on output columns
validForestClasses <- c(41, 43)

lossWithMetadata <- lossWithMetadata$filter(
  ee$Filter$And(
    ee$Filter$inList('nlcd_2001_class', as.list(validForestClasses)),
    ee$Filter$inList('nlcd_2021_class', as.list(validForestClasses))
  )
)

# -----------------------------------------------------------------------------
# 10. Export to Google Drive as CSV
# -----------------------------------------------------------------------------
task <- ee_table_to_drive(
  collection  = lossWithMetadata,
  description = 'forest_loss_polygons',
  folder      = 'Reidy_research',
  fileFormat  = 'CSV',
  timePrefix  = FALSE
)

task$start()
ee_monitoring(task)

# -----------------------------------------------------------------------------
# 11. Optional: Interactive map visualization
# -----------------------------------------------------------------------------
Map$setCenter(-79.862539, 37.829550, 10)

Map$addLayer(
  hansen$updateMask(stableForestMask),
  list(bands = 'treecover2000', min = 0, max = 100,
       palette = c('black', 'green')),
  'Tree Cover 2000 (stable forest only)'
)
Map$addLayer(
  hansen$updateMask(stableForestMask),
  list(bands = 'lossyear', min = 0, max = 23,
       palette = c('yellow', 'red')),
  'Tree Loss Year (stable forest only)'
)
Map$addLayer(
  persistentLoss500m$rename('persistent_loss_500m'),
  list(bands = 'persistent_loss_500m', min = 0, max = 1,
       palette = c('green', 'red')),
  'Persistent Loss at 500 m (MODIS grid)'
)
