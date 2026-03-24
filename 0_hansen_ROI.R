# =============================================================================
# Script 1: Forest Loss Identification via rgee
# Purpose: Identify persistent forest loss polygons (loss without recovery)
#          within the region of interest using Hansen GFC and NLCD data.
# Output:  CSV of loss polygons with area, centroid coordinates, year of loss,
#          and gain flag — written to Google Drive.
# Requires: rgee installed and authenticated (ee_Initialize() must succeed)
# =============================================================================

library(rgee)

# Initialize Earth Engine
ee_Initialize(drive = TRUE)  # drive = TRUE needed for export to Drive

# -----------------------------------------------------------------------------
# 1. Define region of interest
# -----------------------------------------------------------------------------
# 5 km buffer around center coordinate (Virginia/WV region)
roi <- ee$Geometry$Point(c(-79.862539, 37.829550))$buffer(5000)

# -----------------------------------------------------------------------------
# 2. Load and clip Hansen Global Forest Change 2023
# -----------------------------------------------------------------------------
hansen <- ee$Image('UMD/hansen/global_forest_change_2024_v1_12')$clip(roi)

# -----------------------------------------------------------------------------
# 3. Load NLCD 2001 and create deciduous/mixed forest mask
# -----------------------------------------------------------------------------
# Using 2001 as baseline land cover (classes 41 = deciduous, 43 = mixed forest)
nlcd <- ee$ImageCollection("USGS/NLCD_RELEASES/2019_REL/NLCD")
nlcd2001 <- nlcd$filter(ee$Filter$eq('system:index', '2001'))$first()
landcover <- nlcd2001$select('landcover')

forestClasses <- c(41L, 43L)  # integer values required for remap
forestMask <- landcover$remap(
  from = as.list(forestClasses),
  to   = as.list(forestClasses),
  defaultValue = 0
)$neq(0)

# -----------------------------------------------------------------------------
# 4. Identify persistent loss (loss == 1 AND gain == 0)
#    and disturbed-but-recovered areas (loss == 1 AND gain == 1)
#
#    Persistent loss = areas where canopy was lost and DID NOT regrow.
#    These are the study polygons. Recovered areas are excluded to ensure
#    we are studying true disturbance events, not temporary canopy gaps
#    that have since reforested — which would confound phenology comparisons.
# -----------------------------------------------------------------------------
loss <- hansen$select('loss')
gain <- hansen$select('gain')

# Persistent loss: lost forest that did not recover
persistentLoss <- loss$gt(0)$And(gain$eq(0))

# Disturbed-but-recovered: lost forest that regrew (excluded from analysis)
# Retained here for reference/QA only
recoveredLoss <- loss$gt(0)$And(gain$eq(1))

# Apply forest mask to persistent loss layer
persistentLossMasked <- persistentLoss$updateMask(forestMask)

# -----------------------------------------------------------------------------
# 5. Extract loss polygons with metadata
#    Uses lossyear band to preserve year of loss per polygon
# -----------------------------------------------------------------------------
# Build multi-band image: loss year + constant band (needed for reduceToVectors)
multiImage <- hansen$select('lossyear')$
  updateMask(persistentLossMasked)$   # apply persistent loss + forest mask
  addBands(ee$Image$constant(1L)$updateMask(persistentLossMasked))

# Vectorize: convert raster loss patches to polygon features
lossPolygons <- multiImage$reduceToVectors(
  reducer      = ee$Reducer$first(),
  geometry     = roi,
  scale        = 30,
  geometryType = 'polygon',
  maxPixels    = 1e9
)

# Add area (m²) and centroid coordinates to each polygon feature
lossWithMetadata <- lossPolygons$map(ee_utils_pyfunc(function(feature) {
  feature$set(
    'area_sqm',  feature$geometry()$area(maxError = 1),
    'longitude', feature$geometry()$centroid(maxError = 1)$coordinates()$get(0),
    'latitude',  feature$geometry()$centroid(maxError = 1)$coordinates()$get(1),
    'gain',      gain$reduceRegion(
      reducer  = ee$Reducer$max(),
      geometry = feature$geometry(),
      scale    = 30,
      maxPixels = 1e6
    )$get('gain')
  )
}))

# -----------------------------------------------------------------------------
# 6. Export to Google Drive as CSV
# -----------------------------------------------------------------------------
#note this will have to be manually moved to a file withi
task <- ee_table_to_drive(
  collection  = lossWithMetadata,
  description = 'forest_loss_polygons',
  folder      = 'Reidy_research',
  fileFormat  = 'CSV',
  timePrefix  = FALSE
)

task$start()

# Monitor export progress (optional — can also check in GEE Tasks panel)
ee_monitoring(task)

# -----------------------------------------------------------------------------
# 7. Optional: Interactive map visualization (runs in rgee viewer)
# -----------------------------------------------------------------------------
# Uncomment to preview layers in the rgee interactive map

 Map$setCenter(-79.862539, 37.829550, 12)

 Map$addLayer(
   hansen$updateMask(forestMask),
  list(bands = 'treecover2000', min = 0, max = 100,
       palette = c('black', 'green')),
  'Tree Cover 2000'
)
Map$addLayer(
  hansen$updateMask(forestMask),
  list(bands = 'lossyear', min = 0, max = 23,
       palette = c('yellow', 'red')),
  'Tree Loss Year'
)
Map$addLayer(
  hansen$updateMask(forestMask),
  list(bands = 'gain', min = 0, max = 1,
       palette = c('green', 'blue')),
  'Tree Gain'
)
Map$addLayer(
  persistentLossMasked$rename('persistent_loss'),
  list(bands = 'persistent_loss', min = 0, max = 1,
       palette = c('green', 'red')),
  'Persistent Loss (no recovery)'
)
