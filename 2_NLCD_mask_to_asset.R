# Step 2: NLCD Forest Mask — Fractional Cover at MODIS Scale (revised)
# -----------------------------------------------------------------------------
# Order of operations:
#   1. Establish ROI
#   2. Load NLCD 2021, clip to ROI, mask to forest types 41 and 43
#   3. Reproject to MODIS 500m sinusoidal grid
#   4. Compute fractional forest cover MEAN *and* STD DEV per MODIS pixel
#      (the stdDev band is new — it's the within-cell heterogeneity of forest
#      cover, and it's what Step 3's pixel table needs alongside the mean)
#   5. Save as GEE asset (2-band image: forest_fraction_mean, forest_fraction_stdDev)
#   6. Visualize fractional cover
#   7. Export point-sample CSV to Drive
# -----------------------------------------------------------------------------
# Requires: rgee installed and authenticated (ee_Initialize() must succeed)
# =============================================================================

library(rgee)

ee_Initialize()

# -----------------------------------------------------------------------------
# 1. Establish ROI
# -----------------------------------------------------------------------------
roi <- ee$Geometry$Point(c(-79.862539, 37.829550))$buffer(100000)

# -----------------------------------------------------------------------------
# 2. Load NLCD 2021, clip to ROI, mask to forest types 41 and 43
# -----------------------------------------------------------------------------
nlcd <- ee$Image("USGS/NLCD_RELEASES/2021_REL/NLCD/2021")$
  select('landcover')$
  clip(roi)

forestClasses <- c(41L, 43L)

forestMask <- nlcd$remap(
  from         = as.list(forestClasses),
  to           = list(1L, 1L),
  defaultValue = 0
)$rename('forest_mask')

# -----------------------------------------------------------------------------
# 3. Get MODIS 500m sinusoidal projection
# -----------------------------------------------------------------------------
modisProjection <- ee$ImageCollection('MODIS/061/MOD13A1')$
  first()$
  select('NDVI')$
  projection()

# -----------------------------------------------------------------------------
# 4. Compute fractional forest cover mean AND stdDev per MODIS 500m pixel
#    combine(sharedInputs = TRUE) runs both reducers over the same 30m pixels
#    in one pass, producing forest_mask_mean and forest_mask_stdDev bands.
# -----------------------------------------------------------------------------
coverReducer <- ee$Reducer$mean()$combine(
  reducer2      = ee$Reducer$stdDev(),
  sharedInputs  = TRUE
)

forestFractional <- forestMask$
  reduceResolution(
    reducer    = coverReducer,
    bestEffort = TRUE,
    maxPixels  = 1024
  )$
  reproject(
    crs   = modisProjection,
    scale = 500
  )$
  rename(c('forest_fraction_mean', 'forest_fraction_stdDev'))

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

ee_monitoring(task)

# -----------------------------------------------------------------------------
# 6. Visualize fractional cover
# -----------------------------------------------------------------------------
Map$setCenter(-79.862539, 37.829550, 10)

Map$addLayer(
  forestMask$selfMask(),
  list(min = 0, max = 1, palette = c('lightgreen')),
  'NLCD Forest Mask 30m (41 + 43)'
)

Map$addLayer(
  forestFractional$select('forest_fraction_mean')$selfMask(),
  list(min = 0, max = 1, palette = c('white', 'lightgreen', 'darkgreen')),
  'Forest Fraction MEAN at MODIS 500m'
)

Map$addLayer(
  forestFractional$select('forest_fraction_stdDev')$selfMask(),
  list(min = 0, max = 0.5, palette = c('white', 'purple')),
  'Forest Fraction STDDEV at MODIS 500m (edge heterogeneity)'
)

# -----------------------------------------------------------------------------
# 7. Export fractional cover (mean + stdDev) as CSV to Google Drive
# -----------------------------------------------------------------------------
forestFractional <- ee$Image('projects/breidyee/assets/nlcd_forest_fraction_modis')

pixelCentroids <- forestFractional$sample(
  region     = roi,
  scale      = 500,
  projection = modisProjection,
  geometries = TRUE
)

pixelCentroids <- pixelCentroids$map(ee_utils_pyfunc(function(feature) {
  feature$set(
    'longitude', feature$geometry()$coordinates()$get(0),
    'latitude',  feature$geometry()$coordinates()$get(1)
  )
}))

csvTask <- ee_table_to_drive(
  collection  = pixelCentroids,
  description = 'nlcd_forest_fraction_points',
  folder      = 'Reidy_research',
  fileFormat  = 'CSV',
  timePrefix  = FALSE
)

csvTask$start()
ee_monitoring(csvTask)

# -----------------------------------------------------------------------------
# 8. Pull the CSV down locally
# -----------------------------------------------------------------------------
library(googledrive)

# Set this to wherever you want local files saved on YOUR machine
outputDir <- '~/Reidy_research'
if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)

Sys.sleep(15)
fracFile <- drive_ls(path = 'Reidy_research', pattern = 'nlcd_forest_fraction_points')
if (nrow(fracFile) >= 1) {
  drive_download(
    file      = fracFile[1, ],
    path      = file.path(outputDir, 'nlcd_forest_fraction_points.csv'),
    overwrite = TRUE
  )
}
