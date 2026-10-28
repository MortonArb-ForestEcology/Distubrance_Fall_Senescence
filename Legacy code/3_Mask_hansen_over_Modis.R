# Step 3: Attribute Tables — Patch-level AND Pixel-level
# -----------------------------------------------------------------------------
# Produces two CSVs (the deliverables — these are attribute tables, not
# rasters, so nothing here is saved as a GEE image asset):
#
#   A. patch_attribute_table.csv   — one row per Hansen persistent-loss patch
#        patch_uuid, loss_year, area_ha, pixel_count, centroid_lon/lat,
#        elevation_mean/stdDev, aspect_deg, forest_cover_mean/stdDev,
#        meets_forest_threshold (>= 0.75 mean cover — a FLAG, not a filter,
#        so you can decide later which patches to include in analysis)
#
#   B. pixel_attribute_table.csv   — one row per MODIS 500m pixel that
#        touches at least one persistent-loss pixel
#        pixel_uuid, centroid_lon/lat, loss_pixel_count, dominant_lossyear,
#        elevation_mean, aspect_deg, forest_cover_mean/stdDev
#
# Order of operations:
#   1. Establish ROI, load DEM, load NLCD 30m forest mask
#   2. PATCH TABLE: zonal stats (elevation, aspect, cover) over Step 1 polygons
#   3. PIXEL TABLE: zonal stats aggregated to the MODIS grid, then sampled
#   4. Assign pixel_uuid client-side (patch_uuid already exists from Step 1)
#   5. Write both CSVs locally and to Drive
# -----------------------------------------------------------------------------
# Requires: rgee, sf, uuid, googledrive
# =============================================================================

library(rgee)
library(sf)
library(uuid)
library(googledrive)

ee_Initialize(drive = TRUE)

# -----------------------------------------------------------------------------
# 1. Establish ROI, DEM, and 30m forest mask
# -----------------------------------------------------------------------------
roi <- ee$Geometry$Point(c(-79.862539, 37.829550))$buffer(100000)

modisProjection <- ee$ImageCollection('MODIS/061/MOD13A1')$
  first()$select('NDVI')$projection()

# -- DEM: elevation + circular-safe aspect components -------------------------
# Aspect is circular (0 and 360 are the same direction), so a plain mean would
# be wrong near due-north slopes. Instead we average sin/cos components and
# recover the mean angle with atan2 — done once here as image bands so both
# the patch and pixel sections can reuse them.
dem       <- ee$Image('USGS/SRTMGL1_003')$clip(roi)
elevation <- dem$select('elevation')
aspectDeg <- ee$Terrain$aspect(dem)$rename('aspect')
aspectRad <- aspectDeg$multiply(pi / 180)
sinAspect <- aspectRad$sin()$rename('sinAspect')
cosAspect <- aspectRad$cos()$rename('cosAspect')

# -- NLCD binary forest mask at native 30m (classes 41, 43) -------------------
nlcd <- ee$Image("USGS/NLCD_RELEASES/2021_REL/NLCD/2021")$select('landcover')$clip(roi)
forestMask30 <- nlcd$remap(
  from = list(41L, 43L), to = list(1L, 1L), defaultValue = 0
)$rename('forest_cover')

# -----------------------------------------------------------------------------
# 2. PATCH TABLE — zonal stats over each Hansen persistent-loss polygon
# -----------------------------------------------------------------------------
hansenVectors <- ee$FeatureCollection(
  'projects/breidyee/assets/hansen_persistent_loss_vectors'  # carries patch_uuid
)

statsImg1 <- forestMask30$addBands(elevation)  # bands: forest_cover, elevation
reducer1  <- ee$Reducer$mean()$combine(ee$Reducer$stdDev(), sharedInputs = TRUE)

patchStats1 <- statsImg1$reduceRegions(
  collection = hansenVectors,
  reducer    = reducer1,
  scale      = 30,
  tileScale  = 4
)
# -> adds forest_cover_mean, forest_cover_stdDev, elevation_mean, elevation_stdDev

statsImg2 <- sinAspect$addBands(cosAspect)
patchStats2 <- statsImg2$reduceRegions(
  collection = patchStats1,   # chain — geometry/properties carry through
  reducer    = ee$Reducer$mean(),
  scale      = 30,
  tileScale  = 4
)
# -> adds sinAspect, cosAspect (means)
patchTable <- patchStats2$map(ee_utils_pyfunc(function(feature) {
  centroid <- feature$geometry()$centroid(maxError = 1)$coordinates()
  
  # Guard against null sinAspect/cosAspect — can happen for tiny edge-case
  # patches where the DEM aspect band doesn't cleanly resolve. Defaulting
  # to 0/0 gives aspect_deg = 0 for those rather than crashing the export;
  # worth spot-checking flagged patch_uuids afterward if you want to be sure.
  sinSafe <- ee$Number(ee$Algorithms$If(feature$get('sinAspect'), feature$get('sinAspect'), 0))
  cosSafe <- ee$Number(ee$Algorithms$If(feature$get('cosAspect'), feature$get('cosAspect'), 0))
  aspectRadM <- sinSafe$atan2(cosSafe)
  aspectDegM <- aspectRadM$multiply(180 / pi)$add(360)$mod(360)
  
  # Same guard on forest_cover_mean, so a null there doesn't silently break
  # the threshold flag either.
  coverSafe <- ee$Number(ee$Algorithms$If(feature$get('forest_cover_mean'), feature$get('forest_cover_mean'), 0))
  
  feature$set(
    'centroid_lon',           centroid$get(0),
    'centroid_lat',           centroid$get(1),
    'pixel_count',            ee$Number(feature$get('area_ha'))$multiply(10000)$divide(900)$round(),
    'aspect_deg',             aspectDegM,
    'meets_forest_threshold', coverSafe$gte(0.75)
  )
}))

patchTable <- patchTable$select(list(
  'patch_uuid', 'loss_year', 'area_ha', 'pixel_count',
  'centroid_lon', 'centroid_lat',
  'elevation_mean', 'elevation_stdDev', 'aspect_deg',
  'forest_cover_mean', 'forest_cover_stdDev',
  'meets_forest_threshold'
))

patchTask <- ee_table_to_drive(
  collection  = patchTable,
  description = 'patch_attribute_table',
  folder      = 'Reidy_research',
  fileFormat  = 'CSV',
  timePrefix  = FALSE
)
patchTask$start()
cat('Patch attribute table export started.\n')
ee_monitoring(patchTask)

# -----------------------------------------------------------------------------
# 3. PIXEL TABLE — same stats aggregated to the MODIS 500m grid
# -----------------------------------------------------------------------------
hansenImg      <- ee$Image('UMD/hansen/global_forest_change_2024_v1_12')$clip(roi)
persistentLoss <- hansenImg$select('loss')$eq(1)$And(hansenImg$select('gain')$eq(0))
lossyearPix    <- hansenImg$select('lossyear')$updateMask(persistentLoss)

lossPixelCount <- persistentLoss$
  reduceResolution(reducer = ee$Reducer$sum(), bestEffort = TRUE, maxPixels = 1024)$
  reproject(crs = modisProjection, scale = 500)$
  rename('loss_pixel_count')

dominantLossYear <- lossyearPix$
  reduceResolution(reducer = ee$Reducer$mode(), bestEffort = TRUE, maxPixels = 1024)$
  reproject(crs = modisProjection, scale = 500)$
  rename('dominant_lossyear')

coverReducer <- ee$Reducer$mean()$combine(ee$Reducer$stdDev(), sharedInputs = TRUE)
coverPix <- forestMask30$
  reduceResolution(reducer = coverReducer, bestEffort = TRUE, maxPixels = 1024)$
  reproject(crs = modisProjection, scale = 500)$
  rename(c('forest_cover_mean', 'forest_cover_stdDev'))

elevationPix <- elevation$
  reduceResolution(reducer = ee$Reducer$mean(), bestEffort = TRUE, maxPixels = 1024)$
  reproject(crs = modisProjection, scale = 500)$
  rename('elevation_mean')

sinAspectPix <- sinAspect$
  reduceResolution(reducer = ee$Reducer$mean(), bestEffort = TRUE, maxPixels = 1024)$
  reproject(crs = modisProjection, scale = 500)
cosAspectPix <- cosAspect$
  reduceResolution(reducer = ee$Reducer$mean(), bestEffort = TRUE, maxPixels = 1024)$
  reproject(crs = modisProjection, scale = 500)

aspectDegPix <- sinAspectPix$atan2(cosAspectPix)$
  multiply(180 / pi)$add(360)$mod(360)$
  rename('aspect_deg')

pixelStack <- lossPixelCount$
  addBands(dominantLossYear)$
  addBands(coverPix)$
  addBands(elevationPix)$
  addBands(aspectDegPix)

pixelPoints <- pixelStack$sample(
  region     = roi,
  scale      = 500,
  projection = modisProjection,
  geometries = TRUE
)

# Only keep MODIS cells that actually touch persistent loss — drop this
# filter if you want a full background grid for context/mapping.
pixelPoints <- pixelPoints$filter(ee$Filter$gt('loss_pixel_count', 0))

pixelPoints <- pixelPoints$map(ee_utils_pyfunc(function(feature) {
  feature$set(
    'longitude', feature$geometry()$coordinates()$get(0),
    'latitude',  feature$geometry()$coordinates()$get(1)
  )
}))

pixelTask <- ee_table_to_drive(
  collection  = pixelPoints,
  description = 'pixel_attribute_table_raw',
  folder      = 'Reidy_research',
  fileFormat  = 'CSV',
  timePrefix  = FALSE
)
pixelTask$start()
cat('Pixel attribute table (raw, no UUID yet) export started.\n')
ee_monitoring(pixelTask)

# -----------------------------------------------------------------------------
# 4. Download the raw pixel CSV, assign pixel_uuid client-side, re-save
#    (GEE has no UUID generator — same reasoning as Step 1's patch_uuid)
# -----------------------------------------------------------------------------
# Set this to wherever you want local files saved on YOUR machine
outputDir <- '~/Reidy_research'
if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)

Sys.sleep(15)
rawFile <- drive_ls(path = 'Reidy_research', pattern = 'pixel_attribute_table_raw')
stopifnot(nrow(rawFile) >= 1)

localRaw <- file.path(tempdir(), 'pixel_attribute_table_raw.csv')
drive_download(file = rawFile[1, ], path = localRaw, overwrite = TRUE)

pixelDF <- read.csv(localRaw)
pixelDF$pixel_uuid <- vapply(seq_len(nrow(pixelDF)), function(i) UUIDgenerate(), character(1))
stopifnot(!any(duplicated(pixelDF$pixel_uuid)))

localFinal <- file.path(outputDir, 'pixel_attribute_table.csv')
write.csv(pixelDF, localFinal, row.names = FALSE)

drive_upload(
  media     = localFinal,
  path      = 'Reidy_research/',
  name      = 'pixel_attribute_table.csv',
  overwrite = TRUE
)

cat(sprintf('Pixel attribute table finalized: %d rows, saved locally and to Drive.\n', nrow(pixelDF)))

# -----------------------------------------------------------------------------
# 5. Also pull the finished patch table down locally for convenience
# -----------------------------------------------------------------------------
Sys.sleep(15)
patchFile <- drive_ls(path = 'Reidy_research', pattern = 'patch_attribute_table')
if (nrow(patchFile) >= 1) {
  drive_download(
    file      = patchFile[1, ],
    path      = file.path(outputDir, 'patch_attribute_table.csv'),
    overwrite = TRUE
  )
  cat('Patch attribute table also saved locally.\n')
}

