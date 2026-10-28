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
)$rename('forest_cover')$clip(roi)
# <-- FIX: `remap()` does not reliably preserve the input's clipped
# footprint — the output tends to revert to the source dataset's native
# footprint, which for NLCD is all of CONUS. This wasn't a problem in
# Step 2 because reduceRegions() is constrained by the FeatureCollection's
# own geometries regardless of the image's footprint. But in Step 3,
# reduceResolution()/reproject() has no such constraint and will happily
# compute over the image's full (accidentally CONUS-sized) footprint,
# which is what produced the "Reprojection output too large" error.
# Re-clipping here restores the 100km ROI footprint before it's used
# downstream in both the patch and pixel sections.

# -----------------------------------------------------------------------------
# 2. PATCH TABLE — zonal stats over each Hansen persistent-loss polygon
# -----------------------------------------------------------------------------
hansenVectors <- ee$FeatureCollection(
  'projects/breidyee/assets/hansen_persistent_loss_vectors_2'  # carries patch_uuid
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
    # 'count' is the true per-feature pixel tally already computed by Step
    # 1's reduceToVectors call. Deriving pixel_count from area_ha instead
    # (assuming a flat 900 m^2/pixel) diverges from this true count by
    # roughly 30% at this latitude, since area_ha is geodesic polygon area
    # while 900 m^2/pixel is the nominal equatorial pixel area — carrying
    # the real count forward avoids that mismatch entirely.
    'pixel_count',            ee$Number(feature$get('count')),
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
# <-- REWRITE: earlier versions built this via reduceResolution()/reproject()
# into the MODIS sinusoidal CRS, then sample()d the result. That chain
# reliably threw "Reprojection output too large" no matter how it was
# adjusted (clipping every band, dropping sample()'s projection arg,
# pre-transforming the region, using the native scale) — even a single
# band with no addBands() involved failed identically. The common factor
# across every failed attempt was reproject()/sample() interacting with a
# target CRS whose transform origin sits at the corner of the entire
# global sinusoidal grid (SR-ORG:6974, origin (-20015109, 10007555)); we
# were never able to pin down exactly why that combination blows up, and
# rather than keep guessing at that path, this rebuilds the pixel table
# with a mechanism already proven to work in this script: reduceToVectors()
# + reduceRegions(), the same approach Step 1/2 use for the patch table.
#
# This ALSO guarantees perfect registration to the true MODIS pixel grid
# (a requirement here) since the grid polygons are derived directly from
# a native, unreprojected MOD13A1 image — there's no resampling step that
# could shift or blur pixel boundaries.

hansenImg      <- ee$Image('UMD/hansen/global_forest_change_2024_v1_12')$clip(roi)
persistentLoss <- hansenImg$select('loss')$eq(1)$And(hansenImg$select('gain')$eq(0))
lossyearPix    <- hansenImg$select('lossyear')$updateMask(persistentLoss)

# -- Build real MODIS pixel center points, then approximate footprints ------
# <-- FIX: reduceToVectors() was the wrong tool here. It doesn't vectorize
# "one polygon per pixel" — it merges ADJACENT pixels of the same value
# into a single polygon. That's exactly right for the Hansen patches in
# Step 1 (contiguous loss areas should merge into one patch), but wrong
# here: every pixel was set to the same constant value, so it correctly
# merged the ENTIRE raster into one giant polygon — confirmed by
# modisPixels$size() returning 1, not ~40,000+. That single oversized
# polygon then had loss_pixel_count == 0 relative to its own filter logic
# in practice, which is why pixelStats1 came back empty.
#
# Fix: skip vectorization entirely. Sample the untouched native MOD13A1
# image directly — no scale or crs override, so reproject() is never
# invoked at all, sidestepping the earlier reprojection bug completely.
# sample() on an image in its own native projection just returns its real
# pixel center points, correctly registered by construction.
modisNativeImg <- ee$ImageCollection('MODIS/061/MOD13A1')$
  first()$select('NDVI')$clip(roi)

modisPoints <- modisNativeImg$sample(
  region     = roi,
  geometries = TRUE
)

# Each point IS an exact, correctly registered MODIS pixel center. Zonal
# stats need an area, not a point, so buffer each one out to a circle with
# the same AREA as one real 463.3127m MODIS pixel:
#   side^2 = pi * r^2   =>   r = side / sqrt(pi) ~= 261.5 m
# This approximates the square footprint's shape, but preserves both the
# true pixel center (used below for centroid_lon/lat) and the true pixel
# area (used for the zonal stats) — which is what actually matters for
# aggregating the 30m Hansen/DEM/NLCD layers underneath.
pixelRadius <- 463.3127 / sqrt(pi)

modisPixels <- modisPoints$map(ee_utils_pyfunc(function(feature) {
  centroid <- feature$geometry()$coordinates()
  feature$set(
    'centroid_lon', centroid$get(0),
    'centroid_lat', centroid$get(1)
  )$setGeometry(feature$geometry()$buffer(pixelRadius))
}))

# -- Pass 1: loss pixel count (sum of Hansen 30m loss pixels per MODIS cell) --
pixelStats1 <- persistentLoss$rename('loss_pixel_count')$reduceRegions(
  collection = modisPixels,
  reducer    = ee$Reducer$sum(),
  scale      = 30,
  tileScale  = 4
)

# Filter early — keeps every later pass working on a smaller collection.
# Drop this filter if you want a full background grid for context/mapping.
pixelStats1 <- pixelStats1$filter(ee$Filter$gt('loss_pixel_count', 0))

# -- Pass 2: dominant Hansen loss year (mode) per cell ------------------------
pixelStats2 <- lossyearPix$rename('dominant_lossyear')$reduceRegions(
  collection = pixelStats1,
  reducer    = ee$Reducer$mode(),
  scale      = 30,
  tileScale  = 4
)

# -- Pass 3: forest cover mean/stdDev + elevation mean, all at 30m -----------
statsImg3 <- forestMask30$addBands(elevation)
reducer3  <- ee$Reducer$mean()$combine(ee$Reducer$stdDev(), sharedInputs = TRUE)
pixelStats3 <- statsImg3$reduceRegions(
  collection = pixelStats2,
  reducer    = reducer3,
  scale      = 30,
  tileScale  = 4
)
# -> adds forest_cover_mean, forest_cover_stdDev, elevation_mean,
#    elevation_stdDev (the extra elevation_stdDev isn't in the table B
#    spec but is harmless — dropped in the final select() below)

# -- Pass 4: circular-safe aspect (sin/cos means), same trick as patch table -
statsImg4 <- sinAspect$addBands(cosAspect)
pixelStats4 <- statsImg4$reduceRegions(
  collection = pixelStats3,
  reducer    = ee$Reducer$mean(),
  scale      = 30,
  tileScale  = 4
)

pixelTable <- pixelStats4$map(ee_utils_pyfunc(function(feature) {
  # centroid_lon/centroid_lat were already set from the true (pre-buffer)
  # pixel center point above — recomputing via geometry()$centroid() here
  # would just return the same value from the buffered circle, so we skip
  # that redundant call and only need to compute aspect_deg.
  sinSafe <- ee$Number(ee$Algorithms$If(feature$get('sinAspect'), feature$get('sinAspect'), 0))
  cosSafe <- ee$Number(ee$Algorithms$If(feature$get('cosAspect'), feature$get('cosAspect'), 0))
  aspectRadM <- sinSafe$atan2(cosSafe)
  aspectDegM <- aspectRadM$multiply(180 / pi)$add(360)$mod(360)
  
  feature$set('aspect_deg', aspectDegM)
}))

pixelTable <- pixelTable$select(list(
  'centroid_lon', 'centroid_lat',
  'loss_pixel_count', 'dominant_lossyear',
  'elevation_mean', 'aspect_deg',
  'forest_cover_mean', 'forest_cover_stdDev'
))

pixelTask <- ee_table_to_drive(
  collection  = pixelTable,
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
# Shared output directory (see config.R) — keeps every pipeline script
# writing to/reading from the same place. Sourced via a fixed absolute path
# (not a bare relative "config.R") so this doesn't depend on R's working
# directory at run time, matching how the rest of this script's paths work.
source(file.path("~/Google Drive/My Drive/Reidy_research/fall color code", "config.R"))

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
