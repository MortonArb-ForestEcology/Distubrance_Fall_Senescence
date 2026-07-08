# Step 1: Hansen Persistent Loss Asset (revised — adds patch_uuid, single export)
# -----------------------------------------------------------------------------
# Order of operations:
#   1. Establish ROI
#   2. Load Hansen 2024, clip to ROI
#   3. Apply persistent loss filter (loss == 1 AND gain == 0)
#   4. Vectorize at native 30m resolution using labelProperty
#   5. Add area property, filter small polygons, bake in patch_uuid
#   6. Save as GEE asset — ONE export task, same as the original script
#   7. Optional: pull a local copy down afterward for .gpkg/.csv output
#   8. Optional: visualize in interactive map
# -----------------------------------------------------------------------------
# WHY THIS VERSION IS DIFFERENT FROM THE PRIOR "UUID ROUND TRIP" VERSION:
# GEE has no random-UUID generator, but it DOES give every feature a unique,
# permanent 'system:index' the moment reduceToVectors creates it. Instead of
# pulling the collection into R just to stamp a real UUID and re-uploading
# (which meant a slow drive round trip AND a chunked re-upload to dodge the
# 10MB payload limit), we just copy that existing unique id into a permanent
# 'patch_uuid' property server-side, before ever leaving GEE. That means:
#   - ONE export task, exactly like the original script — no chunking,
#     no per-chunk task-queue overhead, no re-upload.
#   - patch_uuid is still unique and permanent for this asset, which is all
#     Steps 3 and 4 actually need it for (a stable join key).
# Trade-off: it's a deterministic id tied to this computation/asset lineage,
# not a globally-random UUID — irrelevant for a within-project join key.
# Requires: rgee installed and authenticated (ee_Initialize() must succeed)
# =============================================================================

library(rgee)

ee_Initialize(drive = TRUE)

# -----------------------------------------------------------------------------
# 1. Establish ROI
# -----------------------------------------------------------------------------
roi <- ee$Geometry$Point(c(-79.862539, 37.829550))$buffer(100000)

# -----------------------------------------------------------------------------
# 2. Load Hansen Global Forest Change 2024, clip to ROI
# -----------------------------------------------------------------------------
hansen <- ee$Image('UMD/hansen/global_forest_change_2024_v1_12')$clip(roi)

loss     <- hansen$select('loss')
gain     <- hansen$select('gain')
lossyear <- hansen$select('lossyear')

# -----------------------------------------------------------------------------
# 3. Apply persistent loss filter
# -----------------------------------------------------------------------------
persistentLoss <- loss$eq(1)$And(gain$eq(0))
lossYearMasked <- lossyear$updateMask(persistentLoss)

# -----------------------------------------------------------------------------
# 4. Vectorize at native 30m resolution
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
# 5. Add area (ha), calendar year, and patch_uuid — all server-side, in the
#    same map() call the original script already used for area/year.
#    feature$id() returns the feature's system:index, assigned uniquely by
#    reduceToVectors — we just copy it into a named, permanent property.
# -----------------------------------------------------------------------------
lossVectors <- lossVectors$map(ee_utils_pyfunc(function(feature) {
  feature$set(
    'area_ha',    feature$geometry()$area(maxError = 1)$divide(10000),
    'loss_year',  ee$Number(feature$get('lossyear'))$add(2000),
    'patch_uuid', ee$String('patch_')$cat(feature$id())
  )
}))

lossVectors <- lossVectors$filter(ee$Filter$gt('area_ha', 0.09))

# -----------------------------------------------------------------------------
# 6. Save as GEE asset — single export task, same as the original script
# -----------------------------------------------------------------------------
assetId <- 'projects/breidyee/assets/hansen_persistent_loss_vectors'

task <- ee$batch$Export$table$toAsset(
  collection  = lossVectors,
  description = 'hansen_persistent_loss_vectors',
  assetId     = assetId
)

task$start()
ee_monitoring(task)

# -----------------------------------------------------------------------------
# 7. Optional: pull a local copy down AFTER the asset exists, purely for
#    convenience (.gpkg for GIS, .csv for a flat attribute table). This is a
#    one-way download — nothing gets re-uploaded, so none of the payload/
#    chunking concerns from the round-trip version apply here.
# -----------------------------------------------------------------------------
library(sf)

# Set this to wherever you want local files saved on YOUR machine — it will
# be created if it doesn't already exist.
outputDir <- '~/Reidy_research'
if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)

lossVectors_ee <- ee$FeatureCollection(assetId)

lossVectors_sf <- ee_as_sf(
  lossVectors_ee,
  via         = 'drive',
  maxFeatures = 100000
)

st_write(lossVectors_sf, file.path(outputDir, 'hansen_persistent_loss_patches.gpkg'),
         delete_dsn = TRUE)

centroids <- st_coordinates(st_centroid(st_geometry(lossVectors_sf)))
lossVectors_csv <- st_drop_geometry(lossVectors_sf)
lossVectors_csv$centroid_lon <- centroids[, 1]
lossVectors_csv$centroid_lat <- centroids[, 2]

write.csv(
  lossVectors_csv,
  file.path(outputDir, 'hansen_persistent_loss_patches.csv'),
  row.names = FALSE)

# -----------------------------------------------------------------------------
# 8. Optional: visualize in interactive map
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
