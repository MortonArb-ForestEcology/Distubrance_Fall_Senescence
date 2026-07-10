# Step 4: NDVI Yearly Series — One Value Per Year, Full MODIS Record
# -----------------------------------------------------------------------------
# Design change from the original Step 4:
#   - OLD: monthly windows, relativized to each patch's disturbance year
#     (months_from_disturbance), one export task per loss year.
#   - NEW: one NDVI value per CALENDAR year (mean of cloud-free 16-day
#     composites falling in JULY), for every patch, across the FULL MODIS
#     record (2000-2024) — regardless of that patch's own loss_year.
#
# WHY:
#   - July is used because peak growing-season NDVI is the most phenologically
#     stable window (least sensitive to green-up/senescence timing), so a
#     single July value per year is a defensible "one value per year" summary.
#   - Pulling the FULL record (not just pre/post relative to disturbance)
#     means every patch has enough years of data to tell whether any single
#     year's NDVI was anomalously high/low relative to that patch's own
#     multi-year baseline, and lets you compare against pixels/patches that
#     were never disturbed at all. Relativizing to loss_year, and deciding
#     how many pre/post years to actually use, happens downstream — this
#     script just gets everything into one long table so you filter later.
#
# Output: one long-format CSV — one row per (patch_uuid, year).
#   patch_uuid, loss_year, year, ndvi_july_mean, ndvi_july_stdDev, ndvi_july_count
# -----------------------------------------------------------------------------
# Requires: rgee installed and authenticated (ee_Initialize() must succeed)
# =============================================================================

library(rgee)

ee_Initialize(drive = TRUE)

# -----------------------------------------------------------------------------
# 1. Establish ROI
# -----------------------------------------------------------------------------
roi <- ee$Geometry$Point(c(-79.862539, 37.829550))$buffer(100000)

# -----------------------------------------------------------------------------
# 2. Load Step 1 asset — canonical patches, carrying patch_uuid + loss_year
#    Deliberately NOT the Step 3 forest-masked raster: we want NDVI for every
#    patch regardless of the forest-cover threshold, so the threshold flag
#    (from the patch attribute table) can be applied as a filter afterward.
# -----------------------------------------------------------------------------
hansenVectors <- ee$FeatureCollection(
  'projects/breidyee/assets/hansen_persistent_loss_vectors'
)$select(list('patch_uuid', 'loss_year'))

# -----------------------------------------------------------------------------
# 3. MODIS MOD13A1 NDVI, QA-masked
#    SummaryQA: 0 = good, 1 = marginal (both kept); 2 = snow/ice, 3 = cloudy
#    (masked out).
# -----------------------------------------------------------------------------
modisProjection <- ee$ImageCollection('MODIS/061/MOD13A1')$
  first()$select('NDVI')$projection()

applyQaMask <- ee_utils_pyfunc(function(image) {
  qa   <- image$select('SummaryQA')
  mask <- qa$lte(1)
  image$updateMask(mask)$
    select('NDVI')$
    multiply(0.0001)$
    copyProperties(image, list('system:time_start'))
})

modisNDVI <- ee$ImageCollection('MODIS/061/MOD13A1')$
  filterBounds(roi)$
  map(applyQaMask)

# -----------------------------------------------------------------------------
# 4. Build one year's July-mean NDVI, reduced over every patch
#    mean/stdDev/count per patch per year — count lets you see, downstream,
#    how many cloud-free composites actually went into a given patch-year
#    (useful for flagging low-confidence years).
# -----------------------------------------------------------------------------
startYear <- 2000   # MOD13A1 begins Feb 2000; July is fully covered from 2000
endYear   <- 2024   # matches the Hansen 2024 vintage used in Step 1

yearlyJulyNDVI <- function(yr) {
  julyMean <- modisNDVI$
    filter(ee$Filter$calendarRange(yr, yr, 'year'))$
    filter(ee$Filter$calendarRange(7, 7, 'month'))$
    mean()$
    rename('ndvi_july')
  
  reducer <- ee$Reducer$mean()$
    combine(ee$Reducer$stdDev(), sharedInputs = TRUE)$
    combine(ee$Reducer$count(),  sharedInputs = TRUE)
  
  fc <- julyMean$reduceRegions(
    collection = hansenVectors,
    reducer    = reducer,
    scale      = 500,
    tileScale  = 4
  )
  
 
  fc$map(ee_utils_pyfunc(function(feature) {
    ndviMean  <- ee$Algorithms$If(feature$get('ndvi_july_mean'),   feature$get('ndvi_july_mean'),   feature$get('mean'))
    ndviSd    <- ee$Algorithms$If(feature$get('ndvi_july_stdDev'), feature$get('ndvi_july_stdDev'), feature$get('stdDev'))
    ndviCount <- ee$Algorithms$If(feature$get('ndvi_july_count'),  feature$get('ndvi_july_count'),  feature$get('count'))
    
    feature$set(
      'year',              yr,
      'ndvi_july_mean',    ndviMean,
      'ndvi_july_stdDev',  ndviSd,
      'ndvi_july_count',   ndviCount
    )
  }))
}

allYearsNDVI <- ee$FeatureCollection(
  lapply(startYear:endYear, yearlyJulyNDVI)
)$flatten()

allYearsNDVI <- allYearsNDVI$select(list(
  'patch_uuid', 'loss_year', 'year',
  'ndvi_july_mean', 'ndvi_july_stdDev', 'ndvi_july_count'
))

# -----------------------------------------------------------------------------
# 5. Export as one long CSV (single batch task — no per-year download loop
#    needed since we're not joining anything client-side)
# -----------------------------------------------------------------------------
task <- ee_table_to_drive(
  collection  = allYearsNDVI,
  description = 'ndvi_yearly_july_by_patch',
  folder      = 'Reidy_research',
  fileFormat  = 'CSV',
  timePrefix  = FALSE
)
task$start()
cat('NDVI yearly (July-mean) export started: ndvi_yearly_july_by_patch\n')
ee_monitoring(task)

# -----------------------------------------------------------------------------
# 6. Pull it down locally
# -----------------------------------------------------------------------------
library(googledrive)

# Set this to wherever you want local files saved on YOUR machine
outputDir <- '~/Reidy_research'
if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)

Sys.sleep(15)
ndviFile <- drive_ls(path = 'Reidy_research', pattern = 'ndvi_yearly_july_by_patch')
if (nrow(ndviFile) >= 1) {
  drive_download(
    file      = ndviFile[1, ],
    path      = file.path(outputDir, 'ndvi_yearly_july_by_patch.csv'),
    overwrite = TRUE
  )
  cat('NDVI yearly series saved locally.\n')
}
