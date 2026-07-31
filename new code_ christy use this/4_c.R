# Step 4c: Regional Reference Phenology Baseline
# -----------------------------------------------------------------------------
# Builds a year-by-year MidGreendown baseline from STABLE, NEVER-DISTURBED,
# high-forest-cover pixels within the ROI. Purpose: patch-level pre/post
# comparisons in Step 5 (Parts D/E) are confounded with ordinary year-to-year
# climate variability -- a warm autumn shifts senescence later for every
# patch, disturbed or not, so a patch's own "pre vs post" comparison can't
# distinguish a real disturbance effect from a year that just happened to
# run late or early everywhere. This reference series lets you subtract that
# common year-effect out: (patch_value - reference_value_that_year) isolates
# what's different about the disturbed patch specifically.
#
# Reference cell selection mirrors the pixel-table method fixed in Step 3:
# sample real MODIS pixel centers directly from the native MOD13A1 image (no
# reduceResolution/reproject — that's the chain that caused repeated
# "Reprojection output too large" errors earlier), then buffer each to an
# equal-area circle. Cells are kept only if they're >=90% covered by stable,
# high-cover, never-disturbed forest, so the reference reflects real mature
# forest phenology rather than any mix of disturbed/regrowth pixels.
#
# Output: reference_phenology_by_year.csv -- one row per year:
#   year, ref_midgreendown_mean, ref_midgreendown_stdDev, ref_midgreendown_count
# -----------------------------------------------------------------------------
# Requires: rgee, googledrive
# =============================================================================

library(rgee)
ee_Initialize(drive = TRUE)

# -----------------------------------------------------------------------------
# 1. Establish ROI (same as all prior steps)
# -----------------------------------------------------------------------------
roi <- ee$Geometry$Point(c(-79.862539, 37.829550))$buffer(100000)

modisProjection <- ee$ImageCollection('MODIS/061/MOD13A1')$
  first()$select('NDVI')$projection()

# -----------------------------------------------------------------------------
# 2. Stable-forest reference mask (30m): never lost, high 2000 tree cover
#    "Never lost" = Hansen lossyear band is exactly 0 (its no-loss sentinel).
#    ">=75% cover in 2000" matches the same threshold used elsewhere in this
#    project (meets_forest_threshold), so the reference is drawn from the
#    same kind of forest the disturbed patches originally were.
# -----------------------------------------------------------------------------
hansenImg     <- ee$Image('UMD/hansen/global_forest_change_2024_v1_12')$clip(roi)
neverLost     <- hansenImg$select('lossyear')$unmask(0)$eq(0)
highCover2000 <- hansenImg$select('treecover2000')$gte(75)
referenceMask30 <- neverLost$And(highCover2000)$rename('stable_forest')

# -----------------------------------------------------------------------------
# 3. Real MODIS pixel centers -> equal-area circles (same method as Step 3's
#    corrected pixel table -- native sample(), no reproject/reduceResolution)
# -----------------------------------------------------------------------------
modisNativeImg <- ee$ImageCollection('MODIS/061/MOD13A1')$
  first()$select('NDVI')$clip(roi)

modisPoints <- modisNativeImg$sample(region = roi, geometries = TRUE)

pixelRadius <- 463.3127 / sqrt(pi)  # equal-area circle for one MODIS pixel

candidateCells <- modisPoints$map(ee_utils_pyfunc(function(feature) {
  feature$setGeometry(feature$geometry()$buffer(pixelRadius))
}))

# -----------------------------------------------------------------------------
# 4. Keep only cells that are overwhelmingly stable, high-cover forest
# -----------------------------------------------------------------------------
cellStats <- referenceMask30$reduceRegions(
  collection = candidateCells,
  reducer    = ee$Reducer$mean(),
  scale      = 30,
  tileScale  = 4
)
# -> adds 'stable_forest' property = fraction of the cell that is stable forest

referenceCells <- cellStats$filter(ee$Filter$gte('stable_forest', 0.9))

# Cap the number of reference cells for a tractable reduceRegion() below --
# 300 cells at 500m each is already a large, spatially spread sample across
# a 100km ROI. Random subsampling (not just $limit(), which would silently
# bias toward however GEE happens to order the collection) keeps it an
# unbiased sample of the qualifying cells.
referenceCells <- referenceCells$randomColumn('rand', seed = 42)$
  sort('rand')$
  limit(300)

referenceGeometry <- referenceCells$geometry()
cat(sprintf('Reference cells selected: %d\n', referenceCells$size()$getInfo()))

# -----------------------------------------------------------------------------
# 5. MCD12Q2 phenology extraction -- same conversion logic as Step 4b
#    (duplicated here rather than sourced, so this script stays self-
#    contained like Step 4b itself is relative to Step 4).
# -----------------------------------------------------------------------------
convertToDoy <- function(image) {
  imgYear   <- ee$Date(image$get('system:time_start'))$get('year')
  yearStart <- ee$Date$fromYMD(imgYear, 1, 1)
  epochDaysAtYearStart <- yearStart$difference(ee$Date('1970-01-01'), 'day')
  
  doyImg <- image$select('MidGreendown_1')$
    subtract(epochDaysAtYearStart)$
    rename('MidGreendown_1_doy')
  
  validMask <- doyImg$gte(1)$And(doyImg$lte(366))
  doyImg$updateMask(validMask)$copyProperties(image, list('system:time_start'))
}

phenologyCollection <- ee$ImageCollection('MODIS/061/MCD12Q2')$
  filterBounds(roi)$
  map(ee_utils_pyfunc(convertToDoy))

# -----------------------------------------------------------------------------
# 6. One aggregate value per year, reduced over ALL reference cells at once
#    (reduceRegion over the union geometry, not reduceRegions per feature --
#    we want a single regional number per year, not per-cell tracking).
# -----------------------------------------------------------------------------
startYear <- 2001
endYear   <- 2022

statsReducer <- ee$Reducer$mean()$
  combine(reducer2 = ee$Reducer$stdDev(), sharedInputs = TRUE)$
  combine(reducer2 = ee$Reducer$count(), sharedInputs = TRUE)

yearlyReference <- function(yr) {
  yrImage <- phenologyCollection$
    filter(ee$Filter$calendarRange(yr, yr, 'year'))$
    first()
  
  stats <- yrImage$select('MidGreendown_1_doy')$reduceRegion(
    reducer   = statsReducer,
    geometry  = referenceGeometry,
    scale     = 463.3127,
    maxPixels = 1e13,
    bestEffort = TRUE
  )
  
  ee$Feature(NULL, list(
    year                    = yr,
    ref_midgreendown_mean   = stats$get('MidGreendown_1_doy_mean'),
    ref_midgreendown_stdDev = stats$get('MidGreendown_1_doy_stdDev'),
    ref_midgreendown_count  = stats$get('MidGreendown_1_doy_count')
  ))
}

allYearsReference <- ee$FeatureCollection(lapply(startYear:endYear, yearlyReference))

# -----------------------------------------------------------------------------
# 7. Export + download
# -----------------------------------------------------------------------------
referenceColumns <- c('year', 'ref_midgreendown_mean',
                      'ref_midgreendown_stdDev', 'ref_midgreendown_count')

task <- ee_table_to_drive(
  collection  = allYearsReference,
  description = 'reference_phenology_by_year',
  folder      = 'Reidy_research',
  fileFormat  = 'CSV',
  selectors   = referenceColumns,
  timePrefix  = FALSE
)
task$start()
cat('Regional reference phenology export started: reference_phenology_by_year\n')
ee_monitoring(task)

library(googledrive)
source(file.path("~/Google Drive/My Drive/Reidy_research/fall color code", "config.R"))

Sys.sleep(15)
refFile <- drive_ls(path = 'Reidy_research', pattern = 'reference_phenology_by_year')
if (nrow(refFile) >= 1) {
  drive_download(
    file      = refFile[1, ],
    path      = file.path(outputDir, 'reference_phenology_by_year.csv'),
    overwrite = TRUE
  )
  cat('Regional reference phenology series saved locally.\n')
}
