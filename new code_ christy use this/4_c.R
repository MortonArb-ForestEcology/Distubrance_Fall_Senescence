# Step 4c: Regional Reference Phenology Baseline
# -----------------------------------------------------------------------------
# Builds a year-by-year MidGreendown baseline from STABLE, NEVER-DISTURBED,
# high-forest-cover pixels within the ROI.
#
# PURPOSE: patch-level pre/post comparisons in Step 5 (Parts D/E) are
# confounded with ordinary year-to-year climate variability -- a warm autumn
# shifts senescence later for every patch, disturbed or not, so a patch's own
# "pre vs post" comparison can't by itself distinguish a real disturbance
# effect from a year that happened to run late or early everywhere. This
# reference series lets Step 5 subtract that common year-effect out:
# (patch_value - reference_value_that_year) isolates what's different about
# the disturbed patch specifically.
#
# SCOPE NOTE: Step 5's Part G mixed model -- response ~ post +
# (1|patch_uuid) + (1|year) -- already controls for common year effects
# using the 615-patch panel itself, with no separate reference sample and
# no sampling noise, by exploiting the fact that disturbance years are
# staggered across patches. This script is therefore best treated as an
# INDEPENDENT CROSS-CHECK on that model rather than the primary climate
# control. If the two disagree, weight the mixed model.
#
# METHOD: reference cell selection mirrors the pixel-table approach that
# was fixed in Step 3 -- sample real MODIS pixel centers directly from the
# native MOD13A1 image (no reduceResolution/reproject, which is the chain
# that caused repeated "Reprojection output too large" failures), then
# buffer each to an equal-area circle. Cells are kept only if they're >=90%
# covered by stable, high-cover, never-disturbed forest, so the reference
# reflects real mature forest phenology rather than a mix of disturbed and
# regrowth pixels.
#
# Output: reference_phenology_by_year.csv -- one row per year:
#   year, ref_midgreendown_mean, ref_midgreendown_stdDev, ref_midgreendown_count
# -----------------------------------------------------------------------------
# Requires: rgee, googledrive
# =============================================================================

library(rgee)
library(googledrive)

ee_Initialize(drive = TRUE)

# Shared output directory. Defined directly rather than sourced from a
# config.R -- the config file referenced by earlier drafts of this pipeline
# doesn't exist at that path, so every source() call failed and left
# outputDir undefined, which then broke the download at the end of the
# script. Set this once, here.
outputDir <- file.path("~/Google Drive/My Drive", "Reidy_research")
if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)

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
# 3. Real MODIS pixel centers -> equal-area circles, then SUBSAMPLE
# -----------------------------------------------------------------------------
modisNativeImg <- ee$ImageCollection('MODIS/061/MOD13A1')$
  first()$select('NDVI')$clip(roi)

# Sample a RANDOM SUBSET of MODIS pixel centers directly, via numPixels --
# do NOT sample all of them and subset afterward.
#
# Why this matters: sampling the full ROI returns ~144,600 points. Taking a
# subsample afterward with randomColumn()$sort()$limit() does NOT save any
# work, because sort() is a blocking operation -- GEE must materialize all
# 144,600 points AND buffer every one of them before it can rank them. That
# full materialization is the actual bottleneck, and it exceeds GEE's
# ~5-minute interactive limit ("Computation timed out") regardless of how
# few cells survive the subsequent limit().
#
# numPixels pushes the subsampling into the sample() call itself, so the
# expensive full-ROI enumeration never happens. seed makes it reproducible.
# No scale/projection override here, so sampling stays on MOD13A1's native
# grid -- the same reason Step 3's corrected pixel table avoids reproject().
CANDIDATE_POOL_SIZE <- 3000L

modisPoints <- modisNativeImg$sample(
  region     = roi,
  numPixels  = CANDIDATE_POOL_SIZE,
  seed       = 42,
  geometries = TRUE
)

pixelRadius <- 463.3127 / sqrt(pi)  # equal-area circle for one MODIS pixel

candidatePool <- modisPoints$map(ee_utils_pyfunc(function(feature) {
  feature$setGeometry(feature$geometry()$buffer(pixelRadius))
}))
# NOTE: this does not bias the reference sample. The original design also
# ended at a random limit(300); a random draw that is then filtered on
# stable-forest cover yields a random draw of the qualifying population
# either way. 3000 is generous headroom given only 300 need to survive.

# -----------------------------------------------------------------------------
# 4. Keep only cells that are overwhelmingly stable, high-cover forest
# -----------------------------------------------------------------------------
cellStats <- referenceMask30$reduceRegions(
  collection = candidatePool,
  reducer    = ee$Reducer$mean(),
  scale      = 30,
  tileScale  = 4
)

# NOTE the property name: 'mean', NOT 'stable_forest'. When reduceRegions()
# runs on a SINGLE-BAND image with a SINGLE-OUTPUT reducer, GEE names the
# output property after the REDUCER, not the band -- $rename() on the image
# has no effect on it. Only multi-band inputs get band-named outputs.
# Filtering on 'stable_forest' (as earlier drafts did) matches zero
# features, which silently cascades: referenceCells empty -> reference
# geometry empty -> every yearly reduceRegion below returns null with
# count 0 -> the exported CSV is completely blank -> Step 5's Part F #1
# renders two empty figure panels, with no error anywhere in the chain.
# Confirmed by running reduceRegions on 5 test cells and inspecting
# propertyNames(): "mean", "NDVI", "system:index".

# Cap at 300 cells for a tractable reduceRegion in section 6. No further
# randomization needed -- modisPoints was already a seeded random draw, so
# survivors of the filter are already an unbiased random sample.
qualifyingCells   <- cellStats$filter(ee$Filter$gte('mean', 0.9))
referenceCells    <- qualifyingCells$limit(300)
referenceGeometry <- referenceCells$geometry()

# Report both counts. n_qualifying is the diagnostic that matters: if it's
# under ~50, the reference baseline is too thin to trust -- loosen the 0.9
# threshold or raise CANDIDATE_POOL_SIZE rather than pushing ahead.
n_qualifying <- qualifyingCells$size()$getInfo()
n_used       <- referenceCells$size()$getInfo()

cat(sprintf('Stable-forest cells passing >=0.9 filter (of %d sampled): %d\n',
            CANDIDATE_POOL_SIZE, n_qualifying))
cat(sprintf('Reference cells actually used: %d\n', n_used))

if (n_qualifying < 50) {
  warning('Only ', n_qualifying, ' cells passed the 0.9 stable-forest threshold. ',
          'Loosen the threshold or raise CANDIDATE_POOL_SIZE before trusting ',
          'this reference baseline.')
}

# -----------------------------------------------------------------------------
# 5. MCD12Q2 phenology extraction -- same conversion logic as Step 4b
#    (duplicated here rather than sourced, so this script stays self-
#    contained, the same way Step 4b is self-contained relative to Step 4).
#
#    MCD12Q2 date bands are days-since-1970-01-01 (epoch days), NOT
#    day-of-year -- treating the raw value as a DOY silently produces
#    nonsense. Convert per image using that image's own year, then mask to
#    a plausible [1, 366] range, which also removes the large fill values
#    carried by pixels with no detected phenological cycle that year.
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
#
#    Unlike reduceRegions, reduceRegion on a single band with a COMBINED
#    reducer DOES use band_reducer naming, so 'MidGreendown_1_doy_mean' etc.
#    are the correct keys here.
# -----------------------------------------------------------------------------
startYear <- 2001   # MCD12Q2 begins with the 2001 growing season
endYear   <- 2022   # most recent finalized vintage; matches Step 4b

statsReducer <- ee$Reducer$mean()$
  combine(reducer2 = ee$Reducer$stdDev(), sharedInputs = TRUE)$
  combine(reducer2 = ee$Reducer$count(), sharedInputs = TRUE)

yearlyReference <- function(yr) {
  yrImage <- phenologyCollection$
    filter(ee$Filter$calendarRange(yr, yr, 'year'))$
    first()
  
  stats <- yrImage$select('MidGreendown_1_doy')$reduceRegion(
    reducer    = statsReducer,
    geometry   = referenceGeometry,
    scale      = 463.3127,
    maxPixels  = 1e13,
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

# -----------------------------------------------------------------------------
# 8. VERIFY the export actually contains data before moving on.
#    The original failure mode here was a completely blank CSV that produced
#    no error at any stage -- the export "succeeded", the download
#    "succeeded", and Step 5 silently drew two empty figures. Check
#    explicitly rather than assuming.
# -----------------------------------------------------------------------------
ref_check <- read.csv(file.path(outputDir, 'reference_phenology_by_year.csv'))
cat('\nExported reference series:\n')
print(ref_check)

n_valid_years <- sum(!is.na(ref_check$ref_midgreendown_mean))
if (n_valid_years == 0) {
  stop('reference_phenology_by_year.csv is EMPTY (all years NA). Do NOT run ',
       "Step 5's Part F #1 against this -- it will silently produce blank ",
       'figures. Re-check the reference-cell counts printed in section 4.')
} else {
  cat(sprintf('\n%d of %d years have a valid reference mean.\n',
              n_valid_years, nrow(ref_check)))
}
