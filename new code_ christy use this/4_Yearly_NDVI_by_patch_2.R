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
# MODIS-resolvable patch filter: one MODIS pixel footprint is 500m x 500m =
# 25 ha. Patches smaller than that are overwhelmingly diluted by surrounding
# (undisturbed) forest at 500m reduceRegions scale, so NDVI "recovery"
# measured on them mostly reflects the neighborhood rather than the patch
# itself (median patch here is ~0.25 ha — three orders of magnitude below
# one MODIS pixel). Restricting to patches >= 25 ha keeps this script's
# output interpretable as patch-level canopy NDVI. This filter is local to
# this script only — Steps 1-3's full patch catalog (all sizes) is
# untouched, so mapping/attribute outputs still cover every patch.
MODIS_PIXEL_AREA_HA <- 25

hansenVectors <- ee$FeatureCollection(
  'projects/breidyee/assets/hansen_persistent_loss_vectors_2'
)$select(list('patch_uuid', 'loss_year', 'area_ha'))$
  filter(ee$Filter$gte('area_ha', MODIS_PIXEL_AREA_HA))$
  select(list('patch_uuid', 'loss_year'))  # drop area_ha once it's done its job as
                                            # a filter — no need to carry it through
                                            # 25 years of reduceRegions chaining

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
    mean()
  
  # ROOT CAUSE OF THE PREVIOUS TWO FAILED ATTEMPTS: when reduceRegions() is
  # given a SINGLE band + a SIMPLE (single-output) reducer, EE drops any
  # custom band name and just uses the reducer's GENERIC output name
  # ('mean', 'stdDev', 'count') for the property — renaming the band via
  # $rename() has no effect on that. That's why both prior fixes produced
  # columns under names we weren't asking for, so select() found nothing
  # and silently dropped them.
  #
  # The reliable fix: Reducer.setOutputs() explicitly forces the output
  # property name, overriding EE's automatic (and here, unhelpful) naming
  # convention entirely — no guessing required.
  # mean + stdDev combined into one reducer/one reduceRegions pass (same
  # pattern already used in Steps 2 and 3) instead of two separate passes —
  # cuts this loop from 3 reduceRegions calls/year to 2. setOutputs() above
  # is still required, for the same reason given in the note above this
  # function.
  meanSdReducer <- ee$Reducer$mean()$combine(
    reducer2     = ee$Reducer$stdDev(),
    sharedInputs = TRUE
  )$setOutputs(list('ndvi_july_mean', 'ndvi_july_stdDev'))
  countReducer  <- ee$Reducer$count()$setOutputs(list('ndvi_july_count'))

  fc <- julyMean$reduceRegions(
    collection = hansenVectors, reducer = meanSdReducer, scale = 500, tileScale = 4
  )
  fc <- julyMean$reduceRegions(
    collection = fc, reducer = countReducer, scale = 500, tileScale = 4
  )
  
  # NOTE: ndvi_july_stdDev/mean can be genuinely UNDEFINED (not just null —
  # the key itself absent) on features where zero valid (non-cloud-masked)
  # pixels intersect the patch at 500m that year. Reducer.count() never has
  # this problem (0 is always a defined answer), which is why count alone
  # kept showing up in exports while mean/stdDev vanished — GEE's CSV
  # exporter infers its column headers from a feature's actual property
  # keys, and if the sampled feature is missing a key, that whole column
  # gets dropped from the export, even for other rows that DO have values.
  #
  # Fix: explicitly re-set each property to itself. get() on a missing
  # property returns null, and set() with that null value still creates
  # the key on the feature — so every feature is now GUARANTEED to carry
  # all three keys (possibly null), and no feature can cause the column to
  # disappear from the export.
  fc$map(ee_utils_pyfunc(function(feature) {
    feature$set(
      'year',              yr,
      'ndvi_july_mean',    feature$get('ndvi_july_mean'),
      'ndvi_july_stdDev',  feature$get('ndvi_july_stdDev'),
      'ndvi_july_count',   feature$get('ndvi_july_count')
    )
  }))
}

allYearsNDVI <- ee$FeatureCollection(
  lapply(startYear:endYear, yearlyJulyNDVI)
)$flatten()

# -----------------------------------------------------------------------------
# 5. Export as one long CSV (single batch task — no per-year download loop
#    needed since we're not joining anything client-side)
#
#    IMPORTANT: columns are specified via `selectors` on the export call
#    itself, NOT via collection$select() beforehand. select() + the
#    exporter's automatic header-inference is what was silently dropping
#    ndvi_july_mean/stdDev in earlier attempts. `selectors` builds the CSV
#    header from this explicit list directly, independent of which
#    properties any individual feature happens to carry.
# -----------------------------------------------------------------------------
ndviColumns <- c('patch_uuid', 'loss_year', 'year',
                 'ndvi_july_mean', 'ndvi_july_stdDev', 'ndvi_july_count')

task <- ee_table_to_drive(
  collection  = allYearsNDVI,
  description = 'ndvi_yearly_july_by_patch',
  folder      = 'Reidy_research',
  fileFormat  = 'CSV',
  selectors   = ndviColumns,
  timePrefix  = FALSE
)
task$start()
cat('NDVI yearly (July-mean) export started: ndvi_yearly_july_by_patch\n')
ee_monitoring(task)

# -----------------------------------------------------------------------------
# 6. Pull it down locally
# -----------------------------------------------------------------------------
library(googledrive)

# Shared output directory (see config.R) — keeps every pipeline script
# writing to/reading from the same place. Sourced via a fixed absolute path
# (not a bare relative "config.R") so this doesn't depend on R's working
# directory at run time, matching how the rest of this script's paths work.
source(file.path("~/Google Drive/My Drive/Reidy_research/fall color code", "config.R"))

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
