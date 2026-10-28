# Step 4: NDVI Time Series Analysis — Pre/Post Disturbance
# -----------------------------------------------------------------------------
# Order of operations:
#   1. Establish ROI
#   2. Load Step 3 asset (Hansen loss masked to >= 75% forest cover)
#   3. Load MODIS MOD13A1 NDVI collection
#   4. Client-side loop over loss years 2001-2024:
#      - Skip years with no disturbed pixels
#      - Extract monthly NDVI for 1 year pre and 3 years post disturbance
#      - Submit one export task per year
#   5. Monitor all tasks
#   6. Download and combine CSVs into one file, upload to Drive
# -----------------------------------------------------------------------------
# MODIS MOD13A1: 500m 16-day NDVI composites. Monthly values are computed
# as the mean of all 16-day composites falling within each calendar month.
# NDVI is scaled by 0.0001 to convert from raw integer to 0-1 range.
# Requires: rgee installed and authenticated (ee_Initialize() must succeed)
# =============================================================================

library(rgee)
library(purrr)
library(googledrive)

ee_Initialize(drive = TRUE)

# -----------------------------------------------------------------------------
# 1. Establish ROI
# -----------------------------------------------------------------------------
roi <- ee$Geometry$Point(c(-79.862539, 37.829550))$buffer(100000)

# -----------------------------------------------------------------------------
# 2. Load Step 3 asset — Hansen loss masked to >= 75% forest cover
#    Raster at 30m where pixel values = Hansen lossyear (1-24).
# -----------------------------------------------------------------------------
hansenMasked <- ee$Image('projects/breidyee/assets/hansen_loss_forest_masked')

# -----------------------------------------------------------------------------
# 3. Load MODIS MOD13A1 NDVI
#    16-day composites at 500m. NDVI scaled by 0.0001.
#    Quality mask applied using SummaryQA band:
#      0 = good data, 1 = marginal data — both retained.
#      2 = snow/ice, 3 = cloudy — masked out.
# -----------------------------------------------------------------------------
modisProjection <- ee$ImageCollection('MODIS/061/MOD13A1')$
  first()$
  select('NDVI')$
  projection()

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
# 4. Client-side loop over loss years
#    For each loss year (Hansen encoding 1-24 = calendar year 2001-2024):
#      - Check if any disturbed pixels exist for that year
#      - If yes, build monthly NDVI time series for pre and post windows
#      - Submit one export task per year to Google Drive
# -----------------------------------------------------------------------------
tasks <- list()

for (rawYear in 1:24) {
  
  calYear <- rawYear + 2000
  
  # -- Check if this loss year has any pixels in the masked asset ------------
  pixelCount <- hansenMasked$eq(rawYear)$
    reduceRegion(
      reducer   = ee$Reducer$sum(),
      geometry  = roi,
      scale     = 500,
      maxPixels = 1e10
    )$getInfo()$lossyear
  
  if (is.null(pixelCount) || pixelCount == 0) {
    cat(sprintf('  Year %d: no pixels, skipping\n', calYear))
    next
  }
  
  cat(sprintf('  Year %d: %.0f pixels found\n', calYear, pixelCount))
  
  # Mask Hansen raster to this loss year only
  yearMask <- hansenMasked$eq(rawYear)$selfMask()
  
  # -- Date windows -----------------------------------------------------------
  preStart  <- sprintf('%d-01-01', calYear - 1)
  preEnd    <- sprintf('%d-12-31', calYear - 1)
  postStart <- sprintf('%d-01-01', calYear + 1)
  postEnd   <- sprintf('%d-12-31', calYear + 3)
  
  # -- Helper: build monthly NDVI feature collection for a date window --------
  #    Each feature = one MODIS pixel for one month with NDVI and metadata.
  buildMonthly <- function(startDate, endDate, period) {
    
    dates <- seq(as.Date(startDate), as.Date(endDate), by = 'month')
    
    ee$FeatureCollection(
      ee$List(
        lapply(dates, function(d) {
          ms <- format(d, '%Y-%m-%d')
          # Last day of month: advance to next month, subtract 1 day
          me <- format(
            as.Date(format(d, '%Y-%m-01')) +
              32 - as.integer(format(
                as.Date(format(d, '%Y-%m-01')) + 32, '%d'
              )),
            '%Y-%m-%d'
          )
          
          # Monthly mean NDVI masked to disturbed pixels for this year
          monthImg <- modisNDVI$
            filterDate(ms, me)$
            mean()$
            updateMask(yearMask)$
            rename('ndvi')
          
          # months_from_disturbance: negative = pre, positive = post
          # 0 = January of the disturbance year
          monthsSinceDist <- as.integer(format(d, '%Y')) * 12L +
            as.integer(format(d, '%m')) -
            (calYear * 12L + 1L)
          
          yr <- as.integer(format(d, '%Y'))
          mo <- as.integer(format(d, '%m'))
          
          sampled <- monthImg$sample(
            region     = roi,
            scale      = 500,
            projection = modisProjection,
            geometries = TRUE
          )
          
          sampled$map(ee_utils_pyfunc(function(feature) {
            feature$set(
              'longitude',               feature$geometry()$coordinates()$get(0),
              'latitude',                feature$geometry()$coordinates()$get(1),
              'loss_year',               calYear,
              'year',                    yr,
              'month',                   mo,
              'period',                  period,
              'months_from_disturbance', monthsSinceDist
            )
          }))
        })
      )
    )$flatten()
  }
  
  # Combine pre and post windows into one feature collection
  yearSeries <- ee$FeatureCollection(list(
    buildMonthly(preStart, preEnd,   'pre'),
    buildMonthly(postStart, postEnd, 'post')
  ))$flatten()
  
  # -- Submit export task for this year ---------------------------------------
  desc <- sprintf('ndvi_time_series_%d', calYear)
  
  task <- ee_table_to_drive(
    collection  = yearSeries,
    description = desc,
    folder      = 'Reidy_research',
    fileFormat  = 'CSV',
    timePrefix  = FALSE
  )
  task$start()
  cat(sprintf('  Started task: %s\n', desc))
  tasks[[as.character(calYear)]] <- task
}

# -----------------------------------------------------------------------------
# 5. Wait for all tasks to complete
#    Polls Drive every 60 seconds until the number of NDVI CSV files
#    stabilises — handles years with no output gracefully.
# -----------------------------------------------------------------------------
cat('Waiting for all tasks to complete...\n')

nExpected <- length(tasks)
prev_count <- -1L
stable_checks <- 0L

repeat {
  found <- nrow(drive_ls(path = 'Reidy_research', pattern = 'ndvi_time_series_'))
  cat(sprintf('[%s] Files in Drive: %d / %d\n', Sys.time(), found, nExpected))
  if (found == prev_count) {
    stable_checks <- stable_checks + 1L
  } else {
    stable_checks <- 0L
  }
  if (found >= nExpected || stable_checks >= 5L) break
  prev_count <- found
  Sys.sleep(60)
}

cat('All tasks complete. Downloading and combining...\n')

# -----------------------------------------------------------------------------
# 6. Download all per-year CSVs and combine into one file
# -----------------------------------------------------------------------------
drive_files <- drive_ls(path = 'Reidy_research', pattern = 'ndvi_time_series_')
cat(sprintf('Files found in Drive: %d\n', nrow(drive_files)))

year_dfs <- lapply(seq_len(nrow(drive_files)), function(i) {
  fname    <- drive_files$name[i]
  out_file <- file.path(tempdir(), fname)
  tryCatch({
    drive_download(file = drive_files[i, ], path = out_file, overwrite = TRUE)
    df <- read.csv(out_file)
    cat(sprintf('  %s: %d rows\n', fname, nrow(df)))
    df
  }, error = function(e) {
    warning(sprintf('  %s failed: %s', fname, e$message))
    NULL
  })
})

year_dfs <- Filter(Negate(is.null), year_dfs)
combined  <- do.call(rbind, year_dfs)
combined  <- combined[!duplicated(combined[, c('longitude', 'latitude',
                                               'loss_year', 'year', 'month')]), ]

cat(sprintf('Combined dataset: %d rows\n', nrow(combined)))

# -----------------------------------------------------------------------------
# 7. Upload combined CSV to Drive
# -----------------------------------------------------------------------------
tmp_path <- file.path(tempdir(), 'ndvi_time_series_combined.csv')
write.csv(combined, tmp_path, row.names = FALSE)

drive_upload(
  media     = tmp_path,
  path      = 'Reidy_research/',
  name      = 'ndvi_time_series_combined.csv',
  overwrite = TRUE
)

cat('Uploaded combined CSV to Google Drive: Reidy_research/ndvi_time_series_combined.csv\n')
