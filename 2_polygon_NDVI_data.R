# =============================================================================
# Script 3: NDVI and Full Phenometrics Extraction via rgee
# Purpose: For each polygon in floss500_dat.csv, extract annual MODIS NDVI
#          (as disturbance severity check only — NOT a productivity measure)
#          and all 8 MCD12Q2 phenometrics (DOY and Date) for years 2001-2023.
#
# NDVI NOTE: NDVI is extracted here solely to verify that Hansen loss polygons
#            correspond to detectable canopy disturbance in MODIS imagery.
#            NDVI is NOT used as a productivity estimate. NDVI ≠ NPP.
#
# Phenometrics extracted from MCD12Q2:
#   Greenup, MidGreenup, Maturity, Peak, MidGreendown, Senescence,
#   Dormancy, EVI_Minimum
#   Raw integer values (days since 1970-01-01) are dropped.
#   Output columns: DOY and Date (YYYY-MM-DD) only.
#
# Input:  floss500_dat.csv (from Script 2)
# Output: modis_pheno_ndvi.csv (exported to Google Drive -> Reidy_research,
#         then manually moved to Reidy_research/Hansen exploration/2026/)
# =============================================================================

library(geojsonio)
library(geojsonsf)
library(rgee)
library(sf)
library(dplyr)

# Initialize Earth Engine
ee_Initialize(drive = TRUE)

# -----------------------------------------------------------------------------
# 1. File paths
# -----------------------------------------------------------------------------
path.google <- "~/Google Drive/My Drive/"
path.dat    <- file.path(path.google, "Reidy_research/Hansen exploration/2026")

# -----------------------------------------------------------------------------
# 2. Read polygon data from Script 2 and build sf object
# -----------------------------------------------------------------------------
floss500 <- read.csv(file.path(path.dat, "floss500_dat.csv"))

geoms <- sf::st_sf(
  Label    = floss500$Label,
  Year     = floss500$Year,
  area_sqm = floss500$area_sqm,
  geometry = geojsonsf::geojson_sf(floss500$.geo)$geometry
)

sf::st_crs(geoms) <- 4326

# -----------------------------------------------------------------------------
# 3. Load polygon asset from GEE
#    NOTE: floss500_tmp was manually uploaded to GEE via the code editor.
#    To re-upload: run st_write below to regenerate the shapefile, then
#    upload all four files (.shp, .shx, .dbf, .prj) from path.dat to
#    projects/ee-breidy/assets/floss500_tmp in the GEE Assets panel.
#
#    sf::st_write(geoms,
#                 file.path(path.dat, "floss500_polygons.shp"),
#                 delete_dsn = TRUE)
# -----------------------------------------------------------------------------
table <- ee$FeatureCollection("projects/ee-breidy/assets/floss500_tmp")

# -----------------------------------------------------------------------------
# 4. Define analysis parameters
# -----------------------------------------------------------------------------
startYear <- 2001
endYear   <- 2023

# All 8 MCD12Q2 phenometric bands (cycle 1 suffix)
phenoBands <- c(
  'Greenup_1',
  'MidGreenup_1',
  'Maturity_1',
  'Peak_1',
  'MidGreendown_1',
  'Senescence_1',
  'Dormancy_1',
  'EVI_Minimum_1'
)

# Clean output names (drop _1 suffix)
phenoNames <- c(
  'Greenup',
  'MidGreenup',
  'Maturity',
  'Peak',
  'MidGreendown',
  'Senescence',
  'Dormancy',
  'EVI_Minimum'
)

# -----------------------------------------------------------------------------
# 5. NDVI extraction function (disturbance severity check only)
#    Uses MOD13Q1 annual median composite at 250m resolution.
#    Purpose: confirm disturbance signal in MODIS — NOT productivity estimate.
# -----------------------------------------------------------------------------
calculateYearNDVI <- function(feature, year) {
  feature <- ee$Feature(feature)
  year    <- ee$Number(year)
  
  startDate <- ee$Date$fromYMD(year, 1, 1)
  endDate   <- ee$Date$fromYMD(year, 12, 31)
  
  modisColl <- ee$ImageCollection('MODIS/061/MOD13Q1')$
    filterDate(startDate, endDate)$
    filterBounds(feature$geometry())
  
  composite <- modisColl$median()
  
  ndvi <- ee$Algorithms$If(
    modisColl$size()$gt(0),
    composite$select('NDVI')$
      reduceRegion(
        reducer   = ee$Reducer$mean(),
        geometry  = feature$geometry(),
        scale     = 250,
        maxPixels = 1e9
      )$get('NDVI'),
    NULL
  )
  
  ndvi <- ee$Algorithms$If(
    ndvi,
    ee$Number(ndvi)$divide(10000),
    NULL
  )
  
  ee$Feature(NULL, list(
    'Label'           = feature$get('Label'),
    'Year'            = year,
    'NDVI'            = ndvi,
    'NDVI_ImageCount' = modisColl$size()
  ))
}

# -----------------------------------------------------------------------------
# 6. Phenometrics extraction function
#    Extracts all 8 MCD12Q2 bands. Returns DOY and Date only (no Raw).
#    Raw values are days since 1970-01-01; valid range 11138-32766.
#    32767 = fill value (excluded). Values < 11138 = invalid (excluded).
# -----------------------------------------------------------------------------
calculateYearPhenology <- function(feature, year) {
  feature <- ee$Feature(feature)
  year    <- ee$Number(year)
  
  startDate <- ee$Date$fromYMD(year, 1, 1)
  endDate   <- ee$Date$fromYMD(year, 12, 31)
  
  phenoColl  <- ee$ImageCollection('MODIS/061/MCD12Q2')$
    filterDate(startDate, endDate)$
    filterBounds(feature$geometry())
  
  phenoImage <- ee$Image(phenoColl$first())
  
  # Process a single phenology band — returns DOY and Date only
  processPhenoBand <- function(bandName) {
    phenoValue <- ee$Algorithms$If(
      phenoColl$size()$gt(0),
      phenoImage$select(bandName)$
        reduceRegion(
          reducer   = ee$Reducer$mean(),
          geometry  = feature$geometry(),
          scale     = 500,
          maxPixels = 1e9
        )$get(bandName),
      NULL
    )
    
    # Filter fill values and out-of-range values
    phenoValue <- ee$Algorithms$If(
      ee$Algorithms$IsEqual(phenoValue, NULL),
      NULL,
      ee$Algorithms$If(
        ee$Number(phenoValue)$gte(32767),
        NULL,
        ee$Algorithms$If(
          ee$Number(phenoValue)$lt(11138),
          NULL,
          phenoValue
        )
      )
    )
    
    # Convert days-since-1970 to calendar date
    calendarDate <- ee$Algorithms$If(
      phenoValue,
      ee$Date('1970-01-01')$advance(ee$Number(phenoValue), 'day'),
      NULL
    )
    
    # DOY (1-365/366)
    dayOfYear <- ee$Algorithms$If(
      calendarDate,
      ee$Date(calendarDate)$getRelative('day', 'year')$add(1),
      NULL
    )
    
    # Date string YYYY-MM-DD
    dateString <- ee$Algorithms$If(
      calendarDate,
      ee$Date(calendarDate)$format('YYYY-MM-dd'),
      NULL
    )
    
    list(doy = dayOfYear, date = dateString)
  }
  
  # Process all 8 bands
  results <- mapply(processPhenoBand, phenoBands, SIMPLIFY = FALSE)
  names(results) <- phenoNames
  
  # Build output feature — DOY and Date only, no Raw
  props <- list(
    'Label' = feature$get('Label'),
    'Year'  = year,
    'Phenology_ImageCount' = phenoColl$size()
  )
  
  for (nm in phenoNames) {
    props[[paste0(nm, '_DOY')]]  <- results[[nm]]$doy
    props[[paste0(nm, '_Date')]] <- results[[nm]]$date
  }
  
  ee$Feature(NULL, props)
}

# -----------------------------------------------------------------------------
# 7. Map over all features and years and export
# -----------------------------------------------------------------------------
yearRange <- ee$List$sequence(startYear, endYear)

combinedCollection <- table$map(ee_utils_pyfunc(function(feature) {
  yearlyData <- yearRange$map(ee_utils_pyfunc(function(year) {
    feature <- ee$Feature(feature)
    year    <- ee$Number(year)
    
    ndviData  <- calculateYearNDVI(feature, year)
    phenoData <- calculateYearPhenology(feature, year)
    
    props <- list(
      'Label'                = feature$get('Label'),
      'Year'                 = year,
      'NDVI'                 = ndviData$get('NDVI'),
      'NDVI_ImageCount'      = ndviData$get('NDVI_ImageCount'),
      'Phenology_ImageCount' = phenoData$get('Phenology_ImageCount')
    )
    
    for (nm in phenoNames) {
      props[[paste0(nm, '_DOY')]]  <- phenoData$get(paste0(nm, '_DOY'))
      props[[paste0(nm, '_Date')]] <- phenoData$get(paste0(nm, '_Date'))
    }
    
    ee$Feature(NULL, props)
  }))
  ee$FeatureCollection(yearlyData)
}))$flatten()

task <- ee_table_to_drive(
  collection  = combinedCollection,
  description = 'modis_pheno_ndvi',
  folder      = 'Reidy_research',
  fileFormat  = 'CSV',
  timePrefix  = FALSE
)
task$start()
ee_monitoring(task)

# -----------------------------------------------------------------------------
# NOTE: After export completes, manually move modis_pheno_ndvi.csv from
# Reidy_research/ into Reidy_research/Hansen exploration/2026/
# before running Script 4.
# -----------------------------------------------------------------------------