# Patch Step 1: find every qualifying forest loss patch inside the Level II boundary
# -----------------------------------------------------------------------------
# Self-contained. Open in RStudio and click "Source".
# Needs: rgee, googledrive, sf, dplyr, ggplot2
#
# Boundary: EPA Level II ecoregions, used only as the outer boundary (no strata)
#   8.4 Ozark/Ouachita-Appalachian Forests   5.3 Atlantic Highlands
#   5.2 Mixed Wood Shield                     8.1 Mixed Wood Plains
#
# What it does:
#   1. Gets the boundary from Earth Engine and lists its Level III regions
#   2. For EACH Level III region, one Earth Engine job:
#        - turns Hansen forest loss (2004-2020, loss with no gain -- same rule
#          as your pilot Step 1) into patches; keeps every patch >= 5 ha so
#          Step 2 can show the full size distribution
#        - measures each patch with NLCD: share that was hardwood (41/43) in
#          every NLCD year, and what the land was in 2021
#      All jobs start at once and Earth Engine works through them.
#   3. Downloads each region's patches (polygons) as it finishes. Regions
#      already downloaded are skipped, so if some fail or you stop R, just
#      Source again and it picks up the rest.
#   4. Combines the regions. A patch QUALIFIES if:
#        - its centre is inside the boundary (each patch counted once, in the
#          region its centre falls in)
#        - >= 75% hardwood in the latest NLCD year at least 2 years before the cut
#        - not "converted": mostly farm/developed/bare in 2021, or still mostly
#          grass 7+ years after the cut (mines, pasture, development)
#   5. Saves the results, prints how many survive each filter, maps them.
#
# Expect the full run to take hours. Progress prints every 2 minutes, and you
# can also watch the Tasks tab in the GEE Code Editor.
# =============================================================================

library(rgee)
library(googledrive)
library(sf)
library(dplyr)
library(ggplot2)

ee_Initialize(drive = TRUE)

# ---- Settings ---------------------------------------------------------------
L2_CODES        <- c("8.4", "5.3", "5.2", "8.1")
TEST_L3         <- NA    # quick test on ONE Level III region (69 = Central Appalachians,
# contains your pilot area). Set to NA for the real run.
LOSS_YEAR_MIN   <- 2004    # >= 3 years of before data (phenology starts 2001)
LOSS_YEAR_MAX   <- 2020    # >= 3 years of after data (phenology ends 2023)
MIN_FIND_HA     <- 5       # smallest patch recorded (not the sampling cutoff)
BUFFER_M        <- 2000    # look past each region's edge so edge patches stay whole
DECID_THRESHOLD <- 0.75    # share of the patch that was hardwood before the cut
BEFORE_LAG      <- 2       # "before" NLCD year must be >= this many years before the cut
GRASS_YEARS     <- 7       # still grass this long after a cut -> converted
EPOCHS          <- c(2001, 2004, 2006, 2008, 2011, 2013, 2016, 2019)  # NLCD 2019 release
AFTER_YEAR      <- 2021

outDir      <- file.path("~/Google Drive/My Drive/Reidy_research", "scaleup_patches")
regDir      <- file.path(outDir, "by_region")
dir.create(regDir, recursive = TRUE, showWarnings = FALSE)
driveFolder <- "Reidy_research"
TEST_RUN    <- !is.na(TEST_L3)
tag         <- if (TEST_RUN) "_TEST" else ""   # test outputs never overwrite real ones

# ---- 1. Boundary ------------------------------------------------------------
ecoAll <- ee$FeatureCollection("EPA/Ecoregions/2013/L3")
ecoFC  <- ecoAll$filter(ee$Filter$inList("na_l2code", as.list(L2_CODES)))

if (ecoFC$size()$getInfo() == 0) {
  codes <- unlist(ecoAll$aggregate_array("na_l2code")$distinct()$getInfo())
  stop("No Level II code matched. Codes in the dataset: ",
       paste(sort(codes), collapse = ", "), call. = FALSE)
}
if (TEST_RUN) {
  ecoFC <- ecoFC$filter(ee$Filter$eq("us_l3code", as.character(TEST_L3)))
  if (ecoFC$size()$getInfo() == 0) stop("TEST_L3 = ", TEST_L3, " isn't inside the chosen Level II regions.")
  cat("\nTEST RUN: only Level III region", TEST_L3, "\n")
}
ecoFC <- ecoFC$map(ee_utils_pyfunc(function(f) f$simplify(maxError = 250)))

eco <- ee_as_sf(ecoFC, via = "drive", maxFeatures = 10000) %>%
  st_transform(5070) %>%
  st_make_valid() %>%
  group_by(na_l2code, na_l2name, us_l3code, us_l3name) %>%
  summarise(.groups = "drop") %>%
  mutate(us_l3code = as.character(us_l3code),
         area_km2  = round(as.numeric(st_area(.)) / 1e6)) %>%
  arrange(na_l2code, us_l3code)

cat("\nLevel III regions inside the boundary:\n")
print(as.data.frame(st_drop_geometry(eco)))
cat(sprintf("Total area: %s km2\n", format(sum(eco$area_km2), big.mark = ",")))
st_write(eco, file.path(outDir, paste0("boundary_ecoregions", tag, ".gpkg")),
         delete_dsn = TRUE, quiet = TRUE)

l3_codes <- eco$us_l3code

# ---- 2. Images used by every region's job ------------------------------------
hansen   <- ee$Image("UMD/hansen/global_forest_change_2024_v1_12")
lossyear <- hansen$select("lossyear")

qualLoss <- hansen$select("loss")$eq(1)$
  And(hansen$select("gain")$eq(0))$
  And(lossyear$gte(LOSS_YEAR_MIN - 2000))$
  And(lossyear$lte(LOSS_YEAR_MAX - 2000))

# Skip tiny clumps before making polygons (keeps the jobs manageable). Clumps
# are built like the patches: same loss year, touching side-to-side.
clumpHa <- lossyear$updateMask(qualLoss)$
  connectedPixelCount(maxSize = 1024L, eightConnected = FALSE)$
  multiply(ee$Image$pixelArea())$divide(10000)
lossBig <- lossyear$updateMask(qualLoss$And(clumpHa$gte(0.8 * MIN_FIND_HA)))

# NLCD: hardwood share in every epoch, and land cover groups in 2021
nlcdCol <- ee$ImageCollection("USGS/NLCD_RELEASES/2019_REL/NLCD")
attrImg <- NULL
for (yr in EPOCHS) {
  b <- nlcdCol$filter(ee$Filter$eq("system:index", as.character(yr)))$first()$
    select("landcover")$
    remap(list(41L, 43L), list(1L, 1L), 0)$
    rename(paste0("decid_", yr))
  attrImg <- if (is.null(attrImg)) b else attrImg$addBands(b)
}
lcAfter <- ee$Image(sprintf("USGS/NLCD_RELEASES/2021_REL/NLCD/%d", AFTER_YEAR))$select("landcover")
group_band <- function(classes, name) {
  lcAfter$remap(as.list(classes), as.list(rep(1L, length(classes))), 0)$rename(name)
}
attrImg <- attrImg$
  addBands(group_band(c(41L, 43L),           "after_decid"))$
  addBands(group_band(42L,                   "after_evergreen"))$
  addBands(group_band(52L,                   "after_shrub"))$
  addBands(group_band(71L,                   "after_grass"))$
  addBands(group_band(c(81L, 82L),           "after_farm"))$
  addBands(group_band(c(21L, 22L, 23L, 24L), "after_developed"))$
  addBands(group_band(31L,                   "after_bare"))

after_cols <- c("after_decid", "after_evergreen", "after_shrub", "after_grass",
                "after_farm", "after_developed", "after_bare")
keep_cols  <- c("patch_id", "run_l3", "loss_year", "area_ha", "lon", "lat",
                paste0("decid_", EPOCHS), after_cols)

find_patches <- function(code) {
  region <- ecoFC$filter(ee$Filter$eq("us_l3code", code))$
    geometry(maxError = 250)$
    dissolve(maxError = 250)$
    buffer(BUFFER_M, 250)
  
  fc <- lossBig$reduceToVectors(
    geometry       = region,
    scale          = lossyear$projection()$nominalScale(),   # Hansen's own grid
    geometryType   = "polygon",
    eightConnected = FALSE,
    labelProperty  = "lossyear",
    maxPixels      = 1e13,
    tileScale      = 8
  )
  
  fc <- fc$map(ee_utils_pyfunc(function(f) {
    cen <- f$geometry()$centroid(maxError = 1)$coordinates()
    f$set(
      "patch_id",  ee$String(paste0("p", code, "_"))$cat(f$id()),
      "run_l3",    code,
      "loss_year", ee$Number(f$get("lossyear"))$add(2000),
      "area_ha",   f$geometry()$area(maxError = 1)$divide(10000),
      "lon",       cen$get(0),
      "lat",       cen$get(1)
    )
  }))$filter(ee$Filter$gte("area_ha", MIN_FIND_HA))
  
  # Several bands + one reducer -> columns are named after the bands.
  attrImg$reduceRegions(collection = fc, reducer = ee$Reducer$mean(),
                        scale = 30, tileScale = 8)$
    select(as.list(keep_cols))
}

# ---- 3. One Earth Engine job per region; download as they finish ------------
local_file <- function(code) file.path(regDir, paste0("patches_L3_", code, tag, ".geojson"))

todo <- l3_codes[!file.exists(vapply(l3_codes, local_file, ""))]
cat(sprintf("\nRegions already downloaded: %d of %d\n", length(l3_codes) - length(todo), length(l3_codes)))

if (length(todo) > 0) {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  tasks <- list()
  for (code in todo) {
    nm <- paste0("patches_L3_", code, tag, "_", stamp)   # unique name -> never picks up an old file
    t  <- ee_table_to_drive(collection = find_patches(code), description = nm,
                            folder = driveFolder, fileFormat = "GeoJSON", timePrefix = FALSE)
    t$start()
    tasks[[code]] <- list(task = t, name = nm)
  }
  cat(sprintf("Started %d Earth Engine jobs. Checking every 2 minutes...\n", length(tasks)))
  
  repeat {
    st <- vapply(tasks, function(x) x$task$status()$state, "")
    tb <- table(st)
    cat(format(Sys.time(), "%H:%M"), " ", paste(names(tb), tb, sep = ": ", collapse = " | "), "\n")
    if (all(st %in% c("COMPLETED", "FAILED", "CANCELLED"))) break
    Sys.sleep(120)
  }
  
  for (code in names(tasks)) {
    s <- tasks[[code]]$task$status()
    if (s$state != "COMPLETED") {
      cat("Region", code, "ended as", s$state, "--", if (!is.null(s$error_message)) s$error_message else "", "\n")
      next
    }
    f <- data.frame()
    for (try_i in 1:5) {                       # Drive can lag a little behind GEE
      f <- drive_ls(path = driveFolder, pattern = tasks[[code]]$name)
      if (nrow(f) > 0) break
      Sys.sleep(30)
    }
    if (nrow(f) == 0) { cat("Region", code, "finished but its file isn't in Drive yet.\n"); next }
    drive_download(as_id(f$id[1]), path = local_file(code), overwrite = TRUE)
  }
}

have    <- file.exists(vapply(l3_codes, local_file, ""))
missing <- l3_codes[!have]
if (length(missing) > 0) {
  warning("Regions not downloaded yet: ", paste(missing, collapse = ", "),
          ". Source the script again to retry just those.", call. = FALSE)
}

# ---- 4. Combine regions and apply the filters (in R) ------------------------
read_region <- function(path) {
  x <- st_read(path, quiet = TRUE)
  if (nrow(x) == 0) return(NULL)
  x$run_l3   <- as.character(x$run_l3)
  x$patch_id <- as.character(x$patch_id)
  x[, keep_cols]
}
poly <- do.call(rbind, lapply(vapply(l3_codes[have], local_file, ""), read_region))
cat(sprintf("\nPatches >= %d ha found (regions + edge buffers): %d\n", MIN_FIND_HA, nrow(poly)))

pts <- st_drop_geometry(poly) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(5070) %>%
  st_join(eco[, c("na_l2code", "na_l2name", "us_l3code", "us_l3name")], join = st_within)

before_epoch <- vapply(pts$loss_year, function(ly) max(EPOCHS[EPOCHS <= ly - BEFORE_LAG]), 0)
decid_mat    <- as.matrix(st_drop_geometry(pts)[, paste0("decid_", EPOCHS)])

pts <- pts %>%
  mutate(
    # Counted once: only from the job of the region its centre falls in
    in_boundary  = !is.na(us_l3code) & us_l3code == run_l3,
    before_epoch = before_epoch,
    decid_before = decid_mat[cbind(seq_len(n()), match(before_epoch, EPOCHS))],
    hardwood_ok  = decid_before >= DECID_THRESHOLD,
    years_after  = AFTER_YEAR - loss_year,
    after_label  = case_when(
      after_farm + after_developed + after_bare >= 0.5   ~ "converted",
      after_grass >= 0.5 & years_after >= GRASS_YEARS     ~ "converted",
      after_decid + after_evergreen + after_shrub >= 0.5  ~ "regrowing",
      TRUE                                                ~ "unclear"
    ),
    qualifies = in_boundary & hardwood_ok & after_label != "converted"
  )

pts_df <- st_drop_geometry(pts) %>% filter(in_boundary)
write.csv(pts_df, file.path(outDir, paste0("patches_all_labeled", tag, ".csv")), row.names = FALSE)

# Polygons of the qualifying patches (Step 3 uploads the sampled ones to GEE)
qual_poly <- poly %>%
  filter(patch_id %in% pts_df$patch_id[pts_df$qualifies]) %>%
  left_join(pts_df %>% select(patch_id, na_l2code, na_l2name, us_l3code, us_l3name,
                              decid_before, after_label),
            by = "patch_id")
st_write(qual_poly, file.path(outDir, paste0("patches_qualifying", tag, ".gpkg")),
         delete_dsn = TRUE, quiet = TRUE)

# ---- 5. Summary + map -------------------------------------------------------
funnel <- pts_df %>%
  group_by(na_l2code, na_l2name) %>%
  summarise(
    patches_5ha_plus = n(),
    hardwood_before  = sum(hardwood_ok),
    converted        = sum(hardwood_ok & after_label == "converted"),
    qualify          = sum(qualifies),
    of_which_regrow  = sum(qualifies & after_label == "regrowing"),
    of_which_unclear = sum(qualifies & after_label == "unclear"),
    .groups = "drop"
  )
cat("\nPatches surviving each filter, by Level II region:\n")
print(as.data.frame(funnel))

cat("\nQualifying patches by Level III region:\n")
print(as.data.frame(pts_df %>% filter(qualifies) %>% count(na_l2code, us_l3code, us_l3name)))
cat(sprintf("\nQualifying patches in total: %d\n", sum(pts_df$qualifies)))

p <- ggplot() +
  geom_sf(data = eco, aes(fill = na_l2name), color = "grey50", linewidth = 0.1, alpha = 0.3) +
  geom_sf(data = filter(pts, in_boundary & qualifies), color = "darkgreen", size = 0.05, alpha = 0.3) +
  labs(title = "Qualifying forest loss patches inside the Level II boundary",
       subtitle = sprintf(">= %d ha, lost %d-%d, >= %.0f%% hardwood before, not converted",
                          MIN_FIND_HA, LOSS_YEAR_MIN, LOSS_YEAR_MAX, 100 * DECID_THRESHOLD),
       fill = "Level II ecoregion") +
  theme_minimal()
print(p)
ggsave(file.path(outDir, paste0("map_qualifying_patches", tag, ".png")), p, width = 11, height = 8, dpi = 200)

cat("\nDone. Files in", outDir, "(names end in _TEST for a test run):\n",
    " patches_all_labeled.csv   (every patch, with its filters)\n",
    " patches_qualifying.gpkg   (polygons of the qualifying patches)\n",
    " boundary_ecoregions.gpkg, map_qualifying_patches.png\n",
    " by_region/                (raw download for each Level III region)\n")
if (TEST_RUN) {
  cat("\nTEST RUN finished. If the numbers look sensible, set TEST_L3 <- NA and run again.\n")
} else if (length(missing) == 0) {
  cat("Next: patch_step2_size_distribution.R\n")
}
