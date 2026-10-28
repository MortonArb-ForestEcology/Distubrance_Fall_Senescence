# Patch Step 1: find every qualifying forest loss patch in the study area
# -----------------------------------------------------------------------------
# Open in RStudio and use Source with Echo so the tables
# and map show. A plain "Source" still saves every file but won't show them.
# Needs: rgee, googledrive, sf, dplyr, ggplot2
# Reads:  scaleup_patches/ecoregions_conus_l3.gpkg   (made by Patch Step 0)
# Saves (names end in _TEST for a test run):
#   patches_all_labeled.csv     every patch in the study area, with its rules
#   patches_qualifying.gpkg     polygons of the qualifying patches
#   map_qualifying_patches.png
#   by_tile_v2/                 raw download for each tile
#
# HOW PATCHES ARE KEPT WHOLE AND COUNTED ONCE
#   Earth Engine can't do the whole eastern US in one job, so the study area
#   is cut into 2-degree squares ("tiles"), one job each. Two rules make the
#   tiles invisible in the result:
#     - Each job looks 10 km past its tile (and past the study-area edge), so
#       a patch near an edge is found WHOLE, not cut off.
#     - A patch found by two neighboring jobs is kept only from the tile its
#       ceenter falls in, so it is counted ONCE.
#   The script then checks that no kept patch is wider than 10 km (a patch
#   that big could have been cut off) and reports any that are.
#   Ecoregions are not used to cut anything. Each whole patch is assigned to
#   the ONE Level III ecoregion that holds most of its area. A patch is in the
#   study if that ecoregion is included. Patches that cross a boundary are
#   marked (crosses_ecoregion, eco_share).
#
# What it does:
#   1. Reads the included ecoregions from Step 0 and lists the tiles
#   2. One Earth Engine job per tile:
#        - Hansen forest loss 2004-2020 (loss with no gain, as in the pilot)
#          turned into patches >= 5 ha
#        - NLCD hardwood share (41/43) in every NLCD year, and 2021 land cover
#        - a loose pre-filter: keeps patches >= 50% hardwood before the loss
#          (the real 75% rule is applied below). This keeps the downloads small
#          in pine country.
#   3. Downloads each tile as it finishes. Finished tiles are skipped on a
#      re-run, so you can stop R and Source again later.
#   4. Keeps each patch once, assigns it to an ecoregion, applies the rules:
#        - >= 75% hardwood in the latest NLCD year at least 2 years before loss
#        - not "converted": mostly farm/developed/bare in 2021, or still mostly
#          grass 7+ years after the loss (mines, pasture, development)
#   5. Saves the results, shows the counts, maps the qualifying patches.
#
# TEST RUN FIRST: TEST_RUN <- TRUE runs only the tile holding the pilot area.
# The full run will take several hours; a progress line shows every 3 minutes.
# If R stops, wait until the Tasks tab shows nothing running, then Source
# again: it downloads what finished and only restarts what's missing.
# =============================================================================

library(rgee)
library(googledrive)
library(sf)
library(dplyr)
library(ggplot2)

ee_Initialize(drive = TRUE)

# ---- Settings ---------------------------------------------------------------
RUN_TAG          <- "v2"    # names the exports; change it if you change any setting below
TEST_RUN         <- FALSE   # TRUE = only the tile with the pilot area. FALSE = everything.
TEST_POINT       <- c(-79.862539, 37.829550)
LOSS_YEAR_MIN    <- 2004    # >= 3 years of before data (phenology starts 2001)
LOSS_YEAR_MAX    <- 2020    # >= 3 years of after data (phenology ends 2023)
MIN_FIND_HA      <- 5       # smallest patch recorded
TILE_DEG         <- 2       # tile size in degrees
BUFFER_M         <- 10000   # each job looks this far past its tile and the study edge
SERVER_DECID_MIN <- 0.5     # loose pre-filter in Earth Engine
DECID_THRESHOLD  <- 0.75    # the real rule: share of the patch hardwood before the loss
BEFORE_LAG       <- 2       # "before" NLCD year must be >= this many years before the loss
GRASS_YEARS      <- 7       # still grass this long after a loss -> converted
EPOCHS           <- c(2001, 2004, 2006, 2008, 2011, 2013, 2016, 2019)  # NLCD 2019 release
AFTER_YEAR       <- 2021

outDir      <- file.path("~/Google Drive/My Drive/Reidy_research", "scaleup_patches")
tileDir     <- file.path(outDir, paste0("by_tile_", RUN_TAG))
dir.create(tileDir, recursive = TRUE, showWarnings = FALSE)
driveFolder <- "Reidy_research"
tag         <- if (TEST_RUN) "_TEST" else ""   # test outputs never overwrite real ones

# ---- 1. Study area and tiles ------------------------------------------------
eco_file <- file.path(outDir, "ecoregions_conus_l3.gpkg")
if (!file.exists(eco_file)) stop("Run patch_step0_study_area.R first.", call. = FALSE)
eco <- st_read(eco_file, quiet = TRUE) %>%
  mutate(us_l3code = as.character(us_l3code), included = as.logical(included))
inc_codes <- eco$us_l3code[eco$included]
length(inc_codes)   # included Level III ecoregions

study_buf <- eco %>% filter(included) %>% st_union() %>% st_buffer(BUFFER_M) %>% st_transform(4326)
bb <- unname(st_bbox(study_buf))   # xmin, ymin, xmax, ymax
tiles <- expand.grid(i = floor(bb[1] / TILE_DEG):floor(bb[3] / TILE_DEG),
                     j = floor(bb[2] / TILE_DEG):floor(bb[4] / TILE_DEG)) %>%
  mutate(lon0 = i * TILE_DEG, lat0 = j * TILE_DEG,
         id   = sprintf("tile_%dW_%dN", -lon0, lat0))
tile_poly <- st_sfc(lapply(seq_len(nrow(tiles)), function(k) {
  st_polygon(list(matrix(c(tiles$lon0[k], tiles$lat0[k],
                           tiles$lon0[k] + TILE_DEG, tiles$lat0[k],
                           tiles$lon0[k] + TILE_DEG, tiles$lat0[k] + TILE_DEG,
                           tiles$lon0[k], tiles$lat0[k] + TILE_DEG,
                           tiles$lon0[k], tiles$lat0[k]), ncol = 2, byrow = TRUE)))
}), crs = 4326)
tiles <- tiles[lengths(st_intersects(tile_poly, study_buf)) > 0, ]

if (TEST_RUN) {   # test run: only the tile holding the pilot area
  tiles <- tiles %>% filter(lon0 == floor(TEST_POINT[1] / TILE_DEG) * TILE_DEG,
                            lat0 == floor(TEST_POINT[2] / TILE_DEG) * TILE_DEG)
}
nrow(tiles)   # tiles to process

# ---- 2. Earth Engine: images and the per-tile job ---------------------------
studyFC  <- ee$FeatureCollection("EPA/Ecoregions/2013/L3")$
  filter(ee$Filter$inList("us_l3code", as.list(inc_codes)))
hansen   <- ee$Image("UMD/hansen/global_forest_change_2024_v1_12")
lossyear <- hansen$select("lossyear")

qualLoss <- hansen$select("loss")$eq(1)$
  And(hansen$select("gain")$eq(0))$
  And(lossyear$gte(LOSS_YEAR_MIN - 2000))$
  And(lossyear$lte(LOSS_YEAR_MAX - 2000))

# Skip tiny clumps before making polygons (keeps the jobs light). Clumps are
# built like the patches: same loss year, touching side-to-side.
clumpHa <- lossyear$updateMask(qualLoss)$
  connectedPixelCount(maxSize = 1024L, eightConnected = FALSE)$
  multiply(ee$Image$pixelArea())$divide(10000)
lossBig <- lossyear$updateMask(qualLoss$And(clumpHa$gte(0.8 * MIN_FIND_HA)))

nlcdCol <- ee$ImageCollection("USGS/NLCD_RELEASES/2019_REL/NLCD")
attrImg <- NULL
for (yr in EPOCHS) {
  b <- nlcdCol$filter(ee$Filter$eq("system:index", as.character(yr)))$first()$
    select("landcover")$remap(list(41L, 43L), list(1L, 1L), 0)$rename(paste0("decid_", yr))
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

before_epoch_of <- function(ly) max(EPOCHS[EPOCHS <= ly - BEFORE_LAG])
beforeBands <- ee$List(as.list(paste0("decid_", sapply(LOSS_YEAR_MIN:LOSS_YEAR_MAX, before_epoch_of))))

after_cols <- c("after_decid", "after_evergreen", "after_shrub", "after_grass",
                "after_farm", "after_developed", "after_bare")
keep_cols  <- c("patch_id", "tile", "loss_year", "area_ha", "lon", "lat",
                paste0("decid_", EPOCHS), after_cols)

find_patches <- function(tile_id, lon0, lat0) {
  rect   <- ee$Geometry$Rectangle(coords = list(lon0, lat0, lon0 + TILE_DEG, lat0 + TILE_DEG),
                                  proj = "EPSG:4326", geodesic = FALSE)
  big    <- rect$buffer(2 * BUFFER_M, 1000)
  local  <- studyFC$filterBounds(big)$geometry(maxError = 250)$intersection(big, 250)
  region <- local$buffer(BUFFER_M, 250)$intersection(rect$buffer(BUFFER_M, 1000), 250)
  
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
      "patch_id",  ee$String(paste0(RUN_TAG, "_", tile_id, "_"))$cat(f$id()),
      "tile",      tile_id,
      "loss_year", ee$Number(f$get("lossyear"))$add(2000),
      "area_ha",   f$geometry()$area(maxError = 1)$divide(10000),
      "lon",       cen$get(0),
      "lat",       cen$get(1)
    )
  }))$filter(ee$Filter$gte("area_ha", MIN_FIND_HA))
  
  fc <- attrImg$reduceRegions(collection = fc, reducer = ee$Reducer$mean(),
                              scale = 30, tileScale = 8)
  
  # Loose hardwood-before pre-filter (patches outside NLCD come out empty and drop)
  fc <- fc$map(ee_utils_pyfunc(function(f) {
    band <- beforeBands$get(ee$Number(f$get("loss_year"))$subtract(LOSS_YEAR_MIN)$toInt())
    f$set("decid_before_srv", f$get(band))
  }))$filter(ee$Filter$gte("decid_before_srv", SERVER_DECID_MIN))
  
  fc$select(as.list(keep_cols))
}

# ---- 3. Run the jobs and download (finished tiles are skipped) --------------
export_name <- function(id) paste0("patches_", RUN_TAG, "_", id)
tile_folder <- function(id) file.path(tileDir, id)
is_done     <- function(id) file.exists(file.path(tile_folder(id), "done.txt"))

ee_states <- function() {   # newest state of each export, if Earth Engine will tell us
  tl <- tryCatch(ee$data$getTaskList(), error = function(e) NULL)
  if (is.null(tl) || length(tl) == 0) {
    return(data.frame(description = character(), state = character(), error = character()))
  }
  d <- data.frame(
    description = vapply(tl, function(x) if (is.null(x$description)) "" else x$description, ""),
    state       = vapply(tl, function(x) if (is.null(x$state)) "" else x$state, ""),
    error       = vapply(tl, function(x) if (is.null(x$error_message)) "" else x$error_message, ""),
    stringsAsFactors = FALSE)
  d[!duplicated(d$description), ]   # the list is newest first
}

download_ready <- function(ids) {
  f <- drive_ls(path = driveFolder, pattern = paste0("patches_", RUN_TAG, "_"))
  for (id in ids) {
    nm  <- export_name(id)
    hit <- f[f$name == paste0(nm, ".geojson") | startsWith(f$name, paste0(nm, "-")), ]
    if (nrow(hit) == 0) next
    dir.create(tile_folder(id), recursive = TRUE, showWarnings = FALSE)
    for (k in seq_len(nrow(hit))) {
      drive_download(as_id(hit$id[k]), path = file.path(tile_folder(id), hit$name[k]),
                     overwrite = TRUE)
    }
    writeLines(hit$name, file.path(tile_folder(id), "done.txt"))
  }
}

todo <- tiles$id[!vapply(tiles$id, is_done, TRUE)]
length(todo)   # tiles not downloaded yet (0 = all done, skips straight to section 4)
failed <- data.frame()

if (length(todo) > 0) {
  download_ready(todo)                       # anything an earlier run finished
  todo <- todo[!vapply(todo, is_done, TRUE)]
  
  st <- ee_states()
  waiting <- st$description[st$state %in% c("READY", "RUNNING")]
  n_started <- 0
  for (id in todo) {
    if (export_name(id) %in% waiting) next   # still running from an earlier session
    k <- match(id, tiles$id)
    t <- ee_table_to_drive(collection = find_patches(id, tiles$lon0[k], tiles$lat0[k]),
                           description = export_name(id), folder = driveFolder,
                           fileFormat = "GeoJSON", timePrefix = FALSE)
    t$start()
    n_started <- n_started + 1
  }
  
  # Progress lines use message() because nothing inside a loop shows on its own
  message(sprintf("Started %d Earth Engine jobs (%d were already running). Checking every 3 minutes...",
                  n_started, length(todo) - n_started))
  
  t0 <- Sys.time()
  repeat {
    Sys.sleep(180)
    download_ready(todo)
    todo <- todo[!vapply(todo, is_done, TRUE)]
    st <- ee_states()
    st <- st[st$description %in% export_name(todo), ]
    failed <- st[st$state %in% c("FAILED", "CANCELLED"), ]
    message(sprintf("%s  tiles left: %d  (running: %d, waiting: %d, failed: %d)",
                    format(Sys.time(), "%H:%M"), length(todo), sum(st$state == "RUNNING"),
                    sum(st$state == "READY"), nrow(failed)))
    if (length(todo) == 0 || length(todo) == nrow(failed)) break
    if (difftime(Sys.time(), t0, units = "hours") > 24) break
  }
}

failed   # tiles that failed (empty = none). Source again to retry them.

have    <- tiles$id[vapply(tiles$id, is_done, TRUE)]
missing <- setdiff(tiles$id, have)
if (length(missing) > 0) {
  warning(length(missing), " tile(s) not downloaded yet: ", paste(missing, collapse = ", "),
          ". Source again to retry just those.", call. = FALSE)
}

# ---- 4. Keep each patch once, assign ecoregions, apply the rules ------------
read_tile <- function(id) {
  files <- list.files(tile_folder(id), pattern = "\\.geojson$", full.names = TRUE)
  parts <- lapply(files, function(path) {
    x <- tryCatch(st_read(path, quiet = TRUE), error = function(e) NULL)   # empty tile
    if (is.null(x) || nrow(x) == 0) return(NULL)
    x$patch_id <- as.character(x$patch_id)
    x$tile     <- as.character(x$tile)
    x[, keep_cols]
  })
  do.call(rbind, parts)
}
poly <- do.call(rbind, lapply(have, read_tile))
poly <- poly[!duplicated(poly$patch_id), ]
nrow(poly)   # patches found in all tiles, incl. the 10 km overlaps

# 4a. Counted once: keep a patch only from the tile its centre falls in
poly <- poly %>%
  left_join(tiles %>% select(tile = id, lon0, lat0), by = "tile") %>%
  filter(lon >= lon0, lon < lon0 + TILE_DEG, lat >= lat0, lat < lat0 + TILE_DEG) %>%
  select(-lon0, -lat0)
nrow(poly)   # patches after removing tile overlaps

# 4b. Kept whole: a patch narrower than the 10 km look-past can't have been cut
poly_alb <- st_make_valid(st_transform(poly, 5070))
ext <- t(sapply(st_geometry(poly_alb), st_bbox))
poly$extent_km        <- round(pmax(ext[, 3] - ext[, 1], ext[, 4] - ext[, 2]) / 1000, 2)
poly$possibly_cut_off <- poly$extent_km * 1000 >= BUFFER_M
sum(poly$possibly_cut_off)   # patches 10 km or more across (could have been cut off at a job edge)

# 4c. One ecoregion per patch: the Level III ecoregion holding most of its area
eco_alb <- st_transform(eco, 5070)
cen     <- st_as_sf(st_drop_geometry(poly)[, c("patch_id", "lon", "lat")],
                    coords = c("lon", "lat"), crs = 4326) %>% st_transform(5070)
hit     <- st_within(cen, eco_alb)
eco_centre <- eco_alb$us_l3code[vapply(hit, function(x) if (length(x)) x[1] else NA_integer_, 1L)]

# Is the whole patch inside the ecoregion its centre is in? (most are)
inside <- rep(FALSE, nrow(poly))
for (code in unique(na.omit(eco_centre))) {
  idx <- which(eco_centre == code)
  inside[idx] <- lengths(st_within(poly_alb[idx, ], eco_alb[eco_alb$us_l3code == code, ])) > 0
}
poly$crosses_ecoregion <- !inside
poly$us_l3code <- eco_centre
poly$eco_share <- ifelse(inside, 1, NA_real_)

# For the ones that cross: measure how much is in each ecoregion, take the largest
if (any(!inside)) {
  ov <- suppressWarnings(st_intersection(poly_alb[!inside, "patch_id"], eco_alb[, "us_l3code"]))
  ov$a <- as.numeric(st_area(ov))
  best <- st_drop_geometry(ov) %>%
    group_by(patch_id) %>%
    mutate(share = a / sum(a)) %>%
    slice_max(a, n = 1, with_ties = FALSE) %>%
    ungroup()
  m  <- match(poly$patch_id, best$patch_id)
  ok <- !is.na(m)
  poly$us_l3code[ok] <- best$us_l3code[m[ok]]
  poly$eco_share[ok] <- round(best$share[m[ok]], 3)
}
sum(poly$crosses_ecoregion)   # patches crossing an ecoregion boundary (each assigned to the one with most of its area)

poly <- poly %>%
  left_join(st_drop_geometry(eco) %>%
              select(us_l3code, us_l3name, na_l2code, na_l2name, included),
            by = "us_l3code") %>%
  mutate(in_boundary = coalesce(included, FALSE)) %>%
  select(-included)

# 4d. Hardwood before + not converted
d <- st_drop_geometry(poly)
before_epoch <- vapply(d$loss_year, before_epoch_of, 0)
decid_mat    <- as.matrix(d[, paste0("decid_", EPOCHS)])

pts_df <- d %>%
  mutate(
    before_epoch = before_epoch,
    decid_before = decid_mat[cbind(seq_len(n()), match(before_epoch, EPOCHS))],
    hardwood_ok  = coalesce(decid_before >= DECID_THRESHOLD, FALSE),
    years_after  = AFTER_YEAR - loss_year,
    after_label  = case_when(
      after_farm + after_developed + after_bare >= 0.5   ~ "converted",
      after_grass >= 0.5 & years_after >= GRASS_YEARS     ~ "converted",
      after_decid + after_evergreen + after_shrub >= 0.5  ~ "regrowing",
      TRUE                                                ~ "unclear"
    ),
    qualifies = in_boundary & hardwood_ok & after_label != "converted"
  ) %>%
  filter(in_boundary)

write.csv(pts_df, file.path(outDir, paste0("patches_all_labeled", tag, ".csv")), row.names = FALSE)

qual_poly <- poly %>%
  filter(patch_id %in% pts_df$patch_id[pts_df$qualifies]) %>%
  select(patch_id, tile, loss_year, area_ha, lon, lat, us_l3code, us_l3name, na_l2code, na_l2name,
         crosses_ecoregion, eco_share, extent_km, possibly_cut_off) %>%
  left_join(pts_df %>% select(patch_id, decid_before, after_label), by = "patch_id")
st_write(qual_poly, file.path(outDir, paste0("patches_qualifying", tag, ".gpkg")),
         delete_dsn = TRUE, quiet = TRUE)

# ---- 5. Summary + map -------------------------------------------------------
# Patches by Level II region (patches_found = >= 5 ha and >= 50% hardwood before)
funnel <- pts_df %>%
  group_by(na_l2code, na_l2name) %>%
  summarise(
    patches_found    = n(),
    hardwood_75      = sum(hardwood_ok),
    converted        = sum(hardwood_ok & after_label == "converted"),
    qualify          = sum(qualifies),
    cross_boundary   = sum(qualifies & crosses_ecoregion),
    .groups = "drop"
  ) %>%
  as.data.frame()
funnel

# Qualifying patches by Level III region
q <- pts_df %>% filter(qualifies)
by_l3 <- q %>% count(na_l2code, us_l3code, us_l3name) %>% as.data.frame()
by_l3

# Qualifying patches in total, how many cross an ecoregion boundary (and the
# median share in the ecoregion they were assigned to), and how many could
# have been cut off at a job edge
q_summary <- data.frame(
  qualifying       = nrow(q),
  cross_boundary   = sum(q$crosses_ecoregion),
  median_eco_share = median(q$eco_share[q$crosses_ecoregion], na.rm = TRUE),
  possibly_cut_off = sum(q$possibly_cut_off))
q_summary

inc  <- eco %>% filter(included)
qpts <- st_as_sf(q, coords = c("lon", "lat"), crs = 4326) %>% st_transform(5070)

ggplot(data=qpts) +
  geom_sf(data=inc, aes(fill=na_l2name), color="grey50", linewidth=0.1, alpha=0.3) + # included ecoregions, colored by Level II
  geom_sf(color="darkgreen", size=0.05, alpha=0.3) + # one dot per qualifying patch
  theme_bw() +
  labs(title="Qualifying forest loss patches in the study area",
       subtitle="Each dot is one patch that meets the Step 1 rules",
       fill="Level II ecoregion")

ggsave(file.path(outDir, paste0("map_qualifying_patches", tag, ".png")), width=12, height=8.5, dpi=200)