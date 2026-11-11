# Patch Step 4: satellite time series for every patch + all nearby uncut forest
# -----------------------------------------------------------------------------
# Open in RStudio and use Source with Echo (Ctrl+Shift+Enter) so the tables
# and figure show. A plain "Source" still saves every file but won't show them.
# Needs: rgee, googledrive, sf, dplyr, tidyr, ggplot2
# Reads:  scaleup_patches/all_patches.csv / .gpkg, patch_modis_pixels.csv   (Patch Step 3)
# Saves (in scaleup_patches/timeseries_full/):
#   patch_year_series.csv   one row per patch per year; use rows with in_analysis = TRUE
#   patch_summary.csv       one row per patch: core pixels, comparison pixels and
#                           how many pixels passed each rule, in_analysis
#   by_part_<VALUES_TAG>/   raw download for each Earth Engine job
#   first_look_patch_minus_control.png
# Earth Engine assets (in assetFolder):
#   full_v3_patches_<k>      patch outlines, part k
#   full_v3_core_<k>         core MODIS pixel centres, part k
#
# What it does:
#   1. Reads every patch from Step 3 and its core MODIS pixels (>= 50% lost)
#   2. Sets up the comparison rules for every MODIS pixel: the share that was
#      hardwood in 2001 (NLCD 41/43), the share ever lost (Hansen v1.13), the
#      share lost in it and its 8 neighbours (about 500 m around it), and
#      elevation. These are worked out inside each patch job, only for the
#      ring around that patch, in small pieces (TILE_SCALE). Making them as one
#      big image for the whole study area failed: far from the MODIS grid's
#      centre line each MODIS block is a slanted shape on the 30 m grids, and
#      Earth Engine couldn't handle pieces that big ("Reprojection output too large").
#   3. Uploads the patch outlines (simplified by 30 m) and core pixel centres,
#      in parts of PART_SIZE patches. Parts already uploaded are reused.
#   4. One Earth Engine job per part. For each patch:
#        - patch value each year = mean of its core pixels
#        - comparison value each year = mean of EVERY MODIS pixel 1-10 km from
#          the patch that was >= 75% hardwood in 2001, <= 5% ever lost, <= 10%
#          lost in it and its neighbours, and within 100 m elevation of the
#          patch. No random draw.
#        - how many pixels in the 1-10 km ring pass each rule
#      Variables: July NDVI (MOD13A1, good/marginal quality composites only)
#      and Greenup, Senescence, MidGreendown, Dormancy as day of year
#      (Jan 1 = 1) plus the quality flag (MCD12Q2, every year Earth Engine has).
#      Values are read straight from each MODIS pixel -- no resampling.
#   5. Downloads each part as it finishes. Finished parts are skipped on a
#      re-run, so you can stop R and Source again later.
#   6. Builds the patch-year table (patch minus comparison each year). A patch
#      is in the analysis if it has at least MIN_CONTROLS comparison pixels.
#
# Changed from the random-sample version: every patch, every qualifying
# comparison pixel, and "loss nearby" = the pixel and its 8 neighbours.
# If R stops, wait until the Tasks tab shows nothing running, then Source
# again. Status lines use message() (red in RStudio; that's normal).
# =============================================================================

library(rgee)
library(googledrive)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)

ee_Initialize(drive = TRUE)

# ---- Settings ---------------------------------------------------------------
RUN_TAG            <- "v3"     # same as Patch Step 1
VALUES_TAG         <- "v3a"    # names the Earth Engine jobs; change it if you change any setting below
CONTROL_MIN_M      <- 1000     # comparison pixels at least this far from the patch
CONTROL_MAX_M      <- 10000    # ... and no farther than this
MIN_CONTROLS       <- 3        # a patch needs at least this many comparison pixels to be used
CORE_FRAC          <- 0.5      # a "core" pixel is at least this share lost (= Step 3's rule)
CTRL_DECID_MIN     <- 0.75     # comparison: hardwood share in 2001
CTRL_LOSS_MAX      <- 0.05     # comparison: share of the pixel ever lost
CTRL_NEAR_LOSS_MAX <- 0.10     # comparison: share lost in the pixel and its 8 neighbours
ELEV_TOL_M         <- 100      # comparison: elevation within this of the patch
NDVI_YEARS         <- 2001:2025
MCD_LAST_YEAR      <- 2025     # newest MCD12Q2 year to try (only years Earth Engine has are used)
PART_SIZE          <- 1000     # patches per upload part and per Earth Engine job
SIMPLIFY_M         <- 30       # patch outlines simplified by this much before upload
TILE_SCALE         <- 16       # Earth Engine works in pieces 1/16 of the usual size (max 16)

outDir      <- file.path("~/Google Drive/My Drive/Reidy_research", "scaleup_patches")
tsDir       <- file.path(outDir, "timeseries_full")
partDir     <- file.path(tsDir, paste0("by_part_", VALUES_TAG))
dir.create(partDir, recursive = TRUE, showWarnings = FALSE)
driveFolder <- "Reidy_research"
assetFolder <- "projects/breidyee/assets"
patch_asset <- function(k) paste0(assetFolder, "/full_", RUN_TAG, "_patches_", k)
core_asset  <- function(k) paste0(assetFolder, "/full_", RUN_TAG, "_core_", k)

# MODIS sinusoidal grid (sphere R = 6371007.181 m; 2400 x 463.3 m pixels per tile)
SINU <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
X0   <- -20015109.354
Y0   <-  10007554.677
PX   <-  463.312716528
# TRUE where a lon/lat (from Earth Engine) falls in the expected MODIS pixel
grid_ok <- function(lon, lat, col, row) {
  p  <- st_transform(st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"), crs = 4326), SINU)
  xy <- st_coordinates(p)
  floor((xy[, 1] - X0) / PX + 0.25) == col & floor((Y0 - xy[, 2]) / PX + 0.25) == row
}

# Earth Engine helpers
gi <- function(x, tries = 4) {   # getInfo, retrying when the servers are busy
  for (k in seq_len(tries)) {
    out <- tryCatch(x$getInfo(), error = function(e) e)
    if (!inherits(out, "error")) return(out)
    if (k == tries || !grepl("timed out|Too many|429|503|Internal", conditionMessage(out))) stop(out)
    message("  Earth Engine is busy; trying again in 60 seconds...")
    Sys.sleep(60)
  }
}
asset_exists <- function(id) !is.null(tryCatch(ee$data$getAsset(id), error = function(e) NULL))
ee_states <- function() {   # newest state of each task, if Earth Engine will tell us
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

# ---- 1. Every patch and its core pixels -------------------------------------
allp <- read.csv(file.path(outDir, "all_patches.csv")) %>%
  arrange(patch_n) %>%
  mutate(part = ceiling(patch_n / PART_SIZE))
nrow(allp)        # patches
max(allp$part)    # Earth Engine jobs (parts of PART_SIZE patches)

core <- read.csv(file.path(outDir, "patch_modis_pixels.csv")) %>%
  filter(patch_id %in% allp$patch_id, loss_frac >= CORE_FRAC) %>%
  left_join(allp %>% select(patch_id, part), by = "patch_id")
cent <- st_as_sf(data.frame(x = X0 + (core$col + 0.5) * PX, y = Y0 - (core$row + 0.5) * PX),
                 coords = c("x", "y"), crs = SINU) %>% st_transform(4326)
core_sf <- st_sf(core %>% select(patch_id, col, row, part), geometry = st_geometry(cent))
nrow(core_sf)                          # core pixels
summary(as.numeric(table(core$patch_id)))   # core pixels per patch

# ---- 2. Comparison rules for every MODIS pixel ------------------------------
# Defined here, worked out inside each patch job (section 4) for that patch's
# ring only. reduceResolution makes the area-weighted 30 m share for each
# MODIS pixel; the jobs ask for the MODIS grid, so that's the grid it uses.
modisProj <- ee$ImageCollection("MODIS/061/MOD13A1")$first()$select("NDVI")$projection()
to_modis  <- function(img) img$reduceResolution(reducer = ee$Reducer$mean(), maxPixels = 65535L)

hansen  <- ee$Image("UMD/hansen/global_forest_change_2025_v1_13")
lossAny <- hansen$select("loss")$unmask(0)
decid01 <- ee$ImageCollection("USGS/NLCD_RELEASES/2019_REL/NLCD")$
  filter(ee$Filter$eq("system:index", "2001"))$first()$
  select("landcover")$remap(list(41L, 43L), list(1L, 1L), 0)

anyLossM  <- to_modis(lossAny)$rename("any_loss")
nearLossM <- anyLossM$reduceNeighborhood(reducer = ee$Reducer$mean(),
                                         kernel  = ee$Kernel$square(radius = 1, units = "pixels"))$
  rename("near_loss")                                                  # the pixel and its 8 neighbours
decidM    <- to_modis(decid01)$rename("decid2001")
elevM     <- to_modis(ee$Image("USGS/SRTMGL1_003")$select("elevation"))$rename("elev")

passHard  <- decidM$gte(CTRL_DECID_MIN)
passLoss  <- anyLossM$lte(CTRL_LOSS_MAX)
passNear  <- nearLossM$lte(CTRL_NEAR_LOSS_MAX)
validCtrl <- passHard$And(passLoss)$And(passNear)

# ---- 3. Upload patch outlines and core pixels, in parts ---------------------
poly <- st_read(file.path(outDir, "all_patches.gpkg"), quiet = TRUE) %>%
  mutate(part = ceiling(patch_n / PART_SIZE)) %>%
  select(patch_id, part) %>%
  st_transform(5070) %>%
  st_simplify(dTolerance = SIMPLIFY_M, preserveTopology = TRUE) %>%
  st_transform(4326) %>%
  st_make_valid()

uploaded_ids <- function(id) {
  if (!asset_exists(id)) return(character(0))
  unlist(gi(ee$FeatureCollection(id)$aggregate_array("patch_id")))
}
for (k in sort(unique(allp$part))) {
  want <- allp$patch_id[allp$part == k]
  if (!setequal(unique(uploaded_ids(patch_asset(k))), want)) {
    message(sprintf("Uploading patch outlines, part %d of %d", k, max(allp$part)))
    sf_as_ee(poly[poly$part == k, "patch_id"], via = "getInfo_to_asset",
             assetId = patch_asset(k), overwrite = TRUE)
  }
  if (!setequal(unique(uploaded_ids(core_asset(k))), want)) {
    message(sprintf("Uploading core pixels, part %d of %d", k, max(allp$part)))
    sf_as_ee(core_sf[core_sf$part == k, c("patch_id", "col", "row")], via = "getInfo_to_asset",
             assetId = core_asset(k), overwrite = TRUE)
  }
}
message("All parts uploaded.")

# Quick check that our pixel grid matches Earth Engine's MODIS grid
chkFeats <- gi(ee$Image$pixelLonLat()$reduceRegions(
  collection = ee$FeatureCollection(core_asset(1))$limit(20),
  reducer    = ee$Reducer$first(),
  crs        = modisProj))$features
chk <- as.data.frame(t(sapply(chkFeats, function(f)
  as.numeric(unlist(f$properties[c("longitude", "latitude", "col", "row")])))))
names(chk) <- c("longitude", "latitude", "col", "row")
if (!all(grid_ok(chk$longitude, chk$latitude, chk$col, chk$row))) {
  stop("MODIS grid check failed: pixels don't line up with Earth Engine's grid. ",
       "Send this message to Claude before going further.", call. = FALSE)
}
message("MODIS grid check passed.")

# ---- 4. The per-patch job ---------------------------------------------------
# Time series image: one band per variable per year
bands  <- list()
qaNDVI <- ee_utils_pyfunc(function(img) {
  img$select("NDVI")$updateMask(img$select("SummaryQA")$lte(1))$multiply(0.0001)
})
mod13 <- ee$ImageCollection("MODIS/061/MOD13A1")
for (yr in NDVI_YEARS) {
  jul <- mod13$filterDate(sprintf("%d-07-01", yr), sprintf("%d-08-01", yr))$map(qaNDVI)
  bands[[length(bands) + 1]] <- jul$mean()$rename(paste0("ndvi_", yr))
}

mcd <- ee$ImageCollection("MODIS/061/MCD12Q2")
try_years <- 2001:MCD_LAST_YEAR
n_imgs <- unlist(gi(ee$List(lapply(try_years, function(y) {
  mcd$filter(ee$Filter$calendarRange(y, y, "year"))$size()
}))))
mcd_years <- try_years[as.numeric(n_imgs) > 0]
range(mcd_years)   # MCD12Q2 years Earth Engine has

pheno_in  <- c("Greenup_1", "Senescence_1", "MidGreendown_1", "Dormancy_1")
pheno_out <- c("greenup", "senescence", "midgreendown", "dormancy")
for (yr in mcd_years) {
  img <- mcd$filter(ee$Filter$calendarRange(yr, yr, "year"))$first()
  # MCD12Q2 dates are days since 1970-01-01; convert to day of year (Jan 1 = 1).
  off <- as.numeric(as.Date(sprintf("%d-01-01", yr)) - as.Date("1970-01-01")) - 1
  doy <- img$select(pheno_in)$subtract(off)
  doy <- doy$updateMask(doy$gte(1)$And(doy$lte(366)))$rename(paste0(pheno_out, "_", yr))
  qa  <- img$select("QA_Overall_1")
  qa  <- qa$updateMask(qa$lte(3))$rename(paste0("qa_", yr))
  bands[[length(bands) + 1]] <- doy
  bands[[length(bands) + 1]] <- qa
}
tsImg <- Reduce(function(a, b) a$addBands(b), bands)

ts_names    <- c(paste0("ndvi_", NDVI_YEARS),
                 paste0(rep(c(pheno_out, "qa"), times = length(mcd_years)), "_",
                        rep(mcd_years, each = 5)))
ctrl_names  <- paste0("ctrl_", ts_names)
share_names <- c("share_hardwood", "share_pixel_loss", "share_near_loss", "share_elevation", "share_controls")
count_names <- c("n_ring", "n_hardwood", "n_pixel_loss", "n_near_loss", "n_elevation", "n_controls")

per_patch <- ee_utils_pyfunc(function(f) {
  core_pts <- ee$FeatureCollection(ee$List(f$get("core")))$geometry()   # the patch's core pixel centres
  ring     <- f$geometry()$buffer(CONTROL_MAX_M, 50)$
    difference(f$geometry()$buffer(CONTROL_MIN_M, 50), 50)               # 1-10 km from the patch
  
  # Unweighted: every pixel counts once (no partial weights for pixels on the ring's edge)
  pElev    <- ee$Number(elevM$reduceRegion(
    reducer = ee$Reducer$mean()$unweighted(), geometry = core_pts, crs = modisProj,
    tileScale = TILE_SCALE)$get("elev"))
  elevOk   <- elevM$subtract(pElev)$abs()$lte(ELEV_TOL_M)
  ctrlMask <- validCtrl$And(elevOk)
  
  # One pass over the ring: comparison means each year, and the share of ring
  # pixels passing each rule
  shares <- ee$Image$cat(list(passHard, passLoss, passNear, elevOk, ctrlMask))$
    rename(as.list(share_names))
  cv <- tsImg$updateMask(ctrlMask)$rename(as.list(ctrl_names))$addBands(shares)$
    reduceRegion(reducer = ee$Reducer$mean()$unweighted(), geometry = ring, crs = modisProj,
                 maxPixels = 1e9, tileScale = TILE_SCALE)
  nRing <- ee$Image$constant(1)$rename("n_ring")$
    reduceRegion(reducer = ee$Reducer$sum()$unweighted(), geometry = ring, crs = modisProj,
                 maxPixels = 1e9)
  pv <- tsImg$reduceRegion(reducer = ee$Reducer$mean()$unweighted(), geometry = core_pts,
                           crs = modisProj, maxPixels = 1e9)
  ee$Feature(NULL, pv)$set(cv)$set(nRing)$
    set("patch_id", f$get("patch_id"))$set("elev_patch", pElev)
})
part_job <- function(k) {
  joined <- ee$Join$saveAll(matchesKey = "core")$apply(
    primary   = ee$FeatureCollection(patch_asset(k)),
    secondary = ee$FeatureCollection(core_asset(k)),
    condition = ee$Filter$equals(leftField = "patch_id", rightField = "patch_id"))
  ee$FeatureCollection(joined)$map(per_patch)
}

# ---- 5. Run the jobs and download (finished parts are skipped) --------------
export_name <- function(k) paste0("full_", VALUES_TAG, "_values_", k)
part_folder <- function(k) file.path(partDir, paste0("part_", k))
is_done     <- function(k) file.exists(file.path(part_folder(k), "done.txt"))
selectors   <- c("patch_id", "elev_patch", "n_ring", share_names, ts_names, ctrl_names)

download_ready <- function(parts) {
  f <- drive_ls(path = driveFolder, pattern = paste0("full_", VALUES_TAG, "_values_"))
  for (k in parts) {
    nm  <- export_name(k)
    hit <- f[f$name == paste0(nm, ".csv") | startsWith(f$name, paste0(nm, "-")), ]
    if (nrow(hit) == 0) next
    dir.create(part_folder(k), recursive = TRUE, showWarnings = FALSE)
    for (j in seq_len(nrow(hit))) {
      drive_download(as_id(hit$id[j]), path = file.path(part_folder(k), hit$name[j]), overwrite = TRUE)
    }
    writeLines(hit$name, file.path(part_folder(k), "done.txt"))
  }
}

all_parts <- sort(unique(allp$part))
todo <- all_parts[!vapply(all_parts, is_done, TRUE)]
length(todo)   # parts not downloaded yet (0 = all done, skips straight to section 6)
failed <- data.frame()

if (length(todo) > 0) {
  download_ready(todo)                       # anything an earlier run finished
  todo <- todo[!vapply(todo, is_done, TRUE)]
  
  st <- ee_states()
  waiting <- st$description[st$state %in% c("READY", "RUNNING")]
  n_started <- 0
  for (k in todo) {
    if (export_name(k) %in% waiting) next   # still running from an earlier session
    t <- ee_table_to_drive(collection = part_job(k), description = export_name(k),
                           folder = driveFolder, fileFormat = "CSV",
                           selectors = selectors, timePrefix = FALSE)
    t$start()
    n_started <- n_started + 1
  }
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
    message(sprintf("%s  parts left: %d  (running: %d, waiting: %d, failed: %d)",
                    format(Sys.time(), "%H:%M"), length(todo), sum(st$state == "RUNNING"),
                    sum(st$state == "READY"), nrow(failed)))
    if (length(todo) == 0 || length(todo) == nrow(failed)) break
    if (difftime(Sys.time(), t0, units = "hours") > 48) break
  }
}

failed   # parts that failed (empty = none). Source again to retry them.

have    <- all_parts[vapply(all_parts, is_done, TRUE)]
missing <- setdiff(all_parts, have)
if (length(missing) > 0) {
  warning(length(missing), " part(s) not downloaded yet: ", paste(missing, collapse = ", "),
          ". Source again to retry just those.", call. = FALSE)
}

# ---- 6. Patch and comparison series (in R) ----------------------------------
read_part <- function(k) {
  files <- list.files(part_folder(k), pattern = "\\.csv$", full.names = TRUE)
  bind_rows(lapply(files, function(p) read.csv(p, colClasses = c(patch_id = "character"))))
}
w <- bind_rows(lapply(have, read_part)) %>% distinct(patch_id, .keep_all = TRUE)
nrow(w)   # patches downloaded

psum <- allp %>%
  select(patch_id, patch_n, part, loss_year, area_ha, na_l2code, na_l2name, us_l3code, us_l3name,
         region, lat, lon, n_core, max_pixel_frac, core_loss_frac) %>%
  left_join(w %>% select(patch_id, elev_patch, n_ring, all_of(share_names)), by = "patch_id") %>%
  mutate(n_hardwood   = round(share_hardwood   * n_ring),   # pixels in the 1-10 km ring passing each rule
         n_pixel_loss = round(share_pixel_loss * n_ring),
         n_near_loss  = round(share_near_loss  * n_ring),
         n_elevation  = round(share_elevation  * n_ring),
         n_controls   = round(share_controls   * n_ring)) %>%
  mutate(downloaded  = patch_id %in% w$patch_id,
         n_controls  = coalesce(n_controls, 0),
         enough_ctrl = n_controls >= MIN_CONTROLS,
         in_analysis = downloaded & n_core > 0 & enough_ctrl,
         status = case_when(!downloaded  ~ "not downloaded yet",
                            n_core == 0  ~ "not used: no pixel >= 50% lost",
                            !enough_ctrl ~ paste("not used: fewer than", MIN_CONTROLS, "comparison pixels"),
                            TRUE         ~ "in analysis"))

long <- w %>%
  select(patch_id, matches("_\\d{4}$")) %>%
  pivot_longer(-patch_id, names_to = c("metric", "year"), names_pattern = "^(.*)_(\\d{4})$") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(year = as.integer(year))

series <- long %>%
  left_join(psum %>% select(patch_id, loss_year, area_ha, na_l2code, us_l3code, n_controls, in_analysis),
            by = "patch_id") %>%
  mutate(years_since_loss  = year - loss_year,
         ndvi_diff         = ndvi - ctrl_ndvi,
         greenup_diff      = greenup - ctrl_greenup,
         senescence_diff   = senescence - ctrl_senescence,
         midgreendown_diff = midgreendown - ctrl_midgreendown,
         dormancy_diff     = dormancy - ctrl_dormancy) %>%
  arrange(patch_id, year)

write.csv(psum,   file.path(tsDir, "patch_summary.csv"),     row.names = FALSE)
write.csv(series, file.path(tsDir, "patch_year_series.csv"), row.names = FALSE)

# ---- 7. Summary + first-look figure -----------------------------------------
table(psum$status)   # every patch: in the analysis or why not

# By Level II region: patches, how many have enough comparison pixels, and
# how many comparison pixels a patch typically has
comp_by_region <- psum %>%
  filter(downloaded) %>%
  group_by(region) %>%
  summarise(patches          = n(),
            in_analysis      = sum(in_analysis),
            pct_in_analysis  = round(100 * mean(in_analysis), 1),
            median_controls  = median(n_controls),
            .groups = "drop") %>%
  as.data.frame()
comp_by_region

summary(psum$n_controls[psum$in_analysis])   # comparison pixels per patch in the analysis

fig <- series %>%
  filter(in_analysis, years_since_loss >= -5, years_since_loss <= 10) %>%
  select(patch_id, years_since_loss,
         `July NDVI` = ndvi_diff, `MidGreendown (days)` = midgreendown_diff) %>%
  pivot_longer(-c(patch_id, years_since_loss), names_to = "metric", values_to = "diff") %>%
  filter(!is.na(diff)) %>%
  group_by(metric, years_since_loss) %>%
  summarise(n = n(), mean = mean(diff), se = sd(diff) / sqrt(n()), .groups = "drop") %>%
  filter(n >= 5)

ggplot(data=fig) +
  facet_wrap(~ metric, ncol=1, scales="free_y") + # one panel per variable
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dashed", color="red") + # year of loss
  geom_ribbon(aes(x=years_since_loss, ymin=mean - 1.96*se, ymax=mean + 1.96*se), fill="darkgreen", alpha=0.2) + # 95% CI
  geom_line(aes(x=years_since_loss, y=mean), color="darkgreen") +
  geom_point(aes(x=years_since_loss, y=mean), color="darkgreen", size=1) +
  theme_bw() +
  labs(title="Disturbed patch minus nearby uncut forest, every patch",
       x="Years since loss (0 = year of loss)", y="Patch minus comparison pixels",
       caption="Shaded band = 95% confidence interval of the mean across patches (not model-adjusted).")

ggsave(file.path(tsDir, "first_look_patch_minus_control.png"), width=8, height=7, dpi=200)