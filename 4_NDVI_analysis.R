# Patch Step 4: satellite time series for the sampled patches + nearby uncut forest
# -----------------------------------------------------------------------------
# Open in RStudio and use Source with Echo so the tables
# and figure show.
# Needs: rgee, googledrive, sf, dplyr, tidyr, ggplot2
# Reads:  scaleup_patches/sample_patches.gpkg   (made by Patch Step 3)
# Saves (in scaleup_patches/timeseries/):
#   patch_year_series.csv     one row per patch per year; use rows with in_analysis = TRUE
#   patch_summary.csv         one row per drawn patch: pixels, controls, in_analysis
#   pixel_values_long.csv     every pixel, every year
#   patch_pixels_checked.csv, controls_checked.csv, pixel_values_wide.csv
#   first_look_patch_minus_control.png
#
# What it does:
#   1. For each sampled patch (and its replacements), finds the MODIS pixels it
#      touches and works out what share of each pixel was lost ("loss_frac").
#      Step 3 already made sure every patch covers >= 50% of at least one pixel.
#      Pixels are named by their place in the global MODIS grid, so the same
#      pixel always gets the same ID.
#   2. Picks candidate CONTROL pixels for each patch: MODIS pixels 1-10 km away,
#      drawn at random (seeded, so it's repeatable).
#   3. In Earth Engine, for every pixel:
#        - elevation, share that was hardwood in 2001 (NLCD 41/43), share of the
#          pixel ever lost (Hansen, any year), share lost within 500 m around it
#        - July NDVI every year 2001-2024 (MOD13A1, good/marginal quality only)
#          and how many July composites went into it
#        - fall timing every year (MCD12Q2): Greenup, Senescence, MidGreendown,
#          Dormancy as day of year (Jan 1 = 1), plus the quality flag
#      Values are read straight from each MODIS pixel -- no resampling.

#### 
#Greenup, Senescence,Dormancy are relics from previous code and can be removed as needed
####
#   4.in R:
#        - keeps controls that were >= 75% hardwood in 2001, almost never lost
#          (<= 5% of the pixel), little loss within 500 m (<= 10%), and within
#          100 m elevation of the patch; up to 10 per patch
#        - patch value each year = mean of its "core" pixels (>= 50% lost)
#        - control value each year = mean of that patch's controls
#        -  For analysis I walk down the random draw order and take the first
#          100 patches that have at least 3 controls. A sampled patch without
#          enough uncut forest nearby is replaced by the next patch in the order,
#          so the 100 stay a random sample of patches that meet every rule.
#   5. Saves the tables and a first-look figure: patch minus controls, by years
#      since the loss, for the 100 patches in the analysis set.
#
# The Earth Engine part runs once and is saved. To change any threshold in
# the Settings and re-run only the R part, just Source again (it reuses the
# saved file). If the sample from Step 3 changes, the Earth Engine part is
# redone automatically. Set FORCE_EE <- TRUE to redo it anyway.

# IMPORTANT! 

#If R stops while Earth Engine is still working, don't start over: copy the
# export name into RESUME_EXPORT and Source again -- it waits for that export.


# Status lines while waiting on Earth Engine use message() 
# =============================================================================

library(rgee)
library(googledrive)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)

ee_Initialize(drive = TRUE)

# ---- Settings ---------------------------------------------------------------
SEED                 <- 42
CONTROL_MIN_M        <- 1000   # controls at least this far from the patch
CONTROL_MAX_M        <- 10000  # ... and no farther than this
CANDIDATES_PER_PATCH <- 150    # candidate control pixels checked per patch
N_CONTROLS           <- 10     # controls kept per patch
MIN_CONTROLS         <- 3      # a patch needs at least this many controls to be used
N_ANALYSIS           <- 100    # patches in the analysis set
CORE_FRAC            <- 0.5    # a "core" pixel is at least this share lost (= Step 3's rule)
CTRL_DECID_MIN       <- 0.75   # control: hardwood share in 2001
CTRL_LOSS_MAX        <- 0.05   # control: share of the pixel ever lost (5% = ~1 ha)
CTRL_NEAR_LOSS_MAX   <- 0.10   # control: share lost within NEAR_M around it
NEAR_M               <- 500    # about one pixel's width: the range of MODIS signal spill-over
ELEV_TOL_M           <- 100    # control: elevation within this of the patch
NDVI_YEARS           <- 2001:2024
MCD_LAST_YEAR        <- 2023   # newest MCD12Q2 year to use (checked against what exists)
FORCE_EE             <- FALSE
RESUME_EXPORT        <- ""     # if an export was already started, put its name here
# (e.g. "sample_pixel_values_20261021_082751") to
# wait for it instead of starting a new one

outDir      <- file.path("~/Google Drive/My Drive/Reidy_research", "scaleup_patches")
tsDir       <- file.path(outDir, "timeseries")
dir.create(tsDir, recursive = TRUE, showWarnings = FALSE)
driveFolder <- "Reidy_research"
cellAsset   <- "projects/breidyee/assets/scaleup_sample_cells"
ts_file     <- file.path(tsDir, "pixel_values_wide.csv")

# MODIS sinusoidal grid (sphere R = 6371007.181 m; 2400 x 463.3 m pixels per tile)
SINU <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
X0   <- -20015109.354
Y0   <-  10007554.677
PX   <-  463.312716528

cell_square <- function(col, row) {
  x0 <- X0 + col * PX;  x1 <- x0 + PX
  y1 <- Y0 - row * PX;  y0 <- y1 - PX
  st_polygon(list(matrix(c(x0, y0, x1, y0, x1, y1, x0, y1, x0, y0), ncol = 2, byrow = TRUE)))
}
cells_sf <- function(df) {   # df needs columns col, row
  st_sf(df, geometry = st_sfc(mapply(cell_square, df$col, df$row, SIMPLIFY = FALSE), crs = SINU))
}
# TRUE where a lon/lat (from Earth Engine) falls in the expected MODIS pixel
grid_ok <- function(lon, lat, col, row) {
  p  <- st_transform(st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"), crs = 4326), SINU)
  xy <- st_coordinates(p)
  floor((xy[, 1] - X0) / PX + 0.25) == col & floor((Y0 - xy[, 2]) / PX + 0.25) == row
}

# ---- 1 & 2. Patch pixels and candidate control pixels (in R) ---------------
smp <- st_read(file.path(outDir, "sample_patches.gpkg"), quiet = TRUE) %>%
  st_make_valid() %>%
  arrange(draw_order)
table(smp$role)   # sampled patches and replacements

smp_sin <- st_transform(smp, SINU)
smp_alb <- st_transform(smp, 5070)   # equal-distance-ish, for the 1-10 km ring

patch_px <- list()
cand     <- list()
set.seed(SEED)
for (i in seq_len(nrow(smp))) {
  g  <- st_geometry(smp_sin)[i]
  bb <- unname(st_bbox(g))   # xmin, ymin, xmax, ymax
  
  # MODIS pixels the patch touches, and the share of each that is inside it
  cr <- expand.grid(col = floor((bb[1] - X0) / PX):floor((bb[3] - X0) / PX),
                    row = floor((Y0 - bb[4]) / PX):floor((Y0 - bb[2]) / PX))
  ov <- suppressWarnings(st_intersection(cells_sf(cr), g))
  if (nrow(ov) > 0) {
    patch_px[[i]] <- data.frame(patch_id = smp$patch_id[i], col = ov$col, row = ov$row,
                                loss_frac = as.numeric(st_area(ov)) / PX^2)
  }
  
  # Candidate controls: pixel centres 1-10 km from the patch. The search box is
  # padded 2x because sinusoidal coordinates stretch distances up to ~1.75x here.
  pad <- 2 * CONTROL_MAX_M
  cr2 <- expand.grid(col = floor((bb[1] - pad - X0) / PX):floor((bb[3] + pad - X0) / PX),
                     row = floor((Y0 - bb[4] - pad) / PX):floor((Y0 - bb[2] + pad) / PX))
  ctr <- st_as_sf(data.frame(cr2, x = X0 + (cr2$col + 0.5) * PX, y = Y0 - (cr2$row + 0.5) * PX),
                  coords = c("x", "y"), crs = SINU)
  d   <- as.numeric(st_distance(st_transform(ctr, 5070), st_geometry(smp_alb)[i]))
  ok  <- cr2[d > CONTROL_MIN_M & d <= CONTROL_MAX_M, ]
  ok  <- ok[order(ok$col, ok$row), ]                    # fixed order before the random draw
  pick <- ok[sample(nrow(ok), min(CANDIDATES_PER_PATCH, nrow(ok))), ]
  cand[[i]] <- data.frame(patch_id = smp$patch_id[i], col = pick$col, row = pick$row,
                          cand_order = seq_len(nrow(pick)))
}

patch_px <- bind_rows(patch_px) %>% mutate(cell_id = sprintf("c%d_%d", col, row))
cand     <- bind_rows(cand)     %>% mutate(cell_id = sprintf("c%d_%d", col, row))

cells <- bind_rows(patch_px %>% select(cell_id, col, row),
                   cand     %>% select(cell_id, col, row)) %>%
  distinct(cell_id, .keep_all = TRUE) %>%
  arrange(cell_id)
cent <- st_as_sf(data.frame(x = X0 + (cells$col + 0.5) * PX, y = Y0 - (cells$row + 0.5) * PX),
                 coords = c("x", "y"), crs = SINU) %>% st_transform(4326)
ll <- st_coordinates(cent)
cells$lon <- ll[, 1]
cells$lat <- ll[, 2]

write.csv(patch_px, file.path(tsDir, "patch_pixels.csv"),        row.names = FALSE)
write.csv(cand,     file.path(tsDir, "control_candidates.csv"),  row.names = FALSE)

nrow(patch_px)                     # MODIS pixels touched by the patches
median(table(patch_px$patch_id))   # median MODIS pixels per patch
nrow(cand)                         # candidate control pixels
nrow(cells)                        # unique pixels to look up in Earth Engine

# ---- 3. Earth Engine: look up every pixel (runs once; result is saved) -------
# Saved results are reused only if they are for exactly these pixels. If the
# sample or the control settings changed, they are redone automatically.
stale <- FALSE
if (file.exists(ts_file) && !FORCE_EE) {
  saved_ids <- read.csv(ts_file)$cell_id
  stale <- !setequal(saved_ids, cells$cell_id)
  if (stale) {
    message(sprintf(paste0("Saved Earth Engine results are for a different set of pixels ",
                           "(%d saved vs %d needed) -- redoing them."),
                    length(saved_ids), nrow(cells)))
  }
}
if (FORCE_EE || stale || !file.exists(ts_file)) {
  if (nzchar(RESUME_EXPORT)) {
    nm   <- RESUME_EXPORT
    task <- NULL
    message("Not starting a new export -- waiting for the one already started: ", nm)
  } else {
    
    # Interactive Earth Engine requests sometimes time out when the servers are
    # busy. gi() retries a few times before giving up.
    gi <- function(x, tries = 4) {
      for (k in seq_len(tries)) {
        out <- tryCatch(x$getInfo(), error = function(e) e)
        if (!inherits(out, "error")) return(out)
        if (k == tries || !grepl("timed out|Too many|429|503|Internal", conditionMessage(out))) stop(out)
        message("  Earth Engine is busy; trying again in 60 seconds...")
        Sys.sleep(60)
      }
    }
    
    # 3a. Upload the pixel squares, in parts small enough for Earth Engine.
    #     A part already uploaded with exactly the same pixels is reused.
    cells_poly <- cells_sf(cells) %>% st_transform(4326)
    parts_idx  <- split(seq_len(nrow(cells_poly)), ceiling(seq_len(nrow(cells_poly)) / 6000))
    part_ids   <- paste0(cellAsset, "_", seq_along(parts_idx))
    for (k in seq_along(parts_idx)) {
      want <- cells_poly$cell_id[parts_idx[[k]]]
      have <- tryCatch(unlist(gi(ee$FeatureCollection(part_ids[k])$aggregate_array("cell_id"))),
                       error = function(e) character(0))
      if (setequal(have, want)) {
        message(sprintf("Pixels part %d of %d: already uploaded", k, length(parts_idx)))
        next
      }
      message(sprintf("Uploading pixels, part %d of %d", k, length(parts_idx)))
      sf_as_ee(cells_poly[parts_idx[[k]], c("cell_id", "col", "row", "lon", "lat")],
               via = "getInfo_to_asset", assetId = part_ids[k], overwrite = TRUE)
    }
    partFCs <- lapply(part_ids, function(a) ee$FeatureCollection(a))
    n_up    <- sum(vapply(partFCs, function(fc) as.numeric(gi(fc$size())), 0))
    if (n_up != nrow(cells)) stop("Uploaded ", n_up, " pixels, expected ", nrow(cells), call. = FALSE)
    cellsFC <- Reduce(function(a, b) a$merge(b), partFCs)
    
    modisProj <- ee$ImageCollection("MODIS/061/MOD13A1")$first()$select("NDVI")$projection()
    centroid  <- ee_utils_pyfunc(function(f) f$setGeometry(f$geometry()$centroid(maxError = 1)))
    
    # 3b. Quick check that our pixel grid matches Earth Engine's MODIS grid
    chkFeats <- gi(ee$Image$pixelLonLat()$reduceRegions(
      collection = partFCs[[1]]$limit(20)$map(centroid),
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
    
    # 3c. Pixel checks at 30 m: loss ever, hardwood 2001, elevation, loss nearby
    hansen  <- ee$Image("UMD/hansen/global_forest_change_2024_v1_12")
    hProj   <- hansen$select("loss")$projection()
    lossAny <- hansen$select("loss")$unmask(0)
    decid01 <- ee$ImageCollection("USGS/NLCD_RELEASES/2019_REL/NLCD")$
      filter(ee$Filter$eq("system:index", "2001"))$first()$
      select("landcover")$remap(list(41L, 43L), list(1L, 1L), 0)
    elev    <- ee$Image("USGS/SRTMGL1_003")$select("elevation")
    
    statImg <- lossAny$rename("any_loss")$
      addBands(decid01$rename("decid2001"))$
      addBands(elev$rename("elev"))
    s1 <- statImg$reduceRegions(collection = cellsFC, reducer = ee$Reducer$mean(),
                                crs = hProj, tileScale = 4)
    s2 <- s1$map(ee_utils_pyfunc(function(f) f$setGeometry(f$geometry()$buffer(NEAR_M, 10))))
    s2 <- lossAny$reduceRegions(collection = s2,
                                reducer = ee$Reducer$mean()$setOutputs(list("near_loss")),
                                crs = hProj, tileScale = 4)
    pts <- s2$map(centroid)
    
    # 3d. Time series image: one band per variable per year
    bands <- list()
    qaNDVI <- ee_utils_pyfunc(function(img) {
      img$select("NDVI")$updateMask(img$select("SummaryQA")$lte(1))$multiply(0.0001)
    })
    mod13 <- ee$ImageCollection("MODIS/061/MOD13A1")
    for (yr in NDVI_YEARS) {
      jul <- mod13$filterDate(sprintf("%d-07-01", yr), sprintf("%d-08-01", yr))$map(qaNDVI)
      bands[[length(bands) + 1]] <- jul$mean()$rename(paste0("ndvi_", yr))
      bands[[length(bands) + 1]] <- jul$count()$rename(paste0("ndvin_", yr))
    }
    
    mcd <- ee$ImageCollection("MODIS/061/MCD12Q2")
    # Which years exist: ask Earth Engine how many images each year has (one call)
    try_years <- 2001:MCD_LAST_YEAR
    n_imgs <- unlist(gi(ee$List(lapply(try_years, function(y) {
      mcd$filter(ee$Filter$calendarRange(y, y, "year"))$size()
    }))))
    mcd_years <- try_years[as.numeric(n_imgs) > 0]
    if (length(mcd_years) == 0) {
      stop("Earth Engine reported no MCD12Q2 images for 2001-", MCD_LAST_YEAR,
           ". Send this message to Claude.", call. = FALSE)
    }
    message("MCD12Q2 years used: ", min(mcd_years), "-", max(mcd_years),
            if (length(mcd_years) < length(try_years))
              paste0(" (missing: ", paste(setdiff(try_years, mcd_years), collapse = ", "), ")"))
    
    pheno_in  <- c("Greenup_1", "Senescence_1", "MidGreendown_1", "Dormancy_1")
    pheno_out <- c("greenup", "senescence", "midgreendown", "dormancy")
    for (yr in mcd_years) {
      img <- mcd$filter(ee$Filter$calendarRange(yr, yr, "year"))$first()
      # MCD12Q2 dates are days since 1970-01-01; convert to day of year (Jan 1 = 1).
      # (The pilot used Jan 1 = 0, so pilot DOYs are 1 lower.)
      off <- as.numeric(as.Date(sprintf("%d-01-01", yr)) - as.Date("1970-01-01")) - 1
      doy <- img$select(pheno_in)$subtract(off)
      doy <- doy$updateMask(doy$gte(1)$And(doy$lte(366)))$rename(paste0(pheno_out, "_", yr))
      qa  <- img$select("QA_Overall_1")
      qa  <- qa$updateMask(qa$lte(3))$rename(paste0("qa_", yr))
      bands[[length(bands) + 1]] <- doy
      bands[[length(bands) + 1]] <- qa
    }
    bands[[length(bands) + 1]] <- ee$Image$pixelLonLat()$rename(c("chk_lon", "chk_lat"))
    tsImg <- Reduce(function(a, b) a$addBands(b), bands)
    
    # 3e. Read every band at each pixel centre, on the MODIS grid itself
    out <- tsImg$reduceRegions(collection = pts, reducer = ee$Reducer$first(),
                               crs = modisProj, tileScale = 4)
    
    ts_cols <- c(paste0(c("ndvi_", "ndvin_"), rep(NDVI_YEARS, each = 2)),
                 paste0(rep(c(pheno_out, "qa"), times = length(mcd_years)), "_",
                        rep(mcd_years, each = 5)),
                 "chk_lon", "chk_lat")
    sel <- c("cell_id", "col", "row", "lon", "lat",
             "any_loss", "decid2001", "elev", "near_loss", ts_cols)
    
    nm   <- paste0("sample_pixel_values_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    task <- ee_table_to_drive(collection = out, description = nm, folder = driveFolder,
                              fileFormat = "CSV", selectors = sel, timePrefix = FALSE)
    task$start()
    message("Earth Engine task started: ", nm, "\n",
            "If R stops before it finishes, DON'T start a new one: put this name in\n",
            "RESUME_EXPORT at the top of the script and Source again.")
  }
  
  # 3f. Wait for the export: check every 2 minutes
  t0 <- Sys.time()
  repeat {
    if (!is.null(task)) {
      st <- task$status()$state
      if (st %in% c("FAILED", "CANCELLED")) {
        stop("Earth Engine task ended as ", st, ": ", task$status()$error_message, call. = FALSE)
      }
    } else {
      st <- "see Tasks tab"
    }
    f <- drive_ls(path = driveFolder, pattern = nm)
    if (nrow(f) > 0) break
    if (difftime(Sys.time(), t0, units = "hours") > 12) {
      stop("Still no file in Drive after 12 hours. Check the Tasks tab in the GEE Code Editor.",
           call. = FALSE)
    }
    message(format(Sys.time(), "%H:%M"), "  export status: ", st, " -- waiting...")
    Sys.sleep(120)
  }
  Sys.sleep(60)                                # let every part of a split export arrive
  f <- drive_ls(path = driveFolder, pattern = nm)
  message(sprintf("Export found in Drive: %d file(s). Downloading...", nrow(f)))
  parts <- lapply(seq_len(nrow(f)), function(k) {
    p <- file.path(tempdir(), f$name[k])
    drive_download(as_id(f$id[k]), path = p, overwrite = TRUE)
    read.csv(p)
  })
  got <- bind_rows(parts)
  if (!setequal(got$cell_id, cells$cell_id)) {
    stop("The downloaded export (", nm, ") is for a different set of pixels than the current ",
         "sample. If RESUME_EXPORT is set, set it back to \"\" and Source again.", call. = FALSE)
  }
  write.csv(got, ts_file, row.names = FALSE)
} else {
  message("Using saved Earth Engine results (they match the current sample): ", ts_file)
}

# ---- 4. Patch and control series (in R) -------------------------------------
w <- read.csv(ts_file)
if (!setequal(w$cell_id, cells$cell_id)) {
  stop("Saved results don't match the current sample. Set FORCE_EE <- TRUE and Source again.",
       call. = FALSE)
}
aligned <- grid_ok(w$chk_lon, w$chk_lat, w$col, w$row)
table(aligned, useNA = "ifany")   # MODIS grid check on every pixel: all should be TRUE

long <- w %>%
  select(cell_id, matches("_\\d{4}$")) %>%
  pivot_longer(-cell_id, names_to = c("metric", "year"), names_pattern = "(.*)_(\\d{4})$") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(year = as.integer(year))
write.csv(long, file.path(tsDir, "pixel_values_long.csv"), row.names = FALSE)

cellstat   <- w %>% select(cell_id, any_loss, decid2001, elev, near_loss)
patch_info <- st_drop_geometry(smp) %>%
  select(patch_id, role, draw_order, loss_year, area_ha, na_l2code, us_l3code)

# Patch pixels; "other_loss" = loss in the pixel from other years / other patches
pp <- patch_px %>%
  left_join(cellstat, by = "cell_id") %>%
  mutate(other_loss = pmax(any_loss - loss_frac, 0))

psum <- pp %>%
  group_by(patch_id) %>%
  summarise(n_pixels      = n(),
            n_core        = sum(loss_frac >= CORE_FRAC),
            max_loss_frac = max(loss_frac),
            elev_patch    = weighted.mean(elev, loss_frac, na.rm = TRUE),
            .groups = "drop")

used <- pp %>% filter(loss_frac >= CORE_FRAC)   # core pixels only

ctrl <- cand %>%
  left_join(cellstat, by = "cell_id") %>%
  left_join(psum %>% select(patch_id, elev_patch), by = "patch_id") %>%
  mutate(passes = coalesce(decid2001 >= CTRL_DECID_MIN &
                             any_loss  <= CTRL_LOSS_MAX &
                             near_loss <= CTRL_NEAR_LOSS_MAX &
                             abs(elev - elev_patch) <= ELEV_TOL_M, FALSE)) %>%
  arrange(patch_id, cand_order) %>%
  group_by(patch_id) %>%
  mutate(kept = passes & cumsum(passes) <= N_CONTROLS) %>%
  ungroup()

psum <- patch_info %>%
  left_join(psum, by = "patch_id") %>%
  left_join(ctrl %>% group_by(patch_id) %>%
              summarise(n_candidates = n(), n_controls = sum(kept), .groups = "drop"),
            by = "patch_id") %>%
  mutate(has_core    = coalesce(n_core > 0, FALSE),
         enough_ctrl = coalesce(n_controls >= MIN_CONTROLS, FALSE),
         eligible    = has_core & enough_ctrl) %>%
  arrange(draw_order) %>%
  mutate(in_analysis = eligible & cumsum(eligible) <= N_ANALYSIS,
         status = case_when(in_analysis  ~ "in analysis",
                            !has_core    ~ "not used: no pixel >= 50% lost",
                            !enough_ctrl ~ paste("not used: fewer than", MIN_CONTROLS, "controls"),
                            TRUE         ~ "not needed (spare replacement)"))

vars   <- c("ndvi", "ndvin", "greenup", "senescence", "midgreendown", "dormancy", "qa")
nan_na <- function(x) ifelse(is.nan(x), NA_real_, x)

patch_series <- suppressWarnings(used %>% select(patch_id, cell_id) %>% inner_join(long, by = "cell_id")) %>%
  group_by(patch_id, year) %>%
  summarise(across(all_of(vars), ~ nan_na(mean(.x, na.rm = TRUE))), n_px = n(), .groups = "drop")

ctrl_series <- suppressWarnings(ctrl %>% filter(kept) %>% select(patch_id, cell_id) %>%
                                  inner_join(long, by = "cell_id")) %>%
  group_by(patch_id, year) %>%
  summarise(across(all_of(vars), ~ nan_na(mean(.x, na.rm = TRUE)), .names = "ctrl_{.col}"),
            n_ctrl = n(), .groups = "drop")

series <- patch_series %>%
  full_join(ctrl_series, by = c("patch_id", "year")) %>%
  left_join(patch_info, by = "patch_id") %>%
  left_join(psum %>% select(patch_id, in_analysis), by = "patch_id") %>%
  mutate(years_since_loss  = year - loss_year,
         ndvi_diff         = ndvi - ctrl_ndvi,
         greenup_diff      = greenup - ctrl_greenup,
         senescence_diff   = senescence - ctrl_senescence,
         midgreendown_diff = midgreendown - ctrl_midgreendown,
         dormancy_diff     = dormancy - ctrl_dormancy) %>%
  arrange(draw_order, year)

write.csv(psum,   file.path(tsDir, "patch_summary.csv"),        row.names = FALSE)
write.csv(ctrl,   file.path(tsDir, "controls_checked.csv"),     row.names = FALSE)
write.csv(pp,     file.path(tsDir, "patch_pixels_checked.csv"), row.names = FALSE)
write.csv(series, file.path(tsDir, "patch_year_series.csv"),    row.names = FALSE)

# ---- 5. Summary + first-look figure -----------------------------------------
a <- psum %>% filter(in_analysis)
nrow(a)          # patches in the analysis set (target N_ANALYSIS)
table(a$role)    # original picks (Sample) and replacements used
if (nrow(a) < N_ANALYSIS) {
  warning("Ran out of replacements: only ", nrow(a), " patches meet every rule. ",
          "Raise N_REPLACEMENTS in Step 3 (the first draws stay the same), re-run Step 3, ",
          "then re-run this script (it redoes the Earth Engine part by itself).", call. = FALSE)
}
table(psum$status, psum$role)   # status of every drawn patch, by original role
table(psum$enough_ctrl)         # drawn patches with at least MIN_CONTROLS controls (TRUE) or not (FALSE)

# Share of candidate control pixels passing each check on its own
ctrl_checks <- round(colMeans(cbind(
  hardwood   = ctrl$decid2001 >= CTRL_DECID_MIN,
  pixel_loss = ctrl$any_loss <= CTRL_LOSS_MAX,
  near_loss  = ctrl$near_loss <= CTRL_NEAR_LOSS_MAX,
  elevation  = abs(ctrl$elev - ctrl$elev_patch) <= ELEV_TOL_M,
  all        = ctrl$passes), na.rm = TRUE), 2)
ctrl_checks

table(a$n_core)         # core pixels per patch in the analysis set
summary(a$n_controls)   # controls kept per patch in the analysis set
table(a$na_l2code)      # analysis set by Level II region

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
  labs(title="Disturbed patch minus nearby uncut forest, by years since the loss",
       x="Years since loss (0 = year of loss)", y="Patch minus controls")

ggsave(file.path(tsDir, "first_look_patch_minus_control.png"), width=8, height=7, dpi=200)