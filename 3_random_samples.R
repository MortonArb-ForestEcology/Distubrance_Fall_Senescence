# Patch Step 3: the full set of patches (no random sample)
# -----------------------------------------------------------------------------
# Run line by line or with source with Echo so the tables and plots show. A plain
# "Source" still saves every file but won't show them.
# Needs: sf, dplyr, ggplot2
# Reads:  scaleup_patches/patches_all_labeled.csv     (Patch Step 1)
#         scaleup_patches/ecoregion_list.csv          (Patch Step 0, which ecoregions are in)
#         scaleup_patches/patches_qualifying.gpkg     (Patch Step 1, polygons)
#         scaleup_patches/boundary_ecoregions.gpkg    (Patch Step 0, for the map)
# Saves:  all_patches.csv / .gpkg       every patch that meets the rules (Step 4 uses these)
#         patch_modis_pixels.csv        every MODIS pixel each patch touches, and the share it covers
#         map_all_patches.png, all_patches_by_year.png
#
# What it does:
#   1. Starts from the qualifying patches (Step 1) at least MIN_PATCH_HA, in
#      the ecoregions Step 0 currently includes
#   2. MODIS RULE: keeps patches that cover at least half of at least one
#      MODIS pixel, worked out exactly from each patch's outline laid over the
#      MODIS grid. The pixels a patch covers at least half of are its "core"
#      pixels; Step 4 reads the patch's values from those.
#   3. Keeps EVERY patch that meets the rules. There is no random sample and no
#      seed: the analysis uses all of them.
#   4. Shows the counts by region, size and year of loss, and maps them
#
# Replaces patch_step3_random_sample.R.
# The MODIS-pixel check is saved, so re-running is fast.
# Next: patch_step4_timeseries.R
# =============================================================================

library(sf)
library(dplyr)
library(ggplot2)

# ---- Settings ---------------------------------------------------------------
MIN_PATCH_HA   <- 25     # size floor (a patch must be >= 10.7 ha to pass the MODIS rule anyway)
MIN_PIXEL_FRAC <- 0.5    # patch must cover at least this share of at least one MODIS pixel
CHUNK          <- 2000   # patches checked against the MODIS grid at a time

outDir  <- file.path("~/Google Drive/My Drive/Reidy_research", "scaleup_patches")
px_file <- file.path(outDir, "patch_modis_pixels.csv")

# MODIS sinusoidal grid (same as Step 4)
SINU <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
X0   <- -20015109.354
Y0   <-  10007554.677
PX   <-  463.312716528

cell_square <- function(col, row) {
  x0 <- X0 + col * PX;  x1 <- x0 + PX
  y1 <- Y0 - row * PX;  y0 <- y1 - PX
  st_polygon(list(matrix(c(x0, y0, x1, y0, x1, y1, x0, y1, x0, y0), ncol = 2, byrow = TRUE)))
}

# Every MODIS pixel a patch touches, and the share of that pixel inside the patch
pixel_fracs <- function(poly_sin) {
  bb <- t(sapply(st_geometry(poly_sin), st_bbox))   # xmin, ymin, xmax, ymax
  cr <- bind_rows(lapply(seq_len(nrow(bb)), function(i) {
    expand.grid(col = floor((bb[i, 1] - X0) / PX):floor((bb[i, 3] - X0) / PX),
                row = floor((Y0 - bb[i, 4]) / PX):floor((Y0 - bb[i, 2]) / PX))
  })) %>% distinct()
  cells <- st_sf(cr, geometry = st_sfc(mapply(cell_square, cr$col, cr$row, SIMPLIFY = FALSE), crs = SINU))
  ov <- suppressWarnings(st_intersection(poly_sin[, "patch_id"], cells))
  data.frame(patch_id = ov$patch_id, col = ov$col, row = ov$row,
             loss_frac = as.numeric(st_area(ov)) / PX^2)
}

# ---- 1. Qualifying patches above the size floor -----------------------------
pt   <- read.csv(file.path(outDir, "patches_all_labeled.csv"))
pool <- pt %>% filter(qualifies, area_ha >= MIN_PATCH_HA)

# Follow Step 0's CURRENT list of included ecoregions. Leaving ecoregions out
# in Step 0 (ALSO_EXCLUDE_L2 / ALSO_EXCLUDE_L3) works here without re-running
# Step 1. Adding ecoregions in Step 0 does need Step 1 re-run.
eco_list  <- read.csv(file.path(outDir, "ecoregion_list.csv"))
inc_codes <- as.character(eco_list$us_l3code[as.logical(eco_list$included)])
left_out  <- pool %>% filter(!(as.character(us_l3code) %in% inc_codes))
table(paste(left_out$us_l3code, left_out$us_l3name))   # patches left out: their ecoregion is no longer in the study area

pool <- pool %>% filter(as.character(us_l3code) %in% inc_codes)
nrow(pool)   # qualifying patches >= MIN_PATCH_HA in the study area

# ---- 2. MODIS rule: covers >= half of at least one MODIS pixel --------------
# A patch smaller than half a pixel (~10.7 ha) can't pass, so skip those.
need <- pool$patch_id[pool$area_ha >= MIN_PIXEL_FRAC * PX^2 / 1e4]

px <- if (file.exists(px_file)) read.csv(px_file) else
  data.frame(patch_id = character(), col = integer(), row = integer(), loss_frac = numeric())
todo <- setdiff(need, px$patch_id)
length(todo)   # patches still to check against the MODIS grid (0 = all saved from before)

if (length(todo) > 0) {
  poly <- st_read(file.path(outDir, "patches_qualifying.gpkg"), quiet = TRUE)
  poly <- poly[poly$patch_id %in% todo, "patch_id"] %>% st_make_valid() %>% st_transform(SINU)
  chunks <- split(seq_len(nrow(poly)), ceiling(seq_len(nrow(poly)) / CHUNK))
  for (k in seq_along(chunks)) {
    px <- bind_rows(px, pixel_fracs(poly[chunks[[k]], ]))
    message(sprintf("  %d of %d checked", min(k * CHUNK, nrow(poly)), nrow(poly)))   # progress inside the loop
  }
  write.csv(px, px_file, row.names = FALSE)
}

px_sum <- px %>%
  group_by(patch_id) %>%
  summarise(max_pixel_frac = max(loss_frac),
            n_core         = sum(loss_frac >= MIN_PIXEL_FRAC),
            core_loss_frac = mean(loss_frac[loss_frac >= MIN_PIXEL_FRAC]),
            .groups = "drop")

pool <- pool %>%
  left_join(px_sum, by = "patch_id") %>%
  mutate(max_pixel_frac = coalesce(max_pixel_frac, 0),
         n_core         = coalesce(n_core, 0L),
         modis_ok       = max_pixel_frac >= MIN_PIXEL_FRAC)

# Patches covering >= half of at least one MODIS pixel, by size
modis_by_size <- pool %>%
  mutate(size = cut(area_ha, c(0, 25, 50, 100, 250, Inf), right = FALSE,
                    labels = c("under 25 ha", "25-50 ha", "50-100 ha", "100-250 ha", "over 250 ha"))) %>%
  group_by(size) %>%
  summarise(patches = n(), pass = sum(modis_ok),
            pct_pass = round(100 * mean(modis_ok), 1), .groups = "drop") %>%
  as.data.frame()
modis_by_size

# ---- 3. The full set ---------------------------------------------------------
# Ordered by tile, then location, so Step 4's batches each cover one area
allp <- pool %>%
  filter(modis_ok) %>%
  mutate(region = paste(na_l2code, tools::toTitleCase(tolower(na_l2name)))) %>%
  arrange(tile, round(lat, 5), round(lon, 5)) %>%
  mutate(patch_n = row_number())
nrow(allp)   # patches in the full set

# Patches by Level II region
by_region <- allp %>%
  group_by(region) %>%
  summarise(patches = n(), median_ha = round(median(area_ha), 1),
            median_loss_year = median(loss_year), .groups = "drop") %>%
  as.data.frame()
by_region

write.csv(allp, file.path(outDir, "all_patches.csv"), row.names = FALSE)

poly_all <- st_read(file.path(outDir, "patches_qualifying.gpkg"), quiet = TRUE) %>%
  inner_join(allp %>% select(patch_id, patch_n, n_core, max_pixel_frac, core_loss_frac), by = "patch_id") %>%
  arrange(patch_n)
if (nrow(poly_all) != nrow(allp)) {
  stop("Found ", nrow(poly_all), " of ", nrow(allp), " patches in patches_qualifying.gpkg. ",
       "Re-run Patch Step 1, then this script.", call. = FALSE)
}
st_write(poly_all, file.path(outDir, "all_patches.gpkg"), delete_dsn = TRUE, quiet = TRUE)

# ---- 4. Map and year of loss -------------------------------------------------
eco    <- st_read(file.path(outDir, "boundary_ecoregions.gpkg"), quiet = TRUE)
all_sf <- st_as_sf(allp, coords = c("lon", "lat"), crs = 4326) %>% st_transform(5070)

ggplot(data=eco) +
  geom_sf(aes(fill=na_l2name), color="grey60", linewidth=0.1, alpha=0.3) + # study area, colored by Level II
  geom_sf(data=all_sf, color="darkgreen", size=0.1, alpha=0.4) + # one dot per patch
  theme_bw() +
  labs(title="Every forest loss patch in the analysis", fill="Level II ecoregion")

ggsave(file.path(outDir, "map_all_patches.png"), width=11, height=8, dpi=200)

ggplot(data=allp) +
  geom_bar(aes(x=factor(loss_year), fill=region)) +
  theme_bw() +
  labs(title="Year of loss, every patch, by region", x="Year of loss",
       y="Number of patches", fill="Level II ecoregion")

ggsave(file.path(outDir, "all_patches_by_year.png"), width=10, height=5, dpi=200)