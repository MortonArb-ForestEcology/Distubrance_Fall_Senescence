# Patch Step 3: draw a repeatable random sample of 100 patches
# -----------------------------------------------------------------------------
#  Run line by line or use
# Source with Echo so the tables and plots show. 
# Needs: sf, dplyr, ggplot2
# Reads:  scaleup_patches/patches_all_labeled.csv     (Patch Step 1)
#         scaleup_patches/ecoregion_list.csv          (Patch Step 0, which ecoregions are in)
#         scaleup_patches/patches_qualifying.gpkg     (Patch Step 1, polygons)
#         scaleup_patches/boundary_ecoregions.gpkg    (Patch Step 0, for the map)
# Saves:  sample_patches.csv / .gpkg     the 100 patches + 100 replacements
#         sample_pool_ordered.csv        every eligible patch, in random order
#         patch_modis_pixel_cover.csv    MODIS-pixel check for each patch
#         map_sample_patches.png, sample_vs_pool_size.png, sample_vs_pool_year.png
#
# What it does:
#   1. Starts from the qualifying patches at least MIN_SAMPLE_HA
#   2. MODIS RULE (applied BEFORE the draw): keeps only patches that cover at
#      least half of at least one MODIS pixel. Works this out exactly from each
#      patch's outline laid over the MODIS grid. Patches that can't meet it are
#      never in the pool, so nothing has to be dropped or flagged afterwards.
#   3. Sorts the eligible patches into a fixed order, then shuffles with SEED.
#      The first N_SAMPLE are the sample; the next N_REPLACEMENTS are numbered
#      replacements. Step 4 uses replacements, in order, for any sampled patch
#      that has too little uncut forest nearby to compare against.
#   4. Maps the sample and checks it looks like the pool
#   5. Saves the sampled patch polygons for Step 4
#
# The MODIS-pixel check is saved, so re-running (e.g. with a new seed) is fast.
# Same input file + same settings = same sample every time.
# =============================================================================

library(sf)
library(dplyr)
library(ggplot2)

# ---- Settings ---------------------------------------------------------------
SEED           <- 42
N_SAMPLE       <- 100
N_REPLACEMENTS <- 100    # extra patches kept in order as replacements (only the first ones needed are used)
MIN_SAMPLE_HA  <- 25     # size floor (a patch must be >= 10.7 ha to pass the MODIS rule anyway)
MIN_PIXEL_FRAC <- 0.5    # patch must cover at least this share of at least one MODIS pixel

outDir     <- file.path("~/Google Drive/My Drive/Reidy_research", "scaleup_patches")
cover_file <- file.path(outDir, "patch_modis_pixel_cover.csv")

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

# For each patch: the largest share of any one MODIS pixel it covers.
pixel_cover <- function(poly_sin) {
  bb <- t(sapply(st_geometry(poly_sin), st_bbox))   # xmin, ymin, xmax, ymax
  cr <- bind_rows(lapply(seq_len(nrow(bb)), function(i) {
    expand.grid(col = floor((bb[i, 1] - X0) / PX):floor((bb[i, 3] - X0) / PX),
                row = floor((Y0 - bb[i, 4]) / PX):floor((Y0 - bb[i, 2]) / PX))
  })) %>% distinct()
  cells <- st_sf(cr, geometry = st_sfc(mapply(cell_square, cr$col, cr$row, SIMPLIFY = FALSE), crs = SINU))
  ov <- suppressWarnings(st_intersection(poly_sin[, "patch_id"], cells))
  ov$frac <- as.numeric(st_area(ov)) / PX^2
  st_drop_geometry(ov) %>%
    group_by(patch_id) %>%
    summarise(max_pixel_frac = max(frac), .groups = "drop")
}

# ---- 1. Qualifying patches above the size floor -----------------------------
pt   <- read.csv(file.path(outDir, "patches_all_labeled.csv"))
pool <- pt %>% filter(qualifies, area_ha >= MIN_SAMPLE_HA)

# Follow Step 0's CURRENT list of included ecoregions. Leaving ecoregions out
# in Step 0 still works without re-running
# Step 1.If ecoregions are added in Step 0 rerun Step 1. 
eco_list  <- read.csv(file.path(outDir, "ecoregion_list.csv"))
inc_codes <- as.character(eco_list$us_l3code[as.logical(eco_list$included)])
left_out  <- pool %>% filter(!(as.character(us_l3code) %in% inc_codes))
table(paste(left_out$us_l3code, left_out$us_l3name))   # patches left out: their ecoregion is no longer in the study area

pool <- pool %>% filter(as.character(us_l3code) %in% inc_codes)
nrow(pool)   # qualifying patches >= MIN_SAMPLE_HA in the study area

# ---- 2. MODIS rule: covers >= half of at least one MODIS pixel --------------
# A patch smaller than half a pixel (~10.7 ha) can't pass, so skip those.
need <- pool$patch_id[pool$area_ha >= MIN_PIXEL_FRAC * PX^2 / 1e4]

cover <- if (file.exists(cover_file)) read.csv(cover_file) else
  data.frame(patch_id = character(), max_pixel_frac = numeric())
todo <- setdiff(need, cover$patch_id)
length(todo)   # patches still to check against the MODIS grid (0 = all saved from before)

if (length(todo) > 0) {
  poly <- st_read(file.path(outDir, "patches_qualifying.gpkg"), quiet = TRUE)
  poly <- poly[poly$patch_id %in% todo, "patch_id"] %>% st_make_valid() %>% st_transform(SINU)
  chunks <- split(seq_len(nrow(poly)), ceiling(seq_len(nrow(poly)) / 2000))
  for (k in seq_along(chunks)) {
    cover <- bind_rows(cover, pixel_cover(poly[chunks[[k]], ]))
    message(sprintf("  %d of %d done", min(k * 2000, nrow(poly)), nrow(poly)))   # progress inside the loop
  }
  cover <- cover %>% distinct(patch_id, .keep_all = TRUE)
  write.csv(cover, cover_file, row.names = FALSE)
}

pool <- pool %>%
  left_join(cover %>% select(patch_id, max_pixel_frac), by = "patch_id") %>%
  mutate(max_pixel_frac = coalesce(max_pixel_frac, 0),
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

pool <- pool %>% filter(modis_ok)
nrow(pool)               # the pool: patches eligible for the sample
table(pool$na_l2name)    # the pool by Level II region
if (nrow(pool) < N_SAMPLE + N_REPLACEMENTS) {
  stop("Only ", nrow(pool), " eligible patches; lower MIN_SAMPLE_HA or N_SAMPLE.", call. = FALSE)
}

# ---- 3. Fixed order -> seeded shuffle ---------------------------------------
# Sort on things that don't change between runs (year of loss, then location),
# so the same seed always gives the same draw.
set.seed(SEED)
pool <- pool %>%
  arrange(loss_year, round(lon, 5), round(lat, 5)) %>%
  mutate(draw_order = sample(n())) %>%
  arrange(draw_order) %>%
  mutate(role = case_when(
    draw_order <= N_SAMPLE                  ~ "Sample",
    draw_order <= N_SAMPLE + N_REPLACEMENTS ~ "Replacement",
    TRUE                                    ~ "Not drawn"
  ))

write.csv(pool, file.path(outDir, "sample_pool_ordered.csv"), row.names = FALSE)
smp <- pool %>% filter(role != "Not drawn")
write.csv(smp, file.path(outDir, "sample_patches.csv"), row.names = FALSE)

table(smp$na_l2name, smp$role)   # sampled patches by Level II region
sort(table(smp$us_l3name[smp$role == "Sample"]), decreasing = TRUE)   # sampled patches by Level III region

# ---- 4. Map + "does the sample look like the pool?" -------------------------
eco <- st_read(file.path(outDir, "boundary_ecoregions.gpkg"), quiet = TRUE)
pool_sf <- st_as_sf(pool, coords = c("lon", "lat"), crs = 4326) %>% st_transform(5070)

ggplot(data=eco) +
  geom_sf(aes(fill=na_l2name), color="grey60", linewidth=0.1, alpha=0.3) + # study area, colored by Level II
  geom_sf(data=filter(pool_sf, role=="Not drawn"), aes(color=role), size=0.15, alpha=0.4) + # every other eligible patch
  geom_sf(data=filter(pool_sf, role=="Replacement"), aes(color=role), size=1.3) +
  geom_sf(data=filter(pool_sf, role=="Sample"), aes(color=role), size=1.6) +
  scale_color_manual(name="Patch", values=c("Sample"="red", "Replacement"="orange", "Not drawn"="grey55")) +
  theme_bw() +
  labs(title="Random sample of forest loss patches", fill="Level II ecoregion")

ggsave(file.path(outDir, "map_sample_patches.png"), width=11, height=8, dpi=200)

check <- bind_rows(
  pool %>% mutate(set = "All eligible"),
  pool %>% filter(role == "Sample") %>% mutate(set = "Sample")
)

ggplot(data=check) +
  geom_boxplot(aes(x=set, y=area_ha, fill=set), outlier.size=0.4) +
  scale_y_log10() +
  scale_fill_manual(values=c("All eligible"="grey70", "Sample"="red"), guide="none") +
  theme_bw() +
  labs(title="Patch size: sample vs. everything it was drawn from",
       subtitle="The boxes should look alike.", x=NULL, y="Patch size (ha, log scale)")

ggsave(file.path(outDir, "sample_vs_pool_size.png"), width=6, height=5, dpi=200)

ggplot(data=check) +
  geom_histogram(aes(x=loss_year, y=after_stat(density), fill=set), binwidth=1, position="identity", alpha=0.5) +
  scale_fill_manual(name="", values=c("All eligible"="grey50", "Sample"="red")) +
  theme_bw() +
  labs(title="Year of loss: sample vs. everything it was drawn from", x="Year of loss", y="Share of patches")

ggsave(file.path(outDir, "sample_vs_pool_year.png"), width=8, height=5, dpi=200)

# ---- 5. Save the sampled patch polygons for Step 4 --------------------------
poly_all <- st_read(file.path(outDir, "patches_qualifying.gpkg"), quiet = TRUE)
smp_poly <- poly_all %>%
  inner_join(smp %>% select(patch_id, role, draw_order, max_pixel_frac), by = "patch_id") %>%
  select(patch_id, role, draw_order, loss_year, area_ha, max_pixel_frac, na_l2code, us_l3code,
         decid_before, after_label, lon, lat) %>%
  arrange(draw_order)
if (nrow(smp_poly) != nrow(smp)) {
  stop("Found ", nrow(smp_poly), " of ", nrow(smp), " sampled patches in patches_qualifying.gpkg. ",
       "Was Patch Step 1 re-run after this sample was drawn?", call. = FALSE)
}
st_write(smp_poly, file.path(outDir, "sample_patches.gpkg"), delete_dsn = TRUE, quiet = TRUE)