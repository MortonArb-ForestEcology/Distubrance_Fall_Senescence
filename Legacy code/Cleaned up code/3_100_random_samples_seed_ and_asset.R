# Patch Step 3: draw a repeatable random sample of 100 patches
# -----------------------------------------------------------------------------
# Open in RStudio and click "Source".
# Needs: rgee, sf, dplyr, ggplot2
# Reads:  scaleup_patches/patches_all_labeled.csv     (Patch Step 1)
#         scaleup_patches/patches_qualifying.gpkg     (Patch Step 1, polygons)
#         scaleup_patches/boundary_ecoregions.gpkg    (Patch Step 1, for the map)
#
# What it does:
#   1. Keeps qualifying patches at least MIN_SAMPLE_HA (choose it from Step 2)
#   2. Sorts them into a fixed order, then shuffles with SEED
#   3. The first N_SAMPLE in the shuffled order are the sample, the next ones
#      are numbered replacements (use them in order if a sampled patch has to
#      be dropped later, for example if there is bad MODIS data). The full order is saved.
#   4. Maps the sampled centroids and checks the sample looks like the pool
#   5. Saves the sampled patch polygons (with replacements) and uploads them
#      to a GEE asset for the NDVI and phenology steps
#
# Same input file + same settings = same sample every time.
# =============================================================================

xlibrary(rgee)
library(sf)
library(dplyr)
library(ggplot2)

ee_Initialize()

# ---- Settings ---------------------------------------------------------------
SEED           <- 42
N_SAMPLE       <- 100
N_REPLACEMENTS <- 20     # extra patches kept in order as replacements
MIN_SAMPLE_HA  <- 25     # choose from Step 2 (25 ha ~ one MODIS pixel and a bit)

outDir      <- file.path("~/Google Drive/My Drive/Reidy_research", "scaleup_patches")
sampleAsset <- "projects/breidyee/assets/scaleup_sample_patches"

# ---- 1. The pool ------------------------------------------------------------
pt   <- read.csv(file.path(outDir, "patches_all_labeled.csv"))
pool <- pt %>% filter(qualifies, area_ha >= MIN_SAMPLE_HA)

cat(sprintf("\nPool: %d qualifying patches >= %g ha\n", nrow(pool), MIN_SAMPLE_HA))
print(table(pool$na_l2name))
if (nrow(pool) < N_SAMPLE + N_REPLACEMENTS) {
  stop("Only ", nrow(pool), " patches in the pool; lower MIN_SAMPLE_HA or N_SAMPLE.", call. = FALSE)
}

# ---- 2 & 3. Fixed order -> seeded shuffle -----------------------------------
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

cat("\nSampled patches by Level II region:\n")
print(table(smp$na_l2name, smp$role))
cat("\nSampled patches by Level III region:\n")
print(sort(table(smp$us_l3name[smp$role == "Sample"]), decreasing = TRUE))

# ---- 4. Map + "does the sample look like the pool?" -------------------------
eco <- st_read(file.path(outDir, "boundary_ecoregions.gpkg"), quiet = TRUE)
pool_sf <- st_as_sf(pool, coords = c("lon", "lat"), crs = 4326) %>% st_transform(5070)

p_map <- ggplot() +
  geom_sf(data = eco, aes(fill = na_l2name), color = "grey60", linewidth = 0.1, alpha = 0.3) +
  geom_sf(data = filter(pool_sf, role == "Not drawn"), color = "grey55", size = 0.15, alpha = 0.4) +
  geom_sf(data = filter(pool_sf, role == "Replacement"), color = "orange", size = 1.3) +
  geom_sf(data = filter(pool_sf, role == "Sample"), color = "red", size = 1.6) +
  labs(title = sprintf("Random sample of %d patches (red), %d replacements (orange)",
                       N_SAMPLE, N_REPLACEMENTS),
       subtitle = sprintf("Grey = all %d qualifying patches >= %g ha. Seed = %d.",
                          nrow(pool), MIN_SAMPLE_HA, SEED),
       fill = "Level II ecoregion") +
  theme_minimal()
print(p_map)
ggsave(file.path(outDir, "map_sample_patches.png"), p_map, width = 11, height = 8, dpi = 200)

check <- bind_rows(
  pool %>% mutate(set = "All qualifying"),
  pool %>% filter(role == "Sample") %>% mutate(set = "Sample")
)
p_chk <- ggplot(check, aes(x = set, y = area_ha, fill = set)) +
  geom_boxplot(outlier.size = 0.4) +
  scale_y_log10() +
  scale_fill_manual(values = c("All qualifying" = "grey70", "Sample" = "red"), guide = "none") +
  labs(title = "Patch size: sample vs. everything it was drawn from",
       subtitle = "The boxes should look alike.", x = NULL, y = "Patch size (ha, log scale)") +
  theme_minimal()
ggsave(file.path(outDir, "sample_vs_pool_size.png"), p_chk, width = 6, height = 5, dpi = 200)

p_yr <- ggplot(check, aes(x = loss_year, fill = set)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, position = "identity", alpha = 0.5) +
  scale_fill_manual(values = c("All qualifying" = "grey50", "Sample" = "red"), name = NULL) +
  labs(title = "Year of loss: sample vs. everything it was drawn from",
       x = "Year of loss", y = "Share of patches") +
  theme_minimal()
ggsave(file.path(outDir, "sample_vs_pool_year.png"), p_yr, width = 8, height = 5, dpi = 200)

# ---- 5. Upload the sampled patch polygons to their own GEE asset ----------
poly <- st_read(file.path(outDir, "patches_qualifying.gpkg"), quiet = TRUE)
smp_poly <- poly %>%
  inner_join(smp %>% select(patch_id, role, draw_order), by = "patch_id") %>%
  select(patch_id, role, draw_order, loss_year, area_ha, na_l2code, us_l3code,
         decid_before, after_label, lon, lat) %>%
  arrange(draw_order)
if (nrow(smp_poly) != nrow(smp)) {
  stop("Found ", nrow(smp_poly), " of ", nrow(smp), " sampled patches in patches_qualifying.gpkg. ",
       "Was Patch Step 1 re-run after this sample was drawn?", call. = FALSE)
}
st_write(smp_poly, file.path(outDir, "sample_patches.gpkg"), delete_dsn = TRUE, quiet = TRUE)

sf_as_ee(st_transform(smp_poly, 4326), via = "getInfo_to_asset",
         assetId = sampleAsset, overwrite = TRUE)

cat("\nDone. Files in", outDir, ":\n",
    " sample_patches.csv / .gpkg (the", N_SAMPLE, "patches + replacements)\n",
    " sample_pool_ordered.csv   (every patch in the pool, in random order)\n",
    " map_sample_patches.png, sample_vs_pool_size.png, sample_vs_pool_year.png\n",
    "Earth Engine asset:", sampleAsset, "\n")
