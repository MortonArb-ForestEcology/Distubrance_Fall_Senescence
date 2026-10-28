# Patch Step 2: look at the size distribution of the qualifying patches
# -----------------------------------------------------------------------------
# Runs entirely in R (no Earth Engine). Open in RStudio and click "Source".
# Needs: dplyr, tidyr, ggplot2
# Reads:  scaleup_patches/patches_all_labeled.csv   (made by Patch Step 1)
#
# Use this to choose MIN_SAMPLE_HA for Step 3. The question to answer:
# how many qualifying patches are left at each size cutoff, and where do they
# sit relative to one MODIS pixel (~21.5 ha)?
# =============================================================================

library(dplyr)
library(ggplot2)

outDir      <- file.path("~/Google Drive/My Drive/Reidy_research", "scaleup_patches")
MODIS_PX_HA <- 21.5                          # one MODIS 500 m pixel
CUTOFFS     <- c(5, 10, 21.5, 25, 50, 100)   # candidate minimum sizes to compare

pt <- read.csv(file.path(outDir, "patches_all_labeled.csv"))
q  <- pt %>% filter(qualifies)
cat(sprintf("\nQualifying patches: %d (of %d found)\n", nrow(q), nrow(pt)))

# ---- 1. Basic numbers -------------------------------------------------------
cat("\nPatch size (ha), all qualifying patches:\n")
print(round(quantile(q$area_ha, c(0, .10, .25, .50, .75, .90, .99, 1)), 1))

# How many patches survive each candidate cutoff, overall and per Level II region
cut_tbl <- bind_rows(lapply(CUTOFFS, function(cmin) {
  q %>%
    filter(area_ha >= cmin) %>%
    group_by(na_l2name) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(min_ha = cmin)
})) %>%
  tidyr::pivot_wider(names_from = na_l2name, values_from = n, values_fill = 0) %>%
  mutate(total = rowSums(across(-min_ha)))
cat("\nPatches left at each minimum size:\n")
print(as.data.frame(cut_tbl))
write.csv(cut_tbl, file.path(outDir, "patch_counts_by_cutoff.csv"), row.names = FALSE)

# Patches in each size class. Each class includes its lower edge, so a patch
# of exactly 25 ha counts in "25-50 ha", not "10-25 ha".
SIZE_BREAKS <- c(5, 10, 25, 50, 100, 250, Inf)
SIZE_LABELS <- c("5-10 ha", "10-25 ha", "25-50 ha", "50-100 ha", "100-250 ha", "over 250 ha")

q <- q %>%
  mutate(size_class = cut(area_ha, breaks = SIZE_BREAKS, labels = SIZE_LABELS, right = FALSE))

size_tbl <- q %>%
  group_by(size_class) %>%
  summarise(patches = n(), total_area_ha = round(sum(area_ha)), .groups = "drop") %>%
  mutate(pct_of_patches = round(100 * patches / sum(patches), 1),
         pct_of_area    = round(100 * total_area_ha / sum(total_area_ha), 1))
cat("\nPatches by size class (all regions):\n")
print(as.data.frame(size_tbl), row.names = FALSE)

size_by_region <- q %>%
  count(na_l2name, size_class) %>%
  tidyr::pivot_wider(names_from = size_class, values_from = n, values_fill = 0)
cat("\nPatches by size class and Level II region:\n")
print(as.data.frame(size_by_region), row.names = FALSE)

write.csv(size_tbl,       file.path(outDir, "patch_size_classes.csv"),           row.names = FALSE)
write.csv(size_by_region, file.path(outDir, "patch_size_classes_by_region.csv"), row.names = FALSE)

cat("\nTypical patch size by Level II region:\n")
print(as.data.frame(q %>%
                      group_by(na_l2name) %>%
                      summarise(patches = n(),
                                median_ha = round(median(area_ha), 1),
                                pct_25ha_plus = round(100 * mean(area_ha >= 25), 1),
                                .groups = "drop")))

# ---- 2. Histogram of patch size (log scale) ---------------------------------
p1 <- ggplot(q, aes(x = area_ha)) +
  geom_histogram(bins = 60, fill = "darkgreen", color = "white", linewidth = 0.1) +
  geom_vline(xintercept = MODIS_PX_HA, linetype = "dashed") +
  annotate("text", x = MODIS_PX_HA, y = Inf, label = " 1 MODIS pixel (21.5 ha)",
           hjust = 0, vjust = 1.5, size = 3.5) +
  scale_x_log10(breaks = c(5, 10, 25, 50, 100, 250, 500, 1000)) +
  labs(title = "Size of qualifying forest loss patches",
       subtitle = sprintf("%d patches, lost 2004-2020, hardwood before, not converted", nrow(q)),
       x = "Patch size (ha, log scale)", y = "Number of patches") +
  theme_minimal()
print(p1)
ggsave(file.path(outDir, "patch_size_histogram.png"), p1, width = 9, height = 5, dpi = 200)

# ---- 3. "How many are at least this big?" curve -----------------------------
surv <- q %>%
  arrange(desc(area_ha)) %>%
  mutate(n_at_least = row_number())
p2 <- ggplot(surv, aes(x = area_ha, y = n_at_least)) +
  geom_step(color = "darkgreen") +
  geom_vline(xintercept = MODIS_PX_HA, linetype = "dashed") +
  geom_hline(yintercept = 100, linetype = "dotted", color = "red") +
  annotate("text", x = min(q$area_ha), y = 100, label = "100 = sample size",
           hjust = 0, vjust = -0.5, size = 3.5, color = "red") +
  scale_x_log10(breaks = c(5, 10, 25, 50, 100, 250, 500, 1000)) +
  scale_y_log10() +
  labs(title = "Number of patches at least a given size",
       x = "Minimum patch size (ha, log scale)", y = "Patches this big or bigger (log scale)") +
  theme_minimal()
print(p2)
ggsave(file.path(outDir, "patch_size_at_least.png"), p2, width = 9, height = 5, dpi = 200)

# ---- 4. Size by region and by loss year -------------------------------------
p3 <- ggplot(q, aes(x = na_l2name, y = area_ha)) +
  geom_boxplot(fill = "darkseagreen", outlier.size = 0.4) +
  geom_hline(yintercept = MODIS_PX_HA, linetype = "dashed") +
  scale_y_log10() +
  scale_x_discrete(labels = function(x) sapply(x, function(s)
    paste(strwrap(tools::toTitleCase(tolower(s)), 18), collapse = "\n"))) +
  labs(title = "Patch size by Level II ecoregion", x = NULL, y = "Patch size (ha, log scale)") +
  theme_minimal()
ggsave(file.path(outDir, "patch_size_by_region.png"), p3, width = 7, height = 5, dpi = 200)

p4 <- ggplot(q, aes(x = factor(loss_year), y = area_ha)) +
  geom_boxplot(fill = "darkseagreen", outlier.size = 0.4) +
  geom_hline(yintercept = MODIS_PX_HA, linetype = "dashed") +
  scale_y_log10() +
  labs(title = "Patch size by year of loss", x = "Year of loss", y = "Patch size (ha, log scale)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(file.path(outDir, "patch_size_by_year.png"), p4, width = 9, height = 5, dpi = 200)

cat("\nDone. Files in", outDir, ":\n",
    " patch_counts_by_cutoff.csv, patch_size_classes.csv, patch_size_classes_by_region.csv\n",
    " patch_size_histogram.png, patch_size_at_least.png\n",
    " patch_size_by_region.png, patch_size_by_year.png\n",
    "Next: choose MIN_SAMPLE_HA, then run patch_step3_random_sample.R\n")
