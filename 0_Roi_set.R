# Patch Step 0: choose the study area and map it
# -----------------------------------------------------------------------------
# Open in RStudio and run line by line  or with e Source
# pacakges: rgee, sf, dplyr, ggplot2; maps
#
# What it does:
#   1. Downloads every EPA Level III ecoregion in the lower 48 states
#   2. Finds the westernmost point of the ORIGINAL study regions (Level II 8.4,
#      5.3, 5.2, 8.1) and draws a north-south line there
#   3. Adds every other Level III ecoregion with >= 50% of its
#      area east of that line.
#
#   4. Shows the list and saves it as map: included ecoregions colored by Level
#      II region, original ones outlined in black, the west line dashed red.
#
# Look at the map BEFORE running Step 1. To drop a whole Level II region, put
# its code in ALSO_EXCLUDE_L2. To add or drop a single Level III ecoregion,
# put its Level III code in ALSO_INCLUDE_L3 / ALSO_EXCLUDE_L3. Then run again.
# =============================================================================

library(rgee)
library(sf)
library(dplyr)
library(ggplot2)

ee_Initialize(drive = TRUE)

# ---- Settings ---------------------------------------------------------------
ORIGINAL_L2     <- c("8.4", "5.3", "5.2", "8.1")
BOUNDARY_RULE   <- "east_of_line" # or "level1"
MIN_SHARE_EAST  <- 0.5            # east_of_line: share of an ecoregion east of the line
LEVEL1_CODES    <- c("5", "8")    # level1: Northern Forests + Eastern Temperate Forests
ALSO_EXCLUDE_L2 <- c("9.5", "9.4", "15.4")   # Level II codes, e.g. "9.4" drops all of it
ALSO_INCLUDE_L3 <- character(0)   # Level III codes, e.g. c("48") to force one in
ALSO_EXCLUDE_L3 <- character(0)   # Level III codes, e.g. c("76") to force one out

outDir <- file.path("~/Google Drive/My Drive/Reidy_research", "scaleup_patches")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. All Level III ecoregions --------------------------------------------
ecoFC <- ee$FeatureCollection("EPA/Ecoregions/2013/L3")$
  map(ee_utils_pyfunc(function(f) f$simplify(maxError = 250)))
eco <- ee_as_sf(ecoFC, via = "drive", maxFeatures = 10000) %>%
  st_transform(5070) %>%
  st_make_valid() %>%
  group_by(us_l3code, us_l3name, na_l2code, na_l2name, na_l1code, na_l1name) %>%
  summarise(.groups = "drop") %>%
  mutate(us_l3code = as.character(us_l3code),
         area_km2  = round(as.numeric(st_area(.)) / 1e6))

# ---- 2. The west line -------------------------------------------------------
original <- eco$na_l2code %in% ORIGINAL_L2
west_lon <- unname(st_bbox(st_transform(st_union(eco[original, ]), 4326))["xmin"])
west_lon   # westernmost point of the original regions (degrees longitude)

# Everything east of the line, as a polygon (dense west edge so it stays a
# true line of longitude after projecting)
lats <- seq(20, 55, by = 0.1)
east_ll <- st_polygon(list(rbind(cbind(west_lon, lats),
                                 cbind(-60, rev(lats)),
                                 c(west_lon, 20))))
east <- st_transform(st_sfc(east_ll, crs = 4326), 5070)

ix <- suppressWarnings(st_intersection(eco[, "us_l3code"], east))
ix$east_km2 <- as.numeric(st_area(ix)) / 1e6
eco <- eco %>%
  left_join(st_drop_geometry(ix) %>% group_by(us_l3code) %>%
              summarise(east_km2 = sum(east_km2), .groups = "drop"),
            by = "us_l3code") %>%
  mutate(share_east = coalesce(east_km2 / area_km2, 0)) %>%
  select(-east_km2)

# ---- 3. Which ecoregions are in ---------------------------------------------
eco <- eco %>%
  mutate(share_east = round(pmin(share_east, 1), 3),
         na_l1code  = as.character(na_l1code),
         by_line    = share_east >= MIN_SHARE_EAST,
         by_level1  = na_l1code %in% LEVEL1_CODES,
         rule_ok    = if (BOUNDARY_RULE == "level1") by_level1 else by_line,
         status = case_when(
           us_l3code %in% ALSO_EXCLUDE_L3            ~ "Not included",
           na_l2code %in% ALSO_EXCLUDE_L2            ~ "Not included",
           na_l2code %in% ORIGINAL_L2                ~ "Original",
           us_l3code %in% ALSO_INCLUDE_L3            ~ "Added",
           rule_ok                                   ~ "Added",
           TRUE                                      ~ "Not included"),
         included = status != "Not included",
         l2_label = paste(na_l2code, tools::toTitleCase(tolower(na_l2name))))

inc <- eco %>% filter(included)

# Included ecoregions
inc_table <- st_drop_geometry(inc) %>%
  arrange(na_l2code, us_l3code) %>%
  select(status, na_l2code, na_l2name, us_l3code, us_l3name, area_km2, share_east) %>%
  as.data.frame()
inc_table

# Ecoregions that straddle the line (share east of it; >= 0.5 is included)
straddle <- eco %>%
  filter(share_east > 0.05, share_east < 0.95, !(na_l2code %in% ORIGINAL_L2)) %>%
  st_drop_geometry() %>%
  select(us_l3code, us_l3name, na_l2code, share_east, status) %>%
  as.data.frame()
straddle

# Ecoregions where the two rules disagree (east of line vs Level I 5 + 8)
differ <- eco %>%
  filter(by_line != by_level1, !(na_l2code %in% ORIGINAL_L2)) %>%
  st_drop_geometry() %>%
  select(us_l3code, us_l3name, na_l1name, share_east, by_line, by_level1, included) %>%
  as.data.frame()
differ

st_write(eco, file.path(outDir, "ecoregions_conus_l3.gpkg"), delete_dsn = TRUE, quiet = TRUE)
st_write(inc, file.path(outDir, "boundary_ecoregions.gpkg"), delete_dsn = TRUE, quiet = TRUE)
write.csv(st_drop_geometry(eco) %>% arrange(desc(included), na_l2code, us_l3code),
          file.path(outDir, "ecoregion_list.csv"), row.names = FALSE)

# ---- 4. Map -----------------------------------------------------------------
excluded  <- eco %>% filter(!included)
original_outline <- st_union(filter(inc, status == "Original"))
states    <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) %>% st_transform(5070)
west_line <- st_transform(st_sfc(st_linestring(cbind(west_lon, seq(24, 50, by = 0.1))), crs = 4326), 5070)

ggplot(data=inc) +
  geom_sf(data=excluded, fill="grey92", color="grey70", linewidth=0.15) + # ecoregions left out
  geom_sf(aes(fill=l2_label), color="grey35", linewidth=0.15) + # included ecoregions, colored by Level II
  geom_sf(data=states, fill=NA, color="white", linewidth=0.25) + # state lines
  geom_sf(data=original_outline, fill=NA, color="black", linewidth=0.7) + # outline of the 4 original Level II regions
  geom_sf(data=west_line, color="red", linetype="dashed", linewidth=0.6) + # the west line
  geom_sf_text(aes(label=us_l3code), size=2.2) + # Level III code on each ecoregion
  scale_fill_brewer(name="Level II ecoregion", palette="Set3") +
  theme_bw() +
  labs(title="Study area: EPA Level III ecoregions",
       subtitle="Colored = included. Black outline = original regions. Red dashed line = their westernmost point.")

ggsave(file.path(outDir, "map_study_ecoregions.png"), width=12, height=8.5, dpi=200)

