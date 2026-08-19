library(ggplot2); library(nlme)
library(tidyverse)

outputDir <- "~/Google Drive/My Drive/Reidy_research/"
phenology_long <- read.csv(file.path(outputDir, 'phenology_yearly_by_patch.csv'))
summary(phenology_long)

patch_attrs <- read.csv(file.path(outputDir, 'patch_attribute_table.csv'))
summary(patch(attrs))

ndvi_long <- ndvi_long %>%
  select(-any_of('area_ha')) %>%
  left_join(patch_attrs %>% select(patch_uuid, area_ha), by = 'patch_uuid') %>%
  filter(area_ha >= 25)

# Drop rows with no valid NDVI (ndvi_july_count == 0 -> ndvi_july_mean is NA,
# per the fix in Step 4 that guarantees the column exists even when blank).
ndvi_long <- ndvi_long %>% filter(!is.na(ndvi_july_mean))

# Relativize to each patch's own loss_year HERE — this is the deliberate
# downstream step Step 4 was designed to leave open, rather than baking a
# fixed pre/post window into the extraction itself.
ndvi_long <- ndvi_long %>%
  mutate(years_since_loss = year - loss_year)
