library(sf)
library(spdep)
library(tidyverse)
library(tmap)

data <- read_csv("D:/Rabies_Analysis/rbs/points_with_wards.csv") %>%
  group_by(ADM2_EN) %>%
  summarise(Incidence = mean(Incidence, na.rm = TRUE)) %>%
  ungroup()

# View(data)
# str(data)

ADM2_EN_sf <- st_read("D:/DATA CENTER/ken_adm_iebc_20191031_shp/ken_admbnda_adm2_iebc_20191031.shp") %>%
  st_transform(4326) %>%
  select(ADM2_EN)


merged_sf <- ADM2_EN_sf %>%
  # Left-join 
  left_join(data, by = "ADM2_EN") %>%     
  #  Replace those NAs with zeros
  mutate(Incidence = replace_na(Incidence, 0))

invalid_idx <- which(!st_is_valid(merged_sf))
if (length(invalid_idx) > 0) {
  message("Found ", length(invalid_idx), " invalid geometries; attempting repair…")
  merged_sf[invalid_idx, ] <- st_make_valid(merged_sf[invalid_idx, ])
  # Double‐check
  stopifnot(all(st_is_valid(merged_sf)))
}

#  Creating spatial weights matrix (row‐standardized)
nb <- poly2nb(merged_sf, queen = TRUE)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

#  Global Moran's I
moran_results <- moran.test(merged_sf$Incidence, 
                            listw = lw, 
                            randomisation = TRUE)

print(moran_results)

#  Local Moran’s I
local_mi <- localmoran(merged_sf$Incidence, lw, zero.policy = TRUE)
local_mi <- bind_cols(data.frame(local_mi,check.names = F), attr(local_mi, "quadr"))
colnames(local_mi)
# View(local_mi)

# Attach local Moran’s I statistics back to the sf object
merged_sf <- merged_sf %>%
  mutate(
    Local_I       = local_mi[, "Ii"],    
    Expectation   = local_mi[, "E.Ii"],   
    Variance      = local_mi[, "Var.Ii"], 
    Z_score       = local_mi[, "Z.Ii"],  
    P_value       = local_mi[, "Pr(z != E(Ii))"],
    quadr         = local_mi[, "mean"]
  )


# View(merged_sf)

# save shp
# st_write(
#   merged_sf,
#   dsn         = "D:/Rabies_Analysis/rbs/ADM2_EN_Mean_LISA.shp",
#   delete_layer = TRUE
# )

# Save csv
# write_csv(
#   merged_sf %>%
#     st_drop_geometry() %>%
#     select(ADM2_EN, Incidence, Z_score, Local_I, P_value, quadr),
#   "D:/Rabies_Analysis/rbs/local_moran_results_ADM2_EN_Mean.csv"
# )
