
# Spatial map of the Dominican Republic study areas ----

----

# 1. load in the data ---- 

library(sf)
library(ggplot2)

df_0 <- read_excel("Data/DR_surveillance/EN/Copy of DR_serosurvey_MASTER_with_MBA_UPDATED.xlsx")
df <- df_0

dr0 <- st_read("Data/Spatial/gadm41_DOM_shp/gadm41_DOM_0.shp")
dr1 <- st_read("Data/Spatial/gadm41_DOM_shp/gadm41_DOM_1.shp")
dr2 <- st_read("Data/Spatial/gadm41_DOM_shp/gadm41_DOM_2.shp")

ggplot(dr0) + geom_sf() + theme_minimal()
ggplot(dr1) + geom_sf() + theme_minimal()
ggplot(dr2) + geom_sf() + theme_minimal()

----
# 2. check the data ---- 

df$gps_latitude[df$gps_latitude == "NA"] <- NA
df$gps_longitude[df$gps_longitude == "NA"] <- NA

df$gps_latitude  <- as.numeric(df$gps_latitude)
df$gps_longitude <- as.numeric(df$gps_longitude)

df_clean <- df[!is.na(df$gps_latitude) & !is.na(df$gps_longitude), ]

points_sf <- st_as_sf(
  df_clean,
  coords = c("gps_longitude", "gps_latitude"), 
  crs = 4326   
)

ggplot() +
  geom_sf(data = dr2, fill = "grey95") +
  geom_sf(data = points_sf, color = "red", size = 1) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

