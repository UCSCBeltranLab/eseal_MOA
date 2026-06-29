
# Setup: load libraries and data -----------------------------------------------

library(here)
library(tidyverse)
library(lubridate)
library(sf)
library(viridis)

##suppress scientific notation
options(scipen = 999)

##turn off spherical geometry
sf_use_s2(FALSE)

# SEAL DENSITY -----------------------------

#### 1) Read spatial data ####

# read merged Picterra seal detection polygons
picterra.output <- st_read(here("IntermediateData", "merged_polygons.gpkg"), quiet = TRUE)

# read Año Nuevo beach polygons
beaches <- st_read(here("IntermediateData", "beaches.gpkg"), quiet = TRUE)

#### 2) Convert to NAD83 / California zone 3 ####

# check original coordinate reference system
st_crs(picterra.output)

# transform seal detections to NAD83 / California zone 3
picterra.output <- st_transform(picterra.output, "EPSG:26943")

# confirm updated CRS
st_crs(picterra.output)

# transform beach polygons to the same CRS
beaches <- st_transform(beaches, "EPSG:26943")

#### 3) Assign age/sex class to seals ####

# standardize age/sex labels and keep only females, males, and pups
uas.data <- picterra.output %>%
  filter(age_sex %in% c("female", "male", "pup"))

#### 4) Assign seals to beaches ####

# convert seal polygons to centroid points
centroids <- st_centroid(uas.data)

# spatially assign each seal centroid to a beach polygon
seals_by_beach <- st_intersection(beaches, centroids)

#### 5) Calculate female density ####

# create vector of survey dates from the seal data
survey_dates <- unique(seals_by_beach$date)

# initialize empty list to store density estimates
density.list <- list()

# loop through all survey dates
for (i in seq_along(survey_dates)) {
  
  # subset female seals for the current survey date
  survey.subset <- seals_by_beach %>%
    filter(date == survey_dates[i],
           age_sex == "female")
  
  # skip dates with no female seals
  if (nrow(survey.subset) == 0) next
  
  # calculate centroids for each female seal
  seal.centroids <- st_centroid(survey.subset)
  
  # create 10 m buffers around each centroid
  seal.buffer <- st_buffer(seal.centroids, 10)
  
  # identify neighboring seals within each buffer
  int <- st_intersects(seal.buffer, seal.centroids)
  
  # calculate 10 m density as the number of neighboring seals
  survey.subset$density <- lengths(int) - 1
  
  # add survey results to output list
  density.list[[i]] <- survey.subset
}

# combine all survey results into one sf object
density.df <- do.call(rbind, density.list)

#### 6) Save density dataframe ####

# calculate centroid coordinates
density.centroids <- st_centroid(density.df)

# extract coordinates
coords <- st_coordinates(density.centroids)

# create final density dataframe
seal_density <- density.centroids %>%
  mutate(lon = coords[,1],
         lat = coords[,2]) %>%
  mutate(date = as.integer(substr(date,1, 4))) %>%
  rename(season = date) %>%
  select(season, age_sex, lat, lon, Beach, density) %>%
  st_drop_geometry()

# save density dataframe as CSV
write.csv(seal_density, here("IntermediateData", "seal_density.csv"), row.names = FALSE, na = "NA")

# Plot two density radius examples --------------------

#### 1. low density example ####

# grab light color from mako palette
light_mako <- viridis(190, option = "mako", direction = -1)[20]

# choose one low-density female from density dataframe
low_point <- density.df %>%
  filter(age_sex == "female",
         density == 4) %>%
  slice(1)

# get original seal polygons from same date
same <- uas.data %>%
  filter(date == low_point$date,
         age_sex == "female")

# find original polygon closest to selected density point
focal_id <- st_nearest_feature(low_point, same)

# pull original focal seal polygon
low_example <- same[focal_id, ]

# create 10 m buffer around focal seal centroid
seal.buffer <- st_buffer(st_centroid(low_example), 10)

# identify neighboring females using centroids
hits <- st_intersects(seal.buffer, st_centroid(same))[[1]]

# keep neighboring seal polygons, excluding focal seal
int.poly <- same[setdiff(hits, focal_id), ]

# save figure
png(here("TablesFigures", "density_example_low.png"),
    width = 2000,
    height = 2000,
    res = 300,
    bg = "transparent")

#show plot
par(bg = NA, mar = c(0, 0, 0, 0))

plot(st_geometry(seal.buffer),
     border = "grey40")

plot(st_geometry(int.poly),
     add = TRUE,
     col = light_mako,
     border = "grey40")

plot(st_geometry(low_example),
     add = TRUE,
     col = light_mako,
     border = "black",
     lwd = 3)

dev.off()

#### 2. high density example ####

# grab dark color from mako palette
dark_mako <- viridis(50, option = "mako", direction = 1)[20]

# choose one high-density female from density dataframe
high_point <- density.df %>%
  filter(age_sex == "female",
         density == 20) %>%
  slice(1)

# get original seal polygons from same date
same <- uas.data %>%
  filter(date == high_point$date,
         age_sex == "female")

# find original polygon closest to selected density point
focal_id <- st_nearest_feature(high_point, same)

# pull original focal seal polygon
high_example <- same[focal_id, ]

# create 10 m buffer around focal seal centroid
seal.buffer <- st_buffer(st_centroid(high_example), 10)

# identify neighboring females using centroids
hits <- st_intersects(seal.buffer, st_centroid(same))[[1]]

# keep neighboring seal polygons, excluding focal seal
int.poly <- same[setdiff(hits, focal_id), ]

# save figure
png(here("TablesFigures", "density_example_high.png"),
    width = 2000,
    height = 2000,
    res = 300,
    bg = "transparent")

#show plot
par(bg = NA, mar = c(0, 0, 0, 0))

plot(st_geometry(seal.buffer),
     border = "grey40")

plot(st_geometry(int.poly),
     add = TRUE,
     col = dark_mako,
     border = "grey40")

plot(st_geometry(high_example),
     add = TRUE,
     col = dark_mako,
     border = "black",
     lwd = 5)

dev.off()


