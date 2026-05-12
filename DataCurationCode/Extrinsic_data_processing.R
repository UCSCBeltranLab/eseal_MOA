
##In this file, I (1) Clean up the wave and tide files imported from NDBC and NOAA for the data processing code
##and (2) Calculate conspecific density using the drone shapefiles and beach polygons

# Setup: load libraries and data -----------------------------------------------

library(here)
library(tidyverse)
library(lubridate)

##suppress scientific notation
options(scipen = 999)

##turn off spherical geometry
sf_use_s2(FALSE)

# WAVE AND TIDE ---------------------------

#all data for tides are available to upload at: https://tidesandcurrents.noaa.gov/waterlevels.html?id=9413450
#all data for wave metrics are available to upload at: https://www.ndbc.noaa.gov/station_history.php?station=46042

##Identify seasonal date range: 
##Calculate earliest and latest dates are in the breeding season
##Used for: selecting dates when uploading the wave and tide data
breeding_season_dates <- metadata %>% 
  group_by(season) %>%
  summarize(earliest_date = min(date, na.rm = TRUE), 
            latest_date = max(date, na.rm = TRUE))

#### 1) Read in and clean tide data ####

#list only CSV files containing "tides" in the filename
tide_files <- list.files(path = here("RawData"),
                         pattern = "tides.*\\.csv$", 
                         full.names = TRUE)


# read and combine all tide CSVs
tide_data_all <- tide_files %>%
  map_dfr(read_csv)

#clean and create date variables
tide_data_all <- tide_data_all %>%
  rename(tide_height = "Verified (ft)") %>%
  rename(time_GMT = "Time (GMT)") %>%
  mutate(Date = as.Date(as.character(Date)),
         MM = month(Date),
         YY = year(Date),
         season = case_when(MM == 12 ~ YY + 1, #December observations belong to the *following* breeding season
                            MM %in% c(1, 2, 3) ~ YY, #Jan–Mar belong to the current year’s season
                            TRUE ~ NA_real_))

#select necessary columns from the tide data 
#Date: 
#Time: 
#Verified: 
#season: 
tide_data_clean <- tide_data_all %>%
  select("Date", "time_GMT", "tide_height", "season") %>%
  mutate(tide_datetime = as.POSIXct(paste(Date, time_GMT),
                                    format = "%Y-%m-%d %H:%M",
                                    tz = "UTC"))

#### 2) Read in and clean wave data ####

#list only CSV files containing "waves" in the filename
wave_files <- list.files(path = here("RawData"),
                         pattern = "waves.*\\.csv$", 
                         full.names = TRUE)

#read and combine only the wave data files
wave_data_all <- wave_files %>%
  map_dfr(~ read_csv(.x, col_types = cols(.default = "c")))

#convert all wave columns to numeric
wave_data_all <- wave_data_all %>%
  mutate(YY = as.numeric(trimws(YY)), #trimws to make sure there are no weird blank spaces and make each column numeric
         MM = as.numeric(trimws(MM)),
         DD = as.numeric(trimws(DD)),
         hh = as.numeric(trimws(hh)),
         WVHT = as.numeric(trimws(WVHT)),
         DPD = as.numeric(trimws(DPD))) %>%
  filter(WVHT != 99, DPD != 99) #99 is used as an NA placeholder in the data, so we can filter these out

#fit wave data to specific breeding seasons
wave_data_all <- wave_data_all %>%
  mutate(season = case_when(MM == 12 ~ YY + 1, #December observations belong to the *following* breeding season
                            MM %in% c(1, 2, 3) ~ YY, #Jan–Mar belong to the current year’s season
                            TRUE ~ NA_real_))

#select necessary columns from the wave data
#YY: 
#MM: 
#DD: 
#hh:
#WVHT: 
#DPD:
#season: 
wave_data_clean <- wave_data_all %>%
  select("YY", "MM", "DD", "hh", "WVHT", "DPD", "season") %>%
  mutate(wave_power = 0.49 * WVHT^2 * DPD,
         wave_datetime = as.POSIXct(sprintf("%04d-%02d-%02d %02d:00:00",
                                            YY, MM, DD, hh),
                                    format = "%Y-%m-%d %H:%M:%S",
                                    tz = "UTC"))

#### 3) Save processed CSVs ####

write.csv(tide_data_clean, here("IntermediateData", "tide_data_processed.csv"), row.names = FALSE, na = "NA")

write.csv(wave_data_clean, here("IntermediateData", "wave_data_processed.csv"), row.names = FALSE, na = "NA")

# CONSPECIFIC DENSITY -----------------------------

#### 1) Read spatial data ####

# read merged Picterra seal detection polygons
picterra.output <- st_read(here("RawData", "merged_polygons.gpkg"), quiet = TRUE)

# read Año Nuevo beach polygons
beaches <- st_read(here("RawData", "beaches.gpkg"), quiet = TRUE)

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
  select(date, age_sex, lat, lon, Beach, density) %>%
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


