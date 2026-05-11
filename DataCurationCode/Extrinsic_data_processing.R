
# Setup: load libraries and data -----------------------------------------------

library(here)
library(tidyverse)
library(lubridate)

##suppress scientific notation
options(scipen = 999)

##turn off spherical geometry
sf_use_s2(FALSE)

# WAVE AND TIDE ---------------------------

##Identify seasonal date range: 
##Calculate earliest and latest dates are in the breeding season
##Used for: selecting dates when uploading the wave and tide data
dates <- metadata %>% 
  group_by(season) %>%
  summarize(earliest_date = min(date, na.rm = TRUE), 
            latest_date = max(date, na.rm = TRUE))

#### 1) Read in and clean tide data ####

#list only CSV files containing "tides" in the filename
tide_files <- list.files(path = "./RawData/", 
                         pattern = "tides.*\\.csv$", 
                         full.names = TRUE)

#read and combine all CSVs
tide_data_all <- tide_files %>%
  lapply(read_csv) %>%
  bind_rows()

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
wave_files <- list.files(path = "./RawData/", 
                         pattern = "waves.*\\.csv$", 
                         full.names = TRUE)

#read and combine only the wave data files
wave_data_all <- wave_files %>%
  map_dfr(~ read_csv(.x, colClasses = "character"))

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

# read in beach shapefile
list.files("..", pattern = "Ano Nuevo.*\\.shp$", recursive = TRUE, full.names = TRUE)

beaches <- st_read("./RawData/ANM map/Ano Nuevo Map Final.shp", quiet = TRUE) %>%
  select(-id)

# read in seal detections
base_dir <- "./RawData/Picterra_outputs_updated"

shps <- list.files(base_dir, recursive = TRUE, pattern = "polygons\\.shp$", full.names = FALSE)

dates <- dirname(shps)
years <- substr(dates, 1, 4)

read_one <- function(rel_shp, d, y) {
  
  shp_path <- file.path(base_dir, rel_shp)
  
  x <- st_read(shp_path, quiet = TRUE) %>%
    mutate(date = d, year = y)
  
  if (is.na(st_crs(x))) st_crs(x) <- 4326
  x
}

picterra.list <- pmap(list(rel_shp = shps, d = dates, y = years), read_one) |> compact()

picterra.output <- do.call(rbind, picterra.list)

### Convert to NAD83 / California zone 3 ###
st_crs(picterra.output)

picterra.output <- st_transform(picterra.output, "EPSG:26943")

st_crs(picterra.output)

beaches <- st_transform(beaches, "EPSG:26943")

#### 2) Assign age/sex class to seals ####
uas.data <- picterra.output %>%
  mutate(age_sex = tolower(as.character(age_sex))) %>%
  filter(age_sex %in% c("female", "male", "pup"))

#### 3) Assign seals to beach ####
centroids <- st_centroid(uas.data)

seals_by_beach <- st_intersection(beaches, centroids)

#### 4) Calculate female density ####
density.df <- data.frame()

for (i in 1:length(dates)) {
  
  survey.subset <- seals_by_beach %>%
    filter(date == dates[i],
           age_sex == "female")
  
  if (nrow(survey.subset) == 0) next
  
  seal.centroids <- st_centroid(survey.subset)
  
  ## set to 10m, but you can modify that to whatever you would like!
  seal.buffer <- st_buffer(seal.centroids, 10)
  
  int <- st_intersects(seal.buffer, seal.centroids)
  
  survey.subset$density <- lengths(int) - 1
  
  density.df <- rbind(density.df, survey.subset)
}

#### 5) Save density dataframe ####

seal_density <- density.df %>%
  mutate(centroid = st_centroid(geometry),
         lon = st_coordinates(centroid)[,1],
         lat = st_coordinates(centroid)[,2]) %>%
  select(date, age_sex, lat, lon, Beach, density) %>%
  st_drop_geometry()

write.csv(seal_density, here("IntermediateData", "seal_density.csv"), row.names = FALSE, na = "NA")

#### 6) Plot two density radius examples ####

###low density version 

# grab a LIGHT color from mako
light_mako <- viridis(190, option = "mako", direction = -1)[20]

# find a seal with 5 neighbors
target_n <- 5
best_diff <- Inf
low_example <- NULL

for (i in 1:nrow(uas.data)) {
  
  focal <- uas.data[i, ]
  
  if (focal$age_sex[1] != "female") next
  
  same <- uas.data %>%
    filter(date == focal$date[1],
           age_sex == "female")
  
  hits <- st_intersects(
    st_buffer(st_centroid(focal), 10),
    st_centroid(same)
  )[[1]]
  
  n_neighbors <- length(hits) - 1
  
  if (n_neighbors == target_n) {
    low_example <- focal
    break
  }
  
  if (abs(n_neighbors - target_n) < best_diff) {
    best_diff <- abs(n_neighbors - target_n)
    low_example <- focal
  }
}

# plot
seal.buffer <- st_buffer(st_centroid(low_example), 10)

same <- uas.data %>%
  filter(date == low_example$date[1],
         age_sex == "female")

hits <- st_intersects(seal.buffer, st_centroid(same))[[1]]

int.poly <- same[hits[-1], ]   # drop focal seal itself

png("./TablesFigures/density_example_low.png",
    width = 2000,
    height = 2000,
    res = 300,
    bg = "transparent")   # <-- makes outside transparent

par(bg = NA, mar = c(0, 0, 0, 0))  # remove outer margins

plot(st_geometry(seal.buffer), border = "grey40")

# surrounding seals
plot(st_geometry(int.poly),
     add = TRUE,
     col = light_mako,
     border = "grey40")

# focal seal
plot(st_geometry(low_example),
     add = TRUE,
     col = light_mako,
     border = "black",
     lwd = 3)

dev.off()

###high density version

dark_mako <- viridis(50, option = "mako", direction = 1)[20]

# find a seal with MANY neighbors (target ~18)
target_n <- 18
best_diff <- Inf
high_example <- NULL

for (i in 1:nrow(uas.data)) {
  
  focal <- uas.data[i, ]
  
  if (focal$age_sex[1] != "female") next
  
  same <- uas.data %>%
    filter(date == focal$date[1],
           age_sex == "female")
  
  hits <- st_intersects(
    st_buffer(st_centroid(focal), 10),
    st_centroid(same)
  )[[1]]
  
  n_neighbors <- length(hits) - 1
  
  if (n_neighbors == target_n) {
    high_example <- focal
    break
  }
  
  if (abs(n_neighbors - target_n) < best_diff) {
    best_diff <- abs(n_neighbors - target_n)
    high_example <- focal
  }
}

# plot
seal.buffer <- st_buffer(st_centroid(high_example), 10)

same <- uas.data %>%
  filter(date == high_example$date[1],
         age_sex == "female")

hits <- st_intersects(seal.buffer, st_centroid(same))[[1]]
int.poly <- same[hits[-1], ]   # drop focal seal

png("./TablesFigures/density_example_high.png",
    width = 2000,
    height = 2000,
    res = 300,
    bg = "transparent")   # outside = transparent

par(bg = NA, mar = c(0, 0, 0, 0))  # remove margins

plot(st_geometry(seal.buffer), border = "grey40")

# surrounding seals
plot(st_geometry(int.poly),
     add = TRUE,
     col = dark_mako,
     border = "grey40")

# focal seal
plot(st_geometry(high_example),
     add = TRUE,
     col = dark_mako,
     border = "black",
     lwd = 5)

dev.off()

