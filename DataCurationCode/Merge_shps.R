library(tidyverse)
library(sf)
library(here)

# read in beach shapefile -----------------------------------------------------

# locate Año Nuevo shapefile
list.files(here(), pattern = "Ano Nuevo.*\\.shp$", recursive = TRUE, full.names = TRUE)

# read beach polygons and drop unnecessary ID column
beaches <- st_read(here("RawData", "ANM map/Ano Nuevo Map Final.shp"), quiet = TRUE) %>%
  select(-id)

# save beaches as a geopackage too
st_write(beaches, here("IntermediateData", "beaches.gpkg"), delete_dsn = TRUE)

# read in Picterra seal detections -------------------------------------------

# define directory containing Picterra outputs
base_dir <- here("RawData", "Picterra_outputs_updated")

# find all polygon shapefiles
shps <- list.files(base_dir,
                   recursive = TRUE,
                   pattern = "polygons\\.shp$",
                   full.names = FALSE)

# extract date folder names from file paths
dates <- dirname(shps)

# extract year from date strings
years <- substr(dates, 1, 4)

# function to read individual shapefiles
read_one <- function(rel_shp, d, y) {
  
  # construct full shapefile path
  shp_path <- file.path(base_dir, rel_shp)
  
  # read shapefile and add date/year columns
  x <- st_read(shp_path, quiet = TRUE) %>%
    mutate(date = d,
           year = y)
  
  # assign WGS84 CRS if missing for standardized CRS
  if (is.na(st_crs(x))) st_crs(x) <- 4326
  
  x
}

# read all shapefiles into a list
picterra.list <- pmap(list(rel_shp = shps,
                           d = dates,
                           y = years),
                      read_one) %>%
  compact()

# combine all shapefiles into one sf object
picterra.output <- do.call(rbind, picterra.list)

# save merged Picterra polygons ----------------------------------------------

# save merged polygons as a geopackage
st_write(picterra.output, here("IntermediateData", "merged_polygons.gpkg"), delete_dsn = TRUE)

