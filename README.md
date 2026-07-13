Description of the data and file structure

Multiple datasets were integrated for this project. We calculated mother-offspring association for 936 adult female seals using observational data (resight_data.csv). Weanling mass measurements were collected for a subset of known-age females each year (wean_mass_data.csv). Yearly UAS drone surveys and a map of the colony (beaches.gpkg) allowed us to calculate average seal density of each breeding location (seal_density.csv). We obtained wave height and dominant period records from the National Data Buoy Center station 46042 (wave_data_processed.csv) and tide height data from National Oceanic and Atmospheric Administration station 9413450 (tide_data_processed.csv), which were used to calculate the number of extreme events during each breeding season. Filtering and processing of these datasets are included for transparency in MOA_data_processing.R. We were then able to examine the drivers and consequences of mother-offspring association (the proportion of observations during lactation where a mother was seen with a pup) using generalized linear mixed models and a linear mixed model (MOA_Analysis_v1.R).

For all files, missing values are shown as "NA".

File: MOA_data_pull.csv

Description: This data file contains all processed variables necessary to reproduce the statistical analyses and figures presented in the manuscript. Each row represents data for an individual mother in a given breeding season.
Variables

    animalID: unique maternal identifier

    season: breeding season year

    age: mother's chronological age (years)

    birth_date: offspring birth date (YYYY-MM-DD)

    pupping_exp: maternal reproductive experience, defined as the number of previous pups produced

    total_resights: total number of post-birth observation days for the mother during the breeding season

    count_1_pup: number of post-birth observation days on which the mother was observed with one pup

    MOA_proportion: mother-offspring association (MOA) proportion, calculated as:

    MOA_proportion = count_1_pup / total_resights

    dominant_area: location where the mother was observed most frequently during the breeding season

    avg_density: average conspecific density in the mother's dominant area (seals in a 10m radius)

    n_extreme_both: number of events within a breeding season where extreme wave and tide conditions occurred simultaneously

    pupID: unique offspring identifier

    Wt_wean_corrected: offspring weaning mass

File: seal_density.csv

Description: This file contains conspecific density calculations for each focal female from the UAS drone imagery. Each row is a focal female seal detected in the drone imagery.
Variables

    season: breeding season year
    age_sex: demographic classification of observed individuals (all female)
    lat: seal latitude coordinate (decimal degrees)
    lon: seal longitude coordinate (decimal degrees)
    Beach: beach location name
    density: conspecific density for each focal female in a 10m radius

File: beaches.gpkg

Description: Geospatial polygon layer defining beach boundaries and colony locations used for defining a mother's primary breeding location and calculating local conspecific density.
File: tide_data_processed.csv

Description: Tide height data used to identify extreme tide conditions during breeding seasons. Each row is an hourly measurement.
Variables

    Date: date of tide observation (YYYY-MM-DD)
    time_GMT: time of tide observation in GMT
    tide_height: tide height (ft)
    season: breeding season year
    tide_datetime: combined date-time variable to match the wave data

File: wave_data_processed.csv

Description: Processed wave condition data used to identify extreme wave events during breeding seasons. Each row is an hourly measurement.
Variables

    YY: year
    MM: month
    DD: day
    hh: hour
    WVHT: significant wave height (m)
    DPD: dominant wave period (s)
    season: breeding season year
    wave_power: calculated wave power (kW/m)
    wave_datetime: combined date-time variable to match the tide data

File: resight_data.csv

Description: 
Variables

    animalID: unique maternal identifier
    season: breeding season year
    date: observation date (YYYY-MM-DD)
    area: observed area or location on the colony
    observer: the observer's name
    withpup: pup sighting information (NA, 0,1,2+)
    age: mother's chronological age (years)
    birth_date: approximate pup birth date
    precision: days between the previous observation and the pup's approximate birth date

File: wean_mass_data.csv

Description: Raw offspring weaning mass measurements.
Variables

    season: breeding season year
    pupID: unique offspring identifier
    animalID: unique maternal identifier
    Wt: measured offspring body mass at weighing (kg)
    weighingdate: date offspring was weighed (YYYY-MM-DD)
    weandate: estimated offspring weaning date (YYYY-MM-DD)

Code/software

All code was run using R statistical software (version 4.3.3).
Code File: MOA_Analysis_v1.R

Description: Analysis script used to reproduce all statistical analyses, model outputs, and manuscript figures using MOA_data_pull.csv.
Code File: MOA_data_processing.R

Description: This script reads in and processes the raw data files (density, tide, wave, weaning mass, and resight observations). The output file from this processing is MOA_data_pull.csv.
Access information

Other publicly accessible locations of the data:

    Wave data can be accessed at: https://www.ndbc.noaa.gov/station_history.php?station=46042.
    Tide data can be accessed at: https://tidesandcurrents.noaa.gov/waterlevels.html?id=9413450.


