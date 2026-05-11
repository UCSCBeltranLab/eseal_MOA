
#Setup: combine and read in raw data ----------------------------------

### Load necessary packages ###

library(here)
library(conflicted)
library(tidyverse)
library(lubridate)
conflicts_prefer(dplyr::lag)

# read in tide data
tide_data_clean <- read_csv(here("IntermediateData", "tide_data_processed.csv"), show_col_types = FALSE)

# read in wave data
wave_data_clean <- read_csv(here("IntermediateData", "wave_data_processed.csv"), show_col_types = FALSE)

# read in seal density data
seal_density_data <- read_csv(here("IntermediateData", "seal_density.csv"), show_col_types = FALSE)

# read in seal density data
weaner_data <- read_csv(here("RawData", "weanling_weight_master.csv"), show_col_types = FALSE)

# read in the raw resight data
# includes animalID, season, date, area, entry date, entry = observer, withpup = pup status of the observed female (0,1,2,etc.)
raw_data <- read_csv(here("RawData", "Aditi 2024 Data Pull 2025_04_17_RAW.csv"), show_col_types = FALSE)

# read in summarized individual seal data
# specific per-animal data: includes animalID, age, season, pup birth date with precision
summarized_data <- read_csv(here("RawData", "Aditi 2024 Data Pull 2025_04_17_SUMMARIZED.csv"), show_col_types = FALSE)

# Resight data processing -----------------------------------

##how many individuals/sightings are in the sample pre-filtering
raw_data %>%
  mutate(date = as.Date(date)) %>%
  summarise(n_seals = n_distinct(animalID), n_sightings = n())

raw_data %>%
  filter(!(withpup %in% c(0, 1))) %>%
  summarize(n_obs = n())

##clean the raw resight data
raw_data <- raw_data %>%
  arrange(animalID) %>%
  group_by(animalID, date) %>%
  slice_min(order_by = entry, n = 1, with_ties = FALSE) %>% #gets rid of double resights
  filter(withpup %in% c(0, 1)) %>% #restrict to only values of withhpup score = 0 or 1 (either was not seen with a pup or was seen with 1 pup)
  ungroup() 

### Calculate previous pupping experience from the raw data ###

pupping_experience <- raw_data %>%
  group_by(animalID, season) %>%
  summarize(had_pup = as.integer(any(withpup == 1)), #any observation of 1 pup counts as pupping experience
            .groups = "drop") %>%
  group_by(animalID) %>%
  arrange(season, .by_group = TRUE) %>%  #make sure season is chronological
  mutate(pupping_exp = lag(cumsum(had_pup == 1), default = 0)) %>% #calculate prior pup experience
  ungroup()

##join with raw data to include previous pupping experience
raw_data <- raw_data %>%
  left_join(pupping_experience, by = c("animalID", "season"))

##create a metadata table combining summarized and raw data
metadata <- summarized_data %>%
  left_join(raw_data, by = c("animalID", "season"))

#filter metadata so it only includes precise birth_dates and removes weird animals
metadata <- metadata %>%
  rename(age = AgeYears) %>%
  rename(birth_date = BirthDate) %>%
  filter(!is.na(birth_date)) %>% # only known birth_date moms
  filter(Precision <= 7) %>% # pup birth_date precision within a week
  filter(age <= 24) # removes animals outside of normal age range

### Resight effort and set up for calculating mother-offspring association ###

metadata <- metadata %>%
  mutate(date = as.Date(date), 
         birth_date = as.Date(birth_date)) %>% #ensure both are dates
  group_by(season, animalID) %>%
  mutate(total_resights = sum(date >= birth_date, na.rm = TRUE), # total post-birth resights
         count_1_pup = sum(withpup == 1 & date >= birth_date, na.rm = TRUE)) %>% # post-birth 1-pup sightings
  filter(total_resights >= 14) %>% # keep only animals with ≥14 observations post birth
  ungroup()

### Assigning each female pupping beach locations ###

# 1) Count sightings for each female in each area 
area_counts <- metadata %>%
  group_by(animalID, season, area) %>%
  summarize(count = n(), .groups = "drop") #count number of times each animal was seen in each area

## 2) for each animal, find the area with the maximum count
harem_assignment <- area_counts %>%
  group_by(animalID, season) %>%
  filter(count == max(count)) %>% #these are assigned areas for each female
  ungroup() %>%
  select(animalID, season, dominant_area = area) #call dominant area to avoid confusion

# Wave and tide processing -----------------------------------

# 1) Rebuild date/time columns from CSV
tide_data_clean <- tide_data_clean %>%
  mutate(Date = as.Date(as.character(Date)),
         tide_datetime = as.POSIXct(paste(Date, time_GMT),
                                    format = "%Y-%m-%d %H:%M",
                                    tz = "UTC"))

wave_data_clean <- wave_data_clean %>%
  mutate(wave_datetime = as.POSIXct(sprintf("%04d-%02d-%02d %02d:00:00",
                                            YY, MM, DD, hh),
                                    format = "%Y-%m-%d %H:%M:%S",
                                    tz = "UTC"))

# 2) Join the data sets by season
tide_wave <- left_join(wave_data_clean,
                       tide_data_clean,
                       by = c("wave_datetime" = "tide_datetime", 
                              "season")) %>%
  filter(!is.na(tide_height))

# 3) Calculate percentile thresholds
quantile(tide_wave$wave_power,
         probs = c(0.9, 0.95, 0.975, 0.99),
         na.rm = TRUE)

quantile(tide_wave$tide_height,
         probs = c(0.9, 0.95, 0.975, 0.99),
         na.rm = TRUE)

# 4) Flag and categorize extreme events based on 90th percentile
extreme_wave_threshold <- 110.8388
extreme_tide_threshold <- 5.03

tide_wave_flagged <- tide_wave %>%
  mutate(extreme_wave = wave_power >= extreme_wave_threshold,
         extreme_tide = tide_height >= extreme_tide_threshold,
         extreme_both = extreme_wave & extreme_tide) %>%
  group_by(season) %>%
  summarize(n_extreme_both = sum(extreme_both),
            n_extreme_tide = sum(extreme_tide),
            n_extreme_wave = sum(extreme_wave))

# 5) Plot extreme tide and wave events per season
ggplot(tide_wave_flagged, aes(x = season, y = n_extreme_both)) +
  geom_line(linewidth = 1.2, color = "#1f78b4") +
  geom_point(color = "darkblue") +
  scale_x_continuous(n.breaks = 28) +
  labs(title = "Extreme Tide and Storm Events (1996–2025)",
       x = "Year", 
       y = "Number of Extreme Events") +
  theme_minimal()

 # Conspecific density data processing -----------------------------------------

# 1) extract year (YYYY) from date string and store as season to match the resight data
seal_density_data$season <- as.integer(substr(seal_density_data$date, 1, 4))

# 2a) Calculate average density per area-season for modeling, then link to assigned harems
area_density <- seal_density_data %>%
  rename(dominant_area = Beach) %>%
  group_by(dominant_area, season) %>%
  summarize(avg_density = mean(density)) %>% #calculate mean density per area in each season 
  ungroup() %>%
  left_join(harem_assignment, by = c("season", "dominant_area")) %>% #join with location info per animal
  filter(!is.na(animalID)) #only keep animals observed in the 2016-2025 dataset to match the drone data

# 2b) Data for plot below (one row per beach x season)
area_density_plot <- area_density %>%
  group_by(dominant_area, season) %>%
  summarize(avg_density = first(avg_density),
            n_animals = n_distinct(animalID), 
            .groups = "drop")

##Quick check: area density per-season at each location
ggplot(area_density_plot, aes(x = dominant_area,
                              y = avg_density,
                              fill = avg_density)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ season, scales = "free_x") +
  scale_fill_gradient(low  = "#DEEBF7",
                      high = "navy") +
  labs(x = "Beach",
       y = "Average Density",
       fill = "Average Density") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3) Clean up the area density data
# Identify animal-season cases with multiple locations (ties)
location_duplicates <- area_density %>%
  group_by(season, animalID) %>%
  filter(n() > 1) %>%
  arrange(animalID, season)

# Break ties for the 17 duplicates by selecting the location with higher average density
area_density <- area_density %>%
  group_by(season, animalID) %>%
  slice_max(avg_density, with_ties = FALSE) %>%
  ungroup()

# Pup weaning mass data processing -------------------------------------------------

##necessary filtering
weaner_data <- weaner_data %>% 
  rename(pupID = animalID) %>%
  filter(!is.na(MomAnimalID)) %>% #only known moms
  filter(!is.na(year)) #only known years

##what does the sample size look like?
n_distinct(weaner_data$MomAnimalID)
n_distinct(weaner_data$year)

##are there duplicates?
wean_dup_cases <- weaner_data %>%
  count(MomAnimalID, year) %>%
  filter(n > 1)

##why are there duplicates (weighed same pup twice? 2 different pups?)
wean_dup_rows <- weaner_data %>%
  semi_join(wean_dup_cases, by = c("MomAnimalID", "year")) %>%
  arrange(MomAnimalID, year)

# 1) identify mother-year cases with multiple pup IDs
multi_pup_cases <- weaner_data %>%
  group_by(MomAnimalID, year) %>%
  summarise(n_pups = n_distinct(pupID),
            .groups = "drop") %>%
  filter(n_pups > 1)

# 2) remove ambiguous multi-pup mother-year cases
weaner_data <- weaner_data %>%
  anti_join(multi_pup_cases, by = c("MomAnimalID", "year"))

# 3) filter wean weight data
weaner_data <- weaner_data %>%
  select(MomAnimalID, pupID, sex, year, Wt, weighingdate, weandate) %>%
  rename(animalID = MomAnimalID,
         season = year) %>%
  mutate(weighingdate = as.Date(weighingdate),
         weandate = as.Date(weandate),
         days_post_wean = as.numeric(weighingdate - weandate)) %>% #calculate the number of days between weaning and weighing
  filter(days_post_wean >= 0) %>% #only animals weighed as actual weans
  filter(days_post_wean <= 71) %>% #excluding wildly wrong wean dates
  filter(Wt <= 170) %>% #no super weaners
  mutate(Wt_wean_corrected = Wt * exp(0.00596 * days_post_wean)) #back-correction for wean mass based on fasting rate

##check the remaining duplicates (the same pup weighed twice post-weaning)
final_wean_duplicates <- weaner_data %>%
  group_by(animalID, season) %>%
  filter(n() > 1)

# 4) for these 12 remaining duplicates select the weight measurement closest to the weaning date
weaner_data <- weaner_data %>%
  group_by(animalID, season) %>%
  slice_min(days_post_wean, n = 1, with_ties = FALSE) %>% #break ties by the weight closer to the weandate
  ungroup()

# Build final modeling data frame ---------------------------------------------------

##make table with all necessary variables for analyses and figures
processed_data <- metadata %>%
  select(animalID, season, age, date, birth_date, pupping_exp, withpup, total_resights, count_1_pup) %>%
  mutate(birth_date = as.Date(birth_date)) %>%
  left_join(area_density %>% select(animalID, season, dominant_area, avg_density), by = c("animalID", "season")) %>%
  left_join(tide_wave_flagged %>% select(season, n_extreme_both), by = "season") %>%
  left_join(weaner_data %>% select(animalID, pupID, season, Wt_wean_corrected), by = c("animalID", "season"))

##write final CSV for all modeling variables
write.csv(processed_data, here("IntermediateData", "MOA_data_pull.csv"), row.names = FALSE, na = "NA")

          