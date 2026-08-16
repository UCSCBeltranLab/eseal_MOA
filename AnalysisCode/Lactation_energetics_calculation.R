library(tidyverse)

# 1. McHuron milk intake curve --------------------

# S3 from McHuron (2026) contains predicted mass-specific milk intake:
# g milk day^-1 g pup^-0.82
s3 <- read_csv(here("IntermediateData", "McHuron_S3_table.csv"), show_col_types =  FALSE)

# Average milk intake predictions for northern elephant seals across lactation
nes_curve <- s3 %>%
  filter(CommonName == "Northern elephant seal") %>%
  group_by(TimeIntoLactation) %>%
  summarise(intake_rate = mean(response, na.rm = TRUE),
            lower_rate = mean(lower_ci, na.rm = TRUE),
            upper_rate = mean(upper_ci, na.rm = TRUE),
            .groups = "drop")

# Calculate lactation duration for each animal in our dataset
lactation_duration <- metadata %>%
  mutate(date = as.Date(date),
         birth_date = as.Date(birth_date)) %>%
  group_by(animalID, season) %>%
  summarise(birth_date = first(birth_date),
            last_date = max(date, na.rm = TRUE),
            lactation_days = as.numeric(last_date - birth_date),
            .groups = "drop") 

# Biological assumptions
birth_mass_kg <- 40 #from Ortiz et al. 1984
mass_gain_kg <- 91.4 #from Ortiz et al. 1984
lactation_days <- mean(lactation_duration$lactation_days) #from our data
delta_MOA <- 0.10 #change in MOA we want to predict for
retention_fraction <- 0.667 # proportion of milk mass retained as pup mass (Ortiz et al. 1984)

# Assume approximately linear pup growth across lactation (Ortiz et al. 1984)
# convert McHuron's mass-specific rate to absolute kg milk/day
nes_curve <- nes_curve %>%
  mutate(pup_mass_kg = birth_mass_kg + mass_gain_kg * TimeIntoLactation,
         pup_mass_g = pup_mass_kg * 1000,
         milk_kg_day = intake_rate * pup_mass_g^0.82 / 1000,
         milk_lower_kg_day = lower_rate * pup_mass_g^0.82 / 1000,
         milk_upper_kg_day = upper_rate * pup_mass_g^0.82 / 1000,
         lactation_day = TimeIntoLactation * lactation_days)

# 2. Total milk intake across lactation --------------------

# Integrate kg milk/day over the average lactation duration from our females
milk_function <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

total_milk_kg <- milk_function(nes_curve$lactation_day, nes_curve$milk_kg_day)
total_milk_lower <- milk_function(nes_curve$lactation_day, nes_curve$milk_lower_kg_day)
total_milk_upper <- milk_function(nes_curve$lactation_day, nes_curve$milk_upper_kg_day)

# 3. Theoretical effect of a 0.10 MOA difference ----------------

# Treat 0.10 lower MOA as the equivalent of losing 10% of milk intake across lactation!
# This is a theoretical benchmark: withpup = 0 does not necessarily mean that the pup actually missed an entire day of nursing.

equivalent_days_lost <- lactation_days * delta_MOA

milk_difference_kg <- total_milk_kg * delta_MOA
milk_difference_lower <- total_milk_lower * delta_MOA
milk_difference_upper <- total_milk_upper * delta_MOA

# Ortiz et al. 1984: pups retained ~66.7% of milk mass received as body mass
predicted_mass_difference <- milk_difference_kg * retention_fraction
predicted_mass_lower <- milk_difference_lower * retention_fraction
predicted_mass_upper <- milk_difference_upper * retention_fraction

c(total_milk_kg = total_milk_kg,
  equivalent_days_lost = equivalent_days_lost,
  milk_difference_kg = milk_difference_kg,
  predicted_mass_difference = predicted_mass_difference,
  predicted_mass_lower = predicted_mass_lower,
  predicted_mass_upper = predicted_mass_upper)

# 4. Compare theoretical vs. observed mass effect -----------------

# Model predicted weaning-mass difference associated with 0.10 MOA
observed_mass_difference <- fixef(mod_wean_mass)["MOA_proportion"] * delta_MOA

# How large is the observed effect relative to the theoretical full-nursing-loss benchmark?
proportion_observed <- observed_mass_difference / predicted_mass_difference
percent_observed <- proportion_observed * 100

c(observed_mass_difference = observed_mass_difference,
  predicted_mass_difference = predicted_mass_difference,
  proportion_observed = proportion_observed,
  percent_observed = percent_observed)

