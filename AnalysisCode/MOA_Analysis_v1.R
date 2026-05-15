# Setup ----------------------------------------------------------

lapply(
  c("tidyverse",
    "scales",
    "here",
    
    # modeling
    "lme4",
    "lmerTest",
    "merTools",
    "DHARMa",
    "performance",
    "boot",
    "emmeans",
    
    # plotting
    "ggthemes",
    "cowplot",
    "viridis",
    "ggtext",
    "ggspatial",
    
    # tables
    "flextable",
    "broom.mixed",
    
    # spatial
    "sf",
    "lubridate",
    "rnaturalearth",
    "rnaturalearthdata",
    "rnaturalearthhires",
    "maptiles",
    "tidyterra"),
  require,
  character.only = TRUE)

library(conflicted)

# Prefer dplyr for common conflicts
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("count", "dplyr")

### Bootstrap setup ###
set.seed(123)
nsim <- 100 # keep 100 simulations

### Read modeling variables csv ###

data_pull <- read_csv(here("IntermediateData", "MOA_data_pull.csv"), show_col_types = FALSE)

#### Calculate mother-offspring association ####

# mother-offspring association = 
# count_1_pup (number of days observed with 1 pup during lactation) / total_resights (total days observed during lactation)

data_pull <- data_pull %>%
  mutate(MOA_proportion = count_1_pup / total_resights) %>%
  mutate(date = as.Date(date),
         birth_date = as.Date(birth_date)) %>%
  filter(!is.na(MOA_proportion), MOA_proportion>0) #ensure valid association values

### Collapse to 1 row per seal-season for modeling ###

model_variables <- data_pull %>%
  distinct(animalID, season, .keep_all = TRUE)
  
##Setting age threshold and experience threshold (see Tables S10 and S11)
age_thresh <- 9
exp_thresh <- 5
  
##Setting thresholds in data
model_variables <- model_variables %>%
  mutate(age_cat = factor(ifelse(age < age_senesce, "Young", "Old"),
                            levels = c("Young", "Old"))) %>%
  mutate(age10 = (age - age_senesce) / 10) %>% #scaled numeric version of age centered at the threshold
  mutate(experience_cat = factor(ifelse(pupping_exp < exp_thresh, "Inexperienced", "Experienced"),
                        levels = c("Inexperienced", "Experienced"))) %>%
  mutate(exp10 = (pupping_exp - exp_thresh) / 10) %>% #scaled numeric version of previous pupping experience centered at the threshold
  mutate(animalID_fct = factor(animalID), #animalID factor
         season_fct = factor(season)) #season factor

#Summary of the data
data_summary <- model_variables %>%
  summarise(
    n_animal_ID = n_distinct(animalID),
    n_season = n_distinct(season),
    min_age = min(age),
    max_age = max(age),
    mean_MOA = mean(MOA_proportion),
    sd_MOA = sd(MOA_proportion),
    min_MOA = min(MOA_proportion),
    max_MOA = max(MOA_proportion),
    q05_MOA = quantile(MOA_proportion, 0.05),
    q95_MOA = quantile(MOA_proportion, 0.95),
    mean_density = mean(avg_density, na.rm = TRUE),
    sd_density = sd(avg_density, na.rm = TRUE),
    min_density = min(avg_density, na.rm = TRUE),
    max_density = max(avg_density, na.rm = TRUE),
    min_extreme = min(n_extreme_both, na.rm = TRUE),
    max_extreme = max(n_extreme_both, na.rm = TRUE))

#### Pearson's correlations for model predictors ####

vars <- model_variables[, c("age", "avg_density", "n_extreme_both", "pupping_exp")]

cor(vars, use = "pairwise.complete.obs", method = "pearson")

# MODELING FRAMEWORK --------------------

# For Models 1a/1b and 2a/2b the response is MOA as a proportion in (0,1] (where successes = the number of observations a mother was seen with 1 pup)
# All models are a binomial with logit link, weighted by the total observations (weights = total_resights)
# We use bobyqa as a robust optimizer to improve model convergence
# Random effects for animalID and season are included

## Model 1a includes age, conspecific density, and extreme wave and tide events using a smaller subset of data (2016-2023)
## Model 1b includes previous pupping experience, conspecific density, and extreme wave and tide events using a smaller subset of data (2016-2023)
## Model 2a is a piecewise segmented regression for age using the full dataset (1996-2025)
## Model 2b is a piecewise segmented regression for previous pupping experience using the full dataset (1996-2025)

#### Model 1a: 2016-2023 subset model with age, density, extreme ####

#make subset dataframe for the 2016-2023 models and figures (makes code for figures easier)
model_variables_2016_2023 <- model_variables %>%
  filter(!is.na(n_extreme_both)) %>% #no NAs
  filter(!is.na(avg_density)) %>% #no NAs
  mutate(season_fct = droplevels(season_fct))

## model for 2016-2023 subset to test all predictors together
## age = maternal age as a linear predictor
## avg_density = conspecific density of each seal's observed location 
## n_extreme_both = the number of per-year extreme wave and tide events
## random effects of animalID and year

# 2) all predictors model
mod_age_2016_2023 <- glmer(MOA_proportion ~ age + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                             weights = total_resights,
                             family = binomial(link = "logit"),
                             control = glmerControl(optimizer = "bobyqa"),
                             data = model_variables_2016_2023); summary(mod_age_2016_2023)

ranef(mod_age_2016_2023) #random effect variance for each animalID and season
exp(fixef(mod_age_2016_2023)) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_age_2016_2023, plot = TRUE) #plot residuals
check_collinearity(mod_age_2016_2023) #check predictor VIFs

#### Model 1b:  2016-2023 subset model with pupping experience, density, extreme ####

## pupping_exp = number of prior breeding seasons observed at least once with a pup
## avg_density = conspecific density of each seal's observed location 
## n_extreme_both = the number of per-year extreme wave and tide events
## random effects of animalID and year

mod_exp_2016_2023 <- glmer(MOA_proportion ~ pupping_exp + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                                 weights = total_resights,
                                 family = binomial(link = "logit"),
                                 control = glmerControl(optimizer = "bobyqa"),
                                 data = model_variables_2016_2023); summary(mod_exp_2016_2023) #model summary

ranef(mod_exp_2016_2023) #random effect variance for each animalID and season
exp(fixef(mod_exp_2016_2023)) #converts fixed-effect log-odds to odds ratios
simulateResiduals(mod_exp_2016_2023, plot = TRUE) #plot residuals
check_collinearity(mod_exp_2016_2023) #check predictor VIFs

#### Model 2a: 1996-2025 full dataset model for age ####

## breakpoint model for age
## age_cat = "Young", "Old" based on age_senesce
## age10 = (age - age_senesce) / 10) scaled numeric version of age centered at senescence threshold
## random effects of animalID and year

mod_age_1996_2025 <- glmer(MOA_proportion ~ age_cat : age10 + (1 | animalID_fct) + (1 | season_fct),
                           weights = total_resights,
                           family = binomial(link = "logit"),
                           control = glmerControl(optimizer = "bobyqa"),
                           data = model_variables); summary(mod_age_1996_2025)

ranef(mod_age_1996_2025) #random effect variance for each animalID and season
exp(fixef(mod_age_1996_2025)) #converts fixed-effect log-odds to odds ratios
resid_age <- simulateResiduals(mod_age_1996_2025, plot = TRUE) #plot residuals
plotResiduals(resid_age, form = model_variables$age_cat) #age cat residuals

#### Model 2b: 1996-2025 full dataset model for pupping experience ####

## breakpoint model for experience
## experience_cat = "inexperienced" or "experienced", based on <5 or >= 5
## exp10 = (pupping_exp - 5) / 10) scaled numeric version of experience centered at the threshold
## random effects of animalID and year

# 2) Fit pupping experience model
mod_exp_1996_2025 <- glmer(MOA_proportion ~ experience_cat : exp10 + (1 | animalID_fct) + (1 | season_fct),
                           weights = total_resights,
                           family = binomial(link = "logit"),
                           control = glmerControl(optimizer = "bobyqa"),
                           data = model_variables); summary(mod_exp_1996_2025)

ranef(mod_exp_1996_2025) #random effect variance for each animalID and season
exp(fixef(mod_exp_1996_2025)) #converts fixed-effect log-odds to odds ratios
resid_exp <- simulateResiduals(mod_exp_1996_2025, plot = TRUE) #plot residuals
plotResiduals(resid_exp, form = model_variables$experience_cat) #experience cat residuals

#### Model 3: Offspring fitness consequences ####

## linear model for weaning mass against mother-offspring association
## age = maternal age
## random effects of animalID and year

##check distribution for normality
ggplot(model_variables, aes(x = Wt_wean_corrected)) +
  geom_histogram(bins = 100, color = "black", fill = "skyblue") +
  theme_classic() +
  labs(title = "Distribution of Weaning Weight",
       x = "Weight (kg)",
       y = "Count")

#fit linear model
mod_wean_mass <- lmerTest::lmer(Wt_wean_corrected ~ MOA_proportion + age + (1 | season_fct) + (1 | animalID_fct),
                                data = model_variables); summary(mod_wean_mass)
fixef(mod_wean_mass)["MOA_proportion"] * 0.1
ranef(mod_wean_mass) #random effect variance for each animalID and season
simulateResiduals(mod_wean_mass, plot = TRUE) #plot residuals
check_collinearity(mod_wean_mass) #check predictor VIFs

# MAIN FIGURES ----------------------------------------------

#### Figure 1a (MOA conceptual figure) ####

# 1) Choose 3 seal mothers from 2024
example_ids <- c("44397", "48257", "51872") #example seal IDs

# 2) Build resight tiles
resight_tiles <- data_pull %>%
  mutate(date = as.Date(date), #convert resight date
         birth_date = as.Date(birth_date), #convert approximate birth date
         withpup = as.character(withpup)) %>% #make withpup character for labels
  filter(season == 2024, #keep 2024 season
         animalID %in% example_ids, #keep example mothers
         !is.na(date), #remove missing resight dates
         !is.na(withpup)) %>% #remove missing pup status
  distinct(animalID, season, date, .keep_all = TRUE) %>% #one tile per seal-date
  mutate(pup_status = case_when(withpup == "0" ~ "Observed with no pup", #blue: observed without pup
                                withpup == "1" ~ "Observed with one pup")) %>% #green: observed with pup
  select(animalID, season, date, age, MOA_proportion, birth_date, pup_status, withpup) #keep plot columns

# 3) Add approximate birth tiles
MOA_plot_df <- bind_rows(resight_tiles, #observed resight tiles
                         resight_tiles %>%
                           filter(!is.na(birth_date)) %>% #keep known birth dates
                           distinct(animalID, season, birth_date, .keep_all = TRUE) %>% #one birth tile per seal
                           transmute(animalID, season,
                                     date = birth_date, #use birth date as tile date
                                     age, MOA_proportion,
                                     pup_status = "Approximate birth", #pink: approximate birth tile
                                     withpup = "1")) %>% #print 1 on birth tile
  mutate(pup_status = factor(pup_status, levels = c("Observed with no pup",
                                                    "Approximate birth",
                                                    "Observed with one pup"))) %>% #legend order: blue, pink, green
  arrange(animalID, date, pup_status) %>% #order rows before removing duplicates
  distinct(animalID, season, date, .keep_all = TRUE) %>% #one final tile per seal-date
  mutate(panel_lab = paste0("**AnimalID = ", animalID, "**", #to make each animal bold in label
                            "<br>Age = ", age,
                            "<br><span style='color:",
                            ifelse(MOA_proportion == 1, "#6F8D41", "#913161"), #match colors to Fig 1b
                            ";'>Association = ",
                            sprintf("%.2f", MOA_proportion),
                            "</span>"))
# 4) Plot resight tiles
conceptual_MOA_plot <- ggplot(MOA_plot_df, aes(date, 1, fill = pup_status)) +
  geom_tile(width = 0.97, height = 2.5, color = "white", linewidth = 0.3) + #daily tiles
  geom_text(data = filter(MOA_plot_df, withpup %in% c("0", "1")), #0/1 labels on observed tiles
            aes(label = withpup),
            size = 5,
            fontface = "bold") +
  facet_grid(panel_lab ~ ., switch = "y") + #one row per mother
  scale_fill_manual(name = "Pup status", #legend with written-out categories
                    values = c("Observed with no pup" = "#8fd0eb",
                               "Approximate birth" = "#f2b6c6",
                               "Observed with one pup" = "#86d37c"),
                    drop = FALSE) +
  scale_x_date(breaks = seq(min(MOA_plot_df$date),
                 max(MOA_plot_df$date),
                 by = "3 days"), #start at first date, breaks every 3 days
               date_labels = "%b %d",
               expand = c(0, 0)) +
  scale_y_continuous(breaks = NULL, expand = c(0, 0)) + #remove y-axis
  labs(x = "Date", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 22) +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_markdown(angle = 0, hjust = 0,
                                             lineheight = 1.3, margin = margin(r = 20)),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_text(face = "bold", hjust = 0.52, margin = margin(t = 15)),
        panel.spacing.y = unit(1, "lines"),
        legend.position = "top",
        legend.title = element_text(face = "bold", size = 19)); conceptual_MOA_plot

ggsave(here("TablesFigures", "Figure1a.png"), conceptual_MOA_plot, height = 7, width = 15, dpi = 1000)

#### Figure 1b (MOA distribution across years) ####

# 1) Format season and categorize MOA values
MOA_data_descriptive <- model_variables %>%
  mutate(prop_cat = case_when(MOA_proportion == 1 ~ "1", #perfect association
                              MOA_proportion < 1 ~ "< 1")) #below perfect association

# 2) Calculate yearly percentages
MOA_percent_data <- MOA_data_descriptive %>%
  group_by(season_fct, prop_cat) %>% #count each category per year
  summarise(n_cat = n(), .groups = "drop") %>%
  group_by(season_fct) %>%
  mutate(total_n = sum(n_cat), #sample size per year
         percent = n_cat / total_n) %>% #category percentage
  ungroup()

# 3) Get yearly sample-size labels
n_labels <- MOA_percent_data %>%
  distinct(season_fct, total_n) #one n label per year

# 4) Plot yearly MOA distribution
MOA_distribution_years <- ggplot(MOA_percent_data,
                                 aes(season_fct, percent, fill = prop_cat)) +
  geom_col(position = "stack") + #stacked yearly percentages
  geom_text(data = n_labels, #sample-size labels
            aes(season_fct, 1.03, label = total_n),
            inherit.aes = FALSE,
            size = 5) +
  scale_y_continuous(labels = percent_format(accuracy = 1), #percent y-axis
                     limits = c(0, 1.05),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(breaks = levels(MOA_data_descriptive$season_fct),
                   drop = FALSE) +
  scale_fill_manual(values = c("1" = "#6F8D41",
                               "< 1" = "#913161")) +
  labs(x = "Year",
       y = "Percentage",
       fill = "Mother-offspring\nassociation values") +
  theme_classic(base_size = 23) +
  theme(axis.title.x = element_text(margin = margin(t = 15)),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.title = element_text(size = 20, face = "bold"),
        legend.margin = margin(r = -25),
        legend.key.width = unit(1, "cm"),
        legend.position = "left",
        legend.justification = "top"); MOA_distribution_years

ggsave(here("TablesFigures", "Figure1b.png"), MOA_distribution_years, height = 7, width = 17, dpi = 1000)

#### Figure 2a (breakpoint age, 1996-2025) ####

# 1) Set plotting colors for young/old
AGECOL <- c(Young = "#92BAEE", Old = "#EB99D2")

# 2) Define x-axis ages and the set of seasons to draw thin season-specific curves
ages <- sort(unique(model_variables$age)) #unique ages present in data
seasons <- levels(model_variables$season_fct) #factor levels for seasons (used in predictions)

# 3) Compute pooled observed MOA_proportion by age (black points + binomial CI whiskers)
observed_data_1996_2025 <- model_variables %>%
  group_by(age) %>% 
  summarise(n_age = n(), #sample size at each age
            n_success = sum(count_1_pup, na.rm = TRUE), #total successes at each age
            n_trials = sum(total_resights, na.rm = TRUE), #total trials/effort at each age
            avg_prop = n_success/n_trials, #pooled observed MOA_proportion at each age
            lwr = binom.test(n_success, n_trials)$conf.int[1], #binomial CI lower
            upr = binom.test(n_success, n_trials)$conf.int[2], #binomial CI upper
            .groups = "drop") #return ungrouped tibble

# 4) Build standardized newdata for each age A (for post-stratification)
# For each A, keep the full covariate mix of model_variables but overwrite age variables as if everyone were age A
nd_by_age <- lapply(ages, \(A) transform(model_variables, #start from full dataset to preserve covariate distribution and weights
                                         age = A, #overwrite age so every row is evaluated at age A
                                         age_cat = factor(ifelse(A < age_senesce, "Young", "Old"), levels = c("Young","Old")), 
                                         age10 = (A - age_senesce)/10)) 

# 5) Post-stratified curve: predict for each age, then compute weighted average using total_resights
# season = NULL gives overall curve, season ="" forces season RE for thin curves
post_curve <- function(fit, season = NULL) {
  vapply(seq_along(nd_by_age), function(a){ #loop over ages
    nd <- nd_by_age[[a]] #newdata for one age
    if(!is.null(season)) nd$season_fct <- factor(season, levels = seasons) #force all rows to a single season level
    p <- predict(fit, newdata = nd, type = "response", re.form = NULL) #predict probability including all random effects
    weighted.mean(p, w = nd$total_resights, na.rm = TRUE) #post-stratified mean, weighted by total resights
  }, numeric(1)) #numeric vector of predictions across ages
}

# 6) Thin season lines: for each season, compute the post-stratified curve
season_lines_1996_2025 <- bind_rows(lapply(seasons, \(S)
                                           tibble(season_fct = S, #season for grouping/plotting
                                                  age = ages, #x-axis ages
                                                  pred = post_curve(mod_age_1996_2025, season = S)))) #predicted curve for each season

# 7) Bootstrap CI for the overall curve using parametric bootstrapping of the fitted GLMM
post_curve_overall <- function(fit) post_curve(fit, season = NULL) #overall curve only

overall_1996_2025 <- tibble(age = ages, #continuous age in years
                            pred = post_curve_overall(mod_age_1996_2025)) %>% #overall post-stratified predictions
  left_join(observed_data_1996_2025 %>% select(age, n_age), by = "age") #add labels for sample size per age

boot_1996_2025 <- bootMer(mod_age_1996_2025, FUN = post_curve_overall, nsim = nsim, #bootstrap curves
                          type = "parametric", use.u = FALSE) #parametric bootstrap simulate ranefs each time

# 8) Turn bootstrap curves into 95% CI bands; define Young/Old segment for coloring
ci_1996_2025 <- overall_1996_2025 %>%
  mutate(lo = apply(boot_1996_2025$t, 2, quantile, 0.025, na.rm = TRUE), #2.5% quantile at each age
         hi = apply(boot_1996_2025$t, 2, quantile, 0.975, na.rm = TRUE), #97.5% quantile at each age
         age_cat = factor(ifelse(age < age_senesce, "Young", "Old"), levels = c("Young","Old")), #segment label for each age
         seg = age_cat) #used to group lines on either side of the threshold in plot

# 9) Add young/old labels to the season-specific lines for correct grouping across the threshold
season_lines_1996_2025 <- season_lines_1996_2025 %>%
  mutate(age_cat = factor(ifelse(age < age_senesce, "Young", "Old"), levels = c("Young","Old")), #segment label
         seg = age_cat) #keep seg consistent

# 10) Make age plot: thin season curves + bootstrap ribbon + weighted mean line + observed avg points for each age + threshold line + sample size labels
plot_age_1996_2025 <- ggplot() +
  geom_line(data = season_lines_1996_2025, #thin season-specific curves for context
            aes(age, pred, group = interaction(season_fct, seg)),
            color = "grey40",
            alpha = 0.2) +
  geom_ribbon(data = ci_1996_2025, #bootstrap CI band for overall curve
              aes(age, ymin = lo, ymax = hi, fill = age_cat, group = seg), #separate ribbons on each side of threshold
              alpha = 0.3) +
  geom_line(data = subset(ci_1996_2025, age < age_senesce),
            aes(age, pred, color = age_cat, group = seg),
            linewidth = 1.5) +
  geom_line(data = subset(ci_1996_2025, age >= age_senesce),
            aes(age, pred, color = age_cat, group = seg),
            linewidth = 1.5,
            linetype = "dashed") +
  geom_pointrange(data = observed_data_1996_2025, #pooled observed points with binomial CI
                  aes(age, avg_prop, ymin = lwr, ymax = upr),
                  color = "black") +
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  annotate("text", 
           x = age_senesce - 0.5, y = 0.8, 
           label = "Age threshold (9 years)",
           angle = 90, vjust = -0.8, size = 3.5) +
  geom_text(data = observed_data_1996_2025, aes(age, 1.01, label = n_age), vjust = -0.5) + #sample size labels per age
  coord_cartesian(ylim = c(0.75, 1.02), clip = "off") +
  scale_color_manual(name = "Age class", values = AGECOL) +
  scale_fill_manual(name = "Age class", values = AGECOL) +
  scale_x_continuous(breaks = ages) +
  theme_few(base_size = 18) +
  labs(x = "Maternal age (years)", y = "Mother-offspring association"); plot_age_1996_2025

#### Figure 2b (age, 2016-2023) ####

# 1) Define x-axis ages used in model
ages_2016_2023 <- sort(unique(model_variables_2016_2023$age)) #unique ages present (x values)

# 2) Compute pooled observed MOA_proportion by age (black points + binomial CIs)
observed_data_age_2016_2023 <- model_variables_2016_2023 %>%
  group_by(age) %>% #group by rows within each age
  summarise(n_age = n(), #sample size at each age for labels
            n_success = sum(count_1_pup, na.rm = TRUE), #total successes
            n_trials = sum(total_resights, na.rm = TRUE), #total effort
            avg_prop = n_success/n_trials, #observed MOA_proportion
            lwr = binom.test(n_success, n_trials)$conf.int[1], #binomial CI lower
            upr = binom.test(n_success, n_trials)$conf.int[2], #binomial CI upper
            .groups = "drop")

# 3) Build standardized dataframe for each age for post-stratification
nd_by_age_2016_2023 <- lapply(ages_2016_2023, \(A) transform(model_variables_2016_2023, #retain covariate distribution
                                                             age = A)) #set all seals to age A

# 4) Post-stratified prediction curve (row predictions -> weighted mean)
post_curve_age_2016_2023 <- function(fit, season = NULL){
  vapply(seq_along(nd_by_age_2016_2023), function(a){
    nd <- nd_by_age_2016_2023[[a]] #dataset for one age
    if(!is.null(season)) nd$season_fct <- factor(season, levels = seasons_2016_2023) #force season level
    p <- predict(fit, newdata = nd, type = "response", re.form = NULL) #predict including all random effects
    weighted.mean(p, w = nd$total_resights, na.rm = TRUE) #total resight weighted mean
  }, numeric(1))
}

# 5) Thin season curves (one post-stratified curve per season)
season_lines_age_2016_2023 <- bind_rows(lapply(seasons_2016_2023, \(S)
                                               tibble(season_fct = S, #season grouping
                                                      age = ages, #x-axis
                                                      pred = post_curve_age_2016_2023(mod_age_2016_2023, season = S)))) #season-specific predictions

# 6) Overall curve + parametric bootstrap CI
post_curve_age_overall_2016_2023 <- function(fit) post_curve_age_2016_2023(fit, season = NULL) #overall curve

overall_age_2016_2023 <- tibble(age = ages_2016_2023,
                                pred = post_curve_age_overall_2016_2023(mod_age_2016_2023)) %>% #mean prediction
  left_join(observed_data_age_2016_2023 %>% select(age, n_age), by = "age") #add sample sizes

boot_age_2016_2023 <- bootMer(mod_age_2016_2023, FUN = post_curve_age_overall_2016_2023,
                              nsim = nsim, type = "parametric", use.u = FALSE) #parametric bootstrap

# 7) Convert bootstrap draws into 95% CI ribbon
ci_age_2016_2023 <- overall_age_2016_2023 %>%
  mutate(lo = apply(boot_age_2016_2023$t, 2, quantile, 0.025, na.rm = TRUE), #lower CI
         hi = apply(boot_age_2016_2023$t, 2, quantile, 0.975, na.rm = TRUE)) #upper CI

# 8) Plot: thin season curves + bootstrap ribbon + mean curve + observed data
plot_age_2016_2023 <- ggplot() +
  geom_line(data = season_lines_age_2016_2023, #thin season-specific curves
            aes(age, pred, group = season_fct),
            color = "grey40",
            alpha = 0.4) +
  geom_ribbon(data = ci_age_2016_2023, #bootstrap CI band
              aes(age, ymin = lo, ymax = hi),
              fill = "#D295E3", alpha = 0.28) +
  geom_line(data = ci_age_2016_2023, #overall predicted curve
            aes(age, pred),
            color = "#D295E3", linewidth = 1.5) +
  geom_pointrange(data = observed_data_age_2016_2023, #observed pooled MOA_proportions
                  aes(age, avg_prop, ymin = lwr, ymax = upr),
                  color = "black") +
  geom_text(data = observed_data_age_2016_2023, #sample size labels
            aes(age, 1.01, label = n_age),
            vjust = -0.5) +
  scale_x_continuous(breaks = ages) +
  coord_cartesian(ylim = c(0.75, 1.02), clip = "off") +
  theme_few(base_size = 18) +
  labs(x = "Maternal age (years)", y = "Mother-offspring association"); plot_age_2016_2023

#### Figure 2d (breakpoint pupping experience, 1996-2025)  ####

# 1) Colors for inexperienced/experienced
EXPCOL <- c(Inexperienced = "#9FD46C", Experienced = "#7C82F1")

# 2) Define x-axis experience and seasons
exp_vals <- sort(unique(model_variables$pupping_exp))
seasons <- levels(model_variables$season_fct)

# 3) Observed MOA_proportions by experience
observed_exp_data <- model_variables %>%
  group_by(pupping_exp) %>% 
  summarise(n_exp = n(),
            n_success = sum(count_1_pup, na.rm = TRUE),
            n_trials = sum(total_resights, na.rm = TRUE),
            avg_prop = n_success/n_trials,
            lwr = binom.test(n_success, n_trials)$conf.int[1],
            upr = binom.test(n_success, n_trials)$conf.int[2],
            .groups = "drop")

# 4) Standardized newdata for each experience E
nd_by_exp <- lapply(exp_vals, \(E) transform(model_variables,
                                             pupping_exp = E,
                                             experience_cat = factor(ifelse(E < 5, "Inexperienced", "Experienced"), levels = c("Inexperienced","Experienced")),
                                             exp10 = (E - 5)/10))

# 5) Post-stratified curve
post_curve_exp <- function(fit, season = NULL) {
  vapply(seq_along(nd_by_exp), function(e){
    nd <- nd_by_exp[[e]]
    if(!is.null(season)) nd$season_fct <- factor(season, levels = seasons)
    p <- predict(fit, newdata = nd, type = "response", re.form = NULL)
    weighted.mean(p, w = nd$total_resights, na.rm = TRUE)
  }, numeric(1))
}

# 6) Thin season lines
season_lines_exp <- bind_rows(lapply(seasons, \(S)
                                     tibble(season_fct = S,
                                            pupping_exp = exp_vals,
                                            pred = post_curve_exp(mod_exp_1996_2025, season = S))))

# 7) Bootstrap CI
post_curve_exp_overall <- function(fit) post_curve_exp(fit, season = NULL)

overall_exp <- tibble(pupping_exp = exp_vals,
                      pred = post_curve_exp_overall(mod_exp_1996_2025)) %>%
  left_join(observed_exp_data %>% select(pupping_exp, n_exp), by = "pupping_exp")

boot_exp <- bootMer(mod_exp_1996_2025, FUN = post_curve_exp_overall, nsim = nsim,
                    type = "parametric", use.u = FALSE)

# 8) CIs
ci_exp <- overall_exp %>%
  mutate(lo = apply(boot_exp$t, 2, quantile, 0.025, na.rm = TRUE),
         hi = apply(boot_exp$t, 2, quantile, 0.975, na.rm = TRUE),
         experience_cat = factor(ifelse(pupping_exp < 5, "Inexperienced", "Experienced"), levels = c("Inexperienced","Experienced")),
         seg = experience_cat)

# 9) Label season lines
season_lines_exp <- season_lines_exp %>%
  mutate(experience_cat = factor(ifelse(pupping_exp < 5, "Inexperienced", "Experienced"), levels = c("Inexperienced","Experienced")),
         seg = experience_cat)

# 10) Plot
plot_exp_1996_2025 <- ggplot() +
  geom_line(data = season_lines_exp,
            aes(pupping_exp, pred, group = interaction(season_fct, seg)),
            color = "grey40",
            alpha = 0.2) +
  geom_ribbon(data = ci_exp,
              aes(pupping_exp, ymin = lo, ymax = hi, fill = experience_cat, group = seg),
              alpha = 0.4) +
  geom_line(data = subset(ci_exp, pupping_exp < 5),
            aes(pupping_exp, pred, color = experience_cat, group = seg),
            linewidth = 1.5) +
  geom_line(data = subset(ci_exp, pupping_exp >= 5),
            aes(pupping_exp, pred, color = experience_cat, group = seg),
            linewidth = 1.5) +
  geom_pointrange(data = observed_exp_data,
                  aes(pupping_exp, avg_prop, ymin = lwr, ymax = upr),
                  color = "black") +
  geom_vline(xintercept = 5 - 0.5, linetype = "dashed") +
  annotate("text",
           x = 5 - 0.5,
           y = 0.77, label = "Experience threshold (5 pups)",
           angle = 90, vjust = -0.8, size = 3.5) +
  coord_cartesian(ylim = c(0.7, 1.02), clip = "off") +
  geom_text(data = observed_exp_data,
            aes(pupping_exp, 1.01, label = n_exp), vjust = -0.5) +
  scale_color_manual(name = "Experience class", values = EXPCOL) +
  scale_fill_manual(name = "Experience class", values = EXPCOL) +
  scale_x_continuous(breaks = exp_vals) +
  theme_few(base_size = 18) +
  labs(x = "Previous pupping experience (number of pups)", y = "Mother-offspring association"); plot_exp_1996_2025


#### Figure 2d (pupping experience, 2016-2023) ####

# 1) Define x-axis experience values and seasons used for season-specific curves
exp_vals_2016_2023 <- sort(unique(model_variables_2016_2023$pupping_exp)) #unique experience values present (x values)
seasons_2016_2023 <- levels(model_variables_2016_2023$season_fct) #season factor levels used in predictions

# 2) Compute pooled observed MOA_proportion by experience (black points + binomial CIs)
observed_data_exp_2016_2023 <- model_variables_2016_2023 %>%
  group_by(pupping_exp) %>% #group by rows within each experience level
  summarise(n_experience = n(), #sample size at each experience value for labels
            n_success = sum(count_1_pup, na.rm = TRUE), #total successes
            n_trials = sum(total_resights, na.rm = TRUE), #total effort
            avg_prop = n_success/n_trials, #observed MOA_proportion
            lwr = binom.test(n_success, n_trials)$conf.int[1], #binomial CI lower
            upr = binom.test(n_success, n_trials)$conf.int[2], #binomial CI upper
            .groups = "drop")

# 3) Build standardized newdata for each experience value for post-stratification
nd_by_experience_2016_2023 <- lapply(exp_vals_2016_2023, \(E) transform(model_variables_2016_2023, #retain covariate distribution
                                                           pupping_exp = E)) #set all seals to experience value E

# 4) Post-stratified prediction curve (row predictions -> weighted mean)
post_curve_exp_2016_2023 <- function(fit, season = NULL){
  vapply(seq_along(nd_by_experience_2016_2023), function(e){
    nd <- nd_by_experience_2016_2023[[e]] #dataset for one experience value
    if(!is.null(season)) nd$season_fct <- factor(season, levels = seasons_2016_2023) #force season level
    p <- predict(fit, newdata = nd, type = "response", re.form = NULL) #predict including all RE
    weighted.mean(p, w = nd$total_resights, na.rm = TRUE) #total resight weighted mean
  }, numeric(1))
}

# 5) Thin season curves (one post-stratified curve per season)
season_lines_exp_2016_2023 <- bind_rows(lapply(seasons_2016_2023, \(S)
                                           tibble(season_fct = S, #season grouping
                                                  pupping_exp = exp_vals_2016_2023, #x-axis
                                                  pred = post_curve_exp_2016_2023(mod_exp_2016_2023, season = S)))) #season-specific predictions

# 6) Overall curve + parametric bootstrap CI
post_curve_exp_overall_2016_2023 <- function(fit) post_curve_exp_2016_2023(fit, season = NULL) #overall curve

overall_exp_2016_2023 <- tibble(pupping_exp = exp_vals_2016_2023,
                            pred = post_curve_exp_overall_2016_2023(mod_exp_2016_2023)) %>% #mean prediction
  left_join(observed_data_exp_2016_2023 %>% select(pupping_exp, n_experience), by = "pupping_exp") #add sample sizes

boot_exp_2016_2023 <- bootMer(mod_exp_2016_2023, FUN = post_curve_exp_overall_2016_2023,
                          nsim = nsim, type = "parametric", use.u = FALSE) #parametric bootstrap

# 7) Convert bootstrap draws into 95% CI ribbon
ci_exp_2016_2023 <- overall_exp_2016_2023 %>%
  mutate(lo = apply(boot_exp_2016_2023$t, 2, quantile, 0.025, na.rm = TRUE), #lower CI
         hi = apply(boot_exp_2016_2023$t, 2, quantile, 0.975, na.rm = TRUE)) #upper CI

# 8) Plot: thin season curves + bootstrap ribbon + mean curve + observed data
plot_exp_2016_2023 <- ggplot() +
  geom_line(data = season_lines_exp_2016_2023, #thin season-specific curves
            aes(pupping_exp, pred, group = season_fct),
            color = "grey40",
            alpha = 0.3) +
  geom_ribbon(data = ci_exp_2016_2023, #bootstrap CI band
              aes(pupping_exp, ymin = lo, ymax = hi),
              fill = "#2BB295", alpha = 0.28) +
  geom_line(data = ci_exp_2016_2023, #overall predicted curve
            aes(pupping_exp, pred),
            color = "#2BB295", linewidth = 1.5) +
  geom_pointrange(data = observed_data_exp_2016_2023, #observed pooled MOA_proportions
                  aes(pupping_exp, avg_prop, ymin = lwr, ymax = upr),
                  color = "black") +
  geom_text(data = observed_data_exp_2016_2023, #sample size labels
            aes(pupping_exp, 1.01, label = n_experience),
            vjust = -0.5) +
  scale_x_continuous(breaks = exp_vals_2016_2023) +
  coord_cartesian(ylim = c(0.65, 1.02), clip = "off") +
  theme_few(base_size = 18) +
  labs(x = "Previous pupping experience (number of pups)", y = "Mother-offspring association"); plot_exp_2016_2023

### Facet Figure 2 Together ###

# 1) Left column: age plots
age_plots <- plot_grid(plot_age_1996_2025 +
                      theme(axis.title.x = element_blank(),
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank()),
                    
                    plot_age_2016_2023,
                    ncol = 1,
                    labels = c("(a)", "(b)"),
                    label_size = 20,
                    label_fontface = "bold",
                    label_x = 0.02,
                    label_y = 0.98,
                    hjust = 0,
                    vjust = 1,
                    align = "v") #stack age figures

# 2) Right column: experience plots
exp_plots <- plot_grid(plot_exp_1996_2025 +
                      theme(axis.title.x = element_blank(),
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            axis.title.y = element_blank()), #remove right-column y-axis title
                    
                    plot_exp_2016_2023 +
                      theme(axis.title.y = element_blank()), #remove right-column y-axis title
                    
                    ncol = 1,
                    labels = c("(c)", "(d)"),
                    label_size = 20,
                    label_fontface = "bold",
                    label_x = -0.01,
                    label_y = 0.98,
                    hjust = 0,
                    vjust = 1,
                    align = "v") #stack experience figures

# 3) Combine columns
plot_age_exp <- plot_grid(age_plots, exp_plots,
                          ncol = 2,
                          align = "hv",
                          axis = "tblr",
                          rel_widths = c(1, 1)) #equal column widths
plot_age_exp

# 4) Save figure
ggsave(here("TablesFigures", "Figure2.png"), plot_age_exp, width = 17, height = 10, dpi = 400)

#### Figure 3a (density conceptual map) ####

# 1) Read beaches geopackage
beaches <- st_read(here("IntermediateData", "beaches.gpkg"), quiet = TRUE)

# 2) Read density data
seal_density <- read_csv(here("IntermediateData", "seal_density.csv"), show_col_types = FALSE)

# 3) Mean density per beach across years
area_density_mean <- seal_density %>%
  rename(dominant_area = Beach) %>%
  group_by(dominant_area) %>%
  summarize(mean_density = mean(density, na.rm = TRUE), .groups = "drop")

# 1) Replace "Beach" with dominant_area
beaches_density <- beaches %>%
  left_join(area_density_mean, by = c("Beach" = "dominant_area"))

# 2) Get extent for zooming and tile download
bbox <- st_bbox(beaches_density)

# 3) Convert bbox to sf object for maptiles
bbox_sf <- st_as_sfc(bbox) %>%
  st_as_sf()

# 4) expand bbox slightly before downloading tiles (helps with map fit)
buffer_x <- 0.01
buffer_y <- 0.0005

bbox_expanded <- bbox
bbox_expanded["xmin"] <- bbox["xmin"] - buffer_x
bbox_expanded["xmax"] <- bbox["xmax"] + buffer_x
bbox_expanded["ymin"] <- bbox["ymin"] - buffer_y
bbox_expanded["ymax"] <- bbox["ymax"] + buffer_y

bbox_sf_expanded <- st_as_sfc(bbox_expanded) %>%
  st_as_sf()

# 5) download tiles using expanded bbox
sat_map <- get_tiles(x = bbox_sf_expanded,
                     provider = "Esri.WorldImagery",
                     crop = TRUE,
                     zoom = 16)

# manually set the color limits so the same ones correspond between Figs 2a and 2b
density_limits <- range(c(beaches_density$mean_density,
                          observed_density_2016_2023$avg_density),
                        na.rm = TRUE) 

# 6) make seal density heatmap
seal_density_map <- ggplot() +
  geom_spatraster_rgb(data = sat_map) +
  geom_sf(data = beaches_density,
          aes(fill = mean_density),
          color = "white",
          linewidth = 0.4,
          alpha = 0.8) +
  scale_fill_viridis_c(option = "mako",
                       direction = -1,
                       limits = density_limits,
                       name = "Conspecific density\n(seals in a 10m radius)") +
  annotation_scale(location = "bl", width_hint = 0.25, text_face = "bold", text_cex = 1.5, text_col = "white") +
  coord_sf(xlim = c(bbox_expanded["xmin"], bbox_expanded["xmax"]),
           ylim = c(bbox_expanded["ymin"], bbox_expanded["ymax"]),
           expand = FALSE) +
  theme_void(base_size = 22) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()); seal_density_map

ggsave(here("TablesFigures", "Figure3a.png"), plot = seal_density_map, width = 18, height = 14, dpi = 300, bg = "transparent")

#### Figure 3b (conspecific density effect) ####

# 1) Define density bins used for observed points
model_variables_2016_2023$density_bin <- cut_number(model_variables_2016_2023$avg_density,
                                              n = 8, labels = FALSE)

# 2) Compute pooled observed MOA_proportion by density bin (black points + binomial CIs)
observed_density_2016_2023 <- model_variables_2016_2023 %>%
  group_by(density_bin) %>%
  summarise(avg_density = mean(avg_density, na.rm = TRUE),
            n_bin = n(),
            n_success = sum(count_1_pup, na.rm = TRUE),
            n_trials = sum(total_resights, na.rm = TRUE),
            avg_prop = n_success / n_trials,
            lwr = binom.test(n_success, n_trials)$conf.int[1],
            upr = binom.test(n_success, n_trials)$conf.int[2],
            .groups = "drop")

# 1) Define x-axis density values across the observed range
density_vals <- seq(min(model_variables_2016_2023$avg_density, na.rm = TRUE), #minimum observed density
                    max(model_variables_2016_2023$avg_density, na.rm = TRUE), #maximum observed density
                    length.out = 100) #number of x-axis values for smooth curve

# 2) Build standardized newdata for each density value D (for post-stratification)
# For each D, keep the full predictors of model_variables_2016_2023 but overwrite avg_density as if everyone had density D
nd_by_density <- lapply(density_vals, \(D) transform(model_variables_2016_2023, #preserve covariate distribution and weights
                                                     avg_density = D)) #overwrite avg_density so every row is evaluated at density D

# 3) Post-stratified curve: predict for each density, then compute weighted average using total_resights
post_curve_density <- function(fit) {
  vapply(seq_along(nd_by_density), function(d){ #loop over density values
    nd <- nd_by_density[[d]] #newdata for one density value
    p <- predict(fit, newdata = nd, type = "response", re.form = NULL) #predict probability including all random effects
    weighted.mean(p, w = nd$total_resights, na.rm = TRUE) #post-stratified mean, weighted by total resights
  }, numeric(1)) #numeric vector of predictions across density values
}

# 4) Bootstrap CI for the overall density curve using parametric bootstrapping of the fitted GLMM
overall_density_2016_2023 <- tibble(avg_density = density_vals, #continuous density values
                                    pred = post_curve_density(mod_binom_2016_2023)) #overall post-stratified predictions

boot_density_2016_2023 <- bootMer(mod_binom_2016_2023, FUN = post_curve_density, nsim = nsim, #bootstrap curves
                                  type = "parametric", use.u = FALSE) #parametric bootstrap = simulate ranef each time

# 5) Turn bootstrap curves into 95% CI bands
ci_density_2016_2023 <- overall_density_2016_2023 %>%
  mutate(lo = apply(boot_density_2016_2023$t, 2, quantile, 0.025, na.rm = TRUE), #2.5% quantile at each density value
         hi = apply(boot_density_2016_2023$t, 2, quantile, 0.975, na.rm = TRUE)) #97.5% quantile at each density value

# 8) Make density plot: bootstrap ribbon + weighted mean line + binned observed data by season
plot_density <- ggplot() +
  geom_ribbon(data = ci_density_2016_2023, #bootstrap CI band for overall curve
              aes(avg_density, ymin = lo, ymax = hi),
              fill = "black",
              alpha = 0.2) +
  geom_line(data = ci_density_2016_2023, #overall post-stratified prediction line
            aes(avg_density, pred),
            color = "black",
            linewidth = 1) +
  geom_pointrange(data = observed_density_2016_2023, #observed binned MOA_proportions by season
                  aes(avg_density, avg_prop, ymin = lwr, ymax = upr, color = avg_density),
                  size = 1.2,
                  linewidth = 1.2) +
  geom_text(data = observed_density_2016_2023, #sample size labels
            aes(avg_density, 1, label = n_bin),
            size = 5,
            vjust = -0.5) +
  scale_color_viridis_c(option = "mako",
                        begin = 0.1,
                        end = 0.8,
                        direction = -1,
                        limits = density_limits, #use same color scale as map
                        name = "Density") +
  coord_cartesian(ylim = c(0.8, 1), clip = "off") + 
  scale_x_continuous(n.breaks = 6) +
  scale_y_continuous(n.breaks = 3) +
  theme_few(base_size = 28) +
  labs(x = "Conspecific density (seals in a 10m radius)",
       y = "Mother-offspring association"); plot_density

ggsave(here("TablesFigures", "Figure3b.png"), plot_density, width = 11, height = 8, dpi = 800)

########### Figure 3c (extreme wave and tide effect) ############

# Set plotting colors for seasons
pal <- c("#E57373","#E6A64C","#E5C76B","#7FBF7B","#2E7D32","#5C7EE5","#6DAEDB","#D77CC8")

# 1) Define x-axis extreme event values and seasons used for plotting
extreme_vals <- sort(unique(model_variables_2016_2023$n_extreme_both)) #unique extreme event counts present in data

# 2) Compute observed MOA_proportion by extreme event count and season (colored points + binomial CI whiskers)
observed_extreme_2016_2023 <- model_variables_2016_2023 %>%
  group_by(n_extreme_both, season_fct) %>%
  summarise(n_obs = n(), #sample size at each extreme event count within season
            n_success = sum(count_1_pup, na.rm = TRUE), #total successes at each extreme event count within season
            n_trials = sum(total_resights, na.rm = TRUE), #total resights at each extreme event count within season
            avg_prop = n_success / n_trials, #pooled observed MOA_proportion at each extreme event count within season
            lwr = binom.test(n_success, n_trials)$conf.int[1], #CI lower
            upr = binom.test(n_success, n_trials)$conf.int[2], #CI upper
            .groups = "drop") 

# make sure years are in order
observed_extreme_2016_2023 <- observed_extreme_2016_2023 %>%
  mutate(season_fct = factor(season_fct, levels = c("2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")))

# 3) Build standardized newdata for each extreme event count E (for post-stratification)
# For each E, keep the full covariate mix of model_variables_2016_2023 but overwrite n_extreme_both as if everyone experienced E events
nd_by_extreme <- lapply(extreme_vals, \(X) transform(model_variables_2016_2023, #start from full dataset to preserve covariate distribution and weights
                                                     n_extreme_both = X)) #overwrite n_extreme_both so every row is evaluated at event count E

# 4) Post-stratified curve: predict for each extreme event count, then compute weighted average using total_resights
post_curve_extreme <- function(fit) {
  vapply(seq_along(nd_by_extreme), function(x){ #loop over extreme event counts
    nd <- nd_by_extreme[[x]] #newdata for one extreme event count
    p <- predict(fit, newdata = nd, type = "response", re.form = NULL) #predict probability including all random effects
    weighted.mean(p, w = nd$total_resights, na.rm = TRUE) #post-stratified mean, weighted by total resights
  }, numeric(1)) #numeric vector of predictions across extreme-event counts
}

# 5) Bootstrap CI for the overall extreme event curve using parametric bootstrapping of the fitted GLMM
overall_extreme_2016_2023 <- tibble(n_extreme_both = extreme_vals, #discrete extreme event counts
                                    pred = post_curve_extreme(mod_age_2016_2023)) #overall post-stratified predictions

boot_extreme_2016_2023 <- bootMer(mod_age_2016_2023, FUN = post_curve_extreme, nsim = nsim, #bootstrap curves
                                  type = "parametric", use.u = FALSE) #parametric bootstrap = simulate ranefs each time

# 6) Turn bootstrap curves into 95% CI bands
ci_extreme_2016_2023 <- overall_extreme_2016_2023 %>%
  mutate(lo = apply(boot_extreme_2016_2023$t, 2, quantile, 0.025, na.rm = TRUE), #2.5% quantile at each extreme event count
         hi = apply(boot_extreme_2016_2023$t, 2, quantile, 0.975, na.rm = TRUE)) #97.5% quantile at each extreme event count

model_variables_2016_2023$season_fct <- factor(model_variables_2016_2023$season_fct,
                                         levels = c("2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023"))

## stagger CIs so they don't overlap
observed_extreme_2016_2023 <- observed_extreme_2016_2023 %>%
  arrange(n_extreme_both) %>%
  mutate(close_group = cumsum(c(TRUE, diff(n_extreme_both) > 1))) %>% #find nearby x values
  group_by(n_extreme_both) %>% #dodge exact duplicates only for CIs and points
  arrange(season_fct, .by_group = TRUE) %>%
  mutate(lab_id = row_number(), 
         x_plot = if (n() == 1) {
           n_extreme_both
           } else {
             n_extreme_both + (lab_id - mean(lab_id)) * 0.9 
             }) %>%
  ungroup() %>%
  group_by(close_group) %>% #labels: spread nearby values only
  arrange(n_extreme_both, .by_group = TRUE) %>%
  mutate(
    lab_close_id = row_number(),
    x_lab = if (n() == 1) {
      x_plot
    } else {
      mean(n_extreme_both) + (lab_close_id - mean(lab_close_id)) * 2
    },
    y_lab = 1
  ) %>%
  ungroup()

plot_n_extreme <- ggplot() +
  geom_ribbon(data = ci_extreme_2016_2023, 
              aes(n_extreme_both, ymin = lo, ymax = hi),
              fill = "black",
              alpha = 0.2) +
  geom_line(data = ci_extreme_2016_2023,
            aes(n_extreme_both, pred),
            color = "black",
            linewidth = 1.2) +
  geom_pointrange(data = observed_extreme_2016_2023,
                  aes(x_plot, avg_prop,
                      ymin = lwr, ymax = upr,
                      color = season_fct),
                  size = 1.1,
                  linewidth = 1.1) +
  geom_text(data = observed_extreme_2016_2023,
            aes(x_lab, y_lab, label = n_obs),
            color = "black",
            size = 5,
            vjust = 0,
            show.legend = FALSE) +
  scale_color_manual(values = pal, name = "Year") +
  coord_cartesian(ylim = c(0.8, 1), clip = "off") +
  scale_x_continuous(n.breaks = 6) +
  scale_y_continuous(n.breaks = 3) +
  theme_few(base_size = 28) +
  labs(x = "Number of per-year extreme wave and tide events", y = "Mother-offspring association"); plot_n_extreme

ggsave(here("TablesFigures", "Figure3c.png"), plot_n_extreme, width = 11, height = 8, dpi = 800)

#### Figure 3b (with raw data points) ####

# Density plot: bootstrap ribbon + weighted mean line + raw points weighted by total resights
plot_density_jitter <- ggplot() +
  geom_jitter(data = model_variables_2016_2023,
              aes(avg_density,
                  MOA_proportion,
                  size = total_resights,
                  color = avg_density),
              alpha = 0.3,
              width = 0.15,
              height = 0.02) +
  geom_ribbon(data = ci_density_2016_2023,
              aes(avg_density, ymin = lo, ymax = hi),
              fill = "black",
              alpha = 0.3) +
  geom_line(data = ci_density_2016_2023,
            aes(avg_density, pred),
            color = "black",
            linewidth = 1.3) +
  scale_color_viridis_c(option = "mako",
                        direction = -1,
                        limits = density_limits,
                        name = "Conspecific density") +
  scale_size_continuous(name = "Observations") +
  coord_cartesian(ylim = c(0, 1.01), clip = "off") +
  theme_few(base_size = 18) +
  labs(x = "Conspecific density (seals in a 10m radius)",
       y = "Mother-offspring association"); plot_density_jitter

ggsave(here("TablesFigures", "Figure3b_raw.png"), plot_density_jitter, width = 10, height = 8, dpi = 800)

#### Figure 3c (with raw data points) ####

# 9) Make extreme events plot with jittered raw points
plot_n_extreme_jitter <- ggplot() +
  geom_jitter(data = model_variables_2016_2023, #raw observations weighted by resights
              aes(n_extreme_both, MOA_proportion, size = total_resights, color = season_fct),
              alpha = 0.35,
              width = 0.15,
              height = 0.02) +
  geom_ribbon(data = ci_extreme_2016_2023, #bootstrap CI band for overall curve
              aes(n_extreme_both, ymin = lo, ymax = hi),
              fill = "black",
              alpha = 0.3) +
  geom_line(data = ci_extreme_2016_2023, #overall post-stratified prediction line
            aes(n_extreme_both, pred),
            color = "black",
            linewidth = 1.3) +
  scale_color_manual(values = pal, name = "Year") +
  scale_size_continuous(name = "Observations") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 5) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  theme_few(base_size = 18) +
  labs(x = "Number of per-year extreme wave and tide events",
       y = "Mother-offspring association"); plot_n_extreme_jitter

ggsave(here("TablesFigures", "Figure3c_raw.png"), plot_n_extreme_jitter, width = 6, height = 8, dpi = 800)

#### Figure 4 (offspring weaning mass vs. MOA) ####

# 1) Predicted effect of mother-offspring association on wean mass
pred_wean_mass <- ggpredict(mod_wean_mass,
                            terms = "MOA_proportion")

# 2) Plot model predictions
plot_wean <- ggplot() +
  geom_jitter(data = model_variables,
              aes(x = (MOA_proportion), y = Wt_wean_corrected),
              size = 1.5,
              color = "#664388",
              alpha = 0.45,
              width = 0.003,
              height = 0) +
  geom_ribbon(data = pred_wean_mass,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = "black",
              alpha = 0.18) +
  geom_line(data = pred_wean_mass,
            aes(x = x, y = predicted),
            color = "black",
            linewidth = 1.3) +
  coord_cartesian(xlim = range(model_variables_wean$MOA_proportion, na.rm = TRUE)) +
  labs(x = "Mother-offspring association",
       y = "Offspring weaning mass (kg)") +
  theme_classic(base_size = 20); plot_wean

ggsave(here("TablesFigures", "Figure4.png"), plot_wean, width = 10, height = 8, dpi = 800)

# SUPPLEMENTARY TABLES AND FIGURES -----------------------------------

#### Tables S1-S5 (model outputs) ####

# 1) model output table function
make_mod_flextable <- function(model, #extract fixed effects from model
                               predictor_labels = NULL, #relabel predictors if we relabel
                               digits = 3, #
                               se_digits = 2,
                               save_path = NULL) {
  
  fixed_tbl <- tidy(model, effects = "fixed") %>% #fixed effects
    mutate(Predictor = if (!is.null(predictor_labels)) {
      recode(term, !!!predictor_labels, .default = term)
      } else {
        term #otherwise allows us to keep original terms
      },
      Estimate = as.character(round(estimate, digits)), #round estimates
      SE = as.character(round(std.error, se_digits)), #round SEs
      Z = as.character(round(statistic, digits)), #round Z value
      "P-value" = case_when(p.value < 0.001 ~ paste0(formatC(p.value, format = "e", digits = 1), " ***"), #format p-values and significance labels
                            p.value < 0.01 ~ paste0(round(p.value, 4), " **"),
                            p.value < 0.05 ~ paste0(round(p.value, 4), " *"),
                            TRUE ~ as.character(round(p.value, 4)))) %>%
    select(Predictor, Estimate, SE, Z, "P-value")
  
  random_tbl <- tidy(model, effects = "ran_pars") %>% #random effects
    transmute(Predictor = paste0("Random effect: ", 
                                 recode(group,
                                        animalID_fct = "AnimalID", 
                                        season_fct = "Year",
                                        Residual = "Residual")),
              Estimate = as.character(round(estimate^2, digits)), #rounded ranef estimate
              SE = as.character(round(estimate, digits)), #rounded ranef SE
              Z = "",
              "P-value" = "")
  
  out_tbl <- bind_rows(fixed_tbl, random_tbl) #keep both fixed and ranef

  #save format as flextable
  
  ft <- out_tbl %>%
    flextable() %>%
    align(align = "center", part = "all") %>%
    autofit()
  
  if (!is.null(save_path)) {
    save_as_docx(ft, path = save_path)
  }
  
  return(ft)
}

### Table S1: Model output from the 2016-2023 model with maternal age. ###
make_mod_flextable(mod_age_2016_2023,
                   predictor_labels = c("(Intercept)" = "Intercept",
                                        "age" = "Maternal age",
                                        "n_extreme_both" = "Number of extreme wave and tide events",
                                        "avg_density" = "Conspecific density"),
                   save_path = here("TablesFigures", "mod_age_2016_2023_output.docx"))

### Table S2. Model output from the 2016-2023 model with previous pupping experience ###
make_mod_flextable(mod_exp_2016_2023,
                   predictor_labels = c("(Intercept)" = "Intercept",
                                        "pupping_exp" = "Previous pupping experience",
                                        "n_extreme_both" = "Number of extreme wave and tide events",
                                        "avg_density" = "Conspecific density"),
                   save_path = here("TablesFigures", "mod_exp_2016_2023_output.docx"))

### Table S3. Model output from the 1996-2025 piecewise segmented regression for maternal age. ###
make_mod_flextable(mod_age_1996_2025,
                   predictor_labels = c("(Intercept)" = "Intercept",
                                        "age_catYoung:age10" = "Maternal age : Pre-threshold (< 9 years)",
                                        "age_catOld:age10" = "Maternal age : Post-threshold (≥ 9 years)"),
                   save_path = here("TablesFigures", "mod_age_1996_2025_output.docx"))

### Table S4. Model output from the 1996-2025 piecewise segmented regression for previous pupping experience. ###
make_mod_flextable(mod_exp_1996_2025,
                   predictor_labels = c("(Intercept)" = "Intercept",
                                        "experience_catInexperienced:exp10" = "Previous pupping experience : Pre-threshold (< 5 previous pups)",
                                        "experience_catExperienced:exp10" = "Previous pupping experience : Post-threshold (≥ 5 previous pups)"),
                   save_path = here("TablesFigures", "mod_exp_1996_2025_output.docx"))

### Table S5. Model output for the effects of maternal age and mother-offspring association on offspring weaning mass. ###
make_mod_flextable(mod_wean_mass,
                   predictor_labels = c("(Intercept)" = "Intercept",
                                        "MOA_proportion" = "Mother-offspring association",
                                        "age" = "Maternal age"),
                   save_path = here("TablesFigures", "mod_wean_mass_output.docx"))

#### Table S6 (age linear vs. quadratic vs. breakpoint, 1996-2025) ####

### threshold ###

mod_age_thresh_1996_2025 <- glmer(MOA_proportion ~ age_cat : age10 + (1 | animalID_fct) + (1 | season_fct),
                                  weights = total_resights,
                                  family = binomial(link= "logit"),
                                  control = glmerControl(optimizer = "bobyqa"),
                                  data = model_variables); summary(mod_age_thresh_1996_2025) #model summary

### linear ###

mod_age_linear_1996_2025 <- glmer(MOA_proportion ~ age + (1 | animalID_fct) + (1 | season_fct),
                                  weights = total_resights,
                                  family = binomial(link= "logit"),
                                  control = glmerControl(optimizer = "bobyqa"),
                                  data = model_variables); summary(mod_age_linear_1996_2025) #model summary

### quadratic ###

mod_age_quad_1996_2025 <- glmer(MOA_proportion ~ age + I(age^2) + (1 | animalID_fct) + (1 | season_fct),
                                weights = total_resights,
                                family = binomial(link= "logit"),
                                control = glmerControl(optimizer = "bobyqa"),
                                data = model_variables); summary(mod_age_quad_1996_2025) #model summary

## name models to compare
mod_age_comparisons_1996_2025 <- list(Linear = mod_age_linear_1996_2025,
                                      Quadratic = mod_age_quad_1996_2025,
                                      Threshold = mod_age_thresh_1996_2025)

# 2) make table with AIC comparisons
aic_table_age_1996_2025 <- tibble(Model = names(mod_age_comparisons_1996_2025),
                                  AIC = sapply(mod_age_comparisons_1996_2025, AIC)) %>%
  mutate(delta_AIC = AIC - min(AIC),
         aic_weights = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))) %>%
  arrange(delta_AIC)

# 3) relabel table for formatting
aic_table_age_1996_2025 <- aic_table_age_1996_2025 %>%
  transmute(Model = Model,
            AIC = round(AIC, 1),
            "ΔAIC" = round(delta_AIC, 2),
            "AIC weight" = round(aic_weights, 3)) %>%
  arrange("ΔAIC")

best_model <- aic_table_age_1996_2025$Model[1] #best model is the one with lowest AIC difference

# 4) make flextable
aic_table_age_1996_2025 <- flextable(aic_table_age_1996_2025) %>%
  align(align = "center", part = "all") %>%
  bold(i = which(aic_table_age_1996_2025$Model == best_model), #bold best model
       part = "body"); aic_table_age_1996_2025

#save final table
save_as_docx(aic_table_age_1996_2025, path = here("TablesFigures", "2016_2023_Age_Predictor_Comparison.docx"))

#### Table S7 (age linear vs. quadratic vs. breakpoint, 2016-2023) ####

### threshold ###

mod_age_thresh_2016_2023 <- glmer(MOA_proportion ~ age_cat : age10 + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                                  weights = total_resights,
                                  family = binomial(link= "logit"),
                                  control = glmerControl(optimizer = "bobyqa"),
                                  data = model_variables); summary(mod_age_thresh_2016_2023) #model summary

### linear ###

mod_age_linear_2016_2023 <- glmer(MOA_proportion ~ age + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                                  weights = total_resights,
                                  family = binomial(link= "logit"),
                                  control = glmerControl(optimizer = "bobyqa"),
                                  data = model_variables); summary(mod_age_linear_2016_2023) #model summary

### quadratic ###

mod_age_quad_2016_2023 <- glmer(MOA_proportion ~ age + I(age^2) + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                                weights = total_resights,
                                family = binomial(link= "logit"),
                                control = glmerControl(optimizer = "bobyqa"),
                                data = model_variables); summary(mod_age_quad_2016_2023) #model summary

# 1) name models to compare
mod_age_comparisons_2016_2023 <- list(Linear = mod_age_linear_2016_2023,
                                      Quadratic = mod_age_quad_2016_2023,
                                      Threshold = mod_age_thresh_2016_2023)

# 2) make table with AIC comparisons
aic_table_age_2016_2023 <- tibble(Model = names(mod_age_comparisons_2016_2023),
                                  AIC = sapply(mod_age_comparisons_2016_2023, AIC)) %>%
  mutate(delta_AIC = AIC - min(AIC),
         aic_weights = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))) %>%
  arrange(delta_AIC)

# 3) relabel table for formatting
aic_table_age_2016_2023 <- aic_table_age_2016_2023 %>%
  transmute(Model = Model,
            AIC = round(AIC, 1),
            "ΔAIC" = round(delta_AIC, 2),
            "AIC weight" = round(aic_weights, 3)) %>%
  arrange("ΔAIC")

best_model <- aic_table_age_2016_2023$Model[1] #best model is the one with lowest AIC difference

# 4) make flextable
aic_table_age_2016_2023 <- flextable(aic_table_age_2016_2023) %>%
  align(align = "center", part = "all") %>%
  bold(i = which(aic_table_age_2016_2023$Model == best_model), #bold lowest AIC model
       part = "body"); aic_table_age_2016_2023

#save final table
save_as_docx(aic_table_age_2016_2023, path = here("TablesFigures", "2016_2023_Age_Predictor_Comparison.docx"))

#### Table S8 (experience linear vs. quadratic vs. breakpoint, 1996-2025) ####

### threshold ###

mod_exp_thresh_1996_2025 <- glmer(MOA_proportion ~ experience_cat : exp10 + (1 | animalID_fct) + (1 | season_fct),
                                  weights = total_resights,
                                  family = binomial(link = "logit"),
                                  control = glmerControl(optimizer = "bobyqa"),
                                  data = model_variables); summary(mod_exp_thresh_1996_2025)


### linear ###

mod_exp_linear_1996_2025 <- glmer(MOA_proportion ~ pupping_exp + (1 | animalID_fct) + (1 | season_fct),
                                  weights = total_resights,
                                  family = binomial(link = "logit"),
                                  control = glmerControl(optimizer = "bobyqa"),
                                  data = model_variables); summary(mod_exp_linear_1996_2025)

### quadratic ###

mod_exp_quad_1996_2025 <- glmer(MOA_proportion ~ exp10 + I(exp10^2) + (1 | animalID_fct) + (1 | season_fct),
                                  weights = total_resights,
                                  family = binomial(link = "logit"),
                                  control = glmerControl(optimizer = "bobyqa"),
                                  data = model_variables); summary(mod_exp_quad_1996_2025)

# 1) name models to compare
mod_exp_comparisons_1996_2025 <- list(Linear = mod_exp_linear_1996_2025,
                                      Quadratic = mod_exp_quad_1996_2025,
                                      Threshold = mod_exp_thresh_1996_2025)

# 2) make table with AIC comparisons
aic_table_exp_1996_2025 <- tibble(Model = names(mod_exp_comparisons_1996_2025),
                                  AIC = sapply(mod_exp_comparisons_1996_2025, AIC),
                                  K = sapply(mod_exp_comparisons_1996_2025, function(m) attr(logLik(m), "df"))) %>%
  mutate(delta_AIC = AIC - min(AIC),
         aic_weights = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))) %>%
  arrange(delta_AIC)

# 3) relabel table for formatting
aic_table_exp_1996_2025 <- aic_table_exp_1996_2025 %>%
  transmute(Model = Model, 
            AIC = round(AIC, 1),
            "ΔAIC" = round(delta_AIC, 2),
            "AIC weight" = round(aic_weights, 3)) %>%
  arrange("ΔAIC")

best_model <- aic_table_exp_1996_2025$Model[1] #best model is the one with lowest AIC difference

# 4) make flextable
aic_table_exp_1996_2025 <- flextable(aic_table_exp_1996_2025) %>%
  align(align = "center", part = "all") %>%
  bold(i = which(aic_table_exp_1996_2025$Model == best_model), #bold lowest AIC model
       part = "body"); aic_table_exp_1996_2025

#save final table
save_as_docx(aic_table_exp_1996_2025, path = here("TablesFigures", "1996_2025_Experience_Predictor_Comparison.docx"))

############ Table S9 (experience linear vs. quadratic vs. breakpoint, 2016-2023) ##################

### threshold ###

mod_exp_thresh_2016_2023 <- glmer(MOA_proportion ~ experience_cat : exp10 + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                                  weights = total_resights,
                                  family = binomial(link= "logit"),
                                  control = glmerControl(optimizer = "bobyqa"),
                                  data = model_variables); summary(mod_exp_thresh_2016_2023) #model summary

### linear ###

mod_exp_linear_2016_2023 <- glmer(MOA_proportion ~ pupping_exp + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                                  weights = total_resights,
                                  family = binomial(link = "logit"),
                                  control = glmerControl(optimizer = "bobyqa"),
                                  data = model_variables); summary(mod_exp_linear_2016_2023)

### quadratic ###

mod_exp_quad_2016_2023 <- glmer(MOA_proportion ~ exp10 + I(exp10^2) + avg_density + n_extreme_both + (1 | animalID_fct) + (1 | season_fct),
                                weights = total_resights,
                                family = binomial(link = "logit"),
                                control = glmerControl(optimizer = "bobyqa"),
                                data = model_variables); summary(mod_exp_quad_2016_2023)

# 1) name models to compare
mod_exp_comparisons_2016_2023 <- list(Linear = mod_exp_linear_2016_2023,
                                      Quadratic = mod_exp_quad_2016_2023,
                                      Threshold = mod_exp_thresh_2016_2023)

# 2) make table with AIC comparisons
aic_table_exp_2016_2023 <- tibble(Model = names(mod_exp_comparisons_2016_2023),
                                  AIC = sapply(mod_exp_comparisons_2016_2023, AIC)) %>%
  mutate(delta_AIC = AIC - min(AIC),
         aic_weights = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))) %>%
  arrange(delta_AIC)

# 3) relabel table for formatting
aic_table_exp_2016_2023 <- aic_table_exp_2016_2023 %>%
  transmute(Model = Model,
            AIC = round(AIC, 1),
            "ΔAIC" = round(delta_AIC, 2),
            "AIC weight" = round(aic_weights, 3)) %>%
  arrange("ΔAIC")

best_model <- aic_table_exp_2016_2023$Model[1] #best model is the one with lowest AIC difference

# 4) make flextable
aic_table_exp_2016_2023 <- flextable(aic_table_exp_2016_2023) %>%
  align(align = "center", part = "all") %>%
  bold(i = which(aic_table_exp_2016_2023$Model == best_model), #bold lowest AIC model
       part = "body"); aic_table_exp_2016_2023

#save final table
save_as_docx(aic_table_exp_2016_2023, path = here("TablesFigures", "2016_2023_Experience_Predictor_Comparison.docx"))

#### Table S10 (age threshold comparison) ####

# 1) test all possible thresholds
age_cutoff <- 5:14

# 2) For each cutoff:
# a) define Young/Old at a, and center age at a
# b) fit the same model

age_threshold_test <- function(a) {
  
  d <- model_variables %>%
    mutate(age_cat = factor(if_else(age >= a, "Old", "Young"),
                            levels = c("Young","Old")), 
           age10 = (age - a) / 10)
  
  m <- glmer(MOA_proportion ~ age_cat : age10 + (1 | animalID_fct) + (1 | season_fct),
             weights = total_resights,
             family = binomial(link="logit"),
             control = glmerControl(optimizer="bobyqa"),
             data = d)
  
  # 3) Return model fit AIC
  tibble(Threshold = a,
         AIC = AIC(m),
         logLik = as.numeric(logLik(m)))
}

# 4) Fit all thresholds automatically
aic_age_threshold_comparison <- map_dfr(age_cutoff, age_threshold_test) %>%
  mutate(delta_AIC = AIC - min(AIC),
         AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))) %>%
  arrange(Threshold)

age_threshold_comparison_tbl <- aic_age_threshold_comparison %>%
  transmute("Threshold ages" = Threshold,
            AIC = round(AIC, 1),
            "ΔAIC" = round(delta_AIC, 2),
            "AIC weight" = round(AIC_weight, 3)) %>%
  arrange("Threshold ages")

# 5) identify best threshold
best_age_threshold <- age_threshold_comparison_tbl %>%
  slice_min(AIC, n = 1) %>%
  pull("Threshold ages")

# 6) make flextable and bold the lowest AIC threshold
age_threshold_comparison_ft <- flextable(age_threshold_comparison_tbl) %>%
  align(align = "center", part = "all") %>%
  bold(i = which(age_threshold_comparison_tbl$"Threshold ages" == best_age_threshold), part = "body"); age_threshold_comparison_ft

save_as_docx(age_threshold_comparison_ft, path = here("TablesFigures", "AIC_Age_Threshold_Comparison.docx"))

#### Table S11 (experience threshold comparison) ####

# 1) test all possible experience thresholds
exp_cutoff <- 1:11

# 2) For each cutoff:
# a) define inexperienced vs. experienced at a, and center experience at a
# b) fit the same model

threshold_test_exp <- function(a) {
  
  d <- model_variables %>%
    mutate(experience_cat = factor(if_else(pupping_exp >= a, "Experienced", "Inexperienced"),
                                   levels = c("Inexperienced", "Experienced")),
           exp10 = (pupping_exp - a) / 10)
  
  m <- glmer(MOA_proportion ~ experience_cat : exp10 + (1 | animalID_fct) + (1 | season_fct),
             weights = total_resights,
             family = binomial(link = "logit"),
             control = glmerControl(optimizer = "bobyqa"),
             data = d)
  # 3) Return model fit AIC
  tibble(Threshold = a,
         AIC = AIC(m),
         logLik = as.numeric(logLik(m)))
}

# 4) Fit all thresholds automatically
aic_experience_threshold_comparison <- map_dfr(exp_cutoff, threshold_test_exp) %>%
  mutate(delta_AIC = AIC - min(AIC),
         AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))) %>%
  arrange(Threshold)

exp_threshold_comparison_tbl <- aic_experience_threshold_comparison %>%
  transmute("Threshold experience levels" = Threshold,
            AIC = round(AIC, 1),
            "ΔAIC" = round(delta_AIC, 2),
            "AIC weight" = round(AIC_weight, 3)) %>%
  arrange("Threshold experience levels")

# 5) identify best threshold
best_threshold <- exp_threshold_comparison_tbl %>%
  slice_min(AIC, n = 1) %>%
  pull("Threshold experience levels")

# 6) bold the lowest AIC threshold
exp_threshold_comparison_ft <- flextable(exp_threshold_comparison_tbl) %>%
  align(align = "center", part = "all") %>%
  bold(i = which(exp_threshold_comparison_tbl$"Threshold experience levels" == best_threshold), part = "body"); exp_threshold_comparison_ft

save_as_docx(exp_threshold_comparison_ft, path = here("TablesFigures", "AIC_Experience_Threshold_Comparison.docx"))

#### Figure S2 (variation in pupping experience within ages) ####

# reproductive experience vs age figure

# most typical recruitment age
recruitment_age <- model_variables %>%
  filter(pupping_exp == 0) %>% #first observed reproductive event
  count(age, sort = TRUE) %>%
  slice(1) %>%
  pull(age)

recruitment_age

#Most common recruitment age is 4
typical_recruitment_age <- 4

plot_variation_exp <- ggplot(model_variables, 
                             aes(x = age,
                                 y = pupping_exp)) +
  geom_jitter(width = 0.3,
              height = 0.2,
              alpha = 0.4,
              size = 2.2,
              color = "#C8548C") +
  geom_segment(aes(x = 4,
                   y = 0,
                   xend = 21,
                   yend = 17),
               color = "#80113C",
               linewidth = 1) +
  geom_vline(xintercept = typical_recruitment_age,
             linetype = "dashed",
             linewidth = 1) +
  annotate("text",
           x = typical_recruitment_age,
           y = 8,
           label = "Typical recruitment age (4 years old)",
           angle = 90,
           vjust = -0.6,
           size = 4) +
  coord_cartesian(xlim = c(3, 21)) +
  coord_cartesian(ylim = c(0, 18)) +
  scale_x_continuous(breaks = 3:21) +
  scale_y_continuous(n.breaks = 10) +
  theme_few(base_size = 16) +
  labs(x = "Maternal age (years)",
       y = "Previous pupping experience (number of pups)"); plot_variation_exp

ggsave(here("TablesFigures", "FigureS2.png"), plot_variation_exp, width = 8, height = 6, dpi = 600)

#### Figure S3 (age threshold significance comparison) ####

#test all possible thresholds
age_cutoff <- 5:14

# For each cutoff:
# 1) define Young/Old at a, and center age at a
# 2) fit the same model
# 3) pull the post-threshold slope term (Old × age10), its SE, and p-value
thr_res_age <- map_dfr(age_cutoff, \(a){
  d <- model_variables %>%
    mutate(age_cat = factor(if_else(age >= a, "Old", "Young"),
                            levels = c("Young","Old")),
           age10   = (age - a) / 10)
  
  m <- glmer(MOA_proportion ~ age_cat : age10 + (1 | animalID_fct) + (1 | season_fct),
             weights = total_resights,
             family = binomial(link="logit"),
             control = glmerControl(optimizer="bobyqa"),
             data = d)
  
  #coefficient table
  tab <- as.data.frame(summary(m)$coefficients) %>%
    rownames_to_column("term")
  
  #interaction term = slope of age10 for Old
  term_name <- "age_catOld:age10"
  
  row <- filter(tab, term == term_name)
  
  tibble(age_cutoff = a,
         coef = row$Estimate,
         se = row$`Std. Error`,
         p = 2 * pnorm(abs(row$`z value`), lower.tail = FALSE))
}) %>%
  mutate(sig = p < 0.05, # flag significance
         logp_cap = pmin(-log10(p), 5)) # point size weighted by p

# Plot: coefficient ± 95% CI; bigger points = smaller p; color = significance
age_threshold_comparison_figure <- ggplot(thr_res_age, aes(age_cutoff, coef)) +
  geom_hline(yintercept = 0, linetype = "dashed") +  # no post-threshold slope
  geom_pointrange(aes(ymin = coef - 1.96 * se,
                      ymax = coef + 1.96 * se,
                      color = sig,
                      size = logp_cap)) +
  scale_color_manual(values = c(`TRUE` = "deepskyblue", `FALSE` = "tomato"),
                     labels = c(`TRUE` = "p < 0.05", `FALSE` = "p > 0.05"),
                     name = "Significance") +
  scale_size_continuous(range = c(0.5, 2.2), guide = "none") +
  scale_x_continuous(breaks = cutoff) +
  labs(x = "Age senescence threshold",
       y = "Estimated old slope (log-odds)") +
  theme_classic(); age_threshold_comparison_figure

ggsave(here("TablesFigures","FigureS3.png"), age_threshold_comparison_figure, width = 8, height = 6, dpi = 600)

###################### Figure S4 (experience threshold significance comparison) ############################

#test all possible thresholds
exp_cutoff <- 1:11

# For each cutoff:
# 1) define Inexperienced/Experienced at e, and center experience at e
# 2) fit the same model
# 3) pull the post-threshold slope term (Experienced × exp10), its SE, and p-value
thr_res_exp <- map_dfr(exp_cutoff, \(e){
  d <- model_variables %>%
    mutate(experience_cat = factor(if_else(pupping_exp >= e, "Experienced", "Inexperienced"),
                            levels = c("Inexperienced","Experienced")),
           exp10 = (pupping_exp - e)/10)
  
  m <- glmer(MOA_proportion ~ experience_cat : exp10 + (1 | animalID_fct) + (1 | season_fct),
             weights = total_resights,
             family = binomial(link="logit"),
             control = glmerControl(optimizer="bobyqa"),
             data = d)
  
  #coefficient table
  tab <- as.data.frame(summary(m)$coefficients) %>%
    rownames_to_column("term")
  
  #interaction term = slope of age10 for Old
  term_name <- "experience_catExperienced:exp10"
  
  row <- filter(tab, term == term_name)
  
  tibble(exp_cutoff = e,
         coef = row$Estimate,
         se = row$`Std. Error`,
         p = 2 * pnorm(abs(row$`z value`), lower.tail = FALSE))
}) %>%
  mutate(sig = p < 0.05, # flag significance
         logp_cap = pmin(-log10(p), 5)) # point size weighted by p

# Plot: coefficient ± 95% CI; bigger points = smaller p; color = significance
experience_threshold_comparison_figure <- ggplot(thr_res_exp, aes(exp_cutoff, coef)) +
  geom_hline(yintercept = 0, linetype = "dashed") +  # no post-threshold slope
  geom_pointrange(aes(ymin = coef - 1.96 * se,
                      ymax = coef + 1.96 * se,
                      color = sig,
                      size = logp_cap)) +
  scale_color_manual(values = c(`TRUE` = "deepskyblue", `FALSE` = "tomato"),
                     labels = c(`TRUE` = "p < 0.05", `FALSE` = "p > 0.05"),
                     name = "Significance") +
  scale_size_continuous(range = c(0.5, 2.2), guide = "none") +
  scale_x_continuous(breaks = exp_cutoff) +
  labs(x = "Experience senescence threshold",
       y = "Estimated experienced slope (log-odds)") +
  theme_classic(); experience_threshold_comparison_figure

ggsave(here("TablesFigures", "FigureS4.png"), experience_threshold_comparison_figure, width = 8, height = 6, dpi = 600)

