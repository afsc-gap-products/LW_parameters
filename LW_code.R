# Estimating LW parameters for EBS & NBS survey species
# Created by: Liz Dawson
# Contact: liz.dawson@noaa.gov
# Created: 2/10/2021
# Modified: 2023-02-07
# Modified by: Liz Dawson, liz.dawson@noaa.gov


# install and load packages -----------------------------------------------

install.packages("FSA")
install.packages("FSAdata")
install.packages("RODBC")
library(FSA)
library(FSAdata)
library(here)
library(tidyverse)
library(lubridate)
library(RODBC)
library(rstudioapi)
library(janitor)
library(readr)
library(broom)

# renv::init()
# renv::restore()

# update warm and cold years each year before running LW regression --------
#defines warm and cold years

# future updates: prompt add warm or cold from last year; 
# # write year-condition to csv, so you don't have to manually add each year to warm or cold year_condition
# # auto-run both warm and cold conditions and produce both outputs

year_condition <- "all"
# year_condition<-"warm"
# year_condition <- "cold"
present<-year(today())
if(year_condition=="warm") {year_select<-c(2003:2005,2014:2020,present)}
if(year_condition=="cold"){year_select<-c(2006:2013, present)}
if(year_condition=="all"){year_select<-c(2003:present)}
year_select

# need to put in survey ID from last cruise (found in RACEBASE.SURVEYS)
# # SURVEY_DEFINITION_ID = 98 for EBS and 143 for NBS
last_survey_id <-  c(589, 591) #2023---date each year, Liz updated with 2023 values already!
  #c(587, 588) 2022  
  #c(585, 584) 2021

cruise <- paste0(202201, ",", 202202) #update each year
vessel_nums <- c(162, 134)  #check each year
#94 = vesteraalen
#134 = NW explorer
# 162 = AK Knight

# read in specimen data ------------------------------------------------------------

#read in data (This needs to change to ROracle or ROBDC script to pull in data directly from Oracle database)
 channel <- odbcConnect(dsn = "AFSC",
                        uid = rstudioapi::showPrompt(title = "Username",
                                                     message = "Oracle Username", default = ""),
                        pwd = rstudioapi::askForPassword("enter password"),
                        believeNRows = FALSE)
# 
# ## checks to see if connection has been established
odbcGetInfo(channel)
# 
# # pollockEBSNBSwarm<-read.csv(file=here::here("data","POLL_EBS_NBS_shelf_warm.csv"))
# 
# #SQL code to get all EBS NBS BS data---really only need YEAR, SPECIES_CODE, LENGTH, WEIGHT from it.
# start <- Sys.time()
oracle_data <- sqlQuery(channel, "   SELECT * FROM RACEBASE.SPECIMEN a
                                                JOIN SAFE.SURVEY b
                                                ON a.CRUISEJOIN =b.CRUISEJOIN
                                                JOIN RACE_DATA.V_CRUISES c
                                                ON a.CRUISEJOIN = c.CRUISEJOIN
                                                WHERE a.REGION = 'BS'
                                                AND a.CRUISE >= 200301
                                                AND c.SURVEY_DEFINITION_ID = 143
                                                OR c.SURVEY_DEFINITION_ID = 98")

# end <- Sys.time() #took ~5 minutes to run
# 
write_csv(oracle_data, path = here("data", paste0("oracle_data_", today(), ".csv"))) #~45 MB, takes a while
# orc_data <- read_csv(file = here::here("data", "oracle_data_2021-02-22.csv"))
orc_data <- read_csv(file = here::here("data", "oracle_data_2023-02-07.csv"))

orc_data

spec_data <- #oracle_data %>% as_tibble %>% 
  orc_data %>%
  clean_names() %>% 
  dplyr::filter(year >= 2003)

head(spec_data)

# read in previous LW and max L parameters --------------------------------

# Survey ID 576 or 577 
# 2021: 585 (EBS), 584 (NBS)
query_command <- paste0(" select * from race_data.cruises where cruise in (", cruise,");")
cruise_info <- sqlQuery(channel, query_command) %>% 
  as_tibble() %>% 
  clean_names()

cruise_id_nums <- cruise_info %>% dplyr::filter(vessel_id %in% vessel_nums) 

# PAY ATTENTION TO THE FOLLOWING LINE: make sure you are getting all unique survey id numbers!
survey_id <- paste0(unique(cruise_id_nums$survey_id)[1], "," , unique(cruise_id_nums$survey_id)[2])

# get previous parameters using previous survey
query_command <- paste0(" select * from RACE_DATA.SURVEY_LENGTH_STANDARDS where SURVEY_ID in (", survey_id,");")

previous_parameters<-sqlQuery(channel, query_command)

write_csv(previous_parameters, path = here("data", paste0("previous_parameters_", today(), ".csv"))) 
# prev_parameters_csv <- read_csv(file = here::here("data", "previous_parameters_2021-03-01.csv")) #2020
prev_parameters_csv <- read_csv(file = here::here("data", "previous_parameters_2023-02-07.csv")) #update with date of previous parameters
prev_params <- prev_parameters_csv %>% janitor::clean_names()

# read in species codes ---------------------------------------------------

#read in species codes
lw2022 <- read.csv(here::here("data","new_species_params_all_2022-02-25.csv")) #update with date of previous parameters file
# species_codes<-lw2022[,"species_code"]
# species_names<-lw2022[,"common_name"]
# name_code <- lw2022 %>% 
#   dplyr::select(common_name, species_code)
 write_csv(name_code, path = here("data", "lw_species_names.csv"))

name_code <- read_csv(here("input", "lw_species_names.csv"))

#survey = 3 for EBS, survey = 7 for NBS, survey = 4 for EBS Slope, survey = 6 for Chukchi Sea
new_spec <- anti_join(name_code, lw2022) %>% 
  mutate(survey = 3) %>% 
  mutate(poly_species_code = case_when(species_code %% 10200 == 0 ~ 18,
                                       species_code %% 10220 == 0 ~ 22)) #when you're adding a new species, 
                                                                      # add polycorder code that corresponds to species number
                                                                      # race_data.race_species_codes has polycorder numbers
# CIA: add feature to pull polycorder code

lw2022 <- full_join(lw2022, new_spec)

# filter data -------------------------------------------------------------

spec_dat <- spec_data %>% 
  dplyr::filter(survey_definition_id %in% c(98,143)) %>% 
  dplyr::filter(species_code %in% name_code$species_code) %>% 
  dplyr::filter(year %in% year_select) %>% 
  dplyr::select(species_code, length, weight) %>% 
  full_join(name_code)

# get l-w regressions -----------------------------------------------------

# group by species code, get l-w regr by species, pull alpha and beta (broom pkg)

get_lw <- function(df, group){ #get l-w relationship from the data frame (df) by species (group)
  
  spec_test <- df %>% dplyr::filter(species_code == group) #filter by species code
  
  #define columns
  weight<-spec_test$weight #set up for nls(), get weight
  length<-spec_test$length #set up for nls(), get length
  
  # test section
  l_w_index <- bind_cols(weight = weight, length = length) %>% #bind matches together
    dplyr::filter(!is.na(weight)) %>% #remove na values
    dplyr::filter(!is.na(length)) %>% 
    dplyr::filter(weight > 0) %>% #remove 0 values
    dplyr::filter(length > 0)
  
  new_weight <- l_w_index$weight
  new_length <- l_w_index$length
  if(length(new_weight) < 1 || length(new_length) < 1) 
  { print(paste("there is no data for this species within this subset of years:", group))
    final <- bind_cols(species_code = group, # return a tibble with spec code, alpha, beta
                       alpha = NA, 
                       beta = NA)
    print(final)
    return(final)
  }

  #run nonlinear least squares regression and print summary report
  nls_LW_fit<-nls(new_weight~alpha*new_length^beta,
                  data=list(new_length, new_weight),
                  start=list(alpha=.001,beta=5),
                  trace=F, 
                  control=list(maxiter = 500))
  
  lw_results <- broom::tidy(nls_LW_fit) #get tidy results
  print(c("species", group)) #print output as fxn runs
  print(lw_results) #print output as fxn runs
  
  alpha_result <- lw_results$estimate[1]*1000 #extract alpha value from summary results
  beta_result <- lw_results$estimate[2]  #extract beta value from summary results
  final <- bind_cols(species_code = group, # return a tibble with spec code, alpha, beta
                     alpha = alpha_result, 
                     beta = beta_result)
  
  return(final) #return final tibble
}

get_params <- map(unique(spec_dat$species_code), get_lw, df = spec_dat) #for each unique species code, use get_lw function, with spec_dat data frame

for(i in 1:length(unique(spec_dat$species_code))) # re-arrange results returned as list into tibble (data frame)
{
  print(i)
  if(i == 1) {param_dat <- tibble()}
  param_dat <- bind_rows(param_dat, get_params[[i]])
}
param_dat <- param_dat %>% arrange(species_code) #final values for each species

dir.create(here("output"), showWarnings = FALSE)
dir.create(here("output", "histograms"), showWarnings = FALSE)
write_csv(param_dat, path = here("output", paste0("species_LW_param_ests_", year_condition, "_", today(),".csv"))) #save results as .csv to output folder with date

# CHECK against previous parameters

# max lengths -------------------------------------------------------------

# * clean and check max length data ---------------------------------------

# if adding new species, do a quick data check 

# step 1: check to see if there are multiple values for max sp length
prev_params %>% 
  dplyr::filter(is.na(maximum_length_source)) %>% 
  group_by(species_code) %>% 
  summarise(sp_len_min = min(maximum_length),
            sp_len_max = max(maximum_length)) %>% 
  mutate(sp_len_check = sp_len_max - sp_len_min) %>% 
  dplyr::filter(sp_len_check > 0)

# step 2: check which are legit and which are arbitrary 999 or 9999 etc.\
prev_params %>% 
  dplyr::filter(is.na(maximum_length_source)) %>% 
  dplyr::filter(maximum_length != 999) %>% 
  dplyr::filter(maximum_length != 9999) %>% 
  left_join(name_code) #names with NA are species estimated in previous parameters, but NOT currently included in our list of species to do.
# # do these numbers make sense, given the species? Or are they incorrect entries?

# set and check max lengths -----------------------------------------------

# pull current max lengths
len_data <- orc_data %>% 
  clean_names() %>% 
  dplyr::filter(survey_definition_id %in% c(98,143)) %>%
  dplyr::filter(species_code %in% name_code$species_code) %>% 
  group_by(species_code) %>% 
  summarise(max_length = max(length, na.rm = T),
            pctile_99 = quantile(length, .999, na.rm = TRUE)) %>% 
  mutate(species_code = as.integer(species_code)) %>% 
  left_join(name_code)

# use previous_parameters to check prev max lengths
check_lens <- prev_params %>% 
  rename(max_len_prev = maximum_length) %>% 
  dplyr::filter(species_code %in% name_code$species_code) %>% 
  dplyr::filter(survey_id %in% last_survey_id) %>% 
  dplyr::select(species_code, max_len_prev) %>% 
  mutate(species_code = as.integer(species_code)) %>% arrange(species_code) %>% 
  distinct()
     
all_lengths <- full_join(x = len_data, y = check_lens, by = "species_code") %>% 
  mutate(flag = if_else(max_length >= pctile_99*1.1, TRUE, FALSE)) # flag = true means problem
  # flag species with max length >10% 99th percentile values

check_spec <- all_lengths %>% 
  filter(flag == TRUE)

write_csv(check_spec, path = here("output", paste0("species_check_length_results", today(), ".csv")))

hist_plot <- spec_data %>%  
  dplyr::filter(species_code %in% check_spec$species_code) 

for(i in unique(hist_plot$species_code))
{
  to_plot <- hist_plot %>% 
    dplyr::filter(species_code == i)
  title_name <- name_code %>% 
    dplyr::filter(species_code == i)
  
  p <- ggplot() +
    geom_histogram(data = to_plot, aes(length)) +
    labs(title = paste("Species:", title_name$common_name, '\n', "Species code:", i))
  
  ggsave(plot = p, filename = paste0("species_code_", i, ".png"), path = here("output", "histograms"))
}

# look at histograms of length distributions by species- check for outliers
# # way to auto-check for outliers?

# select max length
new_lengths <- NULL
length_source <- NULL
for(i in unique(check_spec$species_code))
{
  use <- check_spec %>% dplyr::filter(species_code == i)
  print(paste("species name:", use$common_name))
  print(paste("species code:", use$species_code))
  print(paste("current max length (new):", use$max_length))
  print(paste("previous max length (old):", use$max_len_prev))
  print(paste("length >99.9% all lengths:", round(use$pctile_99, digits = 0)))
  
  length_mod <- readline(prompt = 
                           "Which length do you want to use: n = new, o = old, p = 99.9th percentile, c = custom? ")
  if(length_mod == "n")
  {
    length_used <- use %>% 
      dplyr::select(species_code, max_length)
  }else if(length_mod == "o")
  {
    length_used <- use %>% 
      dplyr::select(species_code, max_len_prev) %>% 
      rename(max_length = max_len_prev)
  }else if(length_mod == "p")
  {
    length_used <- use %>% 
      dplyr::select(species_code, pctile_99) %>% 
      rename(max_length = pctile_99) %>% 
      mutate(max_length = round(max_length, digits = 0))
  }else if(length_mod == "c")
  {
    custom_len <- readline(prompt = "Please enter a max length: ")
    length_used <- bind_cols(species_code = use$species_code, max_length = as.double(custom_len))
  }else{print("Please enter a valid entry: n = new, o = old, p = 99.9th percentile, c = custom")}
  new_lengths <- bind_rows(new_lengths, length_used)
  length_source <- c(length_source, length_mod)
}
new_lengths
length_source

new_lengths <- new_lengths %>% 
  bind_cols(length_source = length_source)

# output everything (all max lengths)

keep_max_lens <- all_lengths %>% 
  dplyr::filter(flag == FALSE) %>% # keep non-flagged species with new data values
  dplyr::select(species_code, max_length) %>% 
  mutate(length_source = "n") %>% 
  bind_rows(new_lengths) %>% 
  arrange(species_code) %>% 
  full_join(name_code)

write_csv(keep_max_lens, path = here("output", paste0("species_max_lens_", today(), ".csv")))

# keep_max_lens <- read_csv(here("output", "species_max_lens_2022-02-25.csv"))

# full output -------------------------------------------------------------

lw2022

keep_max_lens

param_dat

# get sources for previous max length
prev_dat <- prev_params %>% 
  dplyr::select(species_code, maximum_length_source) %>% 
  dplyr::filter(species_code %in% name_code$species_code) %>% 
  drop_na(maximum_length_source) %>% 
  arrange(species_code)
prev_dat <- distinct(prev_dat) #remove duplicate rows

# get lw source label
if(year_condition == "warm") {lw_label = paste("EBS + NBS shelf warm years: ", paste(year_select, collapse = ', '))
}else if(year_condition =="cold"){lw_label = paste("EBS + NBS shelf cold years: ", paste(year_select, collapse = ', '))
}else if(year_condition =="all"){lw_label = paste("EBS + NBS shelf all years: ", paste("2003 to ", (present-1) ))}

full_output <- lw2022 %>% 
  dplyr::select(-alpha, -beta, -max_length, -current_max_len_source, -lw_parameters_source) %>% 
  full_join(param_dat) %>%
  mutate(lw_parameters_source = lw_label) %>% 
  full_join(keep_max_lens) %>% 
  full_join(prev_dat) %>% 
  mutate(current_max_len_source = case_when(length_source == "n" ~ paste("data: max length from current data"),
                          length_source == "p" ~ paste(" data: 99.9th percentile length from current data"),
                          length_source == "o" ~ paste("max length value used in previous years"))) %>% 
  dplyr::select(-maximum_length_source, -length_source) %>% 
  dplyr::relocate(max_length, .before = lw_parameters_source)

write_csv(full_output, path = here("output", paste0("new_species_params_", year_condition, "_", today(), ".csv")))

# future: check all species (instead of core 16), for max lengths
