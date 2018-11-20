## Location Quotient
# CSSS 510 final project
# Shihao Han

library(tidyverse)
library(readr)

df <- read_csv("LTDB_Std_1980_fullcount.csv") # Loading data

# Calculating the proportion of tract race population to tract population

shihao <- function(df) {
  
df <-                            
  df %>%
  arrange(desc(county)) %>%
  mutate(prop_white_tract_80 = NHWHT80 / POP80) %>%
  mutate(prop_black_tract_80 = NHBLK80 / POP80) %>%
  mutate(prop_hisp_tract_80 = HISP80 / POP80) %>%
  mutate(prop_native_tract_80 = NTV80 / POP80) %>%
  mutate(prop_asian_tract_80 = ASIAN80 / POP80)

# Calculating the total population of all counties and total population of races in counties
df$county <- as.factor(df$county) # Factorize each county
sum_pop_80 <-                                   # Calculate total population and population of races in each county
  df %>%                                   
  group_by(county) %>%
  summarize(sum_pop_80 = sum(POP80, na.rm = TRUE),
            sum_white_80 = sum(NHWHT80, na.rm = TRUE),
            sum_black_80 = sum(NHBLK80, na.rm = TRUE),
            sum_hispanic_80 = sum(HISP80, na.rm = TRUE),
            sum_native_80 = sum(NTV80, na.rm = TRUE),
            sum_asian_80 = sum(ASIAN80, na.rm = TRUE))
df <-                                    # merging total population variables with dataset
  merge(df, sum_pop_80, by = "county", all.x = TRUE) 

# calculating proportion of county race population to county population
df <-
  df %>%
  mutate(prop_white_county_80 = sum_white_80 / sum_pop_80,
         prop_black_county_80 = sum_black_80 / sum_pop_80,
         prop_hisp_county_80 = sum_hispanic_80 / sum_pop_80,
         prop_native_county_80 = sum_native_80 / sum_pop_80,
         prop_asian_county_80 = sum_asian_80 / sum_pop_80)

# Get Location Quotient using tract proportion over county proportion
df <- 
  df %>%
  mutate(LQ_white_80 = prop_white_tract_80 / prop_white_county_80,
         LQ_black_80 = prop_black_tract_80 / prop_black_county_80,
         LQ_hispanic_80 = prop_hisp_tract_80 / prop_hisp_county_80,
         LQ_native_80 = prop_native_tract_80 / prop_native_county_80,
         LQ_asian_80 = prop_asian_tract_80 / prop_asian_county_80)

# Return the results
return(df %>%
         select(TRTID10, LQ_white_80, LQ_balck_80, LQ_hispanic_80, LQ_native_80, LQ_asian_80))
}