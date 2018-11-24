library(tidyverse)
#df_10 <- read_csv("data/LTDB_Std_All_fullcount/LTDB_Std_2010_fullcount.csv") 
#df <- read_csv("data/LTDB_Std_1980_fullcount.csv") # read data 
#df <- full_join(df, df_10['tractid'], by = c('TRTID10'='tractid')) %>% 
#  mutate_if(is.numeric, funs(replace(., is.na(.), 0)))

race_typology <- function(df){
columns <- c("TRTID10","POP80", "NHWHT80","NHBLK80","NTV80","ASIAN80","HISP80") # set columns to grab
labels <- c("TRTID10","pop_full", 'white','black','native','asian','latinx') # set labels to use
re_cols <- c('white','black','native','asian','latinx') # subset of just race and ethnicity columns
df <- df[columns] # select the columns
names(df) <- labels # apply the labels
props <- bind_cols(df["TRTID10"], setNames(df[re_cols]/df$pop_full, re_cols)) # make racial proportions

types <- props %>% gather(race, proportion, -TRTID10) %>% #go to long format
  group_by(TRTID10) %>% # group by tract
  summarise(first = last(race,order_by = proportion), # use summarise to extract the top three labels and values
            second = nth(race,-2,order_by = proportion), # maybe this is clumsy but its fast and its after midnight
            third = nth(race,-3,order_by = proportion), 
            highest = last(proportion,order_by = proportion), 
            middle = nth(proportion,-2,order_by = proportion), 
            low = nth(proportion,-3,order_by = proportion)) %>%
  mutate( 
    type = ifelse(highest>.8&first=='white', paste("predominantly", first), NA), # if the first one is white and over 80, it's predominantly that 
    type = ifelse(highest>.5&first!='white', paste("predominantly", first), type), # if the first one is not white and over 50, it's predominantly that 
    type = ifelse(highest=='white'&second>.4, paste(second, first), type), # if the first is white but the second is over 40%, make it second first
    type = ifelse(highest<=.4, "mixed", type), #if nothing it over 40, it's a mixed neighborhood
    type = ifelse(highest>.4&middle>.1&low>.1,paste(first, second, third), type), # if the second and third heighest ar both over 10, it's all three in that order 
    type = ifelse(middle>.1&low<=.1,paste(first, second), type), # if the second is over 20 but the third isn't, we'll just label it with the first two groups,
    type = ifelse(highest<=.8&middle<=.1,paste(first,'mixed'),type), # if the second is below 20 then we'll call label it with the first group and mixed
    type = ifelse(is.na(highest),'empty',type), # if there's an NA it's because the tract proportion divided by zero: an empty tract
    ####################
    ####################
    collapsed_type = ifelse(highest>.8&first=='white', paste("predominantly", first), NA), # if the first one is white and over 80, it's predominantly that 
    collapsed_type = ifelse(highest>.5&first!='white', paste("predominantly", first),  collapsed_type), # if the first one is not white and over 50, it's predominantly that 
    collapsed_type = ifelse(highest=='white'&second>.4, paste(second, first), collapsed_type), # if the first is white but the second is over 40%, make it second first
    collapsed_type = ifelse(highest<=.4, "mixed", collapsed_type), #if nothing it over 40, it's a mixed neighborhood
    collapsed_type = ifelse(highest>.4&middle>.1&low>.1,paste(first, second), collapsed_type), # if the second is over 10, it's first second
    collapsed_type = ifelse(middle>.1&low<=.1,paste(first, second), collapsed_type), # if the second is over 10 but the third isn't, we'll just label it with the first two groups,
    collapsed_type = ifelse(highest<=.8&middle<=.1,paste(first,'mixed'),collapsed_type), # if the second is below 10 then we'll call label it with the first group and mixed
    collapsed_type = ifelse(is.na(highest),'empty',collapsed_type) # if there's an NA it's because the tract proportion divided by zero: an empty tract
  )
return(types %>% select(TRTID10, type, collapsed_type))
}
