## Neighborhood clustering
# CSSS 510 final project
# Hannah Lee

rm(list=ls())

library(dplyr)
library(withr)
library(stringr)
library(readr)
library(tidyverse)
library(tidyr)
library(magrittr)
library(pander)
library(ggplot2)
options(scipen=999)
library(rgdal)
library(spdep)

hannah <- function(df){
#### Read in 1980 decennial census counts ####
pop1980 <- df %>% # there were 59,187 census tracts in 1980
  select(TRTID10:NHBLK80, ASIAN80, HISP80) %>%
  rename(TRACTID = TRTID10) %>%
  mutate(TRACTID = str_pad(TRACTID, width = 11, side = "left", pad = "0")) %>% # add leading zero to get common width of 11
  mutate(StateFIPS = substr(TRACTID, 0, 2),
         CountyFIPS = substr(TRACTID, 3, 5),
         TractFIPS = substr(TRACTID, 6, 11)) %>%
  mutate(prop.nhwht80 = ifelse(POP80 > 0, NHWHT80/POP80,0), # create proportion of each racial group relative to the total population in the census tract; if there are census tracts with zero population, then leave them as zero
         prop.nhblk80 = ifelse(POP80 > 0, NHBLK80/POP80, 0),
         prop.asian80 = ifelse(POP80 > 0, ASIAN80/POP80, 0),
         prop.hisp80 = ifelse(POP80 > 0, HISP80/POP80, 0))

# get the average proportion of each ethnic group for each county
county.avgpop1980 <- pop1980 %>%
  select(TRACTID, StateFIPS, CountyFIPS, prop.nhwht80:prop.hisp80) %>%
  group_by(StateFIPS, CountyFIPS) %>%
  summarise(mean(prop.nhwht80), 
            mean(prop.nhblk80), 
            mean(prop.asian80), 
            mean(prop.hisp80))

pop1980.foranalysis <- pop1980 %>%
  left_join(county.avgpop1980, by = c("StateFIPS", "CountyFIPS")) %>%
  mutate(whtcnctn = ifelse(`mean(prop.nhwht80)` > 0, prop.nhwht80 / `mean(prop.nhwht80)`, 0), # calculate concentration of each racial group for each census tract relative to the average of the county. If the average proportion in the county is zero, then the concentration is zero
         blkcnctn = ifelse(`mean(prop.nhblk80)` > 0, prop.nhblk80 / `mean(prop.nhblk80)`, 0),
         asiancnctn = ifelse(`mean(prop.asian80)` > 0, prop.asian80 / `mean(prop.asian80)`, 0),
         hispcnctn = ifelse(`mean(prop.hisp80)` > 0, prop.hisp80 / `mean(prop.hisp80)`, 0)) %>%
  select(TRACTID, StateFIPS, CountyFIPS, TractFIPS, 
         prop.nhwht80:prop.hisp80,
         whtcnctn:hispcnctn) %>%
  mutate(TRACTID = as.numeric(TRACTID))

#### Read in 2010 census tract shapefiles and the spatial weights ####
if(!file.exists("2010tractshape.Rdata")){
  tract.shape.2010 <- readOGR(dsn = "TractS2010",
                              layer = "tract2010")
  
  save(tract.shape.2010, 
       file = "2010tractshape.Rdata")
}

load("2010tractshape.Rdata") 

tract.shape.clean <- tract.shape.2010[!tract.shape.2010@data$STATEFP10==72,] # filter out puerto rico 
xy <- coordinates(tract.shape.clean) # retrieve coordinates of tracts to calculate spatial distance 
nb.d6 <- dnearneigh(xy, 0, 6, longlat = TRUE) # distance-based spatial influence at 6 kilometers
nb.d6.list <- nb2listw(nb.d6, style = "W", zero.policy = TRUE) # create spatial weights for neighbors list. This is needed to calculate local moran's i

# create data frame to use to calculate local moran's i
df.foranalysis <- tract.shape.clean@data %>% 
  mutate(TRACTID = as.numeric(as.character(TRACTID))) %>%
  left_join(pop1980.foranalysis, by = "TRACTID") %>%
  mutate(prop.nhwht80 = replace_na(prop.nhwht80, 0), # replace the NA's with zeroes for those 2010 census tract shape files that do not have 1980 data
         prop.nhblk80 = replace_na(prop.nhblk80, 0), 
         prop.asian80 = replace_na(prop.asian80, 0),
         prop.hisp80 = replace_na(prop.hisp80, 0),
         whtcnctn = replace_na(whtcnctn, 0), 
         blkcnctn = replace_na(blkcnctn, 0),
         asiancnctn = replace_na(asiancnctn, 0),
         hispcnctn = replace_na(hispcnctn, 0))

# calculate local moran's i
# create a dummy indicating TRUE if the census tract is a focal tract of a cluster and FALSE if the census tract is not a focal tract of a cluster
df.lm <- df.foranalysis %>%
  cbind(localmoran(x = df.foranalysis$prop.nhwht80, 
                   list = nb.d6.list, 
                   na.action = na.pass, 
                   zero.policy = TRUE)) %>%
  select(-c(Ii, E.Ii, Var.Ii)) %>%
  rename(`wh.lm.pr(z>0)` = `Pr(z > 0)`,
         `wh.lm.Z.Ii` = `Z.Ii`) %>%
  mutate(wht.clst.dmmy = ifelse(whtcnctn > 1 & 
           `wh.lm.pr(z>0)`<0.05 & 
           `wh.lm.Z.Ii` > 0, 1, 0),
         wht.clst.dmmy = replace_na(wht.clst.dmmy, 0)) %>%
  cbind(localmoran(x = df.foranalysis$prop.asian80, 
                   list = nb.d6.list, 
                   na.action = na.pass,
                   zero.policy = TRUE)) %>%
  select(-c(Ii, E.Ii, Var.Ii)) %>%
  rename(`as.lm.pr(z>0)` = `Pr(z > 0)`,
         `as.lm.Z.Ii` = `Z.Ii`) %>%
  mutate(as.clst.dmmy = ifelse(asiancnctn > 1 & 
           `as.lm.pr(z>0)`<0.05 & 
           `as.lm.Z.Ii` > 0, 1, 0),
         as.clst.dmmy = replace_na(as.clst.dmmy, 0)) %>%
  cbind(localmoran(x = df.foranalysis$prop.nhblk80, 
                   list = nb.d6.list, 
                   na.action = na.pass,
                   zero.policy = TRUE)) %>%
  select(-c(Ii, E.Ii, Var.Ii)) %>%
  rename(`blk.lm.pr(z>0)` = `Pr(z > 0)`,
         `blk.lm.Z.Ii` = `Z.Ii`) %>%
  mutate(blk.clst.dmmy = ifelse(blkcnctn > 1 & 
           `blk.lm.pr(z>0)`<0.05 & 
           `blk.lm.Z.Ii` > 0, 1, 0),
         blk.clst.dmmy = replace_na(blk.clst.dmmy, 0)) %>%
  cbind(localmoran(x = df.foranalysis$prop.hisp80, 
                   list = nb.d6.list, 
                   na.action = na.pass,
                   zero.policy = TRUE)) %>%
  select(-c(Ii, E.Ii, Var.Ii)) %>%
  rename(`hisp.lm.pr(z>0)` = `Pr(z > 0)`,
         `hisp.lm.Z.Ii` = `Z.Ii`) %>%
  mutate(hisp.clst.dmmy = ifelse(hispcnctn > 1 & 
           `hisp.lm.pr(z>0)`<0.05 & 
           `hisp.lm.Z.Ii` > 0, 1, 0),
         hisp.clst.dmmy = replace_na(hisp.clst.dmmy, 0)) %>%
  select(STATEFP10:TractFIPS, wht.clst.dmmy, blk.clst.dmmy, as.clst.dmmy, hisp.clst.dmmy)

return(df.lm %>% 
         rename(TRTID10 = TRACTID) %>% 
         select(TRTID10, wht.clst.dmmy, blk.clst.dmmy, as.clst.dmmy, hisp.clst.dmmy))
}
