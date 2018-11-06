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

#### Read in 1980 decennial census counts ####
pop1980 <- read_csv("LTDB_Std_1980_fullcount.csv") %>% 
  select(TRTID10:NHBLK80, ASIAN80, HISP80) %>%
  rename(TRACTID = TRTID10) %>%
  mutate(TRACTID = str_pad(TRACTID, width = 11, side = "left", pad = "0")) %>% # add leading zero to get common width of 11
  mutate(StateFIPS = substr(TRACTID, 0, 2),
         CountyFIPS = substr(TRACTID, 3, 5),
         TractFIPS = substr(TRACTID, 6, 11)) %>%
  mutate(prop.nhwht80 = NHWHT80/POP80, # proportion of whites relative to total population in each tract
         prop.nhblk80 = NHBLK80/POP80, # proportion of blacks relative to total population in each tract
         prop.asian80 = ASIAN80/POP80,
         prop.hisp80 = HISP80/POP80)

# get the average proportion of each ethnic group for each county
avg.county.pop1980 <- pop1980 %>%
  select(TRACTID, StateFIPS, CountyFIPS, prop.nhwht80:prop.hisp80) %>%
  group_by(StateFIPS, CountyFIPS) %>%
  summarise(mean(prop.nhwht80), 
            mean(prop.nhblk80), 
            mean(prop.asian80), 
            mean(prop.hisp80))

pop1980.foranalysis <- pop1980 %>%
  left_join(avg.county.pop1980, by = c("StateFIPS", "CountyFIPS")) %>%
  mutate(whitestval = prop.nhwht80 - `mean(prop.nhwht80)`, # calculate standard values that I will use to find the higher than average concnetrations of whites in clusters
         blakstval = prop.nhblk80 - `mean(prop.nhblk80)`,
         asianstval = prop.asian80 - `mean(prop.asian80)`,
         hispstval = prop.hisp80 - `mean(prop.hisp80)`) %>%
  select(TRACTID, StateFIPS, CountyFIPS, TractFIPS, 
         prop.nhwht80:prop.hisp80,
         whitestval:hispstval) %>%
  mutate(TRACTID = as.numeric(TRACTID))

#### Read in 2010 census tract shapefiles ####
if(!file.exists("2010tractshape.Rdata")){
  tract.shape.2010 <- readOGR(dsn = "TractS2010",
                              layer = "tract2010") %>%
  save(tract.shape.2010, file = "2010tractshape.Rdata")
}

load("2010tractshape.Rdata") 

tract.shape.clean <- tract.shape.2010[!tract.shape.2010@data$STATEFP10==72,] # filter out puerto rico 
xy <- coordinates(tract.shape.clean) # retrieve coordinates of tracts to calculate spatial distance 
nb.d5 <- dnearneigh(xy, 0, 5, longlat = TRUE) # distance-based spatial influence at 5 kilometers
nb.d5.list <- nb2listw(nb.d10, style = "W", zero.policy = TRUE) # create spatial weights for neighbors list. This is needed to calculate local moran's i





## test calculating local moran's i with the 2010 shape file king county
king.shape.2010 <- tract.shape.clean[tract.shape.clean@data$COUNTYFP10=="033" & tract.shape.clean@data$STATEFP10=="53", ]
king.xy.2010 <- coordinates(king.shape.2010)
king.nb.d5.2010 <- dnearneigh(king.xy.2010, 0, 5, longlat = TRUE) # distance-based spatial influence at 5 kilometers
king.nb.d5.list <- nb2listw(king.nb.d5.2010, style = "W", zero.policy = TRUE) # create spatial weights for neighbors list. This is needed to calculate local moran's i
king.lm.2010 <- localmoran(king.combined$NHWHT80, king.nb.d5.list, na.action = na.pass) # calculate local moran's i

king.combined <- king.shape.2010@data %>% 
            mutate(TRACTID = as.numeric(as.character(TRACTID))) %>% 
            left_join(pop1980, by = "TRACTID") %>%
            cbind(king.lm.2010) %>%
            mutate(whitestval = NHWHT80 - mean(NHWHT80))
            mutate(whiteclusterdummy = whitestval>0 & `Pr(z > 0)`<0.05)

## MIGHT NEED TO TRY THE TEST ABOVE AGAIN WITH WA (NOT JUST KING COUNTY) AND THEN GROUP BY COUNTY AND CALCULATE THE AVERAGE POPULATION FOR COUNTY
## WILL NEED TO USE PROPORTION OF EACH RACIAL GROUP (I use a measure of local spatial clustering (the local Moran’s I) to characterize the patterns of spatial autocorrelation in the area, based on the proportion of the Asian ethnic group in each tract relative to the county as a whole (Anselin 1995))
##  an ethnic neighborhood is a “hot spot,” composed of a focal census tract with high ethnic concentration compared with the mean concentration in the county and, if present, all similarly high contiguous tracts. 

wa.shape.2010 <- tract.shape.clean[tract.shape.clean@data$STATEFP10=="53", ]
wa.nb.q <- poly2nb(wa.shape.2010, queen = TRUE) # adjacency neighbors list with queen's definition
wa.nb.list <- nb2listw(wa.nb.q, style="W", zero.policy=TRUE) # what i need for local moran's (adjacency spatial weights matrix)
            
# create data frame to use to calculate local moran's i
wa.foranalysis <- wa.shape.2010@data %>% 
              mutate(TRACTID = as.numeric(as.character(TRACTID))) %>%
              left_join(pop1980.foranalysis, by = "TRACTID") %>%
              cbind(localmoran(wa.foranalysis$prop.nhwht80, wa.nb.list, na.action = na.exclude))

# calculate local moran's i for each racial group
wa.nhwht.lm <- localmoran(wa.foranalysis$prop.nhwht80, wa.nb.list, na.action = na.exclude) # calculate local moran's i for white
wa.nhblk.lm <- localmoran(wa.foranalysis$prop.nhblk80, wa.nb.list, na.action = na.pass) # calculate local moran's i for black
wa.asian.lm <- localmoran(wa.foranalysis$prop.asian80, wa.nb.list, na.action = na.pass) # calculate local moran's i for asian
wa.hisp.lm <- localmoran(wa.foranalysis$prop.hisp80, wa.nb.list, na.action = na.pass) # calculate local moran's i for hispanic



# need to now cbind the local moran's i to the data that we have for stvals
wa.lm.combined <- wa.foranalysis %>%
  cbind(wa.nhwht.lm)
              cbind(king.lm.2010)
              mutate(whitestval = NHWHT80 - mean(NHWHT80))
            mutate(whiteclusterdummy = whitestval>0 & `Pr(z > 0)`<0.05)



king.shape.1980 <- tract.shape.clean[tract.shape.clean@data$COUNTY=="King" & tract.shape.clean@data$STATE=="Washington", ]
plot(king.shape.1980)

king.nb.q <- poly2nb(king.shape.1980, queen = TRUE)
king.nb.list <- nb2listw(king.nb.q, style="W", zero.policy=NULL) # what i need for local moran's (adjacency spatial weights matrix)

#### Create spatial weights matrix ####
# Find contiguous neighbors (queen's)
nb_q <- poly2nb(tract.shape.1980, queen = TRUE) ## create neighbors list from polygon list based on regions with contiguous boundaries
attr(nb_q, "region.id")[1:10]
is.symmetric.nb(nb_q)
# Create a neighbors list with spatial weights (adjacency with 1/n = neighbor divided by neighbors adjacent, 0 = not neighbor)
nb_list <- nb2listw(nb_q, style="W", zero.policy=FALSE) # create spatial weights for neighbors list. This is needed to calculate local moran's i
nb_list$style
### use neighbors list for local moran's i? Is the nb_list object an adjacency spatial weights matrix?
xy <- coordinates(king.shape.1980) # Retrieve coordinates of tracts

### how should i save these objects so that i don't have to wait 5-10 mins to run them everytime i close the workspace?
# save all the things that take a long time to rds

### are the following how to distance-based and nearest neighbors spatial weights matrixes?
### should i be using these thresholds?
# Distance-based spatial influence
nb.d10 <- dnearneigh(xy, 0, 5, longlat = TRUE)
nb.d10.list <- nb2listw(nb.d10, style = "W", zero.policy = TRUE)
nb.d25 <- dnearneigh(xy, 0, 25, longlat=TRUE)
# Nearest neighbors spatial influence
nb.k3 <- knn2nb(knearneigh(xy, k=3, RANN=FALSE))
nb.k6 <- knn2nb(knearneigh(xy, k=6, RANN=FALSE))
## Warning in knearneigh(xy, k = 6, RANN = FALSE): k greater than one-third of
## the number of data points

# adjacency spatial weights matrix
king.lm <- king.shape.1980@data$PASIAN80 %>%
  cbind(stval=. - mean(.), localmoran(., king.nb.list)) %>%
  data.frame %>%
  mutate(GISJOIN = king.shape.1980@data$GISJOIN) %>%
  left_join(king.shape.1980@data) %>%
  mutate(asianclusterdummy = stval>0 & `Pr.z...0.`<0.05)

# higher than average and statistically significant
plot(king.shape.1980[king.lm[,"stval"] > 0 & king.lm[,"Pr.z...0."] < .05,])

ggplot(king.lm, aes(stval, group=asianclusterdummy, fill=asianclusterdummy)) +
  geom_density(alpha=.2)

# check the distance spatial weights matrix
king.lm.distance <- king.shape.1980@data$PASIAN80 %>%
  cbind(stval=. - mean(.), localmoran(., nb.d10.list)) %>%
  data.frame %>%
  mutate(GISJOIN = king.shape.1980@data$GISJOIN) %>%
  left_join(king.shape.1980@data)

plot(king.shape.1980[king.lm.distance[,"stval"] > 0 & king.lm.distance[,"Pr.z...0."] < .05,])

# makes more sense to distance spatial weights matrix 
# justify using 5km = distance people would travel to go to a community center or grocery store in their neighborhood
# average household primarily shopped at a store 3.79 miles from home (https://www.ers.usda.gov/amber-waves/2015/august/most-us-households-do-their-main-grocery-shopping-at-supermarkets-and-supercenters-regardless-of-income/)
# 3.79 miles = 6.10 km
