library(raster)
library(feather)
library(tidyverse)
library(stringr)
library(viridis)
library(sf)
library(FedData)
library(mapview)
library(lubridate)
library(feather)
library(tidycensus)
library(USAboundaries)
library(ggthemes)
library(ggplot2)
library(kableExtra)
library(knitr)
library(tmap)
library(tmaptools)
library(RColorBrewer)
library(zoo)
library(nhdplusTools)

city_geom <- get_acs(geography = "metropolitan statistical area/micropolitan statistical area", variables = "B01003_001", year = 2013, geometry = T) %>%
  filter(NAME %in% c("Boulder, CO Metro Area", "Fort Collins, CO Metro Area", "Denver-Aurora-Lakewood, CO Metro Area", "Colorado Springs, CO Metro Area", "Pueblo, CO Metro Area")) %>%
  st_transform(2163)

censu <- read_rds("data/census.rds")

test <- censu %>% 
  st_join(city_geom) %>%
  filter(!is.na(estimate))

c_coli <- st_contains(test, censu) %>%
  unlist()

look <- censu[c_coli, ]

census_fr <- read_rds("data/census_fr.rds")
ft <- look %>%
  mutate(percent_color = (color/population)*100,
         percent_degrees = (degrees/education_25_total)*100,
         percent_renters = (renters/tenure_total)*100) %>%
  select(percent_color, median_income, percent_degrees, percent_renters)


tm_shape(ft) +
  tm_fill(col = 'percent_color') +
  tm_shape(city_geom) +
  tm_borders()


