# This script is an analysis of the relationship between census data publicly 
# available E.coli data in the front range of Colorado. Census data is downloaded 
# using the tidycenbsus package. E.coli data was retrieved by Matt Ross in 2020
# Using the national water quality portal. 
#
# The analysis uses multiple learner regressions to asses is a relationship
# exists between demographic variables (race, income, rental status, and education).
#
# Created by Spencer Rhea (spencer.rhea@duke.edu) ~January 2020. 
# Last Edited 11/26/23

#### Load in Packages ####
library(raster)
library(feather)
library(tidyverse)
library(stringr)
library(ggpubr)
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
library(nhdplusTools)
library(glue)


#### Retrive Census Data ####
# All cenus varibles we tought we might be interested in 
codes <- c(population = "B01003_001",
           median_income = "B19013_001",
           education_25_total = "B15003_001",
           hs_diploma = "B15003_017",
           GED = "B15003_018",
           some_college_1year = "B15003_019",
           some_college = "B15003_020",
           associates = "B15003_021",
           bachelors = "B15003_022",
           masters = "B15003_023",
           professional = "B15003_024",
           doctoral = "B15003_025",
           tenure_total = "B25003_001",
           renters = "B25003_003",
           owners = "B25003_002",
           race_total = "B03002_001",
           white = "B03002_003",
           black = "B03002_004",
           native = "B03002_005",
           asian = "B03002_006",
           hawaiian = "B03002_007",
           other = "B03002_008",
           two_or_more = "B03002_009",
           hispanic = "B03002_012")

# Counties in each metro area are defined by the census bureau 
metro_areas <- list(ft_metro = 'Larimer',
       boulder_metro = 'Boulder',
       denver_metro = c('Denver', 'Arapahoe', 'Jefferson', 'Adams', 'Douglas', 
                  'Broomfield', 'Elbert', 'Park', 'Clear Creek', 'Gilpin'),
       colospgs_metro = 'El Paso',
       pueblo_metro = 'Pueblo')

# Look through each metro area and download data 
for(i in 1:length(metro_areas)) {
    
    # Get all demographic variables without geometries (very  slow to download together)
    block_group_census <- get_acs(geography = "block group", codes, year = 2013, 
                                  state = "CO", county = metro_areas[[i]], geometry = F) 
    # Get census geometries 
    block_group_geom <- get_acs(geography = "block group", 'B01003_001', 
                                year = 2013, state = "CO", county =  metro_areas[[i]], 
                                geometry = T) %>%
        select(GEOID)
    
    # Pivot to wide data frame 
    bg_wide <- block_group_census %>%
        pivot_wider(names_from = c("variable"), values_from = c("estimate", "moe")) %>%
        # calculate people of color (POC) as non-white race catagories in census
        mutate(color = estimate_black+estimate_native+estimate_asian+estimate_hawaiian+estimate_other+estimate_two_or_more+estimate_hispanic,
               # Calculate degree holders as someone with an associates, bachelor, professional, or doctoral degree 
               degrees = estimate_associates+estimate_bachelors+estimate_masters+estimate_professional+estimate_doctoral) %>%
        rename(renters = estimate_renters,
               population = estimate_population,
               tenure_total = estimate_tenure_total,
               education_25_total = estimate_education_25_total,
               median_income = estimate_median_income) %>%
        select(GEOID, NAME, color, population, renters, tenure_total, degrees, education_25_total, median_income)
    
    # Join with geomoetry 
    bg_wide_geom <- inner_join(block_group_geom, bg_wide, by = "GEOID") %>%
        mutate(metro_area = !!names(metro_areas)[i])
    
    # Save file 
    write_rds(bg_wide_geom, file = glue('data_final/census/{m}.rds',
                                        m = names(metro_areas)[i]))
}

#### Load and prep E.coli data ####

# Read in raw e.coli data (downloaded by Matt)
ecoli_raw <- read_feather('data_final/ecoli/Colorado_ecoli.feather')

# Read in inventory of site locations 
inv <- read_feather('data_final/ecoli/wqp_inv.feather') %>%
    select(SiteID=MonitoringLocationIdentifier,
           resultCount,
           parameter,
           lat,
           lon) %>%
    distinct() %>%
    filter(!is.na(lat),
           !is.na(lon)) %>%
    st_as_sf(coords=c('lon','lat'),crs=4326) %>%
    distinct(SiteID, .keep_all = T)

# clean e.coli data
clean_ecoli <-  ecoli_raw %>%
    # select and rename relevent columns 
    dplyr::select(date=ActivityStartDate,
                  parameter=CharacteristicName,
                  units=ResultMeasure.MeasureUnitCode,
                  SiteID=MonitoringLocationIdentifier,
                  org=OrganizationFormalName,
                  org_id=OrganizationIdentifier,
                  time=ActivityStartTime.Time,
                  raw_value=ResultMeasureValue,
                  sample_method=SampleCollectionMethod.MethodName,
                  analytical_method=ResultAnalyticalMethod.MethodName,
                  particle_size=ResultParticleSizeBasisText,
                  date_time=ActivityStartDateTime,
                  media=ActivityMediaName,
                  sample_depth=ActivityDepthHeightMeasure.MeasureValue,
                  sample_depth_unit=ActivityDepthHeightMeasure.MeasureUnitCode,
                  fraction=ResultSampleFractionText,
                  status=ResultStatusIdentifier) %>%
    mutate(units = trimws(units)) %>%
    # Keep only samples that are water samples
    filter(media=='Water') %>%
    filter(units != "ug/l")

# Calculate mean values and filter to select dates 
# CHANGE (comments suggest a change to geometric mean)
ecoli_mean <- clean_ecoli %>%
    # Filter to 5 years before and after census data collections 
    filter(year(date) >= 2005 & year(date) <= 2015) %>%
    filter(parameter == "Escherichia coli") %>%
    group_by(SiteID) %>%
    summarize(mean = mean(raw_value,na.rm=T) %>% round(.,2),
              raw_mean = mean(raw_value,na.rm=T),
              sd = sd(raw_value, na.rm = T),
              count=n()) %>%
    # Only keep sites with more than 9 samples 
    filter(count > 9) %>% 
    # remove less than 1 samples 
    # (Can't remember why we did this, maybe worth removing?)
    mutate(log_mean=log10(mean)) %>%
    filter(log_mean > 0.0000001)

# Set to crs of E.coli data same as census 
ecoli_mean_sf <- inv %>%
    right_join(.,ecoli_mean,by='SiteID') %>%
    st_transform(2163) 



#### Load NHD flowlines #### 
nhd <- st_read('data_final/nhd/nhd_sf.gpkg', "NHDFlowline_Network")

# Set crs to same as census data 
nhd_order <- nhd %>%
    st_transform(2163) %>%
    # Only include somewhat large streams 
    filter(streamorde > 2) 


#remove sites not on river #### 
river_ecoli <- nhd_order %>%
    st_buffer(100) %>%
    st_join(ecoli_mean_sf)

#Filter sites not in the front range ####
census_data <- list.files('data_final/census', full.names = TRUE)

front_range_whole <- map_dfr(census_data, read_rds) %>%
    st_transform(crs = 2163)

# Create a single polygon for front range 
front_range_footprint <- st_union(front_range_whole) 

# Remove sites not on a stream 
river_ecoli_sites <- ecoli_mean_sf %>%
    filter(SiteID %in% !!unique(river_ecoli$SiteID))

# remove sites not in the front range 
fr_ecoli <- st_filter(st_transform(river_ecoli_sites, st_crs(front_range_footprint)), 
                      front_range_footprint) %>%
    st_transform(2163) %>%
    # Create site buffer to join with census data 
    st_buffer(3000)


#### Ecoli linear models ####

# Function to create linear model in each metro area 
census_ecoli <- function(metro_area) {
    
    # read in all metro areas (global) or just one 
    if(metro_area == 'global') {
        fils <- list.files('data_final/census/', full.names = T)
        
        census_raw <-map_dfr(fils, read_rds)
    } else {
        
        metro_path <- glue('data_final/census/{m}.rds',
                           m = metro_area)
        
        census_raw <- read_rds(metro_path) 
    }
    
    census_raw <- st_transform(census_raw, 2163)
    
    # Join census to ecoli data
    colorado_3km_ecoli <- st_join(census_raw, fr_ecoli)
    
    #Removing sites where there is overlap
    grouped_3km <- colorado_3km_ecoli %>%
        filter(!is.na(SiteID)) %>%
        mutate(weighted_income = median_income*population) %>%
        mutate(weighted_mean = mean*count) %>%
        group_by(SiteID, metro_area) %>%
        summarise(color = sum(color, na.rm = T),
                  mean = mean(mean, na.rm = T),
                  weighted_mean = mean(weighted_mean, na.rm = T),
                  mean_n = mean(count, na.rm = T),
                  total = sum(population, na.rm = T),
                  percent_color = (color/total)*100,
                  income = sum(weighted_income, na.rm = T)/sum(total),
                  renters = sum(renters, na.rm = T),
                  tenure_total = sum(tenure_total, na.rm = T),
                  degrees = sum(degrees, na.rm = T),
                  education = sum(education_25_total, na.rm = T),
                  percent_renters = sum(renters, na.rm = T)/sum(tenure_total, na.rm = T),
                  percent_degree = sum(degrees, na.rm = T)/sum(education_25_total, na.rm = T)) %>%
        # Regroup by income
        group_by(income, metro_area) %>%
        summarise(mean = sum(weighted_mean, na.rm = T)/sum(mean_n, na.rm = T),
                  percent_color = sum(color)/sum(total)*100,
                  percent_renters = sum(renters, na.rm = T)/sum(tenure_total, na.rm = T)*100,
                  percent_degrees = sum(degrees, na.rm = T)/sum(education, na.rm = T)*100) 
    
    # Calculate log values (+1 so there are no negative values)
    grouped_3km_lm <- grouped_3km %>%
        mutate(log_color = log(percent_color+1),
               log_renters = log(percent_renters+1),
               log_degrees = log(percent_degrees+1),
               log_mean = log(mean),
               log_income = log(income)) %>%
        select(log_mean, log_color, log_renters, log_degrees, log_income, metro_area)
    
    return(grouped_3km_lm) 
}

#### Model Summaries ####
# Global model 
gloabl_ecoli_census <- census_ecoli(metro_area='global')

# Create each metro areas as a factor 
gloabl_ecoli_census$metro_area <- factor(gloabl_ecoli_census$metro_area)

summary(lm(log_mean ~ log_color+log_renters+log_income+log_degrees+metro_area, gloabl_ecoli_census))

# Fort Collins 
ft_ecoli_census <- census_ecoli(metro_area='ft_metro')
summary(lm(log_mean ~ log_color+log_renters+log_income+log_degrees, ft_ecoli_census))

# Boulder
boulder_ecoli_census <- census_ecoli(metro_area='boulder_metro')
summary(lm(log_mean ~ log_color+log_renters+log_income+log_degrees, boulder_ecoli_census))

# Denver 
denver_ecoli_census <- census_ecoli(metro_area='denver_metro')
summary(lm(log_mean ~ log_color+log_renters+log_income+log_degrees, denver_ecoli_census))

#Colo Spgs
colospgs_ecoli_census <- census_ecoli(metro_area='colospgs_metro')
summary(lm(log_mean ~ log_color+log_renters+log_income+log_degrees, colospgs_ecoli_census))

#Pueblo 
pueblo_metro_ecoli_census <- census_ecoli(metro_area='pueblo_metro')
summary(lm(log_mean ~ log_color+log_renters+log_income+log_degrees, pueblo_metro_ecoli_census))

#### Figures  #### 

## E.coli box plots by metro ####
metro_areas <- c('front_range', 'ft_metro', 'boulder_metro', 'denver_metro', 'colospgs_metro', 'pueblo_metro')

get_city_distri <- function(metro_area) {
    
    if(metro_area == 'front_range') {
        fils <- list.files('data_final/census', full.names = TRUE)
        
        census_raw <- map_dfr(fils, read_rds) %>%
            st_transform(2163)
    } else {
        metro_path <- glue('data_final/census/{m}.rds',
                           m = metro_area)
        
        census_raw <- read_rds(metro_path) %>%
            st_transform(2163)
    }
    
    metro_footprint <- census_raw %>%
        st_union() %>%
        st_as_sf()
    
    colorado_3km_ecoli <- st_filter(fr_ecoli, metro_footprint) %>%
        as_tibble() %>%
        select(-geometry) %>%
        select(mean) %>%
        mutate(metro_area = !!metro_area)
    
}

ft_stats <-tibble()
for(i in 1:length(metro_areas)) {
    area <- get_city_distri(metro_areas[i])
    
    ft_stats <- rbind(ft_stats, area)
}

ft_stats$metro_area <- factor(ft_stats$metro_area,
                       levels = c('boulder_metro', 'ft_metro', 'denver_metro', 'pueblo_metro', 'colospgs_metro', 'front_range'),ordered = TRUE)

# CHANGE (units are incorrect, see comment in manuscript)
save_dist <- ggplot(ft_stats, aes(x = metro_area, y = mean)) +
    geom_boxplot() +
    scale_x_discrete(labels = c('Boulder', 'Fort\nCollins', 'Denver', 'Pueblo', 'Colorado\nSprings', 'Regional')) +
    scale_y_log10() +
    theme_few() +
    geom_hline(yintercept = 235, linetype = 'dashed') +
    labs(x = 'Metro Area', y = 'Mean Ecoli Concentration \nat Sites (#/ml)')

ggsave(filename="figures/ecoli.dist.png", plot=save_dist, device="png",
       height=70, width=110, units="mm", dpi=500)

## Demographics vs E.coli conetrations ####
metro_areas <- list('front_range', 'ft_metro', 'boulder_metro', 'denver_metro', 'colospgs_metro', 'pueblo_metro')

get_city_stats <- function(metro_area) {
    
    if(metro_area == 'front_range') {
        fils <- list.files('data_final/census', full.names = TRUE)
        
        census_raw <- map_dfr(fils, read_rds) %>%
            st_transform(2163)
    } else {
        metro_path <- glue('data_final/census/{m}.rds',
                           m = metro_area)
        
        census_raw <- read_rds(metro_path) %>%
            st_transform(2163)
    }
    
    metro_footprint <- census_raw %>%
        st_union() %>%
        st_as_sf()
    
    colorado_3km_ecoli <- st_filter(fr_ecoli, metro_footprint) 
    
    tibble(ecoli_mean = mean(colorado_3km_ecoli$mean, na.rm = T),
           ecoli_median = median(colorado_3km_ecoli$mean, na.rm = T),
           ecoli_max = max(colorado_3km_ecoli$mean, na.rm = T),
           poc = sum(census_raw$color, na.rm = T)/sum(census_raw$population, na.rm = T),
           income = mean(census_raw$median_income, na.rm = T),
           education = sum(census_raw$degrees, na.rm = T)/sum(census_raw$education_25_total, na.rm = T),
           renters = sum(census_raw$renters, na.rm = T)/sum(census_raw$tenure_total, na.rm = T),
           population = sum(census_raw$tenure_total, na.rm = T),
           density = sum(census_raw$tenure_total, na.rm = T)/(st_area(metro_footprint)/1000000),
           metro_area = metro_area)
    
}

metro_areas <- list('front_range', 'ft_metro', 'boulder_metro', 'denver_metro', 'colospgs_metro', 'pueblo_metro')

ft_stats <- tibble()
for(i in 1:length(metro_areas)) {
    area <- get_city_stats(metro_areas[i])
    
    ft_stats <- rbind(ft_stats, area)
}


gloabl_ecoli_census_non_log <- gloabl_ecoli_census %>%
    mutate(mean = exp(log_mean),
           POC = exp(log_color),
           Renters = exp(log_renters),
           Degrees = exp(log_degrees),
           Income = exp(log_income))

bi_plot <- gloabl_ecoli_census_non_log %>%
    pivot_longer(cols = c('POC', 'Renters', 'Degrees', 'Income')) %>%
    ggplot(aes(mean, value)) +
    geom_point(size = 0.5) +
    facet_wrap(~name, scales = 'free_y') +
    scale_y_log10() +
    scale_x_log10() +
    theme_few() +
    labs(x = 'Mean Ecoli Concentration', y = "Percent of Block Groups or \nMean Median Household Income") +
    geom_vline(xintercept=235, linetype = 'dashed') + 
    theme(axis.title = element_text(size=6),
          axis.text = element_text(size=6), 
          plot.title = element_text(size=6), 
          plot.subtitle = element_text(size=6),
          strip.text.x = element_text(size=6))


ggsave(filename="figures/bi_plot.png", plot=bi_plot, device="png",
       height=2.9, width=2.9, units="in", dpi=300)


## Visualize model outputs #### 
lm_predict <- gloabl_ecoli_census %>%
    select(log_mean, log_color, log_income, log_renters, log_degrees)

lm_multi <- lm(log_mean ~ log_color+log_income+log_renters+log_degrees, lm_predict)


mean(lm_predict$log_color)
mean(lm_predict$log_income)
mean(lm_predict$log_renters)
mean(lm_predict$log_degrees)

x <- 0:100 

#means for all inpputs 
lm_new_data <- tibble(log_color = log(x+1), 
                      log_income = 11.02037,
                      log_renters = 3.236696,
                      log_degrees = 3.754468) 

fake_data <- predict(object = lm_multi,
                     newdata = lm_new_data) 



#Add renters 
lm_new_data_plus_5 <- lm_new_data %>%
    mutate(log_renters = log_renters+log(5))

fake_data_plus_5 <- predict(object = lm_multi,
                            newdata = lm_new_data_plus_5) 

lm_new_data_plus_10 <- lm_new_data %>%
    mutate(log_renters = log_renters+log(10))

fake_data_plus_10 <- predict(object = lm_multi,
                             newdata = lm_new_data_plus_10) 

#Minus renters 
lm_new_data_minus_5 <- lm_new_data %>%
    mutate(log_renters = log_renters-log(5))

fake_data_minus_5 <- predict(object = lm_multi,
                             newdata = lm_new_data_minus_5) 

lm_new_data_minus_10 <- lm_new_data %>%
    mutate(log_renters = log_renters-log(10))

fake_data_minus_10 <- predict(object = lm_multi,
                              newdata = lm_new_data_minus_10) 


fake_data_frame <- lm_new_data %>%
    mutate(prediction_mean = fake_data,
           prediction_plus_5 = fake_data_plus_5,
          # prediction_plus_10 = fake_data_plus_10,
           prediction_minus_5 = fake_data_minus_5) %>%
           #prediction_minus_10 = fake_data_minus_10) %>%
    select(log_color, prediction_mean, prediction_plus_5, prediction_minus_5) %>%
    pivot_longer(cols = c(prediction_mean, prediction_plus_5, prediction_minus_5))

example_lm <- ggplot(fake_data_frame, aes(x = exp(log_color), y = exp(value), colour = name)) +
    geom_line() +
    labs(x = "Percent People of Color", y = "E.coli concentration (#/ml)") +
    scale_colour_manual(name = "Percent Renters \nChange From Mean",
                          breaks = c("prediction_mean", "prediction_plus_5", "prediction_minus_5"), 
                          labels = c("Mean", "Plus %5", "Minus %5"),
                          values = c('#E41A1C', '#377EB8', '#4DAF4A')) +
    theme_few() +
    geom_hline(yintercept = 235, linetype = 'dashed') +
    theme(legend.position = c(0.32, 0.77), legend.background = element_rect(colour = NA, fill = NA),
          axis.title = element_text(size=6),
          axis.text = element_text(size=6), 
          plot.title = element_text(size=6), 
          plot.subtitle = element_text(size=6),
          strip.text.x = element_text(size=6)) 


ggsave(filename="figures/example_lm.png", plot=example_lm, device="png",
       height=2.9, width=2.9, units="in", dpi=300)

## Map of census data #### 

fils <- list.files('data_final/census', full.names = TRUE)

fr_census <- tibble()
for(i in 1:length(fils)) {
    
    site <- read_rds(fils[i])
    
    fr_census <- rbind(fr_census, site)
}

fr_census <- fr_census %>%
    mutate(POC = (color/population)*100,
           Renters = (renters/tenure_total)*100,
           'Degree Holders' = (degrees/education_25_total)*100)  %>%
    rename('Median Income' = median_income)

fr_census_ft <- fr_census %>%
    group_by(metro_area) %>%
    summarise(geometry =  st_union(geometry))

fr_census[is.na(fr_census)] <- 0

city_loc <- tibble(city = c('fort_collins', 'boulder', 'denver', 'colo_spgs', 'pueblo'),
    lat = c('40.585991','40.019313','39.751563','38.833329','38.269811'),
    long = c('-105.075056', '-105.273481', '-104.987354', '-104.822678', '-104.612052')) %>%
    st_as_sf(coords = c('long', 'lat'), crs = 4623)

map_poc <- tm_shape(fr_census) +
    tm_fill(col = 'POC', palette = 'Blues') +
    tm_shape(fr_census_ft) +
    tm_polygons(alpha = 0, border.col = '#000000', border.alpha = 1) +
    tm_shape(city_loc) +
    tm_dots(size=0.1) +
    tm_legend(legend.bg.color = '#FFFFFF', legend.frame = T)

map_renters <- tm_shape(fr_census) +
    tm_fill(col = 'Renters', palette = 'PuRd') +
    tm_shape(fr_census_ft) +
    tm_polygons(alpha = 0, border.col = '#000000', border.alpha = 1) +
    tm_shape(city_loc) +
    tm_dots(size=0.1) +
    tm_legend(legend.bg.color = '#FFFFFF', legend.frame = T)

map_edu <- tm_shape(fr_census) +
    tm_fill(col = 'Degree Holders', palette = 'BuPu') +
    tm_shape(fr_census_ft) +
    tm_polygons(alpha = 0, border.col = '#000000', border.alpha = 1) +
    tm_shape(city_loc) +
    tm_dots(size=0.1) +
    tm_legend(legend.bg.color = '#FFFFFF', legend.frame = T)

map_income <- tm_shape(fr_census) +
    tm_fill(col = 'Median Income', palette = 'BuGn', breaks = c(0, 50000, 100000, 150000, 200000)) +
    tm_shape(fr_census_ft) +
    tm_polygons(alpha = 0, border.col = '#000000', border.alpha = 1) +
    tm_shape(city_loc) +
    tm_dots(size=0.1) +
    tm_legend(legend.bg.color = '#FFFFFF', legend.frame = T)

figure_1 <- tmap_arrange(map_poc, map_income, map_edu, map_renters)

tmap_save(figure_1, filename = 'figures/figure_1.png', 
          width = 5.24, height = 6.55, dpi = 500)

## Site map ####
# Read in all rivers (retrieved using NHDplustools with unknown script, i.e. 
# I cant find it)
paltte<-read_rds('data_out/platte_shp.rds')
poudre<-read_rds('data_out/poudre.rds')[["mainstem"]]
bg<-read_rds('data_out/big_thompson.rds')[["mainstem"]]
vrain<-read_rds('data_out/st_vrain.rds')[["mainstem"]]
boulder<-read_rds('data_out/boulder_creek.rds')[["mainstem"]]
monument<-read_rds('data_out/monument_creek.rds')[["mainstem"]]
fountain<-read_rds('data_out/fountain.rds')[["mainstem"]]
arkanasa<-read_rds('data_out/arkansas.rds')[["mainstem"]]

cherry_creek <- read_rds('data_out/cherry_creek.rds') %>%
    st_make_valid()

cherry_creek_2 <- read_rds('data_out/cherry_creek_2.rds') %>%
    st_make_valid()

clear_creek <- read_rds('data_out/clear_creek.rds') %>%
    filter(comid != 2884976)

clear_creek_2 <- read_rds('data_out/clear_creek_2.rds') %>%
    filter(comid != 2885002,
           comid != 2885048,
           comid != 2885668,
           comid != 2885662)

bear_creek <- read_rds('data_out/brear_creek.rds') %>%
    filter(comid != 188153) %>%
    st_make_valid()

bear_creek_3 <- read_rds('data_out/bear_3.rds') %>%
    filter(comid != '188161') %>%
    st_make_valid()


# Get elevation 
dem <- elevatr::get_elev_raster(locations = st_cast(fr_census_ft, "MULTIPOLYGON"), z = 7)

extent <- extent(st_buffer(fr_census_ft, 1))

extent[3] <- extent[3] + .75
extent[4] <- extent[4] - .75

croped <- crop(dem, extent)

masked <- trim(croped)

png(file="rivers.png", width = 400, height = 500)
tm_shape(masked) +
    tm_raster(palette = 'Greys', title = 'Elevtion (m)') +
    tm_shape(fr_census_ft) +
    tm_polygons(alpha = 0, border.col = '#000000', border.alpha = 1) +
    tm_shape(paltte) +
    tm_lines(col = 'blue') +
    tm_shape(poudre) +
    tm_lines(col = 'blue') +
    tm_shape(bg) +
    tm_lines(col = 'blue') +
    tm_shape(vrain) +
    tm_lines(col = 'blue') +
    tm_shape(boulder) +
    tm_lines(col = 'blue') +
    tm_shape(monument) +
    tm_lines(col = 'blue') +
    tm_shape(fountain) +
    tm_lines(col = 'blue') +
    tm_shape(arkanasa) +
    tm_lines(col = 'blue') +
    tm_shape(city_loc) +
    tm_dots(size=0.3, col = '#BF1515') +
    tm_legend(legend.bg.color = '#FFFFFF', legend.frame = T, 
              legend.position = c(0.01, 0.01), legend.text.size = 1)
dev.off()

## Readlining comparison maps #### 

# Denver 
# get denver census data 
den_census <- read_rds('data_final/census/denver_metro.rds')  %>%
    mutate(POC = (color/population)*100,
           Renters = (renters/tenure_total)*100,
           'Degree Holders' = (degrees/education_25_total)*100)  %>%
    rename('Median Income' = median_income)

den_footprint <- den_census %>%
    summarise(geometry =  st_union(geometry)) 

# get denver e.coli data 
den_ecoli <- st_filter(st_transform(river_ecoli_sites, 4269), den_footprint)

den_census_c <- den_census %>%
    mutate(keep = grepl('Denver County', NAME)) %>%
    filter(keep == TRUE) %>%
    mutate(remove_ = ifelse(GEOID %in% c(080010084021, 080050071013, 080050071014, 080050071012, 080050071011,
                                           080010084022, 080010084012, 080050071061, 080050071031, 080010084011,
                                           080590120584), TRUE, FALSE)) %>%
    filter(remove_ == FALSE)

bigft <- st_buffer(den_census_c, .1) %>%
    st_union() %>%
    st_as_sf()

den_ecoli <- st_filter(st_transform(river_ecoli_sites, 4269), bigft) %>%
    arrange(mean) %>%
    filter(mean <= 3000) %>%
    mutate('Ecoli Conc.'=mean)



denver_redlining_map <- tm_shape(den_census, bbox=bb(bigft)) +
    tm_fill(col = 'POC', palette = 'Reds', colorNA = NULL, title = '% POC') +
    tm_shape(paltte) +
    tm_lines(vol = 'Blue') +
    tm_shape(cherry_creek, size = 1) +
    tm_lines(vol = 'Blue') +
    tm_shape(cherry_creek_2) +
    tm_lines(vol = 'Blue') +
    tm_shape(clear_creek) +
    tm_lines(vol = 'Blue') +
    tm_shape(clear_creek_2) +
    tm_lines(vol = 'Blue') +
    tm_shape(bear_creek) +
    # tm_lines(vol = 'Blue') +
    # tm_shape(bear_creek_2) +
    tm_lines(vol = 'Blue') +
    tm_shape(bear_creek_3) +
    tm_lines(vol = 'Blue') +
    tm_shape(den_ecoli) +
    tm_symbols(col = 'mean', size = 0.75, 
               palette = 'Blues', 
               border.col = '#000000',
               border.lwd = 1, 
               border.alpha = 1) +
    tm_legend(legend.bg.color = '#FFFFFF', frame = TRUE) +
    tm_scale_bar() +
    tm_compass(position = c(0.02,0.82))

tmap_save(denver_redlining_map, filename = 'figures/denver_redlining.png', 
          width = 3.11, height = 2.86, dpi = 1500)

# Pueblo
puebloe_census <- read_rds('data_final/census/pueblo_metro.rds')  %>%
    mutate(POC = (color/population)*100,
           Renters = (renters/tenure_total)*100,
           'Degree Holders' = (degrees/education_25_total)*100)  %>%
    rename('Median Income' = median_income)

pueblo_footprint <- puebloe_census %>%
    summarise(geometry =  st_union(geometry)) 

pueblo_ecoli <- st_filter(st_transform(river_ecoli_sites, 4269), pueblo_footprint)

pueb_ecoli <- st_filter(st_transform(river_ecoli_sites, 4269), pueblo_ecoli) %>%
    arrange(mean) %>%
    filter(mean <= 1500) %>%
    mutate('Ecoli Conc.'=mean) %>%
    filter(SiteID != '21COL001_WQX-7520',
           SiteID != 'USGS-07109500')

fountain_creek <- read_rds('data_out/fountain.rds')[["mainstem"]]
arkanasa <- read_rds('data_out/arkansas.rds')[["mainstem"]]

pueblo_readlining <- tm_shape(puebloe_census, bbox=bb(pueb_ecoli, ext = 1.2)) +
    tm_fill(col = 'POC', palette = 'Reds', colorNA = NULL, title = '% POC') +
    tm_shape(arkanasa) +
    tm_lines(vol = 'Blue') +
    tm_shape(fountain_creek) +
    tm_lines(vol = 'Blue') +
    tm_shape(pueb_ecoli) +
    tm_symbols(col = 'Ecoli Conc.', size = 0.65, palette = 'Blues', border.col = '#000000',
               border.lwd = 1, border.alpha = 1) +
    tm_legend(legend.bg.color = '#FFFFFF', frame = TRUE, title.size = 0.8, legend.hist.size = 0.5, 
              panel.label.size = 0.5, main.title.size = 0.7, legend.text.size = 0.5) +
    tm_scale_bar(position = c(0.05,0.01)) +
    tm_compass(position = c(0.02,0.83))

tmap_save(pueblo_readlining, filename = 'figures/pueblo_readlining.png', 
          width = 3.11, height = 2.86, dpi = 1500)
# Colo Spgs 

colopgs_census <- read_rds('data_final/census/colospgs_metro.rds')  %>%
    mutate(POC = (color/population)*100,
           Renters = (renters/tenure_total)*100,
           'Degree Holders' = (degrees/education_25_total)*100)  %>%
    rename('Median Income' = median_income)

colospgs_footprint <- colopgs_census %>%
    summarise(geometry =  st_union(geometry)) 

colospgs_ecoli <- st_filter(st_transform(river_ecoli_sites, 4269), colospgs_footprint)

colospgs_ecoli <- st_filter(st_transform(river_ecoli_sites, 4269), colospgs_ecoli) %>%
    arrange(mean) %>%
   # filter(mean <= 1500) %>%
    mutate('Ecoli Conc.'=mean) %>%
    filter(SiteID != 'USGS-383347104373401',
           SiteID != 'USGS-07106000',
           SiteID != 'USGS-383854104413601')

fountain_creek <- read_rds('data_out/fountain.rds')[["mainstem"]]
monument <- read_rds('data_out/monument_creek.rds')[["mainstem"]]

png(file="colospgs_red.png", width = 300, height = 300)
tm_shape(colopgs_census, bbox=bb(colospgs_ecoli, ext = 1.2)) +
    tm_fill(col = 'POC', palette = 'Reds', colorNA = NULL, title = '% POC') +
    tm_shape(monument) +
    tm_lines(vol = 'Blue') +
    tm_shape(fountain_creek) +
    tm_lines(vol = 'Blue') +
    tm_shape(colospgs_ecoli) +
    tm_symbols(col = 'Ecoli Conc.', size = 0.75, palette = 'Blues', border.col = '#000000',
               border.lwd = 1, border.alpha = 1) +
    tm_legend(legend.bg.color = '#FFFFFF', frame = TRUE) +
    tm_scale_bar() +
    tm_compass(position = c(0.02,0.85))
dev.off()


## Landcover along river plus e.coli #### 


poudre <- read_rds("data_out/poudre.rds")

poudre[["full"]] <- poudre[["full"]] %>%
     filter(name %in% c('forest', 'agriculture', 'developed'))

 
poudre_landcover <- ggplot(poudre[["full"]], aes(x = pathlength, y = value, colour = name)) +
    # geom_vline(xintercept = poudre[["wwtp"]], color = "black", linetype = 'dashed') +
    geom_line() +
    geom_point(aes(x = pathlength, y = mean)) +
    scale_color_manual(values = c("#d95f02", "#e7298a", "#66a61e", "#7570b3", "#000000", "#1b9e77"), 
                       breaks = c("herbaceous", "shrub", "forest", "wetland", "agriculture", "developed"),
                       name = "Landcover") +
    scale_x_continuous(breaks = c(0, 100, 200)) +
    scale_y_continuous(sec.axis = sec_axis(trans = ~.*50,  name = NULL)) +
    theme_few() +
    theme(legend.position = 'bottom') +
    labs(x = '', y = "Watershed landcover (%)") +
    theme(axis.title.y.right = element_blank(),
          axis.text.y.right = element_blank(), 
          axis.ticks.y.right = element_blank(),
          plot.title = element_text(size=6), 
          plot.subtitle = element_text(size=6),
          strip.text.x = element_text(size=6))

ggsave('figures/poudre_lc.png', poudre_landcover, width = 2, height = 2, 
       dpi = 1000, units = 'in')

south_platte <- read_rds("data_out/south_platte.rds")

sp_sweage_codes <- c(281.122, 281.122, 188.124, 188.124, 178.222, 178.222, 178.222, 178.222, 178.222, 178.222, 157.964, 157.964, 157.964, 157.964, 157.964, 157.964)

south_platte <- south_platte %>%
    filter(name %in% c('forest', 'agriculture', 'developed'))

paltte_landcover <- ggplot(south_platte, aes(x = pathlength, y = value, colour = name)) +
    # geom_vline(xintercept = sp_sweage_codes, linetype = 'dashed') +
    geom_line() +
    geom_point(aes(x = pathlength, y = mean)) +
    scale_color_manual(values = c("#d95f02", "#e7298a", "#66a61e", "#7570b3", "#000000", "#1b9e77"), 
                       breaks = c("herbaceous", "shrub", "forest", "wetland", "agriculture", "developed"),
                       name = "Landcover") +
    scale_y_continuous(sec.axis = sec_axis(trans = ~.*50,  name = "E.coli (#/ml)")) +
    theme_few() +
    theme(legend.position = 'none') +
    labs(x = "Length on River (km)", y = "Watershed landcover (%)") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          plot.title = element_text(size=6), 
          plot.subtitle = element_text(size=6),
          strip.text.x = element_text(size=6))

ggsave('figures/paltte_lc.png', paltte_landcover, width = 2, height = 2, 
       dpi = 1000, units = 'in')

fountaun <- read_rds("data_out/fountain.rds")

fountaun[["full"]] <- fountaun[["full"]] %>%
    filter(name %in% c('forest', 'agriculture', 'developed'))

fountain_lancvoer <- ggplot(fountaun[["full"]], aes(x = pathlength, y = value, colour = name)) +
    # geom_vline(xintercept = fountaun[["wwtp"]], color = "black", linetype = 'dashed') +
    geom_line() +
    geom_point(aes(x = pathlength, y = mean)) +
    scale_color_manual(values = c("#d95f02", "#e7298a", "#66a61e", "#7570b3", "#000000", "#1b9e77"), 
                       breaks = c("herbaceous", "shrub", "forest", "wetland", "agriculture", "developed"),
                       name = "Landcover") +
    scale_y_continuous(sec.axis = sec_axis(trans = ~.*50,  name = "E.coli (#/ml)")) +
    theme_few() +
    labs(x = '') +
    #theme(legend.position = 'bottom') +
    theme(legend.position = 'none') +
    theme(axis.title.y.left = element_blank(),
          axis.text.y.left = element_blank(), 
          axis.ticks.y.left = element_blank(),
          plot.title = element_text(size=6), 
          plot.subtitle = element_text(size=6),
          strip.text.x = element_text(size=6))


ggsave('figures/fountain_lc.png', fountain_lancvoer, width = 2, height = 2, 
       dpi = 1000, units = 'in')

all_landcover <- ggarrange(poudre_landcover, paltte_landcover, fountain_lancvoer, 
                           nrow = 1, legend = 'bottom', common.legend = TRUE)
ggsave('figures/landcover_3.png', all_landcover, width = 6, height = 3, 
       dpi = 1000, units = 'in')
