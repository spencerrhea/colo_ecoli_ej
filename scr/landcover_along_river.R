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
library(nhdplusTools)

### NLCD Data ###
nlcd <- raster("data/nlcd/NLCD_2016.img") 

nc <- us_states(states = "NC") %>%
  st_transform(5070)

nlcd_nc <- crop(nlcd, nc)

nlcd_object <- nlcd_nc

codes <- nlcd_object@data@attributes[[1]][["NLCD.Land.Cover.Class"]]

legend <- as_data_frame(codes) %>%
  mutate(id = row_number()) %>%
  rename(type = value)  %>%
  mutate(land = ifelse(type == "", NA, type)) %>%
  filter(!is.na(land)) %>%
  select(-land) %>%
  mutate(type = as.character(type)) %>%
  mutate(id = id-1)
####

### WWTP Data ###
wastewater <- st_read(dsn = "data/wastewater/CWA_summaries_060314.gdb") %>%
  filter(CWP_STATE == "NC")

waste_water_full <- wastewater %>%
  st_as_sf(coords = c("coords.x1", "coords.x2"), crs = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs") %>%
  filter(CWP_STATE == "NC") 

sewage <- waste_water_full %>%
  filter(SIC_CODES == 4952) %>%
  filter(CWP_PERMIT_STATUS_DESC != "Terminated") %>%
  mutate(type = "plant") %>%
  filter(CWP_PERMIT_STATUS_DESC != "Expired") %>%
  st_transform(2163) 
###


define_land_by_catch <- function(catch, row) {
  
  #transform to the nlcd crs 
  whole_catch <- catch[row,] %>%
    st_transform(5070)
  
  #pull the basin codes
  fetid <- pull(whole_catch, featureid)
  
  #get etent of basin
  catch_extent <- extent(whole_catch)
  
  #crops nlcd to basin 
  nlcd_crop <- crop(nlcd_object, catch_extent) 
  
  nlcd_mask <- mask(nlcd_crop, whole_catch) 
  
  nlcd_trim <- trim(nlcd_mask) 
  
  #gets the numver of pixles in sub basin 
  values <- as.data.frame(getValues(nlcd_trim)) %>%
    rename(id = 1) %>%
    group_by(id) %>%
    summarise(num = n()) %>%
    filter(id != 0) %>%
    full_join(legend) %>%
    filter(!is.na(id))%>%
    mutate(featureid = fetid) %>%
    select(-id) 
  
  #Gets the total number of pixles in a basin
  total <- sum(values$num, na.rm = T)
  
  #pivot the table to wide
  fin <- values %>%
    pivot_wider(names_from = type, values_from = num) %>%
    mutate(Total = total)
}

#This function calculates the percent of landocover for all stream segmants on the main stem of the river 
percecnt_along <- function(river_comid, flow_lines = flowline, watersheds = landcover_by_catchment) {
  
  #gets all sub basins that contribute to a river segmant 
  up_stream <- get_UT(flow_lines, river_comid)
  
  #filter flowlines to only ones contributing to a point on the main stem 
  up_segments <- flow_lines %>%
    filter(comid %in% up_stream)
  
  #filter catchment to only the basin contributing to a point on the main stem 
  catch <- watersheds %>%
    filter(featureid %in% pull(up_segments, comid)) %>%
    pivot_longer(cols = c("Evergreen Forest", "Shrub/Scrub", "Unclassified", "Open Water", "Perennial Snow/Ice", "Developed, Open Space", "Developed, Low Intensity", "Developed, Medium Intensity", "Developed, High Intensity", "Barren Land", "Deciduous Forest", "Mixed Forest", "Herbaceuous", "Hay/Pasture", "Cultivated Crops", "Woody Wetlands", "Emergent Herbaceuous Wetlands")) %>%
    group_by(name) %>%
    summarise(num = sum(value, na.rm = T)) %>%
    mutate(num = as.numeric(num)) 
  
  #gets total number of pixles for each catchment 
  total <- sum(catch$num) 
  
  #gets percent of landcover for each catchment 
  final <- catch %>%
    mutate(Total = total) %>%
    mutate(percent = (num/Total)*100) %>%
    mutate(percent = round(percent, 5)) %>%
    select(name, percent) %>%
    mutate(comid = river_comid) %>%
    pivot_wider(names_from = name, values_from = percent)
}

#This function return the landcover percentages for all subasins 
river_longitudinal <- function(lat, long) {
  
  start <- st_sfc(st_point(c(long, lat)), crs = 4269)
  
  #Gets river segment of where the point is 
  start_comid <- discover_nhdplus_id(start)
  
  #get flowlines for river network
  river_flowline <- navigate_nldi(list(featureSource = "comid", 
                                       featureID = start_comid), 
                                  mode = "upstreamTributaries", 
                                  data_source = "")
  
  #get catchments, waterbodies, and flowline 
  river_subset <- subset_nhdplus(comids = river_flowline$nhdplus_comid,
                                 output_file = "data/nhd/subset.gpkg",
                                 nhdplus_data = "download", 
                                 return_data = TRUE, 
                                 overwrite = T,
                                 simplified = TRUE,
                                 flowline_only = F)
  
  #Pull out all dataframes 
  river_flowline <- river_subset$NHDFlowline_Network
  river_catchment <- river_subset$CatchmentSP
  river_waterbody <- river_subset$NHDWaterbody
  
  #prep dataframe for nhdpulstools 
  river_prep <- prepare_nhdplus(river_flowline, min_network_size = 10, min_path_length = 10) %>%
    select(ID = COMID, toID = toCOMID, length = LENGTHKM)
  
  #get the mainstem flowline codes
  river_main <- get_UM(river_flowline, start_comid)
  
  #get dataframe with only mainstem 
  river <- river_flowline %>%
    filter(comid %in% river_main)
  
  #Get pathlength of river segments 
  river_pathlength <- get_pathlength(river_prep)
  
  #join pathlengths with dataframe 
  river_length <- river_flowline %>%
    rename(ID = comid) %>%
    right_join(river_pathlength, by = "ID")
  
  #join mainstem to pethlengths 
  riv <- river_length %>%
    filter(ID %in% river_main) %>%
    st_transform(4269) %>%
    select(-pathlength.x) %>%
    mutate(pathlength = pathlength.y) %>%
    select(-pathlength.y)
  
  #add id column 
  river_catchment_id <- river_catchment %>%
    mutate(row = row_number())
  
  whole_catch <- river_catchment_id
  
  # Us function to get landcover for all catchments 
  river_landcover_by_catchment <- define_land_by_catch(river_catchment_id, 1)
  
  land_length <- as.numeric(length(river_catchment_id$id))
  
  seqens <- 1:(land_length-1)
  for(i in seqens) {
    x <- i+1
    shed <- define_land_by_catch(river_catchment_id, x)
    
    river_landcover_by_catchment <- rbind(river_landcover_by_catchment, shed)
  }
  
  landcover_by_catchment <- river_landcover_by_catchment
  
  flowline <- river_flowline
  
  #get percentage for all river segments on the mainstem 
  river_comids <- pull(riv, ID)
  
  river_percents <- percecnt_along(river_comids[1], flowline, landcover_by_catchment)
  
  comids_length <- as.numeric(length(river_comids)[1])
  
  seq <- 1:(comids_length-1)
  for(i in seq) {
    new <- percecnt_along(river_comids[i+1], flowline, landcover_by_catchment) 
    
    river_percents <- rbind(new, river_percents)
  }
  
  #define all landcover types 
  
  river_percent_along_river <- river_percents %>%
    rename(barren = 2,
           crops = 3,
           deciduous = 4,
           developed_high = 5,
           developed_low = 6,
           developed_mid = 7,
           developed_open = 8,
           wetland_herbaceous = 9,
           evergreen = 10,
           hay_pasture = 11,
           herbaceous = 12,
           forest_mixed = 13,
           water = 14,
           ice = 15,
           shrub = 16,
           wetland_wood = 18) %>%
    select(-Unclassified)
  
  
  #summarise landcover types 
  river_landcover_sumz <- river_percent_along_river %>%
    mutate(developed = developed_high+developed_low+developed_mid+developed_open,
           forest = forest_mixed+deciduous+evergreen,
           wetland = wetland_herbaceous+wetland_wood,
           agriculture = crops+hay_pasture) %>%
    select(-developed_high, -developed_mid, -developed_low, -forest_mixed, -deciduous, -evergreen, -wetland_herbaceous, -wetland_wood, -crops, -hay_pasture, -developed_open) 
  
  #get laength of river 
  max_length <- max(riv$pathlength)
  
  #change length so the outlet is the largest number 
  riv_trans <- riv %>%
    mutate(pathlength = max_length-pathlength) %>%
    rename(comid = ID) %>%
    st_transform(2163) 
  
  #buffer sewage to join to river 
  wwtp_buff <- sewage %>%
    st_buffer(300)
  
  #get all wwtp on the river 
  river_wwtp <- st_join(riv_trans, wwtp_buff) %>%
    filter(!is.na(CWP_NAME))
  
  #get only codes 
  river_sweage_codes <- river_wwtp %>%
    group_by(REGISTRY_ID) %>%
    summarise(length = mean(pathlength, na.rm = T)) %>%
    pull(length)
  
  river_sweage_id <- river_wwtp %>%
    group_by(REGISTRY_ID) %>%
    summarise(length = mean(pathlength, na.rm = T)) %>%
    pull(REGISTRY_ID)
  
  river_sweage_names <- river_wwtp %>%
    filter(REGISTRY_ID %in% river_sweage_id) %>%
    pull(CWP_NAME)
  
  #join ecoli points to rivers
  #riv_jumps <- st_join(riv_trans, ecoli_mean_sf_100m)
  
  #make dataframe with landcover and ecoli points 
  riv_jumps_fil <- riv_trans %>%
    as.data.frame() %>%
    full_join(river_landcover_sumz, by = "comid") %>%
    select(-geometry) %>%
    filter(!is.na(pathlength)) %>%
    select(pathlength, herbaceous, shrub, developed, forest, wetland, agriculture) %>%
   # mutate(mean = mean/50) %>%
    pivot_longer(cols = c(herbaceous, shrub, developed, forest, wetland, agriculture))
  
  #make a list with both landcover dataframe and wwtp locations 
  fin <- list(full = riv_jumps_fil,
              wwtp = river_sweage_codes,
              wwtp_names = river_sweage_names,
              nlcd_raw = river_percent_along_river,
              flowlines = river_flowline, 
              catchment = river_catchment,
              waterbodies = river_waterbody, 
              mainstem = riv)
  
  return(fin)
  
}

ellerbe <- river_longitudinal(36.059943, -78.812995)

ggplot(ellerbe[["full"]], aes(x = pathlength, y = value, colour = name)) +
  geom_vline(xintercept = ellerbe[["wwtp"]], color = "gray") +
  geom_line() +
  scale_color_manual(values = c("#252525", "#e7298a", "#66a61e", "#7570b3", "#d95f02", "#1b9e77"), 
                     breaks = c("herbaceous", "shrub", "forest", "wetland", "agriculture", "developed"),
                     name = "Landcover") +
  theme_few() +
  labs(x = "Length on River", y = "Watershed landcover (%)")

mapview(ellerbe[["waterbodies"]])
