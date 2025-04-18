############################################################################################
# Extracting IMERG Precipitation Values to GPS movement data                               #
# Jayden Skelly, 4/16/2025                                                                 #
# README                                                                                   #
# This script assumes your final version data downloaded from                              #
# https://arthurhouhttps.pps.eosdis.nasa.gov/gpmdata/ is stored locally in                 #
# a YEAR_MO directory for each year-month. The script calculates the total precip          #
# rate at 30 min in the duration of each fix rate, in this case 1 hour. It then            #
# extracts these values onto the data frame.                                               #
# To speed the process, all tiffs are placed in one spat raster in memory.                 #
# Second, I extract wind U and V from a downloaded .nc file from                           # 
# https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download    #
# I limit the extent for faster computation when requesting the data download.             #
# In this example I use pronghorn gps data                                                 #
# from Montana, filtered to winter (Dec 01 - March 31), for 700 individuals.               #
############################################################################################

#### Packages
require(terra)
require(sf)
require(dplyr)
require(ggplot2)
require(lubridate)
require(openair)
require(circular)
require(MoveTools)
### Directory to Analysis folder
setwd("C:/Users/jaydenskelly/Documents/OneDrive - Montana State University/Classes/RemoteS_525/Project/Analysis")

### Filepaths
windpath <- list.files(file.path("./input/"), full.names = TRUE)[1] #Sole .nc file with total range dates, for every hourly time step, stored locally
IMERG_head <- "C:/Users/jaydenskelly/Documents/IMERG_pull_outsideOD" #Where IMERG data is stored
data <- readRDS("C:/Users/jaydenskelly/Documents/OneDrive - Montana State University/Data/PH_GPSLocsClean.rds") #animal movement data
colnames(data)
head(data)
data <- sample_frac(data,0.25)
data <- filter(data, Year > 2020 )
### Prep data by subsetting winter data, adding spatial components
data <- data %>% filter(Season == "Winter")
data <- st_as_sf(data, coords=geometry, dim="XY", crs=4326)
#Note. Ensure you have timestamp in GMT....
#data <- mutate(data, DateTimeGMT = strptime(timestamp, format="%Y-%m-%d %H:%M") %>% as.POSIXct(tz ="UTC") )
# already GMT in "DateTimeGMT"

###  Calculate and extract total precip rate in previous hour for each timestamp (due to hourly fix rate)
#prep data, round to nearest hour, get timestamp into correct format to match to IMERG file naming
data <- data %>% mutate(round_DateTimeGMT = round_date(DateTimeGMT, unit = "30 minute"))
data <- data %>% mutate(year_month = format(as.POSIXct(round_DateTimeGMT), format = "%Y_%m"))

#add start time and end time for first 30 mins in hour period, but formatted 1 hour lag
data <- data %>%
  mutate(
    round_DateTimeGMT_1 = as.POSIXct(round_DateTimeGMT),
    S_time = format(round_DateTimeGMT - minutes(60), format = "%Y%m%d-S%H%M%S"),
    E_time = format(round_DateTimeGMT - minutes(30) - seconds(1), format = "-E%H%M%S")
  ) %>%
  mutate(
    round_DateTimeGMT_1 = paste0(S_time, E_time)
  ) %>%
  select(-S_time, -E_time)
#add start time and end time for second 30 mins in hour period, but formatted 1 hour lag 
data <- data %>%
  mutate(
    round_DateTimeGMT_2 = as.POSIXct(round_DateTimeGMT),
    S_time = format(round_DateTimeGMT - minutes(30), format = "%Y%m%d-S%H%M%S"),
    E_time = format(round_DateTimeGMT - seconds(1), format = "-E%H%M%S")
  ) %>%
  mutate(
    round_DateTimeGMT_2 = paste0(S_time, E_time)
  ) %>%
  select(-S_time, -E_time)
#arrange for viewing
data <- data %>%  arrange(round_DateTimeGMT)
#iterate through each unique hourly period in data, add two files from last 30min in from memory, extract and then sum
#then store these values as a list for each group (day-hour) and add back to main dataframe.
grouped_data <- data %>% group_by(round_DateTimeGMT)
extracted_values <- list()

for (name in unique(grouped_data$round_DateTimeGMT)) {
  group <- grouped_data %>% filter(round_DateTimeGMT == name)
  print(paste("Group:", name))
  value <- NULL
  try({
    # Find matching rasters to import
    files1 <- list.files(path = paste0(IMERG_head, "./imerg_late_30min_12022024/"), pattern = group$round_DateTimeGMT_1, full.names = TRUE, recursive = TRUE)
    raster1 <- rast(files1[1])
    files2 <- list.files(path = paste0(IMERG_head, "./imerg_late_30min_12022024/"), pattern = group$round_DateTimeGMT_2, full.names = TRUE, recursive = TRUE)
    raster2 <- rast(files2[1])
    # Extract values
    values1 <- terra::extract(raster1, st_transform(group, crs = st_crs(raster1)), ID = FALSE)[, 1]
    values2 <- terra::extract(raster2, st_transform(group, crs = st_crs(raster1)), ID = FALSE)[, 1]
    # Sum values
    value <- values1 + values2
  }, silent = TRUE)
  print(value)
  extracted_values[[as.character(name)]] <- value
}

#Finally take the extracted data in each group and add them back to the main dataframe
#As the == is a pass between two lubridate values they are matched, but in different timezones.
for (name in names(extracted_values)) {
  name_posix <- as.POSIXct(as.numeric(name, origin="1970-01-01", tz="GMT"))

  data$prcp_target[grouped_data$round_DateTimeGMT == name_posix] <- extracted_values[[name]]
}

###Save RDS File as backup
saveRDS(data, "data.rds")

### Add wind speed u and v components to dataframe
wind <- rast(windpath)
wind
#subset into u and v components as individual rasters
windu <- wind["u10"] 
windv <- wind["v10"]
#iterate through groups, for each group add wind u and v as columns



grouped_data <- data %>% group_by(round_DateTimeGMT)
extracted_values_u <- list()
extracted_values_v <- list()
for (name in unique(grouped_data$round_DateTimeGMT)) {
  group <- grouped_data %>% filter(round_DateTimeGMT == name)
  print(paste("Group:", name))
  value <- NULL
  try({
    unix <- as.character(name)
    wind1 <- windu[[paste0('u10_valid_time=',unix)]]
    wind2 <- windv[[paste0('v10_valid_time=',unix)]]
    wind1
    # Extract values
    valuesu <- terra::extract(wind1, st_transform(group, crs = st_crs(wind1)), ID = FALSE)[, 1]
    valuesv <- terra::extract(wind2, st_transform(group, crs = st_crs(wind1)), ID = FALSE)[, 1]

  }, silent = TRUE)
  print(valuesu)
  print(valuesv)
  extracted_values_u[[as.character(name)]] <- valuesu
  extracted_values_v[[as.character(name)]] <- valuesv
}
#Finally take the extracted data in each group and add them back to the main dataframe
#As the == is a pass between two lubridate values they are matched, but in different timezones.
#begin with u component
for (name in names(extracted_values_u)) {
  name_posix <- as.POSIXct(as.numeric(name, origin="1970-01-01", tz="GMT"))
  
  data$wind_u_target[grouped_data$round_DateTimeGMT == name_posix] <- extracted_values_u[[name]]
}
#do this for the v component
for (name in names(extracted_values_v)) {
  name_posix <- as.POSIXct(as.numeric(name, origin="1970-01-01", tz="GMT"))
  
  data$wind_v_target[grouped_data$round_DateTimeGMT == name_posix] <- extracted_values_v[[name]]
}
#calculate windspeed and direction with u and v components
#adapted from https://sgichuki.github.io/Atmo/ 
#Create two functions to calculate the wind direction and wind speed from u,v 
# axis. Sonic anemometer gives wind speed outputs as +ve or -ve speeds along the
# U axis, V axis and W(vertical) axis 
windDir <-function(u,v){
  (270-atan2(u,v)*180/pi)%%360 
}
windSpd <-function(u,v){
  sqrt(u^2+v^2)
}
#apply to data
data$winddirec_target <-windDir(data$wind_u_target,data$wind_v_target)
data$windspeed_target <-windSpd(data$wind_u_target,data$wind_v_target)

###Save RDS File as backup
saveRDS(data, "data.rds")

### Finally, add step lengths
#use MoveTools, adapted from Becker Movement Analysis Class, Fall 2024 Montana State.
proj_new <- "EPSG:5072"  # this is a good one for Wyoming/Montana/Idaho (Albers Equal Area)
st_crs(data)
data <- st_transform(data, proj_new)   #this is the function to reproject an sf object
st_crs(data)
st_crs(data)$proj4string
rm(proj_new)

# You must order your database by id and date
data <- data %>% 
  arrange(AnimalID,DateTimeLocal)

# Calculate bursts for your dataset
# A burst vector is based on relocations of one or multiple animals given a max time between intervals
# It is based on two points connected in time, where row i is connected with row i+1.
# When calculating movement parameters (e.g., step length, turning angle, speed), only sequential points 
# in the same burst vector are considered.
# Example with 1 hr relocation data and multiple ids: 
# set Tmax to 3 hours, and thus two points are connected (in the same burst)
# even if there is a missed fix in between them (i.e., a 6 hr step).
data$burst <- CalcBurst(data=data, id = TRUE, id_name="AnimalID", 
                        date_name="DateTimeLocal", Tmax = 3600*3) #set Tmax to a little more than double the fix rate
length(unique(data$burst))

data <- CalcMovParams(data=data, id_name = "AnimalID", 
                      date_name = "DateTimeLocal", burst=data$burst)

