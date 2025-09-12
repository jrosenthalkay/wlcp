
##
## clean GHS and berkeley earth, forming pattern scaling
##
library(readxl)
library(dplyr)
library(countrycode)
library(sf)
library(raster)
library(exactextractr)
library(data.table)
library(rnaturalearth)
library(ggplot2)

# paths
mypath <- "/Users/jordanrosenthalkay/Dropbox/Research/WLCC_replication/"

## read in data ## 

# population data
ghs <- file.path(mypath,'/data/raw/GHS_POP_E2015_GLOBE_R2023A_4326_30ss_V1_0.tif')

# Read GHS population raster
pop <- raster(ghs)

# temp data
berktemp <- read.csv(file.path(mypath,'data/raw/berktemp_1750-2019.csv'))
berktemp$berk_id <- 1:nrow(berktemp)
berktemp <- berktemp[!is.na(berktemp$temp_2000),]

berktemp.geo <- berktemp[,c('latitude','longitude','berk_id')]

# crosswalk between g-econ and berkeley temperature IDs
# Convert to sf objects
berktemp.sf <- st_as_sf(berktemp.geo, coords = c("longitude", "latitude"), crs = 4326)

# Transform to projected CRS for buffer in meters 
berktemp.sf.proj <- st_transform(berktemp.sf, 3857)
berktemp.buff <- st_buffer(berktemp.sf.proj, dist=111000) %>%
  st_transform(4326)

# Extract population sums within buffers
pop_sums <- exact_extract(pop, berktemp.buff, 'sum')
berktemp$population <- pop_sums

# get countries
world <- ne_countries(scale = "medium", returnclass = "sf")

# get countries
nearest <- st_nearest_feature(berktemp.sf, world)
berktemp$iso3 <- world$iso_a3_eh[nearest]
berktemp$country_name <- world$name[nearest]

# fill missing countries
all_countries <- world$iso_a3_eh
countries_with_temp <- unique(berktemp$iso3)
missing_countries <- setdiff(all_countries, countries_with_temp)

sort(unique(missing_countries))

# For each missing country, find the nearest berktemp point and duplicate it
new_rows <- list()
for(country in missing_countries) {
  # Get the country's geometry
  country_geom <- world[world$iso_a3_eh == country,]
  
  # Find the nearest berktemp point to this country's geometry
  nearest_idx <- st_nearest_feature(country_geom, berktemp.sf)
  
  # Get that row from berktemp
  new_row <- berktemp[nearest_idx,]
  
  # Update the country identifiers
  new_row$iso3 <- country
  new_row$country_name <- world$name[world$iso_a3_eh == country]
  
  new_rows[[country]] <- new_row
}

# Combine all new rows into a single dataframe
new_berktemp <- do.call(rbind, new_rows)

# Append to original berktemp
berktemp <- berktemp %>% rbind(new_berktemp)

# First, reshape from wide to long format
long_df <- melt(
  setDT(berktemp), 
  id.vars = c("berk_id", "population", "iso3", "country_name","longitude","latitude"),
  measure.vars = patterns("^temp_"),
  variable.name = "year",
  value.name = "temperature"
)

# Clean up year column to just keep the year number
long_df[, year := as.numeric(gsub("temp_", "", year))]

# Calculate population-weighted country averages by year
country_avg <- long_df[, .(
  weighted_temp = weighted.mean(temperature, w = population),
  weighted_lon = weighted.mean(longitude, w = population),
  weighted_lat = weighted.mean(latitude, w = population)
), by = .(country_name, iso3, year)]

write.csv(country_avg,file.path(mypath,'data/int/country_avg_temp_timeseries_ghs_popweight.csv'),row.names=F)
