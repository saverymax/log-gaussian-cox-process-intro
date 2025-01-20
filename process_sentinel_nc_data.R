# Script to process sentinel netcdf files for later use in a spatial-temporal model
library(ncdf4)
library(ncmeta)
library(stars)
library(terra)
library(ggplot2)
library(tidyterra)
library(lubridate)
library(gganimate)

#TODO: 1-14-2025. Oh, this isn't easy
nl_map <- st_read(file.path("D:", "data", "maps", "netherlands_bestuurlijkegrenzen_2021", "clean_nl_boundary.gpkg"))
plot(nl_map)
# Transform to coords
nl_map_prj <- st_transform(nl_map, crs="EPSG:4326")
plot(nl_map_prj)
nl_extent <- ext(nl_map_prj)
plot(nl_extent)
sentinel_directory <- file.path("D:", "data", "air_quality_data", 
                                "sentinel_download", "small_subset_download")
file.exists(sentinel_directory)
sentinel_files <- list.files(sentinel_directory, pattern='*.nc', full.names=TRUE) 
sentinel_files

# First just test everything loading 1 file.
sentinel_data <- nc_open(sentinel_files[2])
names(sentinel_data$var)
no2 <- ncvar_get(sentinel_data, var="PRODUCT/nitrogendioxide_tropospheric_column")
print(no2)
lat <- ncvar_get(sentinel_data, "PRODUCT/latitude")
lon <- ncvar_get(sentinel_data, "PRODUCT/longitude")
time_var <- ncvar_get(sentinel_data, "PRODUCT/time_utc") 
dim(time_var)
date <- as.Date(as.POSIXct(time_var, origin = "1970-01-01")) 
date
# Check that they match?
ymd(date[1])==ymd(date[length(date)])
nc_close(sentinel_data)

no2_spatraster <- rast(t(no2), crs="EPSG:4326")
# Set the extent of the SpatRaster
lon_values <- lon[!is.na(lon)]
lat_values <- lat[!is.na(lat)]
# Work with full data, not subsetted
#ext(no2_spatraster) <- ext(min(lon), max(lon), min(lat), max(lat))
# For subset
ext(no2_spatraster) <- ext(min(lon_values), max(lon_values), min(lat_values), max(lat_values))
# There is some intersection problem going on. Not all the files actually 
# have NL in them, since I used a big bounding box in the download
terra::intersect(no2_spatraster, nl_extent)
# And then once we crop, we aren't gaurenteed to get the data.
no2_spatraster <- terra::crop(no2_spatraster, nl_map_prj)

no2_spatraster <- extend(no2_spatraster, nl_extent)
no2_spatraster
all(is.na(values(no2_spatraster)))
# Convert the SpatRaster to a stars object
no2_stars <- st_as_stars(no2_spatraster)

# Plot the stars object
plot(no2_stars, 
     main = "Tropospheric NO2 Column", 
     col = viridis::viridis(100)) 

no2_stars
ggplot() +
  geom_stars(data=no2_stars) +
  #geom_spatraster(data=no2_spatraster) +
  geom_sf(data = nl_map_prj, color=alpha("black",0.9),linewidth=1.5, fill='transparent') +
  scale_fill_viridis_c(begin=0.2, end=1,alpha=0.7, na.value="transparent", name = "NO2\n(mol/m²)") + 
  coord_sf() +  # Use coord_sf for map projection
  labs(title = "Tropospheric NO2 Column",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal()

# Loop to process the files
spatraster_list <- list()
for(i in seq_along(sentinel_files)) {
  print(i)
  print(length(spatraster_list))
  sentinel_data <- nc_open(sentinel_files[i])
  # Get the data
  no2 <- ncvar_get(sentinel_data, var="PRODUCT/nitrogendioxide_tropospheric_column")
  no2_spatraster <- rast(t(no2), crs = "EPSG:4326")
  # Use the nl extent to set the raster extent
  # Use this line with the full data, and ther other one for 
  # data subset using the tool.
  #ext(no2_spatraster) <- ext(min(lon), max(lon), min(lat), max(lat))
  # Need to extract non-na values
  lon_values <- lon[!is.na(lon)]
  lat_values <- lat[!is.na(lat)]
  print(res(no2_spatraster))
  # Then set the extent based on those
  ext(no2_spatraster) <- ext(min(lon_values), max(lon_values), min(lat_values), max(lat_values))
  print(res(no2_spatraster))
  # There are a a few oddly large resolution images, so 
  # arbitrarily filter them out.
  if (res(no2_spatraster)[1] > 0.5){
    next
    print("Skipping due to res")
  }
  # Crop based on NL boundary
  no2_spatraster <- crop(no2_spatraster, nl_map_prj)
  print(res(no2_spatraster))
  # But this leaves the extents still not perfectly matching
  # Force them here. This might add new pixels, in which case
  # they will be NA
  no2_spatraster <- extend(no2_spatraster, nl_extent, fill=NA)
  # Not sure if I need the above extend if I do the below ext()
  # This one is just for extra force
  ext(no2_spatraster) <- nl_extent
  print(res(no2_spatraster))
  # Finally, force the crs just in case
  crs(no2_spatraster) <- crs(nl_extent)
  # There were still some extent problems 
  # so resample is the final step.
  if(i != 1){
    no2_spatraster <- resample(no2_spatraster, spatraster_list[[1]])
  }
  # Skip over any full NA rasters
  if (!all(is.na(values(no2_spatraster)))){
    print("Map and data intersect")
    # Process the time
    time_var <- ncvar_get(sentinel_data, "PRODUCT/time_utc") 
    date_vec <- as.Date(as.POSIXct(time_var, origin = "1970-01-01")) 
    stopifnot(ymd(date[1])==ymd(date[length(date)]))
    date <- date_vec[1]
    nc_close(sentinel_data)
    
    values_vec <- values(no2_spatraster)
    # Get the indices of NA values
    na_indices <- which(is.na(values_vec))
    values(no2_spatraster)[na_indices] <- mean(values_vec, na.rm=T)
    
    # Filter based on resolution
        
    p <- ggplot() +
      geom_spatraster(data=no2_spatraster) +
      geom_sf(data = nl_map_prj, color=alpha("black",0.9),linewidth=1.5, fill='transparent') +
      scale_fill_viridis_c(begin=0.2, end=1,alpha=0.7, na.value="transparent", name = "NO2\n(mol/m²)") + 
      coord_sf() +  # Use coord_sf for map projection
      labs(title = "Tropospheric NO2 Column",
           x = "Longitude",
           y = "Latitude") +
      theme_minimal()
    print(p)
    # Add the date as an attribute to the SpatRaster
    attr(no2_spatraster, "time") <- date
    spatraster_list[[i]] <- no2_spatraster
  }else{
    print(paste("Skipping file", sentinel_files[i]))
  }
}

# List gets populated with null values for every iteration where something is skipped
spatraster_list_filtered <- Filter(Negate(is.null), spatraster_list)

ggplot() +
  geom_spatraster(data=spatraster_list_filtered[[1]]) +
  geom_sf(data = nl_map_prj, color=alpha("black",0.9),linewidth=1.5, fill='transparent') +
  scale_fill_viridis_c(begin=0.2, end=1,alpha=0.7, na.value="transparent", name = "NO2\n(mol/m²)") + 
  coord_sf() +  # Use coord_sf for map projection
  labs(title = "Tropospheric NO2 Column",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal()

# Check inheritance, resolutionm, time attribute
for (i in seq_along(spatraster_list_filtered)){
  print(i)
  stopifnot(class(spatraster_list_filtered[[i]])=="SpatRaster")
  print(res(spatraster_list_filtered[[i]]))
  #print(st_as_stars(spatraster_list[[i]]))
  print(attributes(spatraster_list_filtered[[i]])$time)
  print(inherits(attributes(spatraster_list_filtered[[i]])$time, "Date"))
}

no2_spatraster_cube <- rast(spatraster_list_filtered)
no2_spatraster_cube
# Convert to a stars cube
no2_stars_cube <- st_as_stars(no2_spatraster_cube)
print(no2_stars_cube)

time_values <- lapply(spatraster_list, function(t) attributes(t)$time)
time_vector <- do.call(c, time_values)
time_vector

# Set the time attribute 
no2_stars_cube <- st_set_dimensions(no2_stars_cube, 3, 
                                    values = time_vector, 
                                    names = "time")

time_attr <- as.Date(st_get_dimension_values(no2_stars_cube, "time"))
time_attr
duplicated_days <- time_attr[duplicated(time_attr)]
print(duplicated_days) 
# Test out the duplicated day indexing
dup_indices <- which(time_attr == "2023-01-20")
dup_indices
dup_indices[1:(length(dup_indices)-1)]
dim(no2_stars_cube)
no2_stars_cube[,,, -dup_indices[1:(length(dup_indices)-1)]] 
dim(no2_stars_cube[,,, -dup_indices[1:(length(dup_indices)-1)]])[3]
st_get_dimension_values(no2_stars_cube[,,, -dup_indices[1:(length(dup_indices)-1)]],"time")
date_occurences <- st_get_dimension_values(no2_stars_cube[,,, -dup_indices[1:(length(dup_indices)-1)]],"time")
length(which(date_occurences=="2023-01-20"))

# If there are duplicated days, handle them
if (length(duplicated_days) > 0) {
  for (d in seq_along(duplicated_days)) {
    day <- duplicated_days[d]
    # Find indices of rasters with the duplicated day
    dup_indices <- which(time_attr == day)
    print(day)
    print(dup_indices)
    # Code to replace duplicates with average
    # Average the rasters with the duplicated day
    #avg_raster <- mean(no2_stars_cube[,,, dup_indices])
    # Replace the duplicated rasters with the averaged raster in the cube
    #no2_stars_cube[,,, dup_indices[1]] <- avg_raster 
    # But instead, let's just remove the duplicates
    print(-dup_indices[1:(length(dup_indices)-1)])
    #print(no2_stars_cube[,,, -dup_indices[1:(length(dup_indices)-1)]] )
    # Gets the elements you want to remove and does not return them 
    # with -dup_indices
    dim_before <- dim(no2_stars_cube)[3]
    print("dim before")
    print(dim_before)
    no2_stars_cube <- no2_stars_cube[,,, -dup_indices[1:(length(dup_indices)-1)]] 
    # Need to reset time dimension every replacement
    time_attr <- as.Date(st_get_dimension_values(no2_stars_cube, "time"))
    no2_stars_cube <- st_set_dimensions(no2_stars_cube, 3, values = time_attr, names = "time")
    dim_after <- dim(no2_stars_cube)[3]
    print("dim after")
    print(dim_after)
    stopifnot(dim_after < dim_before)
    date_occurences <- st_get_dimension_values(no2_stars_cube,"time")
    print(date_occurences)
    stopifnot(length(which(date_occurences==day))==1)
  }
  
  # Update the time attribute after removing duplicates
  #time_attr <- as.Date(st_get_dimension_values(no2_stars_cube, "time"))
  #no2_stars_cube <- st_set_dimensions(no2_stars_cube, 3, values = time_attr, names = "time")
}

plot(no2_stars_cube)
# Create the animated plot
p <- ggplot() +
  geom_stars(data = no2_stars_cube) +
  geom_sf(data = nl_map_prj, color=alpha("black",0.9),linewidth=1.5, fill='transparent') +
  scale_fill_viridis_c(option = "viridis", name = "NO2 (mol/m²)") +
  coord_sf() +
  #labs(title = "Tropospheric NO2 Column",  # Dynamic title
  labs(title = "Tropospheric NO2 Column - {time}",  # Dynamic title
       x = "Longitude",
       y = "Latitude") +
  transition_states(
    time,
    transition_length = 2,
    state_length = 1
  ) +
  enter_fade() + 
  exit_shrink() +
  ease_aes('linear')
animate(p, renderer = gifski_renderer())

save_dir <- file.path("D:", "data", "air_quality_data",  
                      "sentinel_download",  "processed_no2_nl.nc")

# I guess this works to set the crs of the whole thing
#no2_stars_cube <- st_set_crs(no2_stars_cube, st_crs(nl_map_prj)) 

raster_2_save <- rast(no2_stars_cube)
plot(raster_2_save)
writeCDF(raster_2_save, save_dir, overwrite=T)
# gdalcubes is another option.
#library(gdalcubes)
#write_ncdf(
#  raster_2_save,
#  fname = save_dir,
#  overwrite = T,
#)

