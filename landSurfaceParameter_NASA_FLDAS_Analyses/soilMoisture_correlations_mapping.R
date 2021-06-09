## install packages needed
##install.packages("raster")
##install.packages("ncdf4")
##install.packages("fields")
##install.packages("maptools")
##install.packages("sf")


### load required packages
library(raster)
library(ncdf4)
library(fields)
library(maps)
library(zoo)
library(maptools)
library(sf)

#nc files from https://disc.gsfc.nasa.gov/datasets/FLDAS_NOAH01_C_GL_M_001/summary?keywords=FLDAS and in Box

#extract_nc() function defined in functions_created.R script used to obtain list of values for interested variable
#get_cor_mat() function defined in functions_created.R script used to create matrix of correlations between two variables of interest


############## use function for netcdf file for every month#####################

#create vector of file names - easiest if file names follow a pattern
files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(1982:2020, each = 12), 
               c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"), ".001.nc.SUB.nc4", sep = "")
files <- c(files, paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(2021, 3), 
                        c("01", "02"), ".001.nc.SUB.nc4", sep = ""))

## getting lat/lon coords
f <- nc_open("ges_disc/FLDAS_NOAH01_C_GL_M.A198201.001.nc.SUB.nc4")
lon <- ncvar_get(f, "X")
lat <- ncvar_get(f, "Y")
lat_rng <- c(-1.5, 0.5)
lon_rng <- c(36.5, 38.5)
lat_ind <- which(lat >= lat_rng[1] & lat <= lat_rng[2])
lon_ind <- which(lon >= lon_rng[1] & lon <= lon_rng[2])
start_lat <- min(lat_ind)
start_lon <- min(lon_ind)
count_lat <- length(lat_ind)
count_lon <- length(lon_ind)

lon_x <- ncvar_get(f, "X", start = start_lon, count = count_lon) #if need
lat_y <- ncvar_get(f, "Y", start = start_lat, count = count_lat) #if need

nc_close(f)

### for smc
smc_list <- extract_nc(files, "SoilMoi00_10cm_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
### for rainfall flux
rain_list <- extract_nc(files, "Rainf_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
### for evapotranspiration
evap_list <- extract_nc(files, "Evap_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
### for surface runoff
surfrun_list <- extract_nc(files, "Qs_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
### subsurface runoff
subrun_list <- extract_nc(files, "Qsb_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
###
tavg_list <- extract_nc(files, "Tair_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))

### smc_10_40, 40_100, 100_200
smc_1040_list <- extract_nc(files, "SoilMoi10_40cm_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
smc_40100_list <- extract_nc(files, "SoilMoi40_100cm_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
smc_100200_list <- extract_nc(files, "SoilMoi100_200cm_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))

## get correlation matrix
cor_mat_smc_rain <- get_cor_mat(rain_list, smc_list) #use variables of interest as arguments

##############
#plot map
############

#Make color range
color_range <- colorRampPalette(c("red","orange", "yellow", "green"))
my_colors <- color_range(10000)

shape <- readShapeSpatial("tana_outline.shp") #Upper Tana River Basin shape


#Plot image
image.plot(lon_x, lat_y, cor_mat_smc_rain, main = "Soil Moisture (0-10cm) and Rainfall Flux Correlation", xlab = "Longitude", ylab = "Latitude", horizontal = TRUE, zlim = c(-1, 1), col = rainbow(200000))
maps::map(add = T)
plot(shape, lwd = 5, border = "black", col = NA, add = TRUE)



