## install packages needed
##install.packages("raster")
##install.packages("ncdf4")
##install.packages("fields")
##install.packages("rgdal")


### load required packages
library(raster)
library(ncdf4)
library(fields)
library(maps)
library(zoo)
library(rgdal)
library(sp)

#nc files from https://disc.gsfc.nasa.gov/datasets/FLDAS_NOAH01_C_GL_M_001/summary?keywords=FLDAS and in Box

##########Some useful functions####################3
# nc_open() opens the netcdf file
# if f <- nc_open("file.nc4"), the output of f will show the variables,dimensions,structure of netcdf file
# ncvar_get() extracts values from the netcdf file as an array
# args ncvar_get: 
######## nc = netcdf object; 
######## varid = name of variable to read data from; 
######## start = vector of indices for where to start reading values (order is X,Y,Z,T) - match to dimensions of netcdf file; start at 1 if unspecified
######## count = vector of intigers for ocunt of values to read along each dimension (default is all entries; also indicated by -1)
# nc_close closes the netcdf file


#########################################################################
########## IF want to make map for a single netcdf file#######################
#########################################################################

######Extract necessary information from NC FILES ##############

##open nc file

f <- nc_open("ges_disc/FLDAS_NOAH01_C_GL_M.A202103.001.nc.SUB.nc4")
#f <- nc_open("ges_disc/FLDAS_NOAH01_C_GL_M.A198212.001.nc.SUB.nc4")

## extract lat/long if needed (Region for East Africa: 30, -10, 55, 10); (Region for Upper Tana: 36.5, -1.5, 38.5, 0.5)
lon <- ncvar_get(f, "X")
lat <- ncvar_get(f, "Y")

lat_rng <- c(-10, 10)
lon_rng <- c(30, 55)
lat_ind <- which(lat >= lat_rng[1] & lat <= lat_rng[2])
lon_ind <- which(lon >= lon_rng[1] & lon <= lon_rng[2])
start_lat <- min(lat_ind)
start_lon <- min(lon_ind)
count_lat <- length(lat_ind)
count_lon <- length(lon_ind)


#create data frame of soil moisture (or whichever variable)
soilmois_2103_df <- ncvar_get(f, "SoilMoi00_10cm_tavg", start = c(1, 1, 1), count = c(-1, -1, 1))
soilmois_2103_df_eafrica <- ncvar_get(f, "SoilMoi40_100cm_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
rain_2103_df_ea <- ncvar_get(f, "Rainf_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
evap_2103_df_eaa <- ncvar_get(f, "Evap_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))


#create image
lon_x <- ncvar_get(f, "X", start = start_lon, count = count_lon) #if need
lat_y <- ncvar_get(f, "Y", start = start_lat, count = count_lat) #if need

image.plot(lon_x, lat_y, soilmois_8212_df_eafrica, main = "Soil Moisture 40-100cm (m3 m-3) March 2021", xlab= "Longitude", ylab = "Latitude", col = my_colors, horizontal = TRUE)
map(add = T)

image.plot(lon_x, lat_y, evap_198201_df_eaa, main = "Evapotranspiration (kg m-2 s-1) March 2021", xlab= "Longitude", ylab = "Latitude", col = my_colors, horizontal = TRUE)
map(add = T)
plot(shape, border = "red", lwd = 2, add = TRUE)

########
image.plot(lon_x, lat_y, soilmois_8212_df, main = "Soil Moisture 40-100cm (m3 m-3) March 2021", xlab= "Longitude", ylab = "Latitude", col = my_colors, horizontal = TRUE)
raster(soilmois_8212_df)

#############

image.plot(lon_x, lat_y, rain_1982_df_ea, main = "Rainfall Flux (kg m-2 s-1) March 2021", xlab = "Longitude", ylab = "Latitude", col = my_colors, horizontal = TRUE)
map(add = T)
plot(shape, border = "red", lwd = 2, add = TRUE)

color_range <- colorRampPalette(c("orange", "yellow", "green", "darkcyan", "blue"))
my_colors <- color_range(100)

# close nc file
nc_close(f)


##########################################################################################
#FOR CREATING CSV of MEANS FOR EACH VARIABLE FOR THE REGION --- January 1982 - March 2021
##########################################################################################

#############################################################
################# Extract mean values using function ########
#############################################################

#### create empty vectors for variable means interested in###

total_evap <- numeric(0) #total evapotranspiration (kg m-2 s-1)
spec_hum <- numeric(0) #specific humidity (kg kg-1)
surf_runoff <- numeric(0) #surface runoff (kg m-2 s-1)
sub_runoff <- numeric(0) #subsurface runoff (kg m-2 s-1)
rain_flux <- numeric(0) #rainfall flux (kg m-2 s-1)
smc00_10 <- numeric(0) #soil moisture content 0 - 10 cm underground (m^3 m-3)
smc10_40 <- numeric(0) #soil moisture content 10 - 40 cm underground (m^3 m-3)
smc40_100 <- numeric(0) #soil moisture content 40 - 100 cm underground (m^3 m-3)
smc100_200 <- numeric(0) #soil moisture content 100 - 200 cm underground (m^3 m-3)
tas <- numeric(0) #surface air temperature (K)


###Write the function###
get_means <- function(x, start = c(1, 1, 1), count = c(-1, -1, 1)) {
  #Args:
  #x: vector of netcdf files
  #Output:
  #vector of monthly means for the entire region for each variable specified in function (same as empty vectors above)
  for (i in seq_along(x)) { #the function will go through all files input as the vector x
    f <- nc_open(x[i]) 
    #extract all values for variables interested in: 
    evap_df <- ncvar_get(f, "Evap_tavg", start = start, count = count) 
    hum_df <- ncvar_get(f, "Qair_f_tavg", start = start, count = count)
    surf_df <- ncvar_get(f, "Qs_tavg", start = start, count = count)
    sub_df <- ncvar_get(f, "Qsb_tavg", start = start, count = count)
    rain_df <- ncvar_get(f, "Rainf_f_tavg", start = start, count = count)
    smc0010_df <- ncvar_get(f, "SoilMoi00_10cm_tavg", start = start, count = count)
    smc1040_df <- ncvar_get(f, "SoilMoi10_40cm_tavg", start = start, count = count)
    smc40100_df <- ncvar_get(f, "SoilMoi40_100cm_tavg", start = start, count = count)
    smc100200_df <- ncvar_get(f, "SoilMoi100_200cm_tavg", start = start, count = count)
    tas_df <- ncvar_get(f, "Tair_f_tavg", start = start, count = count)
    #get mean value for each variable:
    mean_evap_df <- mean(evap_df, na.rm = T)
    mean_hum_df <- mean(hum_df, na.rm = T)
    mean_surf_df <- mean(surf_df, na.rm = T)
    mean_sub_df <- mean(sub_df, na.rm = T)
    mean_rain_df <- mean(rain_df, na.rm = T)
    mean_smc0010_df <- mean(smc0010_df, na.rm = T)
    mean_smc1040_df <- mean(smc1040_df, na.rm = T)
    mean_smc40100_df <- mean(smc40100_df, na.rm = T)
    mean_smc100200_df <- mean(smc100200_df, na.rm = T)
    mean_tas_df <- mean(tas_df, na.rm = T)
    #add on mean values to vectors (created outside of the function) with means for each variable:
    total_evap <<- c(total_evap, mean_evap_df)
    spec_hum <<- c(spec_hum, mean_hum_df)
    surf_runoff <<- c(surf_runoff, mean_surf_df)
    sub_runoff <<- c(sub_runoff, mean_sub_df)
    rain_flux <<- c(rain_flux, mean_rain_df)
    smc00_10 <<- c(smc00_10, mean_smc0010_df)
    smc10_40 <<- c(smc10_40, mean_smc1040_df)
    smc40_100 <<- c(smc40_100, mean_smc40100_df)
    smc100_200 <<- c(smc100_200, mean_smc100200_df)
    tas <<- c(tas, mean_tas_df)
    nc_close(f)
  }  
}

########### get lat/lon for tana############
lat_rng <- c(-1.5, 0.5)
lon_rng <- c(36.5, 38.5)
lat_ind <- which(lat >= lat_rng[1] & lat <= lat_rng[2])
lon_ind <- which(lon >= lon_rng[1] & lon <= lon_rng[2])
start_lat <- min(lat_ind)
start_lon <- min(lon_ind)
count_lat <- length(lat_ind)
count_lon <- length(lon_ind)


############## use function for netcdf file for every month#####################

#create vector of file names - easiest if file names follow a pattern
files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(1982:2020, each = 12), 
               c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"), ".001.nc.SUB.nc4", sep = "")
files <- c(files, paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(2021, 3), 
                        c("01", "02", "03"), ".001.nc.SUB.nc4", sep = ""))

#create vectors with means
get_means(files, start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))

##############create dataframe of date, year, month, variable#################

#Create date vector by month-year

#year, month vectors
year <- c(rep(1982:2020, each = 12), rep(2021, 3))
month <- c(rep(1:12, 39), 1, 2, 3)

#date vector - by first of each month
dates <- as.Date(paste(rep(1982:2020, each = 12), "-", 1:12, "-", 1, sep = ""))
dates <- c(dates, paste(rep(2021, 3), "-", c(1, 2, 3), "-", 1, sep = ""))

#for date by days(or weeks, months, years) between two dates
#dates <- seq(as.Date("2020-01-01"), as.Date("2021-01-01"), by = "days") 

# create data frame
ges_disc_df <- cbind(dates, year, month, 
                     total_evap, spec_hum, surf_runoff, sub_runoff,rain_flux, 
                     smc00_10, smc10_40, smc40_100, smc100_200, tas)
colnames(ges_disc_df) <- c("Date", "Year", "Month", "total evapotranspiration (kg m-2 s-1)", 
                           "specific humidity (kg kg-1)", "surface runoff (kg m-2 s-1)", 
                           "subsurface runoff (kg m-2 s-1)", "rainfall flux (kg m-2 s-1)", 
                           "soil moisture content 0 - 10 cm underground (m^3 m-3)", 
                           "soil moisture content 10 - 40 cm underground (m^3 m-3)", 
                           "soil moisture content 40 - 100 cm underground (m^3 m-3)", 
                           "soil moisture content 100 - 200 cm underground (m^3 m-3)", 
                           "surface air temperature (K)")

## save csv
write.csv(ges_disc_df, "ges_disc_df.csv")
