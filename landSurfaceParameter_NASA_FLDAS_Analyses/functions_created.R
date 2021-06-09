### load required packages
library(raster)
library(ncdf4)
library(fields)


##########################################################################
################# Extract mean values using function to create csv########
##########################################################################

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


###Write the function### --- function also in script functions_created.R script
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


# function to get matrices of selected variable
extract_nc <- function(x, y, start = c(1, 1, 1), count = c(-1, -1, 1)) {
  # Args:
  # x: vector of nc files
  # y: name of variable
  # start = c(lon, lat, time)
  # count = c(lon, lat, time)
  # Return:
  # list of matrices for variable
  mat_vals <- list(numeric(0))
  for (i in seq_along(x)) {
    f <- nc_open(x[i])
    mat_vals[[i]] <- ncvar_get(f, y, start = start, count = count)
    nc_close(f)
  }
  mat_vals
}

# function to get correlations
get_cor_mat <- function(x, y) {
  # Args: 
  # x: list for first varible (obtained through extract_nc function)
  # y: list for second variable
  # Return:
  # Matrix of correlations
  list_vecs_x <- list(numeric(0))
  list_vecs_y <- list(numeric(0))
  for (i in seq_along(x[[1]])) {
    list_vecs_x[[i]] <- sapply(x, function(x) x[i])
    list_vecs_y[[i]] <- sapply(y, function(y) y[i])
  }
  cor_mat <- matrix(nrow = nrow(x[[1]]), ncol = ncol(x[[1]]))
  for (j in seq_along(list_vecs_x)) {
    cor_mat[j] <- cor(list_vecs_x[[j]], list_vecs_y[[j]])
  }
  cor_mat
}

#Function to get average values for variable
get_avg_mat <- function(x) {
  # Args:
  # x: list for variable
  # Return:
  # Matrix of average values
  list_vecs_x <- list(numeric(0))
  for (i in seq_along(x[[1]])) {
    list_vecs_x[[i]] <- sapply(x, function(x) x[i])
  }
  avg_mat <- matrix(nrow = nrow(x[[1]]), ncol = ncol(x[[1]]))
  for (j in seq_along(list_vecs_x)) {
    avg_mat[j] <- mean(list_vecs_x[[j]])
  }
  avg_mat
}

list_to_df <- function(x) {
  #Args:
  #x: list for variable
  #Return:
  #Dataframe with values, lat/long, date
  df_x <- data.frame(Year = rep(rep(1982:2020, each = 12), 400), Month = rep(1:12, 15600), Vals = NA, Longitude = rep(seq(from = 36.55, to = 38.45, by = 0.1), each = 9360), Latitude = rep(rep(seq(from = -1.45, to = 0.45, by = 0.1), each = 468), 20))
  num <- 468*(0:399) + 1
  for (i in seq_along(x[[1]])) {
    ind <- seq(num[i], num[i] + 467)
    df_x$Vals[ind] <- sapply(x, function(x) x[i])
  }
  df_x
}