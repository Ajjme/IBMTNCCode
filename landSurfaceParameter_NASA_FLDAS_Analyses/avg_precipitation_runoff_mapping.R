#create vector of file names - easiest if file names follow a pattern
#### Rainy Season months: 3, 4, 5, 10, 11, 12, 1
#### Dry Season months: 2, 6, 7, 8, 9

#files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(2000:2004, each = 7), c("01", "03", "04", "05", "10", "11", "12"), ".001.nc.SUB.nc4", sep = "")

#nc files from https://disc.gsfc.nasa.gov/datasets/FLDAS_NOAH01_C_GL_M_001/summary?keywords=FLDAS and in Box

###
#extract_nc() function defined in functions_created.R script used to obtain list of values for interested variable
#get_avg_mat() function defined in functions_created.R script used to obtain matrix of average values for interested variable over period of time

####Average values for rainfall over Upper Tana (January 1982 - March 2021)
avg_rain_list <- get_avg_mat(rain_list) #rain list obtained using extract_nc() function defined in smc_correlations.R script
avg_rain_list_mmhr <- avg_rain_list*3600
avg_rain_list_mmday <- avg_rain_list*86400

#Create Color Range
color_range <- colorRampPalette(c("cornsilk","darkseagreen1", "cadetblue1", "cyan2", "deepskyblue2", "blue", "dodgerblue4"))
my_colors <- color_range(10000)

##### plot image overall average rainfall
image.plot(lon_x, lat_y, avg_rain_list_mmday, col = my_colors, main = "1982 - 2021 Average Precipitation Rate (mm/day)", xlab = "Longitude", ylab = "Latitude", zlim = c(0, 10))
plot(shape, lwd = 5, border = "black", add = TRUE)


#################################
#####RAINY SEASON##########
######################################3

#1985-1989
files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(1985:1989, each = 7), 
               c("01", "03", "04", "05", "10", "11", "12"), ".001.nc.SUB.nc4", sep = "")
rainy_8589 <- extract_nc(files, "Rainf_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
avg_rainy_8589 <- get_avg_mat(rainy_8589) * 86400
image.plot(lon_x, lat_y, avg_rainy_8589, col = my_colors, main = "1985 - 1989 Average Rainy Season Precipitation Rate (mm/day)", xlab = "Longitude", ylab = "Latitude", zlim = c(0, 10), horizontal = TRUE)
plot(shape, lwd = 5, border = "black", add = TRUE)

#1990-1994
files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(1990:1994, each = 7), 
               c("01", "03", "04", "05", "10", "11", "12"), ".001.nc.SUB.nc4", sep = "")
rainy_9094 <- extract_nc(files, "Rainf_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
avg_rainy_9094 <- get_avg_mat(rainy_9094) * 86400
image.plot(lon_x, lat_y, avg_rainy_9094, col = my_colors, main = "1990 - 1994 Average Rainy Season Precipitation Rate (mm/day)", xlab = "Longitude", ylab = "Latitude", zlim = c(0, 10))
plot(shape, lwd = 5, border = "black", add = TRUE)

#1995-1999
files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(1995:1999, each = 7), 
               c("01", "03", "04", "05", "10", "11", "12"), ".001.nc.SUB.nc4", sep = "")
rainy_9599 <- extract_nc(files, "Rainf_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
avg_rainy_9599 <- get_avg_mat(rainy_9599) * 86400
image.plot(lon_x, lat_y, avg_rainy_9599, col = my_colors, main = "1995 - 1999 Average Rainy Season Precipitation Rate (mm/day)", xlab = "Longitude", ylab = "Latitude", zlim = c(0, 10))
plot(shape, lwd = 5, border = "black", add = TRUE)

#2000-2004
files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(2000:2004, each = 7), 
               c("01", "03", "04", "05", "10", "11", "12"), ".001.nc.SUB.nc4", sep = "")
rainy_0004 <- extract_nc(files, "Rainf_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
avg_rainy_0004 <- get_avg_mat(rainy_0004) * 86400
image.plot(lon_x, lat_y, avg_rainy_0004, col = my_colors, main = "2000 - 2004 Average Rainy Season Precipitation Rate (mm/day)", xlab = "Longitude", ylab = "Latitude", zlim = c(0, 10))
plot(shape, lwd = 5, border = "black", add = TRUE)

#2005-2009
files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(2005:2009, each = 7), 
               c("01", "03", "04", "05", "10", "11", "12"), ".001.nc.SUB.nc4", sep = "")
rainy_0509 <- extract_nc(files, "Rainf_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
avg_rainy_0509 <- get_avg_mat(rainy_0509) * 86400
image.plot(lon_x, lat_y, avg_rainy_0509, col = my_colors, main = "2005 - 2009 Average Rainy Season Precipitation Rate (mm/day)", xlab = "Longitude", ylab = "Latitude", zlim = c(0, 10))
plot(shape, lwd = 5, border = "black", add = TRUE)

#2010-2014
files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(2010:2014, each = 7), 
               c("01", "03", "04", "05", "10", "11", "12"), ".001.nc.SUB.nc4", sep = "")
rainy_1014 <- extract_nc(files, "Rainf_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
avg_rainy_1014 <- get_avg_mat(rainy_1014) * 86400
image.plot(lon_x, lat_y, avg_rainy_1014, col = my_colors, main = "2010 - 2014 Average Rainy Season Precipitation Rate (mm/day)", xlab = "Longitude", ylab = "Latitude", zlim = c(0, 10))
plot(shape, lwd = 5, border = "black", add = TRUE)

#2015-2019
files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(2015:2019, each = 7), 
               c("01", "03", "04", "05", "10", "11", "12"), ".001.nc.SUB.nc4", sep = "")
rainy_1519 <- extract_nc(files, "Rainf_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
avg_rainy_1519 <- get_avg_mat(rainy_1519) * 86400
image.plot(lon_x, lat_y, avg_rainy_1519, col = my_colors, main = "2015 - 2019 Average Rainy Season Precipitation Rate (mm/day)", xlab = "Longitude", ylab = "Latitude", zlim = c(0, 10))
plot(shape, lwd = 5, border = "black", add = TRUE)

##########################
#Runoff
##########################

files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(2015:2020, each = 12), 
               c("01", "02", "03", "04", "05","06","07", "08","09", "10", "11", "12"), ".001.nc.SUB.nc4", sep = "")

#2015 - 2020
runoff1520 <- extract_nc(files, "Qs_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
avg_runoff1520 <- get_avg_mat(runoff1520)*86400

r_avg_runoff <- raster(avg_runoff1520)
extent(r_avg_runoff) <- c(36.5, 38.5, -1.5, 0.5)
projection(r_avg_runoff) <- CRS("+proj=longlat +datum=WGS84")
ggplot() + geom_raster(data = r_avg_runoff)

plot(r_avg_runoff)
image.plot(lon_x, lat_y, avg_runoff1520, col = my_colors, main = "2015-2020 Average Surface Runoff (mm/day)", xlab = "Longitude", ylab = "Latitude", zlim = c(0, 2.5), horizontal = TRUE)
plot(shape, lwd = 5, border = "black", add = TRUE)

#2000-2005
files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(2000:2005, each = 12), 
               c("01", "02", "03", "04", "05","06","07", "08","09", "10", "11", "12"), ".001.nc.SUB.nc4", sep = "")
runoff0005 <- extract_nc(files, "Qs_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
avg_runoff0005 <- get_avg_mat(runoff0005)*86400
image.plot(lon_x, lat_y, avg_runoff0005, col = my_colors, main = "2000-2005 Average Surface Runoff (mm/day)", xlab = "Longitude", ylab = "Latitude", zlim = c(0, 2.5))
plot(shape, lwd = 5, border = "black", add = TRUE)

################################################################
#####DRY SEASON##########
################################################################

#1985-1989
files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(1985:1989, each = 5), 
               c("02", "06", "07", "08", "09"), ".001.nc.SUB.nc4", sep = "")
dry_8589 <- extract_nc(files, "Rainf_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
avg_dry_8589 <- get_avg_mat(dry_8589) * 86400
image.plot(lon_x, lat_y, avg_dry_8589, col = my_colors, main = "1985 - 1989 Average Dry Season Precipitation Rate (mm/day)", xlab = "Longitude", ylab = "Latitude", zlim = c(0, 10))
plot(shape, lwd = 5, border = "black", add = TRUE)

#1990-1994
files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(1990:1994, each = 5), 
               c("02", "06", "07", "08", "09"), ".001.nc.SUB.nc4", sep = "")
dry_9094 <- extract_nc(files, "Rainf_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
avg_dry_9094 <- get_avg_mat(dry_9094) * 86400
image.plot(lon_x, lat_y, avg_dry_9094, col = my_colors, main = "1990 - 1994 Average Dry Season Precipitation Rate (mm/day)", xlab = "Longitude", ylab = "Latitude", zlim = c(0, 10))
plot(shape, lwd = 5, border = "black", add = TRUE)

#1995-1999
files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(1995:1999, each = 5), 
               c("02", "06", "07", "08", "09"), ".001.nc.SUB.nc4", sep = "")
dry_9599 <- extract_nc(files, "Rainf_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
avg_dry_9599 <- get_avg_mat(dry_9599) * 86400
image.plot(lon_x, lat_y, avg_dry_9599, col = my_colors, main = "1995 - 1999 Average Dry Season Precipitation Rate (mm/day)", xlab = "Longitude", ylab = "Latitude", zlim = c(0, 10))
plot(shape, lwd = 5, border = "black", add = TRUE)

#2000-2004
files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(2000:2004, each = 5), 
               c("02", "06", "07", "08", "09"), ".001.nc.SUB.nc4", sep = "")
dry_0004 <- extract_nc(files, "Rainf_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
avg_dry_0004 <- get_avg_mat(dry_0004) * 86400
image.plot(lon_x, lat_y, avg_dry_0004, col = my_colors, main = "2000 - 2004 Average Dry Season Precipitation Rate (mm/day)", xlab = "Longitude", ylab = "Latitude", zlim = c(0, 10))
plot(shape, lwd = 5, border = "black", add = TRUE)

#2005-2009
files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(2005:2009, each = 5), 
               c("02", "06", "07", "08", "09"), ".001.nc.SUB.nc4", sep = "")
dry_0509 <- extract_nc(files, "Rainf_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
avg_dry_0509 <- get_avg_mat(dry_0509) * 86400
image.plot(lon_x, lat_y, avg_dry_0509, col = my_colors, main = "2005 - 2009 Average Dry Season Precipitation Rate (mm/day)", xlab = "Longitude", ylab = "Latitude", zlim = c(0, 10))
plot(shape, lwd = 5, border = "black", add = TRUE)

#2010-2014
files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(2010:2014, each = 5), 
               c("02", "06", "07", "08", "09"), ".001.nc.SUB.nc4", sep = "")
dry_1014 <- extract_nc(files, "Rainf_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
avg_dry_1014 <- get_avg_mat(dry_1014) * 86400
image.plot(lon_x, lat_y, avg_dry_1014, col = my_colors, main = "2010 - 2014 Average Dry Season Precipitation Rate (mm/day)", xlab = "Longitude", ylab = "Latitude", zlim = c(0, 10))
plot(shape, lwd = 5, border = "black", add = TRUE)

#2015-2019
files <- paste("ges_disc/FLDAS_NOAH01_C_GL_M.A", rep(2015:2019, each = 5), 
               c("02", "06", "07", "08", "09"), ".001.nc.SUB.nc4", sep = "")
dry_1519 <- extract_nc(files, "Rainf_f_tavg", start = c(start_lon, start_lat, 1), count = c(count_lon, count_lat, 1))
avg_dry_1519 <- get_avg_mat(dry_1519) * 86400
image.plot(lon_x, lat_y, avg_dry_1519, col = my_colors, main = "2015 - 2019 Average Dry Season Precipitation Rate (mm/day)", xlab = "Longitude", ylab = "Latitude", zlim = c(0, 10))
plot(shape, lwd = 5, border = "black", add = TRUE)


color_range <- colorRampPalette(c("navy", "dodgerblue4", "dodgerblue3",  "cornflowerblue", "deepskyblue2", "cadetblue1", "skyblue2", "lightskyblue2", "lightblue1", "lightcyan2", "lightcyan1", "lightgoldenrodyellow", "white", "papayawhip", "peachpuff", "rosybrown1", "lightsalmon", "darksalmon", "lightcoral", "indianred2", "brown1", "brown2", "brown3", "brown4", "firebrick4"))
my_colors <- color_range(10000)
